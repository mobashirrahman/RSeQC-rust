//! Port of `read_distribution.py`: summarize read distribution across
//! genomic annotation categories (CDS/UTR/intron/upstream/downstream).
//! Contract: see `compatibility/commands.yaml` entry
//! `read_distribution.py`; algorithm ported from
//! `oracle/upstream-src/scripts/read_distribution.py` directly (it builds
//! its own `RegionModel`/`DistributionCounts` rather than delegating to a
//! qcmodule.SAM method).
//!
//! `overlaps_point`: upstream calls `Intersecter.find(position, position)`
//! -- a literal zero-width query. Despite the docstring's "whether a 1-bp
//! interval at position overlaps a region", bx-python's actual overlap
//! test for a zero-width query is `interval.start < position <
//! interval.end` -- STRICT on both ends, not half-open containment.
//! Confirmed empirically via `PYTHONPATH=oracle/upstream-src/src
//! oracle/venv/bin/python3 -c "from bx.intervals.intersection import
//! Intersecter, Interval; ..."`: `Intersecter().add_interval(Interval(100,
//! 200)).find(100, 100)` returns `[]` (the start boundary itself is
//! excluded), while `find(101, 101)`..`find(199, 199)` all match. A
//! genuinely surprising consequence: a width-1 region (`Interval(100,
//! 101)`) can NEVER match any point query, since no integer is strictly
//! between 100 and 101 -- reproduced here exactly, not "fixed", since this
//! command's whole precedence cascade depends on it.
//!
//! SAM-text and CRAM input (DIV-0002/0004) are supported at the CLI layer via
//! `rseqc_formats::open_alignments`; this module's own functions were
//! unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`).

use std::collections::HashMap;
use std::io;

use noodles_bam as bam;
use noodles_sam as sam;
use rseqc_formats::bed;
use rseqc_formats::cigar::fetch_exon_blocks;
use rseqc_formats::interval::{self, Bed3};

/// Per-chromosome, pre-merged (via `interval::union_bed3`/`subtract_bed3`)
/// non-overlapping region set, queryable for point containment.
struct RegionSet {
    by_chrom: HashMap<String, Vec<(i64, i64)>>,
    bases: i64,
}

impl RegionSet {
    fn new(intervals: Vec<Bed3>) -> Self {
        let bases: i64 = intervals.iter().map(|(_, s, e)| e - s).sum();
        let mut by_chrom: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        // Matches upstream's `build_interval_trees`, which keys on
        // `str(interval[0]).upper()` -- the BED file's chrom names
        // aren't uppercased anywhere else, only here and at the
        // per-read query site in `count_read_distribution`.
        for (chrom, s, e) in intervals {
            by_chrom.entry(chrom.to_uppercase()).or_default().push((s, e));
        }
        for ivs in by_chrom.values_mut() {
            ivs.sort();
        }
        Self { by_chrom, bases }
    }

    fn overlaps_point(&self, chrom: &str, position: i64) -> bool {
        let Some(ivs) = self.by_chrom.get(chrom) else { return false };
        // Intervals here are already merged (non-overlapping, sorted), so
        // binary search for the interval starting at or before `position`.
        match ivs.binary_search_by(|&(s, _)| s.cmp(&position)) {
            // `position == s` exactly: never overlaps (strict `position >
            // s` below always fails here), regardless of interval width --
            // see module docs.
            Ok(_) => false,
            Err(0) => false,
            Err(idx) => {
                let (s, e) = ivs[idx - 1];
                position > s && position < e
            }
        }
    }
}

pub struct RegionModel {
    cds_exon: RegionSet,
    intron: RegionSet,
    utr_5: RegionSet,
    utr_3: RegionSet,
    upstream_1kb: RegionSet,
    upstream_5kb: RegionSet,
    upstream_10kb: RegionSet,
    downstream_1kb: RegionSet,
    downstream_5kb: RegionSet,
    downstream_10kb: RegionSet,
}

fn purify(regions: Vec<Bed3>, cds_exon: &[Bed3], utr_5: &[Bed3], utr_3: &[Bed3], intron: &[Bed3]) -> Vec<Bed3> {
    let r = interval::subtract_bed3(&regions, cds_exon);
    let r = interval::subtract_bed3(&r, utr_5);
    let r = interval::subtract_bed3(&r, utr_3);
    interval::subtract_bed3(&r, intron)
}

/// Builds the annotation region model from a BED12 gene-model file.
/// Re-opens the file once per extraction call, matching upstream's
/// `ParseBED` (each `get*` method reads through and seeks back to 0).
pub fn process_gene_model(bed_path: &std::path::Path) -> io::Result<RegionModel> {
    let read = || -> io::Result<std::io::BufReader<std::fs::File>> {
        Ok(std::io::BufReader::new(std::fs::File::open(bed_path)?))
    };

    let utr_3 = interval::union_bed3(&bed::get_utr(read()?, 3)?);
    let utr_5 = interval::union_bed3(&bed::get_utr(read()?, 5)?);
    let cds_exon = interval::union_bed3(&bed::get_cds_exon(read()?)?);
    let (intron_raw, _skipped) = bed::get_intron(read()?)?;
    let intron = interval::union_bed3(&intron_raw);

    let utr_5 = interval::subtract_bed3(&utr_5, &cds_exon);
    let utr_3 = interval::subtract_bed3(&utr_3, &cds_exon);
    let intron = interval::subtract_bed3(&intron, &cds_exon);
    let intron = interval::subtract_bed3(&intron, &utr_5);
    let intron = interval::subtract_bed3(&intron, &utr_3);

    let upstream_1kb = interval::union_bed3(&bed::get_intergenic(read()?, "up", 1_000)?);
    let upstream_5kb = interval::union_bed3(&bed::get_intergenic(read()?, "up", 5_000)?);
    let upstream_10kb = interval::union_bed3(&bed::get_intergenic(read()?, "up", 10_000)?);
    let downstream_1kb = interval::union_bed3(&bed::get_intergenic(read()?, "down", 1_000)?);
    let downstream_5kb = interval::union_bed3(&bed::get_intergenic(read()?, "down", 5_000)?);
    let downstream_10kb = interval::union_bed3(&bed::get_intergenic(read()?, "down", 10_000)?);

    let upstream_1kb = purify(upstream_1kb, &cds_exon, &utr_5, &utr_3, &intron);
    let upstream_5kb = purify(upstream_5kb, &cds_exon, &utr_5, &utr_3, &intron);
    let upstream_10kb = purify(upstream_10kb, &cds_exon, &utr_5, &utr_3, &intron);
    let downstream_1kb = purify(downstream_1kb, &cds_exon, &utr_5, &utr_3, &intron);
    let downstream_5kb = purify(downstream_5kb, &cds_exon, &utr_5, &utr_3, &intron);
    let downstream_10kb = purify(downstream_10kb, &cds_exon, &utr_5, &utr_3, &intron);

    Ok(RegionModel {
        cds_exon: RegionSet::new(cds_exon),
        intron: RegionSet::new(intron),
        utr_5: RegionSet::new(utr_5),
        utr_3: RegionSet::new(utr_3),
        upstream_1kb: RegionSet::new(upstream_1kb),
        upstream_5kb: RegionSet::new(upstream_5kb),
        upstream_10kb: RegionSet::new(upstream_10kb),
        downstream_1kb: RegionSet::new(downstream_1kb),
        downstream_5kb: RegionSet::new(downstream_5kb),
        downstream_10kb: RegionSet::new(downstream_10kb),
    })
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct DistributionCounts {
    pub total_reads: u64,
    pub total_tags: u64,
    pub unassigned_tags: u64,
    pub cds_exon: u64,
    pub intron: u64,
    pub utr_5: u64,
    pub utr_3: u64,
    pub upstream_1kb: u64,
    pub upstream_5kb: u64,
    pub upstream_10kb: u64,
    pub downstream_1kb: u64,
    pub downstream_5kb: u64,
    pub downstream_10kb: u64,
    pub qc_failed_reads: u64,
    pub duplicate_reads: u64,
    pub secondary_reads: u64,
    pub unmapped_reads: u64,
}

/// Matches `assign_tag()`'s exact category precedence order.
fn assign_tag(chrom: &str, midpoint: i64, model: &RegionModel, counts: &mut DistributionCounts) {
    if model.cds_exon.overlaps_point(chrom, midpoint) {
        counts.cds_exon += 1;
        return;
    }

    let in_5utr = model.utr_5.overlaps_point(chrom, midpoint);
    let in_3utr = model.utr_3.overlaps_point(chrom, midpoint);

    if in_5utr && !in_3utr {
        counts.utr_5 += 1;
        return;
    }
    if in_3utr && !in_5utr {
        counts.utr_3 += 1;
        return;
    }
    if in_5utr && in_3utr {
        counts.unassigned_tags += 1;
        return;
    }

    if model.intron.overlaps_point(chrom, midpoint) {
        counts.intron += 1;
        return;
    }

    let up_10kb = model.upstream_10kb.overlaps_point(chrom, midpoint);
    let down_10kb = model.downstream_10kb.overlaps_point(chrom, midpoint);

    if up_10kb && down_10kb {
        counts.unassigned_tags += 1;
        return;
    }

    if model.upstream_1kb.overlaps_point(chrom, midpoint) {
        counts.upstream_1kb += 1;
        counts.upstream_5kb += 1;
        counts.upstream_10kb += 1;
        return;
    }
    if model.upstream_5kb.overlaps_point(chrom, midpoint) {
        counts.upstream_5kb += 1;
        counts.upstream_10kb += 1;
        return;
    }
    if up_10kb {
        counts.upstream_10kb += 1;
        return;
    }

    if model.downstream_1kb.overlaps_point(chrom, midpoint) {
        counts.downstream_1kb += 1;
        counts.downstream_5kb += 1;
        counts.downstream_10kb += 1;
        return;
    }
    if model.downstream_5kb.overlaps_point(chrom, midpoint) {
        counts.downstream_5kb += 1;
        counts.downstream_10kb += 1;
        return;
    }
    if down_10kb {
        counts.downstream_10kb += 1;
        return;
    }

    counts.unassigned_tags += 1;
}

pub fn count_read_distribution<I>(records: I, header: &sam::Header, model: &RegionModel) -> io::Result<DistributionCounts>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut counts = DistributionCounts::default();

    for result in records {
        let record = result?;
        let flags = record.flags();

        if flags.is_qc_fail() {
            counts.qc_failed_reads += 1;
            continue;
        }
        if flags.is_duplicate() {
            counts.duplicate_reads += 1;
            continue;
        }
        if flags.is_secondary() {
            counts.secondary_reads += 1;
            continue;
        }
        if flags.is_unmapped() {
            counts.unmapped_reads += 1;
            continue;
        }

        counts.total_reads += 1;

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom.to_string().to_uppercase();

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let read_start = (pos.get() - 1) as i64;

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let exons = fetch_exon_blocks(read_start as usize, ops);
        counts.total_tags += exons.len() as u64;

        for (start, end) in exons {
            let (start, end) = (start as i64, end as i64);
            let midpoint = start + (end - start) / 2;
            assign_tag(&chrom, midpoint, model, &mut counts);
        }
    }

    Ok(counts)
}

fn tags_per_kb(tag_count: u64, base_count: i64) -> f64 {
    (tag_count as f64) * 1000.0 / (base_count as f64 + 1.0)
}

fn print_row(label: &str, bases: i64, tags: u64) -> String {
    format!("{:<20}{:<20}{:<20}{:<18.2}", label, bases, tags, tags_per_kb(tags, bases))
}

pub fn render_report(model: &RegionModel, counts: &DistributionCounts) -> String {
    let mut lines = Vec::new();
    lines.push(format!("{:<30}{}", "Total Reads", counts.total_reads));
    lines.push(format!("{:<30}{}", "Total Tags", counts.total_tags));
    lines.push(format!(
        "{:<30}{}",
        "Total Assigned Tags",
        counts.total_tags - counts.unassigned_tags
    ));
    lines.push("=".repeat(69));
    lines.push(format!("{:<20}{:<20}{:<20}{:<20}", "Group", "Total_bases", "Tag_count", "Tags/Kb"));

    lines.push(print_row("CDS_Exons", model.cds_exon.bases, counts.cds_exon));
    lines.push(print_row("5'UTR_Exons", model.utr_5.bases, counts.utr_5));
    lines.push(print_row("3'UTR_Exons", model.utr_3.bases, counts.utr_3));
    lines.push(print_row("Introns", model.intron.bases, counts.intron));
    lines.push(print_row("TSS_up_1kb", model.upstream_1kb.bases, counts.upstream_1kb));
    lines.push(print_row("TSS_up_5kb", model.upstream_5kb.bases, counts.upstream_5kb));
    lines.push(print_row("TSS_up_10kb", model.upstream_10kb.bases, counts.upstream_10kb));
    lines.push(print_row("TES_down_1kb", model.downstream_1kb.bases, counts.downstream_1kb));
    lines.push(print_row("TES_down_5kb", model.downstream_5kb.bases, counts.downstream_5kb));
    lines.push(print_row("TES_down_10kb", model.downstream_10kb.bases, counts.downstream_10kb));

    lines.push("=".repeat(69));
    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region_set(intervals: Vec<Bed3>) -> RegionSet {
        RegionSet::new(interval::union_bed3(&intervals))
    }

    #[test]
    fn overlaps_point_strict_both_ends() {
        // Queries use the uppercased form, matching how
        // count_read_distribution() actually calls this (see the
        // uppercasing regression test below for why).
        //
        // Verified against real bx-python: `Intersecter().add_interval(
        // Interval(100, 200)).find(position, position)` for each position
        // below, via `PYTHONPATH=oracle/upstream-src/src oracle/venv/
        // bin/python3 -c`. The start boundary itself (100) does NOT
        // overlap -- bx-python's zero-width point query is strict on
        // both ends, not half-open containment.
        let rs = region_set(vec![("chr1".to_string(), 100, 200)]);
        assert!(!rs.overlaps_point("CHR1", 99));
        assert!(!rs.overlaps_point("CHR1", 100));
        assert!(rs.overlaps_point("CHR1", 101));
        assert!(rs.overlaps_point("CHR1", 199));
        assert!(!rs.overlaps_point("CHR1", 200));
        assert!(!rs.overlaps_point("CHR2", 150));
    }

    #[test]
    fn overlaps_point_width_one_interval_never_matches() {
        // A genuinely surprising bx-python consequence: no integer is
        // strictly between 100 and 101, so a single-base region can never
        // match any point query. Verified via the same real bx-python
        // probe as above.
        let rs = region_set(vec![("chr1".to_string(), 100, 101)]);
        assert!(!rs.overlaps_point("CHR1", 99));
        assert!(!rs.overlaps_point("CHR1", 100));
        assert!(!rs.overlaps_point("CHR1", 101));
    }

    #[test]
    fn region_set_uppercases_chrom_matching_query_side() {
        // Regression test: build_interval_trees() uppercases chrom names
        // upstream, and so does the per-read query in
        // count_read_distribution() -- a BED file's lowercase "chr1" must
        // still match a query for "CHR1" (both sides get uppercased).
        let rs = region_set(vec![("chr1".to_string(), 100, 200)]);
        assert!(rs.overlaps_point("CHR1", 150));
    }

    #[test]
    fn assign_tag_precedence_cds_wins_over_everything() {
        let model = RegionModel {
            cds_exon: region_set(vec![("chr1".to_string(), 0, 100)]),
            intron: region_set(vec![("chr1".to_string(), 0, 100)]),
            utr_5: region_set(vec![("chr1".to_string(), 0, 100)]),
            utr_3: region_set(vec![]),
            upstream_1kb: region_set(vec![]),
            upstream_5kb: region_set(vec![]),
            upstream_10kb: region_set(vec![]),
            downstream_1kb: region_set(vec![]),
            downstream_5kb: region_set(vec![]),
            downstream_10kb: region_set(vec![]),
        };
        let mut counts = DistributionCounts::default();
        assign_tag("CHR1", 50, &model, &mut counts);
        assert_eq!(counts.cds_exon, 1);
        assert_eq!(counts.intron, 0);
        assert_eq!(counts.utr_5, 0);
    }

    #[test]
    fn assign_tag_ambiguous_utr_is_unassigned() {
        let model = RegionModel {
            cds_exon: region_set(vec![]),
            intron: region_set(vec![]),
            utr_5: region_set(vec![("chr1".to_string(), 0, 100)]),
            utr_3: region_set(vec![("chr1".to_string(), 0, 100)]),
            upstream_1kb: region_set(vec![]),
            upstream_5kb: region_set(vec![]),
            upstream_10kb: region_set(vec![]),
            downstream_1kb: region_set(vec![]),
            downstream_5kb: region_set(vec![]),
            downstream_10kb: region_set(vec![]),
        };
        let mut counts = DistributionCounts::default();
        assign_tag("CHR1", 50, &model, &mut counts);
        assert_eq!(counts.unassigned_tags, 1);
    }

    #[test]
    fn print_row_matches_python_format_spec() {
        // Independently derived from oracle/upstream-src/scripts/
        // read_distribution.py's print_row f-string via python3 -c.
        let row = print_row("CDS_Exons", 1000, 50);
        assert_eq!(row, format!("{:<20}{:<20}{:<20}{:<18.2}", "CDS_Exons", 1000, 50, 49.95));
    }
}
