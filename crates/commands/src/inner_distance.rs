//! Port of `inner_distance.py`: estimate the inner distance between
//! paired-end RNA-seq reads. Contract: see `compatibility/commands.yaml`
//! entry `inner_distance.py`; algorithm ported from `mRNA_inner_distance()`
//! in `oracle/upstream-src/src/qcmodule/SAM.py` (lines 3578-3750). The most
//! involved command ported so far -- see module-level notes on each
//! preserved upstream quirk below.
//!
//! SAM-text input (DIV-0002/0004) is supported at the CLI layer via
//! `rseqc_formats::open_alignments`; this module's own functions were
//! unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`). Note `open_alignments`
//! decodes the whole file eagerly, so the `sample_size` early-exit no
//! longer avoids reading records past the cap for large inputs -- a
//! documented, accepted tradeoff (see `open_alignments`'s own doc
//! comment), not a correctness gap: the computed result is identical
//! either way.

use std::collections::{HashMap, HashSet};
use std::io;

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::cigar::op::Kind};
use rseqc_formats::bed;
use rseqc_formats::cigar::{fetch_exon_blocks, fetch_intron_blocks};
use rseqc_formats::interval::{Bed3, MergedRegions};

/// Point-overlap query over RAW (non-merged) named transcript intervals --
/// distinct from `MergedRegions`, which discards interval identity/names.
/// Ported from the `transcript_ranges` dict of `Intersecter`s built in
/// `mRNA_inner_distance()`.
///
/// **Preserves a real upstream bug**: the building loop (lines 3612-3618)
/// only calls `.add_interval(...)` in the `else` branch (chromosome
/// already seen before), never on a chromosome's first occurrence -- so
/// the FIRST transcript per chromosome is silently dropped from the index.
struct TranscriptIndex {
    by_chrom: HashMap<String, Vec<(i64, i64, String)>>,
}

impl TranscriptIndex {
    fn build(ranges: Vec<bed::TranscriptRange>) -> Self {
        let mut by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
        let mut seen: HashSet<String> = HashSet::new();
        for tr in ranges {
            let chrom = tr.chrom.to_uppercase();
            if seen.contains(&chrom) {
                by_chrom.entry(chrom).or_default().push((tr.tx_start, tr.tx_end, tr.name));
            } else {
                seen.insert(chrom);
            }
        }
        Self { by_chrom }
    }

    /// Distinct transcript names of every stored interval overlapping the
    /// half-open `[start, end)` query.
    fn names_at(&self, chrom: &str, start: i64, end: i64) -> HashSet<String> {
        let mut out = HashSet::new();
        if let Some(ivs) = self.by_chrom.get(chrom) {
            for (s, e, name) in ivs {
                if *s < end && start < *e {
                    out.insert(name.clone());
                }
            }
        }
        out
    }
}

pub struct InnerDistanceModel {
    exon_regions: MergedRegions,
    exon_chroms: HashSet<String>,
    transcripts: TranscriptIndex,
}

pub fn build_model(bed_path: &std::path::Path) -> io::Result<InnerDistanceModel> {
    let open = || -> io::Result<std::io::BufReader<std::fs::File>> {
        Ok(std::io::BufReader::new(std::fs::File::open(bed_path)?))
    };

    let raw_exons = bed::get_exon(open()?)?;
    // Upstream: `ref_exons.append([exn[0].upper(), exn[1], exn[2]])`.
    let exons: Vec<Bed3> = raw_exons.into_iter().map(|(c, s, e)| (c.to_uppercase(), s, e)).collect();
    let exon_chroms: HashSet<String> = exons.iter().map(|(c, _, _)| c.clone()).collect();
    let exon_regions = MergedRegions::new(&exons);

    let (transcript_ranges, _skipped) = bed::get_transcript_ranges(open()?)?;
    let transcripts = TranscriptIndex::build(transcript_ranges);

    Ok(InnerDistanceModel { exon_regions, exon_chroms, transcripts })
}

/// One `.inner_distance.txt` row, pre-formatted (the "unknownChromosome"
/// branch upstream is missing its trailing newline -- a real upstream
/// inconsistency, preserved verbatim rather than normalized: see the
/// `line` field's contents rather than joining rows with an added `\n`).
#[derive(Debug)]
pub struct DistanceRecord {
    pub line: String,
    pub histogram_value: Option<i64>,
}

#[derive(Debug, Default)]
pub struct DistanceResult {
    pub pair_num: u64,
    pub records: Vec<DistanceRecord>,
    /// `true` when the input iterator ran out on its own (upstream's
    /// `except StopIteration`, which prints "Done" to stderr); `false`
    /// when `sample_size` was reached first (upstream's `break`, no
    /// "Done"). See `crates/cli/src/bin/inner_distance.rs`.
    pub loop_exhausted_naturally: bool,
}

pub fn compute_inner_distance<I>(
    records: I,
    header: &sam::Header,
    model: &InnerDistanceModel,
    q_cut: u8,
    sample_size: u64,
) -> io::Result<DistanceResult>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut pair_num = 0u64;
    let mut out = Vec::new();
    let mut loop_exhausted_naturally = true;

    for result in records {
        if pair_num >= sample_size {
            loop_exhausted_naturally = false;
            break;
        }
        let record = result?;
        let flags = record.flags();

        if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            continue;
        }
        if !flags.is_segmented() {
            continue;
        }
        if flags.is_mate_unmapped() {
            continue;
        }
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < q_cut {
            continue;
        }

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        // qlen = pysam's query_alignment_length: sum of M/I/=/X (excludes clips).
        let qlen: i64 = ops
            .iter()
            .filter(|op| matches!(op.kind(), Kind::Match | Kind::Insertion | Kind::SequenceMatch | Kind::SequenceMismatch))
            .map(|op| op.len() as i64)
            .sum();

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let read1_start = (pos.get() - 1) as i64;
        let Some(mate_pos) = record.mate_alignment_start().transpose()? else { continue };
        let read2_start = (mate_pos.get() - 1) as i64;

        if read2_start < read1_start {
            continue; // BAM is coordinate-sorted; mate already processed.
        }
        if read2_start == read1_start && flags.is_first_segment() {
            continue; // Upstream sets inner_distance=0 then `continue`s -- not counted at all.
        }

        pair_num += 1;

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some(mate_ref_id) = record.mate_reference_sequence_id().transpose()? else { continue };
        let qname = record.name().map(|n| n.to_string()).unwrap_or_default();

        if ref_id != mate_ref_id {
            out.push(DistanceRecord {
                line: format!("{qname}\tNA\tsameChrom=No\n"),
                histogram_value: None,
            });
            continue;
        }

        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string().to_uppercase();

        let intron_blocks = fetch_intron_blocks(read1_start as usize, ops.iter().copied());
        let splice_intron_size: i64 = intron_blocks.iter().map(|&(s, e)| (e - s) as i64).sum();
        let read1_end = read1_start + qlen + splice_intron_size;

        let inner_distance: i64 = if read2_start >= read1_end {
            read2_start - read1_end
        } else {
            // Overlapping mates: count this read's own exon-derived 1-based
            // positions [es+1, ee] that fall in (read2_start, read1_end].
            let exon_blocks = fetch_exon_blocks(read1_start as usize, ops.iter().copied());
            let mut count = 0i64;
            for (es, ee) in exon_blocks {
                let (es, ee) = (es as i64, ee as i64);
                let lo = (es + 1).max(read2_start + 1);
                let hi = ee.min(read1_end);
                if hi >= lo {
                    count += hi - lo + 1;
                }
            }
            -count
        };

        let read1_names = model.transcripts.names_at(&chrom, read1_end - 1, read1_end);
        let read2_names = model.transcripts.names_at(&chrom, read2_start, read2_start + 1);
        let same_transcript = read1_names.intersection(&read2_names).next().is_some();

        if !same_transcript {
            out.push(DistanceRecord {
                line: format!("{qname}\t{inner_distance}\tsameTranscript=No,dist=genomic\n"),
                histogram_value: Some(inner_distance),
            });
            continue;
        }

        if inner_distance > 0 {
            if model.exon_chroms.contains(&chrom) {
                let size = model.exon_regions.overlap_length(&chrom, read1_end, read2_start);
                if size == inner_distance {
                    out.push(DistanceRecord {
                        line: format!("{qname}\t{size}\tsameTranscript=Yes,sameExon=Yes,dist=mRNA\n"),
                        histogram_value: Some(size),
                    });
                } else if size > 0 && size < inner_distance {
                    out.push(DistanceRecord {
                        line: format!("{qname}\t{size}\tsameTranscript=Yes,sameExon=No,dist=mRNA\n"),
                        histogram_value: Some(size),
                    });
                } else {
                    out.push(DistanceRecord {
                        line: format!("{qname}\t{inner_distance}\tsameTranscript=Yes,nonExonic=Yes,dist=genomic\n"),
                        histogram_value: Some(inner_distance),
                    });
                }
            } else {
                out.push(DistanceRecord {
                    // Upstream's literal string has no trailing '\n' here -- preserved.
                    line: format!("{qname}\t{inner_distance}\tunknownChromosome,dist=genomic"),
                    histogram_value: Some(inner_distance),
                });
            }
        } else {
            out.push(DistanceRecord {
                line: format!("{qname}\t{inner_distance}\treadPairOverlap\n"),
                histogram_value: Some(inner_distance),
            });
        }
    }

    Ok(DistanceResult { pair_num, records: out, loop_exhausted_naturally })
}

pub fn render_distance_file(result: &DistanceResult) -> String {
    result.records.iter().map(|r| r.line.as_str()).collect()
}

/// `window_left_bound = range(low_bound, up_bound, step)` (Python's
/// exclusive-end `range`), and per-bucket counts using the same flipped
/// half-open convention as upstream's `Intersecter.find(st, st+step)`
/// against `(v-1, v)` singleton intervals: `st < v <= st+step`, NOT the
/// more obvious `st <= v < st+step`.
pub fn histogram_buckets(values: &[i64], low_bound: i64, up_bound: i64, step: i64) -> Vec<(i64, i64, u64)> {
    let mut st = low_bound;
    let mut buckets = Vec::new();
    while st < up_bound {
        let count = values.iter().filter(|&&v| v > st && v <= st + step).count() as u64;
        buckets.push((st, st + step, count));
        st += step;
    }
    buckets
}

pub fn render_freq_table(buckets: &[(i64, i64, u64)]) -> String {
    buckets.iter().map(|(a, b, c)| format!("{a}\t{b}\t{c}\n")).collect()
}

pub fn render_r_script(buckets: &[(i64, i64, u64)], step: i64, out_prefix: &str) -> String {
    let sizes: Vec<String> = buckets.iter().map(|(st, _, _)| (*st as f64 + step as f64 / 2.0).to_string()).collect();
    let counts: Vec<String> = buckets.iter().map(|(_, _, c)| c.to_string()).collect();

    let mut lines = Vec::new();
    lines.push(format!("out_file = '{out_prefix}'"));
    lines.push(format!("pdf('{out_prefix}.inner_distance_plot.pdf')"));
    lines.push(format!("fragsize=rep(c({}),times=c({}))", sizes.join(","), counts.join(",")));
    lines.push("frag_sd = sd(fragsize)".to_string());
    lines.push("frag_mean = mean(fragsize)".to_string());
    lines.push("frag_median = median(fragsize)".to_string());
    lines.push("write(x=c(\"Name\",\"Mean\",\"Median\",\"sd\"), sep=\"\t\", file=stdout(),ncolumns=4)".to_string());
    lines.push("write(c(out_file,frag_mean,frag_median,frag_sd),sep=\"\t\", file=stdout(),ncolumns=4)".to_string());
    lines.push(format!(
        "hist(fragsize,probability=T,breaks={},xlab=\"mRNA insert size (bp)\",main=paste(c(\"Mean=\",frag_mean,\";\",\"SD=\",frag_sd),collapse=\"\"),border=\"blue\")",
        buckets.len()
    ));
    lines.push(format!("lines(density(fragsize,bw={}),col='red')", 2 * step));
    lines.push("dev.off()".to_string());
    // Trailing newline: upstream's plain `print(...)` calls each add
    // their own trailing newline, including the final `dev.off()`.
    format!("{}\n", lines.join("\n"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transcript_index_drops_first_transcript_per_chromosome() {
        let ranges = vec![
            bed::TranscriptRange { chrom: "chr1".into(), tx_start: 0, tx_end: 100, strand: "+".into(), name: "first".into() },
            bed::TranscriptRange { chrom: "chr1".into(), tx_start: 200, tx_end: 300, strand: "+".into(), name: "second".into() },
        ];
        let idx = TranscriptIndex::build(ranges);
        // "first" (chr1's first occurrence) is dropped; only "second" is queryable.
        assert!(idx.names_at("CHR1", 50, 51).is_empty());
        assert_eq!(idx.names_at("CHR1", 250, 251), HashSet::from(["second".to_string()]));
    }

    #[test]
    fn histogram_bucket_boundary_is_right_inclusive() {
        // Bucket [0,5): a value of exactly 5 belongs to THIS bucket (v <= st+step),
        // not the next one, and a value of exactly 0 does NOT belong here (v > st).
        let buckets = histogram_buckets(&[0, 5, 3], 0, 10, 5);
        assert_eq!(buckets[0], (0, 5, 2)); // values 5 and 3 (0 excluded)
        assert_eq!(buckets[1], (5, 10, 0)); // value 5 already counted in bucket 0
    }

    #[test]
    fn render_r_script_exact_text() {
        // Independently derived via python3 -c from oracle/upstream-src/
        // src/qcmodule/SAM.py lines 3734-3746.
        let buckets = vec![(-250, -245, 0u64), (-245, -240, 0), (-240, -235, 0)];
        let output = render_r_script(&buckets, 5, "testprefix");
        let expected = "out_file = 'testprefix'\n\
pdf('testprefix.inner_distance_plot.pdf')\n\
fragsize=rep(c(-247.5,-242.5,-237.5),times=c(0,0,0))\n\
frag_sd = sd(fragsize)\n\
frag_mean = mean(fragsize)\n\
frag_median = median(fragsize)\n\
write(x=c(\"Name\",\"Mean\",\"Median\",\"sd\"), sep=\"\t\", file=stdout(),ncolumns=4)\n\
write(c(out_file,frag_mean,frag_median,frag_sd),sep=\"\t\", file=stdout(),ncolumns=4)\n\
hist(fragsize,probability=T,breaks=3,xlab=\"mRNA insert size (bp)\",main=paste(c(\"Mean=\",frag_mean,\";\",\"SD=\",frag_sd),collapse=\"\"),border=\"blue\")\n\
lines(density(fragsize,bw=10),col='red')\n\
dev.off()\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn render_distance_file_preserves_missing_newline_on_unknown_chromosome() {
        let result = DistanceResult {
            pair_num: 2,
            records: vec![
                DistanceRecord { line: "r1\t100\tsameTranscript=No,dist=genomic\n".to_string(), histogram_value: Some(100) },
                DistanceRecord { line: "r2\t50\tunknownChromosome,dist=genomic".to_string(), histogram_value: Some(50) },
            ],
            loop_exhausted_naturally: true,
        };
        let text = render_distance_file(&result);
        assert_eq!(text, "r1\t100\tsameTranscript=No,dist=genomic\nr2\t50\tunknownChromosome,dist=genomic");
    }
}
