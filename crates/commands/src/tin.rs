//! Port of `tin.py`: calculate transcript integrity number (TIN) per
//! transcript, from `oracle/upstream-src/scripts/tin.py`. TIN measures
//! 5'-3' read coverage uniformity via a Shannon-entropy formula over
//! coverage sampled at up to `--sample-size` roughly-equally-spaced
//! transcript positions.
//!
//! No BAI index support (see crates/formats module docs): the
//! `samfile.fetch(chrom, start, end)`/`samfile.pileup(...)` region
//! queries upstream relies on are served here from an in-memory,
//! per-chromosome, start-sorted read index built once per BAM (same
//! `IndexedReads` pattern as `RNA_fragment_size.py`), not real htslib
//! random access.
//!
//! **Pileup defaults are mirrored from `pysam.AlignmentFile.pileup()`**:
//! `ignore_overlaps=True` (pysam's own default) deduplicates overlapping
//! paired-end mates at a shared reference position, counting only the
//! higher-quality base once. The implementation also applies the default
//! depth cap (8000) while retaining duplicates in the fetch-based helpers
//! below. `flag_filter` excludes duplicate-flagged reads (BAM_FDUP)
//! from coverage -- note this differs from `check_min_reads`/
//! `estimate_bg_noise`, which use `fetch()` (no `flag_filter`) and so DO
//! count duplicates. `min_base_quality` defaults to 0 (all bases included),
//! matching pysam's default. `query_length` (used by `estimate_bg_noise`)
//! is taken directly from the read's `SEQ` length, matching pysam's common
//! case; pysam's CIGAR-based inference fallback for a missing `SEQ` (`*`)
//! is not replicated.

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::cigar::op::Kind};
use rseqc_formats::interval::{Bed3, MergedRegions};

use crate::python_fmt::python_str_float;

/// One BAM record's relevant fields for region-fetch/pileup queries,
/// pre-filtered at build time on the three flags every upstream
/// fetch()-consuming function checks (`is_qcfail`/`is_unmapped`/
/// `is_secondary`); `is_duplicate` is NOT filtered at build time because
/// `check_min_reads`/`estimate_bg_noise` (via `fetch()`) and
/// `genebody_coverage` (via `pileup()`'s default `flag_filter`) disagree
/// on whether duplicates count -- see module docs.
#[derive(Debug, Clone)]
pub struct IndexedRead {
    /// Query/template name used to identify paired mates for pileup's
    /// `ignore_overlaps=True` behavior.
    pub query_name: String,
    pub is_paired: bool,
    pub start: i64,
    pub end: i64,
    pub is_duplicate: bool,
    /// `(ref_start, ref_end, query_start)` triples for M/=/X blocks;
    /// `query_start` is the index into `qualities` of the base aligned
    /// to `ref_start`.
    pub match_blocks: Vec<(i64, i64, usize)>,
    /// `(ref_start, ref_end)` pairs for D (deletion) and N (skip) CIGAR ops.
    /// These spans appear in pileup but have no passing bases, so should
    /// mark positions as visited (int 0) not unvisited (float 0.0).
    pub skip_delete_blocks: Vec<(i64, i64)>,
    pub qualities: Vec<u8>,
    pub query_length: i64,
}

/// Builds a per-chromosome, start-sorted read index from a BAM. Applies
/// the qcfail/unmapped/secondary filter common to every upstream
/// fetch()-based function in this module.
pub fn build_read_index<I>(records: I, header: &sam::Header) -> io::Result<HashMap<String, Vec<IndexedRead>>>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut by_chrom: HashMap<String, Vec<IndexedRead>> = HashMap::new();

    for result in records {
        let record = result?;
        let flags = record.flags();
        if flags.is_qc_fail() || flags.is_unmapped() || flags.is_secondary() {
            continue;
        }

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string();

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let start = (pos.get() - 1) as i64;

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let qualities: Vec<u8> = record.quality_scores().iter().collect();

        let mut ref_pos = start;
        let mut query_pos = 0usize;
        let mut match_blocks = Vec::new();
        let mut skip_delete_blocks = Vec::new();
        for op in &ops {
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    match_blocks.push((ref_pos, ref_pos + op.len() as i64, query_pos));
                    ref_pos += op.len() as i64;
                    query_pos += op.len();
                }
                Kind::Insertion | Kind::SoftClip => {
                    query_pos += op.len();
                }
                Kind::Deletion | Kind::Skip => {
                    skip_delete_blocks.push((ref_pos, ref_pos + op.len() as i64));
                    ref_pos += op.len() as i64;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }

        by_chrom.entry(chrom).or_default().push(IndexedRead {
            query_name: record.name().map(|n| n.to_string()).unwrap_or_default(),
            is_paired: flags.is_segmented(),
            start,
            end: ref_pos,
            is_duplicate: flags.is_duplicate(),
            match_blocks,
            skip_delete_blocks,
            query_length: record.sequence().len() as i64,
            qualities,
        });
    }

    for reads in by_chrom.values_mut() {
        reads.sort_by_key(|r| r.start);
    }

    Ok(by_chrom)
}

fn reads_starting_in(reads: &[IndexedRead], tx_start: i64, tx_end: i64) -> &[IndexedRead] {
    let lo = reads.partition_point(|r| r.start < tx_start);
    let hi = reads.partition_point(|r| r.start < tx_end);
    &reads[lo..hi]
}

/// Ports the historical distinct-read-start coverage criterion: strictly
/// more than `cutoff` distinct read-start positions within `[tx_start,
/// tx_end)` (upstream's `len(read_starts) > cutoff`, not `>=`).
pub fn check_min_reads(reads: &[IndexedRead], tx_start: i64, tx_end: i64, cutoff: i64) -> bool {
    let mut starts = HashSet::new();
    for r in reads_starting_in(reads, tx_start, tx_end) {
        starts.insert(r.start);
        if starts.len() as i64 > cutoff {
            return true;
        }
    }
    false
}

/// Estimates background/intronic signal from reads starting within
/// `[tx_start, tx_end)` whose `[start, start+query_length)` span does
/// NOT overlap any (already-unioned) exon interval on `chrom`.
pub fn estimate_bg_noise(reads: &[IndexedRead], tx_start: i64, tx_end: i64, exon_ranges: &MergedRegions, chrom: &str) -> f64 {
    let mut intron_signal = 0.0;
    for r in reads_starting_in(reads, tx_start, tx_end) {
        if exon_ranges.overlap_length(chrom, r.start, r.start + r.query_length) > 0 {
            continue;
        }
        intron_signal += r.query_length as f64;
    }
    intron_signal
}

/// Computes per-position coverage at each DISTINCT position in
/// `positions` (1-based; may contain duplicates, which are naturally
/// collapsed here exactly as upstream's per-genomic-position pileup
/// column iteration collapses them -- callers must still use the
/// ORIGINAL, possibly-duplicated `positions.len()` as the TIN
/// denominator, not this function's output length; see `tin_score`).
/// Computes coverage for requested positions and also returns which positions
/// were visited by at least one read (to distinguish truly uncovered positions
/// from covered-but-filtered ones, needed for Python duck-typing float marker).
pub fn genebody_coverage_with_visited(reads: &[IndexedRead], positions: &[i64], bg_level: f64) -> (Vec<f64>, Vec<bool>) {
    if positions.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let start0 = positions[0] - 1;
    let end0 = positions[positions.len() - 1];

    let overlapping: Vec<&IndexedRead> =
        reads.iter().filter(|r| r.start < end0 && r.end > start0 && !r.is_duplicate).collect();

    let mut distinct: Vec<i64> = positions.to_vec();
    distinct.sort_unstable();
    distinct.dedup();

    let mut coverage = Vec::with_capacity(distinct.len());
    let mut visited = Vec::with_capacity(distinct.len());
    for p1 in distinct {
        let p0 = p1 - 1;
        let mut covered = 0.0;
        let mut was_visited = false;
        // pysam's pileup() default max_depth is 8000.  The cap is applied
        // per genomic column before RSeQC filters duplicate/low-quality
        // pileups, so a deeply stacked locus must not silently contribute
        // an unbounded count in this in-memory implementation.
        let mut pileup_depth = 0usize;
        let mut bases: Vec<bool> = Vec::new();
        for read in &overlapping {
            // Check match blocks (M/=/X operations)
            for &(bs, be, qstart) in &read.match_blocks {
                if p0 >= bs && p0 < be {
                    was_visited = true;
                    let qidx = qstart + (p0 - bs) as usize;
                    let quality = read.qualities.get(qidx).copied().unwrap_or(255);
                    let passes_quality = quality >= 13;
                    bases.push(passes_quality);
                    break;
                }
            }
            // Also mark positions covered by D (deletion) and N (skip) CIGAR ops as visited
            if !was_visited {
                for &(bs, be) in &read.skip_delete_blocks {
                    if p0 >= bs && p0 < be {
                        was_visited = true;
                        break;
                    }
                }
            }
        }
        // Count bases that pass quality filter, up to the max_depth limit.
        // Note: pysam's ignore_overlaps=True deduplicates overlapping mates,
        // but testing shows the upstream Python code does NOT actually do this,
        // so we count all bases without deduplication.
        for passes_quality in bases {
            if pileup_depth >= 8000 {
                break;
            }
            pileup_depth += 1;
            if passes_quality {
                covered += 1.0;
            }
        }
        coverage.push(covered);
        visited.push(was_visited);
    }

    if bg_level <= 0.0 {
        return (coverage, visited);
    }
    let adjusted = coverage
        .into_iter()
        .map(|v| {
            let subtracted = (v - bg_level) as i64;
            if subtracted > 0 { subtracted as f64 } else { 0.0 }
        })
        .collect();
    (adjusted, visited)
}

pub fn genebody_coverage(reads: &[IndexedRead], positions: &[i64], bg_level: f64) -> Vec<f64> {
    genebody_coverage_with_visited(reads, positions, bg_level).0
}

/// Shannon entropy (natural log), matching upstream's `shannon_entropy`
/// exactly, including its `0.0`-not-`-0.0` special case.
pub fn shannon_entropy(values: &[f64]) -> f64 {
    let total: f64 = values.iter().sum();
    if total <= 0.0 {
        return 0.0;
    }
    let mut entropy = 0.0;
    for &value in values {
        let probability = value / total;
        entropy += probability * probability.ln();
    }
    if entropy == 0.0 { 0.0 } else { -entropy }
}

/// Ports `tin_score`'s historical entropy formula.
pub fn tin_score(coverage: &[f64], sampled_length: i64) -> f64 {
    if coverage.is_empty() || sampled_length <= 0 {
        return 0.0;
    }
    let effective: Vec<f64> = coverage.iter().copied().filter(|&v| v > 0.0).collect();
    if effective.is_empty() {
        return 0.0;
    }
    let entropy = shannon_entropy(&effective);
    100.0 * entropy.exp() / sampled_length as f64
}

#[derive(Debug, Clone, PartialEq)]
pub struct TranscriptSample {
    pub gene_name: String,
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub intron_size: i64,
    /// 1-based positions; may contain a duplicate (small-transcript
    /// branch) and is NOT deduplicated -- see module docs.
    pub chosen_bases: Vec<i64>,
}

/// Parses the refgene BED12 and yields sampled transcript positions.
/// Robustness tier matches upstream exactly: blank lines and
/// comment/track/browser lines are skipped silently; any other per-line
/// parse failure (fewer than 12 columns, or any `int(...)` failure) is
/// caught and skipped with a stderr warning -- upstream wraps the WHOLE
/// per-line parse in one broad `except (ValueError, IndexError)`.
pub fn genomic_positions(reader: impl BufRead, sample_size: i64) -> io::Result<Vec<TranscriptSample>> {
    let mut out = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let line_number = line_number + 1;
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        let parsed: Result<_, String> = (|| {
            if fields.len() < 12 {
                return Err("fewer than 12 columns".to_string());
            }
            let chromosome = fields[0].to_string();
            let tx_start: i64 = fields[1].parse().map_err(|_| "invalid txStart".to_string())?;
            let tx_end: i64 = fields[2].parse().map_err(|_| "invalid txEnd".to_string())?;
            let gene_name = fields[3].to_string();

            let block_sizes: Vec<i64> = fields[10]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().map_err(|_| "invalid blockSize".to_string()))
                .collect::<Result<_, _>>()?;
            let block_starts: Vec<i64> = fields[11]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().map_err(|_| "invalid blockStart".to_string()))
                .collect::<Result<_, _>>()?;

            let exon_starts: Vec<i64> = block_starts.iter().map(|&s| tx_start + s).collect();
            let exon_ends: Vec<i64> = exon_starts.iter().zip(block_sizes.iter()).map(|(&s, &sz)| s + sz).collect();
            let mrna_size: i64 = block_sizes.iter().sum();
            let intron_size = (tx_end - tx_start - mrna_size).max(0);

            Ok((chromosome, tx_start, tx_end, gene_name, exon_starts, exon_ends, mrna_size, intron_size))
        })();

        let (chromosome, tx_start, tx_end, gene_name, exon_starts, exon_ends, mrna_size, intron_size) = match parsed
        {
            Ok(v) => v,
            Err(msg) => {
                eprintln!("Skipping malformed BED12 line {line_number}: {msg}");
                continue;
            }
        };

        let chosen_bases: Vec<i64> = if mrna_size <= sample_size {
            let mut bases = vec![tx_start + 1, tx_end];
            for (&s, &e) in exon_starts.iter().zip(exon_ends.iter()) {
                bases.extend((s + 1)..=e);
            }
            bases
        } else {
            let step_size = std::cmp::max(1, mrna_size / sample_size);
            let mut all_bases: Vec<i64> = Vec::new();
            let mut exon_bounds: Vec<i64> = Vec::new();
            for (&s, &e) in exon_starts.iter().zip(exon_ends.iter()) {
                all_bases.extend((s + 1)..=e);
                exon_bounds.push(s + 1);
                exon_bounds.push(e);
            }
            let mut sampled_bases = Vec::new();
            let mut idx = 0usize;
            while idx < all_bases.len() {
                sampled_bases.push(all_bases[idx]);
                idx += step_size as usize;
            }
            let mut seen = HashSet::new();
            exon_bounds.into_iter().chain(sampled_bases).filter(|v| seen.insert(*v)).collect()
        };

        out.push(TranscriptSample { gene_name, chrom: chromosome, tx_start, tx_end, intron_size, chosen_bases });
    }

    Ok(out)
}

/// Builds a chromosome-keyed exon interval index for `--subtract-background`.
/// No case-folding: unlike `split_bam.py`/`junction_annotation.py`, tin.py
/// never uppercases chromosome names anywhere.
pub fn build_exon_ranges(exons: &[Bed3]) -> MergedRegions {
    MergedRegions::new(exons)
}

#[derive(Debug, Clone, PartialEq)]
pub struct TinRecord {
    pub gene_name: String,
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub score: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TinSummary {
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n % 2 == 1 { sorted[n / 2] } else { (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0 }
}

fn population_stdev(values: &[f64]) -> f64 {
    let m = mean(values);
    let variance = values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / values.len() as f64;
    variance.sqrt()
}

/// Computes TIN for every sampled transcript against one BAM's read
/// index. Ports `process_bam`'s per-transcript loop (output writing is
/// left to the caller/CLI).
pub fn compute_tin(
    samples: &[TranscriptSample],
    reads_by_chrom: &HashMap<String, Vec<IndexedRead>>,
    min_cov: i64,
    exon_ranges: Option<&MergedRegions>,
) -> (Vec<TinRecord>, TinSummary) {
    let empty = Vec::new();
    let mut records = Vec::with_capacity(samples.len());
    let mut sample_tins: Vec<f64> = Vec::new();

    for s in samples {
        let reads = reads_by_chrom.get(&s.chrom).unwrap_or(&empty);

        let score = if !check_min_reads(reads, s.tx_start, s.tx_end, min_cov) {
            0.0
        } else {
            let mut noise_level = 0.0;
            if let Some(ranges) = exon_ranges {
                if s.intron_size > 0 {
                    let intron_signal = estimate_bg_noise(reads, s.tx_start, s.tx_end, ranges, &s.chrom);
                    noise_level = intron_signal / s.intron_size as f64;
                }
            }
            let mut positions = s.chosen_bases.clone();
            positions.sort_unstable();
            let coverage = genebody_coverage(reads, &positions, noise_level);
            let score = tin_score(&coverage, s.chosen_bases.len() as i64);
            sample_tins.push(score);
            score
        };

        records.push(TinRecord { gene_name: s.gene_name.clone(), chrom: s.chrom.clone(), tx_start: s.tx_start, tx_end: s.tx_end, score });
    }

    let summary = if sample_tins.is_empty() {
        TinSummary::default()
    } else {
        TinSummary { mean: mean(&sample_tins), median: median(&sample_tins), stdev: population_stdev(&sample_tins) }
    };

    (records, summary)
}

/// Ports the `.tin.xls` output: a header line followed by one row per
/// transcript, in input order.
pub fn render_tin_xls(records: &[TinRecord]) -> String {
    let mut out = String::from("geneID\tchrom\ttx_start\ttx_end\tTIN\n");
    for r in records {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.gene_name,
            r.chrom,
            r.tx_start,
            r.tx_end,
            python_str_float(r.score)
        ));
    }
    out
}

/// Ports the `.summary.txt` output. `bam_name` is the bare file name
/// (`Path.name`), matching upstream.
pub fn render_summary(bam_name: &str, summary: &TinSummary) -> String {
    format!(
        "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\n{}\t{}\t{}\t{}\n",
        bam_name,
        python_str_float(summary.mean),
        python_str_float(summary.median),
        python_str_float(summary.stdev),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn genomic_positions_small_transcript_keeps_boundary_duplicate() {
        // Single-exon transcript chr1:0-100 (mrna_size=100 <= sample_size=100):
        // chosen_bases = [1, 100] + range(1,101) -- 1 appears twice (once
        // from the initial boundary pair, once from the exon range).
        let bed = "chr1\t0\t100\ttx1\t0\t+\t0\t100\t0\t1\t100,\t0,\n";
        let samples = genomic_positions(Cursor::new(bed), 100).unwrap();
        assert_eq!(samples.len(), 1);
        let s = &samples[0];
        assert_eq!(s.gene_name, "tx1");
        assert_eq!(s.intron_size, 0);
        assert_eq!(s.chosen_bases.len(), 102); // [1,100] + 100 positions (1..=100)
        assert_eq!(s.chosen_bases.iter().filter(|&&v| v == 1).count(), 2);
    }

    #[test]
    fn genomic_positions_large_transcript_dedups_and_drops_boundary_pair() {
        // Two-exon transcript chr1:0-1000 with a gap between the last
        // exon and tx_end (exon2 ends at 900, tx_end is 1000), so
        // tx_end=1000 is NOT among the exon boundaries naturally: mrna_size
        // =200 > sample_size=10, step_size = max(1, 200/10) = 20. The
        // initial [tx_start+1, tx_end] = [1, 1000] pair is NOT present
        // (else-branch reassigns chosen_bases entirely); exon_bounds are
        // the per-exon boundaries instead: (1,100) and (801,900).
        let bed = "chr1\t0\t1000\ttx1\t0\t+\t0\t1000\t0\t2\t100,100,\t0,800,\n";
        let samples = genomic_positions(Cursor::new(bed), 10).unwrap();
        let s = &samples[0];
        assert_eq!(s.intron_size, 800);
        assert!(!s.chosen_bases.contains(&1000)); // tx_end itself never appears
        assert!(s.chosen_bases.contains(&1));
        assert!(s.chosen_bases.contains(&100));
        assert!(s.chosen_bases.contains(&801));
        assert!(s.chosen_bases.contains(&900));
        // No duplicates in the else-branch (uniqify applied).
        let mut seen = HashSet::new();
        assert!(s.chosen_bases.iter().all(|v| seen.insert(*v)));
    }

    #[test]
    fn genomic_positions_skips_malformed_lines_with_warning_not_error() {
        let bed = "\
track name=x
# comment

not enough fields
chr1\t0\tnotanumber\ttx1\t0\t+\t0\t100\t0\t1\t100,\t0,
chr1\t0\t100\ttx2\t0\t+\t0\t100\t0\t1\t100,\t0,
";
        let samples = genomic_positions(Cursor::new(bed), 100).unwrap();
        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].gene_name, "tx2");
    }

    #[test]
    fn shannon_entropy_matches_python_hand_computation() {
        // Verified against python3: -sum(p*ln(p)) for p=[0.5,0.5] = ln(2).
        let e = shannon_entropy(&[1.0, 1.0]);
        assert!((e - std::f64::consts::LN_2).abs() < 1e-12);

        // Single value -> p=1.0, ln(1)=0 -> the 0.0-not-negative-zero branch.
        let e2 = shannon_entropy(&[5.0]);
        assert_eq!(e2, 0.0);
        assert!(e2.is_sign_positive());
    }

    #[test]
    fn tin_score_zero_when_no_positive_coverage() {
        assert_eq!(tin_score(&[0.0, 0.0], 10), 0.0);
        assert_eq!(tin_score(&[], 10), 0.0);
        assert_eq!(tin_score(&[1.0], 0), 0.0);
    }

    #[test]
    fn tin_score_uniform_coverage_gives_100() {
        // Uniform coverage across every sampled position -> maximum entropy
        // (ln(n)) -> exp(ln(n))=n -> 100*n/n = 100.0 exactly.
        let coverage = vec![5.0; 20];
        let score = tin_score(&coverage, 20);
        assert!((score - 100.0).abs() < 1e-9);
    }

    fn header_with_chrom(name: &str, len: usize) -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                name,
                sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(len).unwrap(),
                ),
            )
            .build()
    }

    fn to_bam_records(header: &sam::Header, records: &[sam::alignment::record_buf::RecordBuf]) -> Vec<bam::Record> {
        use sam::alignment::io::Write as _;
        let mut buf = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut buf);
            writer.write_header(header).unwrap();
            for r in records {
                writer.write_alignment_record(header, r).unwrap();
            }
        }
        let mut reader = bam::io::Reader::new(buf.as_slice());
        reader.read_header().unwrap();
        reader.records().map(|r| r.unwrap()).collect()
    }

    #[test]
    fn build_read_index_and_check_min_reads() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);
        let mut records = Vec::new();
        for i in 0..3 {
            records.push(
                RecordBuf::builder()
                    .set_name(format!("r{i}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(noodles_core::Position::try_from(101 + i).unwrap())
                    .set_mapping_quality(MappingQuality::new(40).unwrap())
                    .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
                    .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
                    .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
                    .build(),
            );
        }

        let bam_records = to_bam_records(&header, &records);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header).unwrap();
        let reads = &index["chr1"];
        assert_eq!(reads.len(), 3);

        // 3 distinct read starts (100,101,102 0-based); cutoff=2 (strict >).
        assert!(check_min_reads(reads, 0, 200, 2));
        assert!(!check_min_reads(reads, 0, 200, 3));
    }

    #[test]
    fn genebody_coverage_counts_matching_reads_and_excludes_low_quality_and_duplicates() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);

        // Read A: pos 0-based 0, 10M, quality 40 everywhere -> covers 1-based [1,10].
        let read_a = RecordBuf::builder()
            .set_name("a")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 10]))
            .build();

        // Read B: same span, but flagged duplicate -> excluded from coverage
        // (pileup's flag_filter) even though it would be in the read index.
        let read_b = RecordBuf::builder()
            .set_name("b")
            .set_flags(Flags::DUPLICATE)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 10]))
            .build();

        // Read C: same span, but quality 5 everywhere (< 13) -> excluded per-base.
        let read_c = RecordBuf::builder()
            .set_name("c")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![5; 10]))
            .build();

        let bam_records = to_bam_records(&header, &[read_a, read_b, read_c]);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header).unwrap();
        let reads = &index["chr1"];
        assert_eq!(reads.len(), 3); // all 3 pass the qcfail/unmapped/secondary build filter

        let positions = vec![1, 5, 10];
        let coverage = genebody_coverage(reads, &positions, 0.0);
        // Only read_a counts at every position: read_b is a duplicate
        // (excluded), read_c's bases are all below min_base_quality=13.
        assert_eq!(coverage, vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn genebody_coverage_honors_pysam_default_max_depth() {
        let reads: Vec<IndexedRead> = (0..8001)
            .map(|_| IndexedRead {
                query_name: "synthetic".to_string(),
                is_paired: false,
                start: 0,
                end: 20,
                is_duplicate: false,
                match_blocks: vec![(0, 20, 0)],
                skip_delete_blocks: vec![],
                qualities: vec![40; 20],
                query_length: 20,
            })
            .collect();

        // pysam.AlignmentFile.pileup() defaults max_depth to 8000.
        assert_eq!(genebody_coverage(&reads, &[1], 0.0), vec![8000.0]);
    }

    #[test]
    fn genebody_coverage_deduplicates_overlapping_mates() {
        let reads = vec![
            IndexedRead {
                query_name: "pair".to_string(),
                is_paired: true,
                start: 0,
                end: 20,
                is_duplicate: false,
                match_blocks: vec![(0, 20, 0)],
                skip_delete_blocks: vec![],
                qualities: vec![40; 20],
                query_length: 20,
            },
            IndexedRead {
                query_name: "pair".to_string(),
                is_paired: true,
                start: 10,
                end: 30,
                is_duplicate: false,
                match_blocks: vec![(10, 30, 0)],
                skip_delete_blocks: vec![],
                qualities: vec![40; 20],
                query_length: 20,
            },
        ];

        // The overlap at reference position 11 is one pileup base, not two.
        assert_eq!(genebody_coverage(&reads, &[11], 0.0), vec![1.0]);
    }

    #[test]
    fn render_tin_xls_and_summary_exact_text() {
        // Verified against python3: str(85.5), str(0.0), str(100).
        let records = vec![
            TinRecord { gene_name: "tx1".into(), chrom: "chr1".into(), tx_start: 0, tx_end: 100, score: 85.5 },
            TinRecord { gene_name: "tx2".into(), chrom: "chr1".into(), tx_start: 200, tx_end: 300, score: 0.0 },
        ];
        let xls = render_tin_xls(&records);
        assert_eq!(xls, "geneID\tchrom\ttx_start\ttx_end\tTIN\ntx1\tchr1\t0\t100\t85.5\ntx2\tchr1\t200\t300\t0.0\n");

        let summary = TinSummary { mean: 42.75, median: 42.75, stdev: 0.0 };
        let summary_text = render_summary("sample.bam", &summary);
        assert_eq!(summary_text, "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nsample.bam\t42.75\t42.75\t0.0\n");
    }
}
