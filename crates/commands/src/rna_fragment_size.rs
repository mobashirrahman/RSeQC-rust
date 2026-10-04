//! Port of `RNA_fragment_size.py`: calculate fragment-size statistics for
//! each BED12 transcript/gene record. Contract: see
//! `compatibility/commands.yaml` entry `RNA_fragment_size.py`; algorithm
//! ported directly from the script itself (it does not use
//! `qcmodule.SAM`/`qcmodule.BED` at all -- `parse_bed12_record`,
//! `overlap_length2`, `fragment_size`, `format_result` in
//! `oracle/upstream-src/scripts/RNA_fragment_size.py`).
//!
//! Unlike every other command ported so far: a malformed BED12 line
//! (fewer than 12 columns, or `blockSizes`/`blockStarts` count mismatch)
//! is NOT skipped -- upstream's `parse_bed12_record` raises `ValueError`,
//! uncaught within `fragment_size()`, propagating out through `main()`'s
//! `except (OSError, ValueError, RuntimeError)` handler to a hard exit(1).
//! This port returns an `io::Error` from `parse_bed12_line` for the same
//! reason, deliberately NOT silently skipping the bad line.
//!
//! Requires pysam's `fetch(chrom, start, end)` region-query semantics,
//! which normally uses a `.bai` index. No BAI parsing exists yet, so this
//! builds an in-memory per-chromosome index from a full sequential BAM
//! scan instead -- behaviorally equivalent (same overlap semantics as
//! indexed fetch: `read.reference_start < end && start < read.reference_end`).
//!
//! The query side is NOT a full scan per BED line. An earlier version of this
//! file filtered every read on the chromosome for every transcript, which made
//! a run cost O(reads x transcripts); measured against upstream on a 2.2M-read
//! alignment with a 3000-transcript model it came out 1.47x SLOWER than the
//! Python it replaces, because pysam answers the same query from the BAI in
//! O(log R + hits). `ChromReads` below restores that shape: entries stay in BAM
//! file order (which is what determines the output's last bits -- see
//! `mean_median_std`), and a parallel `by_start` array of indices sorted by
//! `ref_start` narrows each query to the only window that can contain a hit.
//! Memory cost is 4 bytes per retained read.

use std::collections::HashMap;
use std::io;

use noodles_bam as bam;
use noodles_sam as sam;
use rseqc_formats::cigar::reference_span;

use crate::python_fmt::python_str_float;

struct ReadEntry {
    ref_start: i64,
    ref_end: i64,
    mate_start: i64,
    query_length: i64,
    mapq: u8,
}

/// One chromosome's reads, in BAM file order, plus a position-sorted index
/// over them.
///
/// `entries` MUST stay in file order: `mean_median_std` sums in insertion
/// order, and `numpy_mean`/`numpy_std`'s pairwise summation is order-sensitive
/// in its last bits, so reordering reads changes the printed mean and std.
/// `by_start` therefore holds *indices* into `entries`, not entries.
struct ChromReads {
    entries: Vec<ReadEntry>,
    /// `entries` indices, ascending by `entries[i].ref_start`.
    by_start: Vec<u32>,
    /// Largest `ref_end - ref_start` in `entries`. Bounds how far left of
    /// `start` a hit can begin, which is what makes the lower window bound
    /// exact rather than a guess.
    max_span: i64,
}

impl ChromReads {
    /// Half-open `[lo, hi)` range into `by_start` holding every entry that
    /// could possibly satisfy the overlap predicate for `[start, end)`.
    ///
    /// `lo >= hi` is a legitimate answer (an empty or left-of-everything
    /// query) and callers must not slice with it unchecked.
    fn window(&self, start: i64, end: i64) -> (usize, usize) {
        let lo_bound = start - self.max_span;
        let lo = self.by_start.partition_point(|&i| self.entries[i as usize].ref_start <= lo_bound);
        let hi = self.by_start.partition_point(|&i| self.entries[i as usize].ref_start < end);
        (lo, hi)
    }
}

/// In-memory stand-in for pysam's indexed `fetch()`, built once from a
/// full sequential BAM scan and queried per BED line.
pub struct IndexedReads {
    by_chrom: HashMap<String, ChromReads>,
}

impl IndexedReads {
    pub fn build<I>(records: I, header: &sam::Header) -> io::Result<Self>
    where
        I: IntoIterator<Item = io::Result<bam::Record>>,
    {
        let mut by_chrom: HashMap<String, Vec<ReadEntry>> = HashMap::new();

        for result in records {
            let record = result?;
            let flags = record.flags();

            // Matches fragment_size()'s per-read filters exactly: is_paired,
            // !is_read2, !mate_is_unmapped, !qcfail, !duplicate,
            // !secondary. Upstream never checks `is_unmapped` here, and
            // (unlike this port's earlier assumption) pysam's indexed
            // `fetch(chrom, start, end)` DOES return "placed" unmapped
            // reads -- reads flagged unmapped whose RNAME/POS were copied
            // from a mapped mate so they sort next to it (htslib's
            // `bam_index` bins these using their POS the same as mapped
            // records). Filtering `is_unmapped()` here dropped exactly
            // those records, undercounting fragments whenever a mate pair
            // has one mapped and one placed-unmapped end.
            if !flags.is_segmented()
                || flags.is_last_segment()
                || flags.is_mate_unmapped()
                || flags.is_qc_fail()
                || flags.is_duplicate()
                || flags.is_secondary()
            {
                continue;
            }

            let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);

            let Some(ref_id) = record.reference_sequence_id().transpose()? else {
                continue;
            };
            let Some((chrom, _)) = header.reference_sequences().get_index(ref_id) else {
                continue;
            };

            let Some(pos) = record.alignment_start().transpose()? else {
                continue;
            };
            let ref_start = (pos.get() - 1) as i64;

            let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
            let (_, span_end) = reference_span(ref_start as usize, ops);
            // A placed-unmapped read has no CIGAR, so `reference_span`
            // returns a zero-width span (`span_end == ref_start`). htslib
            // still registers a 1 bp footprint for such records in the
            // BAM index, and pysam's region `fetch()` returns them
            // accordingly, so widen the span here to match.
            let ref_end = span_end.max(ref_start as usize + 1);

            let Some(mate_pos) = record.mate_alignment_start().transpose()? else {
                continue;
            };
            let mate_start = (mate_pos.get() - 1) as i64;

            let query_length = record.sequence().len() as i64;

            by_chrom.entry(chrom.to_string()).or_default().push(ReadEntry {
                ref_start,
                ref_end: ref_end as i64,
                mate_start,
                query_length,
                mapq,
            });
        }

        let by_chrom = by_chrom
            .into_iter()
            .map(|(chrom, entries)| {
                let max_span =
                    entries.iter().map(|e| e.ref_end - e.ref_start).max().unwrap_or(1).max(1);
                let mut by_start: Vec<u32> = (0..entries.len() as u32).collect();
                by_start.sort_unstable_by_key(|&i| entries[i as usize].ref_start);
                (chrom, ChromReads { entries, by_start, max_span })
            })
            .collect();

        Ok(Self { by_chrom })
    }

    /// Reads overlapping `[start, end)`, in BAM file order.
    ///
    /// A hit must satisfy `ref_start < end` and `ref_end > start`. Since
    /// `ref_end <= ref_start + max_span`, the second implies
    /// `ref_start > start - max_span`, so both bounds are exact: no entry
    /// outside the window can satisfy the predicate. `max_span` is the
    /// longest read on the chromosome, typically a few hundred bases, so the
    /// window is the transcript plus a short left margin rather than the whole
    /// chromosome.
    fn fetch(&self, chrom: &str, start: i64, end: i64) -> Vec<&ReadEntry> {
        let Some(c) = self.by_chrom.get(chrom) else {
            return Vec::new();
        };
        let (lo, hi) = c.window(start, end);
        if lo >= hi {
            return Vec::new();
        }
        let mut window: Vec<u32> = c.by_start[lo..hi].to_vec();
        // Restore file order so the caller's accumulation order -- and
        // therefore the printed mean/std -- matches the unindexed scan.
        window.sort_unstable();
        window
            .into_iter()
            .map(|i| &c.entries[i as usize])
            .filter(|r| r.ref_start < end && start < r.ref_end)
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BedRecord {
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub name: String,
    /// `[exon_start+1, exon_end+1]` per exon, matching upstream's literal
    /// (not "corrected") 1-based conversion.
    pub exon_ranges: Vec<(i64, i64)>,
}

/// Parses one BED12 line. Returns `Err` (propagate, don't skip) on
/// malformed input, matching upstream's uncaught-`ValueError` crash.
pub fn parse_bed12_line(line: &str) -> io::Result<BedRecord> {
    let bad = |msg: String| io::Error::new(io::ErrorKind::InvalidData, msg);

    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 12 {
        return Err(bad(format!("BED line has {} columns; expected 12", fields.len())));
    }

    let chrom = fields[0].to_string();
    let tx_start: i64 = fields[1].parse().map_err(|_| bad("invalid tx_start".into()))?;
    let tx_end: i64 = fields[2].parse().map_err(|_| bad("invalid tx_end".into()))?;
    let name = fields[3].to_string();

    let block_count: usize = fields[9].parse().map_err(|_| bad("invalid blockCount".into()))?;
    let block_sizes: Vec<i64> = fields[10]
        .trim_end_matches(',')
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().map_err(|_| bad("invalid blockSizes".into())))
        .collect::<io::Result<_>>()?;
    let relative_starts: Vec<i64> = fields[11]
        .trim_end_matches(',')
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().map_err(|_| bad("invalid blockStarts".into())))
        .collect::<io::Result<_>>()?;

    if block_sizes.len() != block_count {
        return Err(bad(format!(
            "blockCount={block_count}, but found {} block sizes",
            block_sizes.len()
        )));
    }
    if relative_starts.len() != block_count {
        return Err(bad(format!(
            "blockCount={block_count}, but found {} block starts",
            relative_starts.len()
        )));
    }

    let exon_ranges = relative_starts
        .iter()
        .zip(block_sizes.iter())
        .map(|(&rel_start, &size)| {
            let exon_start = tx_start + rel_start;
            let exon_end = exon_start + size;
            (exon_start + 1, exon_end + 1)
        })
        .collect();

    Ok(BedRecord { chrom, tx_start, tx_end, name, exon_ranges })
}

fn overlap_inclusive(a: (i64, i64), b: (i64, i64)) -> i64 {
    let lo = a.0.max(b.0);
    let hi = a.1.min(b.1);
    (hi - lo + 1).max(0)
}

#[derive(Debug, Clone, PartialEq)]
pub struct FragmentStats {
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub name: String,
    pub count: usize,
    /// `None` when `count < ncut`: upstream leaves mean/median/std as the
    /// Python int `0` (renders as "0"), not a float, in that case.
    pub stats: Option<(f64, f64, f64)>,
}

pub fn compute_fragment_sizes(bed: &BedRecord, reads: &IndexedReads, qcut: u8, ncut: usize) -> FragmentStats {
    let mut fragment_sizes: Vec<i64> = Vec::new();

    for read in reads.fetch(&bed.chrom, bed.tx_start, bed.tx_end) {
        if read.mapq < qcut {
            continue;
        }

        let (mut read_start, mut mate_start) = (read.ref_start, read.mate_start);
        if read_start > mate_start {
            std::mem::swap(&mut read_start, &mut mate_start);
        }
        if read_start < bed.tx_start || mate_start > bed.tx_end {
            continue;
        }

        let mapped_range = (read_start + 1, mate_start);
        let overlap: i64 = bed.exon_ranges.iter().map(|&r| overlap_inclusive(r, mapped_range)).sum();
        fragment_sizes.push(overlap + read.query_length);
    }

    let count = fragment_sizes.len();
    let stats = if count < ncut {
        None
    } else {
        Some(mean_median_std(&fragment_sizes))
    };

    FragmentStats {
        chrom: bed.chrom.clone(),
        tx_start: bed.tx_start,
        tx_end: bed.tx_end,
        name: bed.name.clone(),
        count,
        stats,
    }
}

/// NumPy-equivalent mean/median/std (population std, ddof=0). Mean and std
/// use `numpy_mean`/`numpy_std` (pairwise summation, matching `np.mean`/
/// `np.std` bit-for-bit) over `values` in their *original* (insertion)
/// order -- `np.mean`/`np.std` never sort their input, and pairwise
/// summation's last-bit result depends on element order. Only the median
/// needs a sorted copy.
fn mean_median_std(values: &[i64]) -> (f64, f64, f64) {
    let floats: Vec<f64> = values.iter().map(|&v| v as f64).collect();
    let mean = crate::python_fmt::numpy_mean(&floats);
    let std = crate::python_fmt::numpy_std(&floats, 0);

    let mut sorted = floats;
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = sorted.len() / 2;
    let median = if sorted.len().is_multiple_of(2) {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    };

    (mean, median, std)
}

pub fn format_result(s: &FragmentStats) -> String {
    let (mean_s, median_s, std_s) = match s.stats {
        None => ("0".to_string(), "0".to_string(), "0".to_string()),
        Some((mean, median, std)) => (python_str_float(mean), python_str_float(median), python_str_float(std)),
    };
    format!(
        "{}\t{}\t{}\t{}\t{}\t{mean_s}\t{median_s}\t{std_s}",
        s.chrom, s.tx_start, s.tx_end, s.name, s.count
    )
}

pub const HEADER: &str = "chrom\ttx_start\ttx_end\tsymbol\tfrag_count\tfrag_mean\tfrag_median\tfrag_std";

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
    use noodles_sam::{
        alignment::{
            io::Write as _,
            record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind},
            record_buf::{Cigar, RecordBuf},
        },
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).unwrap()),
            )
            .build()
    }

    fn to_bam_records(header: &sam::Header, records: &[RecordBuf]) -> Vec<bam::Record> {
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
    fn parses_bed12_line() {
        let bed = parse_bed12_line("chr1\t100\t300\tgeneA\t0\t+\t100\t300\t0\t2\t50,50,\t0,150,\n").unwrap();
        assert_eq!(bed.chrom, "chr1");
        assert_eq!(bed.tx_start, 100);
        assert_eq!(bed.tx_end, 300);
        assert_eq!(bed.name, "geneA");
        // exon1: start=100+0=100,end=150 -> (101,151); exon2: start=100+150=250,end=300 -> (251,301)
        assert_eq!(bed.exon_ranges, vec![(101, 151), (251, 301)]);
    }

    #[test]
    fn rejects_short_line() {
        let err = parse_bed12_line("chr1\t100\t300\n").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn rejects_block_count_mismatch() {
        let err = parse_bed12_line("chr1\t100\t300\tgeneA\t0\t+\t100\t300\t0\t2\t50,\t0,\n").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn below_ncut_renders_integer_zero() {
        let bed = BedRecord {
            chrom: "chr1".into(),
            tx_start: 0,
            tx_end: 1000,
            name: "geneA".into(),
            exon_ranges: vec![(1, 1000)],
        };
        let reads = IndexedReads::build(Vec::<io::Result<bam::Record>>::new(), &test_header()).unwrap();
        let stats = compute_fragment_sizes(&bed, &reads, 30, 3);
        assert_eq!(stats.count, 0);
        assert_eq!(format_result(&stats), "chr1\t0\t1000\tgeneA\t0\t0\t0\t0");
    }

    #[test]
    fn computes_fragment_size_from_overlapping_pair() {
        let header = test_header();

        // Read1 at 0-based pos 10, CIGAR 20M, mate at 0-based pos 40.
        let mk = |pos: usize, mate_pos: usize| {
            RecordBuf::builder()
                .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::new(pos).unwrap())
                .set_mate_alignment_start(Position::new(mate_pos).unwrap())
                .set_mapping_quality(MappingQuality::new(40).unwrap())
                .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 20)]))
                .set_sequence(noodles_sam::alignment::record_buf::Sequence::from(vec![b'A'; 20]))
                .build()
        };
        // pos is 1-based in noodles Position; use 11 (0-based 10) and 41 (0-based 40).
        let record = mk(11, 41);

        let bam_records = to_bam_records(&header, &[record]);
        let reads = IndexedReads::build(bam_records.into_iter().map(Ok), &header).unwrap();

        let bed = BedRecord {
            chrom: "chr1".into(),
            tx_start: 0,
            tx_end: 1000,
            name: "geneA".into(),
            exon_ranges: vec![(1, 1000)],
        };

        let stats = compute_fragment_sizes(&bed, &reads, 30, 1);
        assert_eq!(stats.count, 1);
        // mapped_range = (10+1, 40) = (11,40); exon (1,1000) fully contains
        // it: overlap = 40-11+1 = 30; + query_length 20 = 50.
        let (mean, median, std) = stats.stats.unwrap();
        assert_eq!(mean, 50.0);
        assert_eq!(median, 50.0);
        assert_eq!(std, 0.0);
    }

    #[test]
    fn counts_placed_unmapped_mate_like_upstream() {
        // Regression for the count divergence upstream showed for
        // chr1_g3 (py count 6 vs rust count 5): a "placed unmapped" read
        // -- FLAG has UNMAPPED set, but RNAME/POS were copied from its
        // mapped mate so it sorts next to it. Upstream's per-read filter
        // (RNA_fragment_size.py `fragment_size()`, ~oracle/upstream-src/
        // scripts/RNA_fragment_size.py:322-338) never checks `is_unmapped`
        // at all, and pysam's indexed `fetch()` DOES return such reads, so
        // this record must count as read1 of the pair -- not be silently
        // dropped for being flagged unmapped.
        let header = test_header();

        let unmapped_read1 = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::UNMAPPED)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(101).unwrap())
            .set_mate_alignment_start(Position::new(101).unwrap())
            .set_mapping_quality(MappingQuality::new(0).unwrap())
            .set_sequence(noodles_sam::alignment::record_buf::Sequence::from(vec![b'A'; 20]))
            .build();

        let bam_records = to_bam_records(&header, &[unmapped_read1]);
        let reads = IndexedReads::build(bam_records.into_iter().map(Ok), &header).unwrap();

        let bed = BedRecord {
            chrom: "chr1".into(),
            tx_start: 0,
            tx_end: 1000,
            name: "geneA".into(),
            exon_ranges: vec![(1, 1000)],
        };

        let stats = compute_fragment_sizes(&bed, &reads, 0, 1);
        assert_eq!(stats.count, 1, "placed-unmapped read1 must still be counted");
    }

    #[test]
    fn mean_and_std_use_numpy_pairwise_summation() {
        // Regression for the last-digit divergence upstream showed on
        // chr2 (std 113.64928508354112 vs 113.6492850835411): np.mean/
        // np.std sum with numpy's pairwise-blocked kernel, not naive
        // left-to-right summation, and the two differ in the last bit for
        // long-enough inputs -- concretely, in the squared-deviation sum
        // inside `np.std`, which numpy_std/numpy_sum's >128-element
        // recursive-halving path reproduces bit-for-bit but a naive
        // running sum does not. Use >128 elements, mirroring a real
        // per-transcript fragment_sizes vector, and check the whole
        // pipeline (mean_median_std) matches the shared numpy_mean/
        // numpy_std helpers exactly.
        let values: Vec<i64> = (0..200).map(|i| 100 + (i % 37)).collect();
        let floats: Vec<f64> = values.iter().map(|&v| v as f64).collect();

        let (mean, _median, std) = mean_median_std(&values);

        assert_eq!(mean, crate::python_fmt::numpy_mean(&floats));
        assert_eq!(std, crate::python_fmt::numpy_std(&floats, 0));

        // The squared-deviation sum inside std is where naive and
        // pairwise summation actually diverge for this input (mean of
        // integers sums exactly regardless of order, but the squared
        // deviations are non-integer floats): confirm numpy_std is doing
        // real work here, not coincidentally matching a naive sum.
        let squares: Vec<f64> = floats.iter().map(|&v| (v - mean) * (v - mean)).collect();
        let naive_variance = squares.iter().sum::<f64>() / floats.len() as f64;
        let naive_std = naive_variance.sqrt();
        assert_ne!(
            std.to_bits(),
            naive_std.to_bits(),
            "expected pairwise and naive summation to differ in the last bit for this input"
        );
    }

    /// The literal predicate `fetch` replaced: every read on the chromosome,
    /// in file order. Kept as the oracle for the indexed query below.
    fn brute_force_fetch<'a>(
        reads: &'a IndexedReads,
        chrom: &str,
        start: i64,
        end: i64,
    ) -> Vec<&'a ReadEntry> {
        reads
            .by_chrom
            .get(chrom)
            .into_iter()
            .flat_map(|c| c.entries.iter())
            .filter(move |r| r.ref_start < end && start < r.ref_end)
            .collect()
    }

    /// Build reads scattered over `chrom` with pseudo-random positions and
    /// spans, deliberately NOT in coordinate order, so the indexed query is
    /// tested against input it cannot shortcut by assuming sortedness.
    fn scattered_reads(header: &sam::Header, n: usize, seed: u64) -> IndexedReads {
        let mut state = seed;
        let mut next = move || {
            // xorshift64*, so the fixture is reproducible without a dependency
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let records: Vec<RecordBuf> = (0..n)
            .map(|_| {
                let pos = (next() % 9_000) as usize + 1;
                let span = (next() % 300) as usize + 1;
                RecordBuf::builder()
                    .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::new(pos).unwrap())
                    .set_mate_alignment_start(Position::new(pos + 1).unwrap())
                    .set_mapping_quality(MappingQuality::new(40).unwrap())
                    .set_cigar(Cigar::from(vec![Op::new(Kind::Match, span)]))
                    .set_sequence(noodles_sam::alignment::record_buf::Sequence::from(vec![
                        b'A';
                        span
                    ]))
                    .build()
            })
            .collect();
        let bam_records = to_bam_records(header, &records);
        IndexedReads::build(bam_records.into_iter().map(Ok), header).unwrap()
    }

    #[test]
    fn indexed_fetch_matches_the_full_scan_it_replaced() {
        // The window bounds in `fetch` are derived from `max_span`, so they are
        // only correct if no hit can start further left than `start - max_span`.
        // Sweep query windows across the whole chromosome, including ones that
        // start mid-read and ones entirely inside a read's span, and require
        // identical results -- same entries, same order -- every time.
        let header = test_header();
        let reads = scattered_reads(&header, 2_000, 0x5eed_1234);

        for start in (0..9_000).step_by(37) {
            for width in [1, 2, 17, 250, 1_000, 9_000] {
                let end = start + width;
                let indexed: Vec<_> =
                    reads.fetch("chr1", start, end).into_iter().map(|r| (r.ref_start, r.ref_end)).collect();
                let brute: Vec<_> =
                    brute_force_fetch(&reads, "chr1", start, end).into_iter().map(|r| (r.ref_start, r.ref_end)).collect();
                assert_eq!(
                    indexed, brute,
                    "indexed fetch disagreed with the full scan for [{start}, {end})"
                );
            }
        }
    }

    #[test]
    fn indexed_fetch_preserves_file_order_for_unsorted_positions() {
        // `numpy_mean`/`numpy_std` sum in insertion order and their pairwise
        // kernel is order-sensitive in the last bits, so `fetch` must return
        // reads in BAM file order even though it searches a position-sorted
        // index. Built unsorted on purpose: a stable result here cannot be an
        // accident of coordinate-sorted input.
        let header = test_header();
        let reads = scattered_reads(&header, 2_000, 0xabcd_0001);

        let indexed: Vec<i64> = reads.fetch("chr1", 0, 9_000).iter().map(|r| r.ref_start).collect();
        let brute: Vec<i64> = brute_force_fetch(&reads, "chr1", 0, 9_000).iter().map(|r| r.ref_start).collect();

        assert_eq!(indexed, brute, "file order must be preserved");
        assert!(
            indexed.windows(2).any(|w| w[0] > w[1]),
            "fixture must actually be unsorted, or this test proves nothing"
        );
    }

    #[test]
    fn indexed_fetch_handles_a_window_left_of_every_read() {
        let header = test_header();
        let reads = scattered_reads(&header, 500, 0x1111_2222);
        // `lo >= hi` must return empty rather than panic on the slice range.
        assert!(reads.fetch("chr1", -500, -400).is_empty());
        assert!(reads.fetch("chrMissing", 0, 9_000).is_empty());
    }

    #[test]
    fn query_window_is_bounded_by_the_query_not_the_chromosome() {
        // Guards the defect this index was added to fix. `fetch` used to filter
        // every read on the chromosome for every transcript, making a run cost
        // O(reads x transcripts); measured against upstream on a 2.2M-read
        // alignment with a 3000-transcript model that came out 1.47x SLOWER
        // than the Python it replaces, because pysam answers the same query
        // from the BAI in O(log R + hits).
        //
        // The invariant is that the window handed to the overlap filter is
        // sized by the QUERY, not by chromosome depth. A full scan scores
        // depth/width of the chromosome -- 1.0 here -- so the bound below
        // fails loudly if the window ever degenerates back to a scan. Work is
        // counted rather than timed, so this cannot flake on a loaded machine.
        let header = test_header();
        let queries: Vec<(i64, i64)> = (0..9_000).step_by(900).map(|s| (s, s + 400)).collect();

        for depth in [500_usize, 4_000, 32_000] {
            let reads = scattered_reads(&header, depth, 0x9999_0001);
            let c = &reads.by_chrom["chr1"];
            let per_query: usize =
                queries.iter().map(|&(s, e)| { let (lo, hi) = c.window(s, e); hi - lo }).sum();
            let fraction = per_query as f64 / (queries.len() * depth) as f64;
            assert!(
                fraction < 0.25,
                "at depth {depth} the window covered {fraction:.3} of the chromosome \
                 per query; a full scan would be 100%"
            );
        }

        // And the fraction must not creep up as the chromosome deepens, which
        // is what would happen if the lower bound were a guess rather than
        // derived from `max_span`.
        let fraction_at = |depth: usize| -> f64 {
            let reads = scattered_reads(&header, depth, 0x9999_0001);
            let c = &reads.by_chrom["chr1"];
            let per_query: usize =
                queries.iter().map(|&(s, e)| { let (lo, hi) = c.window(s, e); hi - lo }).sum();
            per_query as f64 / (queries.len() * depth) as f64
        };
        let (shallow, deep) = (fraction_at(500), fraction_at(8_000));
        assert!(
            deep <= shallow * 1.5,
            "window/chromosome fraction grew with depth: {shallow:.3} at 500 reads -> \
             {deep:.3} at 8000 reads"
        );
    }
}
