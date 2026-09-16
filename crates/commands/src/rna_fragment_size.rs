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
//! indexed fetch: `read.reference_start < end && start < read.reference_end`),
//! just without the index's algorithmic speedup. A genuinely large BAM
//! will use more memory and time than upstream here; that's a disclosed
//! performance characteristic, not a correctness gap, and is a natural
//! target for the optimization pass (PORTING_PLAN Step 9), not now.

use std::collections::HashMap;
use std::io;

use noodles_bam as bam;
use noodles_sam as sam;
use rseqc_formats::cigar::reference_span;

struct ReadEntry {
    ref_start: i64,
    ref_end: i64,
    mate_start: i64,
    query_length: i64,
    mapq: u8,
}

/// In-memory stand-in for pysam's indexed `fetch()`, built once from a
/// full sequential scan and queried per BED line.
pub struct IndexedReads {
    by_chrom: HashMap<String, Vec<ReadEntry>>,
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

            // Matches fragment_size()'s per-read filters (is_paired,
            // !is_read2, !mate_is_unmapped, !qcfail, !duplicate,
            // !secondary) plus the implicit "fetch() never returns
            // unmapped reads" behavior.
            if !flags.is_segmented()
                || flags.is_last_segment()
                || flags.is_mate_unmapped()
                || flags.is_qc_fail()
                || flags.is_duplicate()
                || flags.is_secondary()
                || flags.is_unmapped()
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
            let (_, ref_end) = reference_span(ref_start as usize, ops);

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

        Ok(Self { by_chrom })
    }

    fn fetch(&self, chrom: &str, start: i64, end: i64) -> impl Iterator<Item = &ReadEntry> {
        self.by_chrom
            .get(chrom)
            .into_iter()
            .flatten()
            .filter(move |r| r.ref_start < end && start < r.ref_end)
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

/// NumPy-equivalent mean/median/std (population std, ddof=0).
fn mean_median_std(values: &[i64]) -> (f64, f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().map(|&v| v as f64).sum::<f64>() / n;

    let mut sorted: Vec<f64> = values.iter().map(|&v| v as f64).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = sorted.len() / 2;
    let median = if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    };

    let variance = values.iter().map(|&v| (v as f64 - mean).powi(2)).sum::<f64>() / n;
    let std = variance.sqrt();

    (mean, median, std)
}

pub fn format_result(s: &FragmentStats) -> String {
    let (mean_s, median_s, std_s) = match s.stats {
        None => ("0".to_string(), "0".to_string(), "0".to_string()),
        Some((mean, median, std)) => (mean.to_string(), median.to_string(), std.to_string()),
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
        let reads = IndexedReads::build(Vec::<io::Result<bam::Record>>::new().into_iter(), &test_header()).unwrap();
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
}
