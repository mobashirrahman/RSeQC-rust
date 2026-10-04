//! Port of `bam_stat.py`: summarize mapping statistics for a BAM file.
//! Contract: see `compatibility/commands.yaml` entry `bam_stat.py`; algorithm
//! ported from `ParseBAM.stat()` in `oracle/upstream-src/src/qcmodule/SAM.py`.
//! Known divergences: see `compatibility/divergences.yaml` DIV-0001
//! (splice detection simplified to a direct CIGAR-Skip check). DIV-0002
//! (BAM-only, no SAM-text support) is CLOSED for this command -- the
//! CLI opens input via `rseqc_formats::open_alignments`, which
//! dispatches on the `.bam`/`.sam` extension and converts SAM-text
//! records to genuine `bam::Record`s before they ever reach this
//! module, so `compute_stats` itself needed no changes at all.

use std::io;

use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind as CigarOpKind;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct BamStatCounts {
    pub total: u64,
    pub qc_fail: u64,
    pub duplicate: u64,
    pub non_primary: u64,
    pub unmapped: u64,
    pub multi_hit: u64,
    pub uniq_hit: u64,
    pub read1: u64,
    pub read2: u64,
    pub forward: u64,
    pub reverse: u64,
    pub non_splice: u64,
    pub splice: u64,
    pub proper_pair: u64,
    pub proper_pair_diff_chrom: u64,
}

/// Computes bam_stat counts for a sequence of BAM records, matching the
/// branch order and category semantics of upstream's `ParseBAM.stat()`.
/// Pure computation: no I/O, no printing.
pub fn compute_stats<I>(records: I, q_cut: u8) -> io::Result<BamStatCounts>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut counts = BamStatCounts::default();

    for result in records {
        let record = result?;
        counts.total += 1;

        let flags = record.flags();

        if flags.is_qc_fail() {
            counts.qc_fail += 1;
            continue;
        }
        if flags.is_duplicate() {
            counts.duplicate += 1;
            continue;
        }
        if flags.is_secondary() {
            counts.non_primary += 1;
            continue;
        }
        if flags.is_unmapped() {
            counts.unmapped += 1;
            continue;
        }

        // pysam's raw `.mapq` reports a missing MAPQ as byte value 255;
        // noodles represents "missing" as `None`. 255 always compares >=
        // any realistic q_cut, so treating `None` as 255 here reproduces
        // upstream's comparison behavior without a special case.
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);

        if mapq < q_cut {
            counts.multi_hit += 1;
            continue;
        }
        counts.uniq_hit += 1;

        if flags.is_first_segment() {
            counts.read1 += 1;
        }
        if flags.is_last_segment() {
            counts.read2 += 1;
        }
        if flags.is_reverse_complemented() {
            counts.reverse += 1;
        } else {
            counts.forward += 1;
        }

        let has_skip = record
            .cigar()
            .iter()
            .any(|op| matches!(op, Ok(op) if op.kind() == CigarOpKind::Skip));
        if has_skip {
            counts.splice += 1;
        } else {
            counts.non_splice += 1;
        }

        if flags.is_properly_segmented() {
            counts.proper_pair += 1;

            let ref_id = record.reference_sequence_id().transpose()?;
            let mate_ref_id = record.mate_reference_sequence_id().transpose()?;
            if let (Some(r), Some(m)) = (ref_id, mate_ref_id) {
                if r != m {
                    counts.proper_pair_diff_chrom += 1;
                }
            }
        }
    }

    Ok(counts)
}

/// Runs the `bam_stat.py` CLI body over an already-opened record stream
/// (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// progress lines, same report bytes in the same order, same error
/// propagation -- except the record source is a caller-supplied iterator
/// and stdout/stderr are caller-supplied sinks instead of the process
/// globals. The standalone binary delegates to this (passing the process
/// streams); `rseqc_multi` passes one record broadcast plus per-command
/// stream files. `compute_stats` itself is untouched.
pub fn run_bam_stat<I>(
    records: I,
    q_cut: u8,
    stdout: &mut dyn io::Write,
    stderr: &mut dyn io::Write,
) -> io::Result<()>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    // Upstream: `if self.bam_format: print("Load BAM file ... ",
    // end=' ') else: print("Load SAM file ... ", end=' ')` --
    // `self.bam_format` comes from trying `pysam.Samfile(path, 'rb')`
    // FIRST and only falling back to `'r'` (bam_format=False) if that
    // raises. Confirmed via a live diff against real upstream: htslib's
    // `'rb'` open is lenient about actual content and succeeds for a
    // genuine plain-text SAM file too (it auto-detects format,
    // effectively ignoring the 'b' mode hint) -- so `bam_format` is
    // `True`, and "Load BAM file" prints, EVEN for `.sam` input. The
    // "Load SAM file" branch is practically dead code for any valid
    // input, not something this port needs a format check to trigger.
    // The literal's own trailing space plus `end=' '` gives two spaces
    // before "Done".
    write!(stderr, "Load BAM file ...  ")?;
    let counts = compute_stats(records, q_cut)?;
    writeln!(stderr, "Done")?;
    print_report(&counts, stdout);
    Ok(())
}

fn print_report(c: &BamStatCounts, out: &mut dyn io::Write) {
    writeln!(out).unwrap();
    writeln!(out, "#==================================================").unwrap();
    writeln!(out, "#All numbers are READ count").unwrap();
    writeln!(out, "#==================================================").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "{:<40}{}", "Total records:", c.total).unwrap();
    writeln!(out).unwrap();
    writeln!(out, "{:<40}{}", "QC failed:", c.qc_fail).unwrap();
    writeln!(out, "{:<40}{}", "Optical/PCR duplicate:", c.duplicate).unwrap();
    writeln!(out, "{:<40}{}", "Non primary hits", c.non_primary).unwrap();
    writeln!(out, "{:<40}{}", "Unmapped reads:", c.unmapped).unwrap();
    writeln!(out, "{:<40}{}", "mapq < mapq_cut (non-unique):", c.multi_hit).unwrap();
    writeln!(out).unwrap();
    writeln!(out, "{:<40}{}", "mapq >= mapq_cut (unique):", c.uniq_hit).unwrap();
    writeln!(out, "{:<40}{}", "Read-1:", c.read1).unwrap();
    writeln!(out, "{:<40}{}", "Read-2:", c.read2).unwrap();
    writeln!(out, "{:<40}{}", "Reads map to '+':", c.forward).unwrap();
    writeln!(out, "{:<40}{}", "Reads map to '-':", c.reverse).unwrap();
    writeln!(out, "{:<40}{}", "Non-splice reads:", c.non_splice).unwrap();
    writeln!(out, "{:<40}{}", "Splice reads:", c.splice).unwrap();
    writeln!(out, "{:<40}{}", "Reads mapped in proper pairs:", c.proper_pair).unwrap();
    writeln!(
        out,
        "{:<40}{}",
        "Proper-paired reads map to different chrom:", c.proper_pair_diff_chrom
    )
    .unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        self as sam,
        alignment::{
            io::Write as _,
            record::{Flags, MappingQuality, cigar::Op},
            record_buf::{Cigar, RecordBuf},
        },
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .add_reference_sequence(
                "chr2",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build()
    }

    /// Round-trips hand-built records through a real BAM byte stream so the
    /// test exercises the same decode path `compute_stats` runs against in
    /// production, not just in-memory builder state.
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
    fn matches_hand_computed_counts() {
        let header = test_header();

        let normal = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 50)]))
            .build();

        let qc_fail = RecordBuf::builder().set_flags(Flags::QC_FAIL).build();
        let duplicate = RecordBuf::builder().set_flags(Flags::DUPLICATE).build();
        let secondary = RecordBuf::builder().set_flags(Flags::SECONDARY).build();
        let unmapped = RecordBuf::builder().set_flags(Flags::UNMAPPED).build();

        let spliced = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(CigarOpKind::Match, 20),
                Op::new(CigarOpKind::Skip, 100),
                Op::new(CigarOpKind::Match, 20),
            ]))
            .build();

        let low_mapq = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(10).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 30)]))
            .build();

        let proper_pair_diff_chrom = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_mate_reference_sequence_id(1)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 30)]))
            .build();

        let records = vec![
            normal,
            qc_fail,
            duplicate,
            secondary,
            unmapped,
            spliced,
            low_mapq,
            proper_pair_diff_chrom,
        ];

        let bam_records = to_bam_records(&header, &records);
        let counts = compute_stats(bam_records.into_iter().map(Ok), 30).unwrap();

        // Hand-computed from the fixture above, independent of any Python run:
        // 8 records in; qc_fail/duplicate/non_primary/unmapped/low_mapq each
        // "continue" out (5), leaving normal+spliced+proper_pair_diff_chrom (3)
        // in the uniq_hit branch.
        let expected = BamStatCounts {
            total: 8,
            qc_fail: 1,
            duplicate: 1,
            non_primary: 1,
            unmapped: 1,
            multi_hit: 1,
            uniq_hit: 3,
            read1: 1,
            read2: 0,
            forward: 3,
            reverse: 0,
            non_splice: 2,
            splice: 1,
            proper_pair: 1,
            proper_pair_diff_chrom: 1,
        };

        assert_eq!(counts, expected);
    }
}
