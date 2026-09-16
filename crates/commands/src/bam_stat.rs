//! Port of `bam_stat.py`: summarize mapping statistics for a BAM file.
//! Contract: see `compatibility/commands.yaml` entry `bam_stat.py`; algorithm
//! ported from `ParseBAM.stat()` in `oracle/upstream-src/src/qcmodule/SAM.py`.
//! Known divergences: see `compatibility/divergences.yaml` DIV-0001, DIV-0002
//! (splice detection simplified to a direct CIGAR-Skip check; BAM only, no
//! SAM-text support yet).

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
