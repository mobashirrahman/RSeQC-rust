//! Port of `split_paired_bam.py`: split a paired-end BAM into read-1,
//! read-2, and unmapped output BAMs. Contract: see
//! `compatibility/commands.yaml` entry `split_paired_bam.py`; algorithm
//! ported from `split_paired_bam()`/`to_single_end_alignment()`/
//! `single_end_flag()` in `oracle/upstream-src/scripts/split_paired_bam.py`
//! (lines 187-280).
//!
//! `.bai` index generation (DIV-0006, closed) is implemented at the CLI
//! layer via `rseqc_formats::write_bai_index`; this module's own
//! `split_paired_bam` function is unaffected -- it only ever writes the
//! plain BAM records.

use std::io;

use noodles_bam as bam;
use noodles_sam::{
    self as sam, alignment::io::Write as _, alignment::record::Flags,
    alignment::record_buf::RecordBuf,
};

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SplitCounts {
    pub total: u64,
    pub unmapped: u64,
    pub read1: u64,
    pub read2: u64,
}

/// Rebuilds the upstream "historical single-end flag": only
/// reverse/secondary/qc-fail/duplicate survive from the original record: all
/// other bits (paired-ness, mate info, proper-pair, segment number) are
/// dropped, matching `single_end_flag()` in the Python.
fn single_end_flags(old: Flags) -> Flags {
    let mut bits: u16 = 0;
    if old.is_reverse_complemented() {
        bits |= 0x0010;
    }
    if old.is_secondary() {
        bits |= 0x0100;
    }
    if old.is_qc_fail() {
        bits |= 0x0200;
    }
    if old.is_duplicate() {
        bits |= 0x0400;
    }
    Flags::from_bits_truncate(bits)
}

/// Splits `records` (in file order) across three BAM writers using the same
/// header/template as the input. Pure logic over generic writers so it's
/// testable with in-memory sinks; no file handling here.
pub fn split_paired_bam<I, W1, W2, W3>(
    records: I,
    header: &sam::Header,
    read1_out: &mut bam::io::Writer<W1>,
    read2_out: &mut bam::io::Writer<W2>,
    unmapped_out: &mut bam::io::Writer<W3>,
) -> io::Result<SplitCounts>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
    W1: io::Write,
    W2: io::Write,
    W3: io::Write,
{
    let mut counts = SplitCounts::default();

    for result in records {
        let old_record = result?;
        counts.total += 1;

        let old_flags = old_record.flags();

        if old_flags.is_unmapped() {
            // Written as-is, unmodified (upstream: `unmapped_output.write(old_alignment)`).
            unmapped_out.write_alignment_record(header, &old_record)?;
            counts.unmapped += 1;
            continue;
        }

        let mut new_record = RecordBuf::try_from_alignment_record(header, &old_record)?;
        *new_record.flags_mut() = single_end_flags(old_flags);
        // Upstream builds a *fresh* `pysam.AlignedSegment(header)` in
        // `to_single_end_alignment()` and never sets mate reference id, mate
        // position, or template length on it, so they keep pysam's
        // defaults: RNEXT `*` (None), PNEXT 0 (None here), TLEN 0. Our
        // `try_from_alignment_record` above instead copies those fields
        // from the paired input record, so they must be reset explicitly.
        *new_record.mate_reference_sequence_id_mut() = None;
        *new_record.mate_alignment_start_mut() = None;
        *new_record.template_length_mut() = 0;

        if old_flags.is_first_segment() {
            read1_out.write_alignment_record(header, &new_record)?;
            counts.read1 += 1;
        } else {
            // Preserve exact upstream behavior: any mapped record not
            // marked read1 goes to read2_out (no separate is_read2 branch).
            read2_out.write_alignment_record(header, &new_record)?;
            counts.read2 += 1;
        }
    }

    Ok(counts)
}

/// Ports `print_report()`'s exact stdout format (upstream: `print(f"{'Total
/// records:':<55}{counts.total}")` etc, one call per line) -- each label is
/// left-padded to 55 characters before the count, matching
/// `split_bam.rs::render_report`'s identical convention.
pub fn render_report(read1_bam: &str, read2_bam: &str, unmapped_bam: &str, counts: &SplitCounts) -> String {
    let total_label = "Total records:".to_string();
    let r1_label = format!("{read1_bam} (Read 1):");
    let r2_label = format!("{read2_bam} (Read 2):");
    let unmap_label = format!("{unmapped_bam} (Unmapped):");
    format!(
        "{total_label:<55}{}\n{r1_label:<55}{}\n{r2_label:<55}{}\n{unmap_label:<55}{}\n",
        counts.total, counts.read1, counts.read2, counts.unmapped,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        alignment::record::{MappingQuality, cigar::Op, cigar::op::Kind},
        alignment::record_buf::Cigar,
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
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
    fn splits_and_rewrites_flags() {
        let header = test_header();

        let unmapped = RecordBuf::builder().set_flags(Flags::UNMAPPED).build();

        let read1_rev_dup = RecordBuf::builder()
            .set_flags(
                Flags::SEGMENTED
                    | Flags::FIRST_SEGMENT
                    | Flags::REVERSE_COMPLEMENTED
                    | Flags::DUPLICATE,
            )
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let read2_plain = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let input = vec![unmapped, read1_rev_dup, read2_plain];
        let bam_records = to_bam_records(&header, &input);

        let mut r1_buf = Vec::new();
        let mut r2_buf = Vec::new();
        let mut unmap_buf = Vec::new();
        let mut r1_writer = bam::io::Writer::new(&mut r1_buf);
        let mut r2_writer = bam::io::Writer::new(&mut r2_buf);
        let mut unmap_writer = bam::io::Writer::new(&mut unmap_buf);
        r1_writer.write_header(&header).unwrap();
        r2_writer.write_header(&header).unwrap();
        unmap_writer.write_header(&header).unwrap();

        let counts = split_paired_bam(
            bam_records.into_iter().map(Ok),
            &header,
            &mut r1_writer,
            &mut r2_writer,
            &mut unmap_writer,
        )
        .unwrap();

        drop(r1_writer);
        drop(r2_writer);
        drop(unmap_writer);

        assert_eq!(
            counts,
            SplitCounts {
                total: 3,
                unmapped: 1,
                read1: 1,
                read2: 1,
            }
        );

        // Spot-check the rebuilt single-end flag bits on the R1 record:
        // reverse (0x0010) + duplicate (0x0400) should survive; segmented/
        // first-segment should NOT (upstream drops everything except
        // reverse/secondary/qc-fail/duplicate).
        let mut r1_reader = bam::io::Reader::new(r1_buf.as_slice());
        r1_reader.read_header().unwrap();
        let r1_record = r1_reader.records().next().unwrap().unwrap();
        let expected_flags = Flags::from_bits_truncate(0x0010 | 0x0400);
        assert_eq!(r1_record.flags(), expected_flags);
    }

    #[test]
    fn clears_mate_fields_like_fresh_pysam_record() {
        // Upstream's to_single_end_alignment() builds a brand-new
        // pysam.AlignedSegment(header) and never touches mate reference id,
        // mate position, or template length, so those keep pysam's
        // defaults (RNEXT `*`, PNEXT 0, TLEN 0). This guards against
        // silently carrying over the paired input record's mate fields.
        let header = test_header();

        let mut mated_read1 = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .build();
        *mated_read1.mate_reference_sequence_id_mut() = Some(0);
        *mated_read1.mate_alignment_start_mut() =
            Some(noodles_core::Position::try_from(100).unwrap());
        *mated_read1.template_length_mut() = 250;

        let bam_records = to_bam_records(&header, std::slice::from_ref(&mated_read1));

        let mut r1_buf = Vec::new();
        let mut r2_buf = Vec::new();
        let mut unmap_buf = Vec::new();
        let mut r1_writer = bam::io::Writer::new(&mut r1_buf);
        let mut r2_writer = bam::io::Writer::new(&mut r2_buf);
        let mut unmap_writer = bam::io::Writer::new(&mut unmap_buf);
        r1_writer.write_header(&header).unwrap();
        r2_writer.write_header(&header).unwrap();
        unmap_writer.write_header(&header).unwrap();

        split_paired_bam(
            bam_records.into_iter().map(Ok),
            &header,
            &mut r1_writer,
            &mut r2_writer,
            &mut unmap_writer,
        )
        .unwrap();

        drop(r1_writer);
        drop(r2_writer);
        drop(unmap_writer);

        let mut r1_reader = bam::io::Reader::new(r1_buf.as_slice());
        r1_reader.read_header().unwrap();
        let r1_record = r1_reader.records().next().unwrap().unwrap();

        assert_eq!(
            r1_record.mate_reference_sequence_id().transpose().unwrap(),
            None
        );
        assert_eq!(
            r1_record.mate_alignment_start().transpose().unwrap(),
            None
        );
        assert_eq!(r1_record.template_length(), 0);
    }

    #[test]
    fn render_report_exact_text() {
        // Verified via python3 -c against the literal upstream f-string
        // expressions in print_report().
        let counts = SplitCounts { total: 100, unmapped: 5, read1: 40, read2: 55 };
        let report = render_report("out.R1.bam", "out.R2.bam", "out.unmap.bam", &counts);
        let expected = "\
Total records:                                         100
out.R1.bam (Read 1):                                   40
out.R2.bam (Read 2):                                   55
out.unmap.bam (Unmapped):                              5
";
        assert_eq!(report, expected);
    }
}
