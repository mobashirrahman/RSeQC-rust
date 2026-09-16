//! Port of `split_bam.py`: split a BAM file into exon-overlapping ("in"),
//! non-overlapping ("ex"), and junk (QC-fail/unmapped) BAM files, by
//! testing read-start and mate-start positions against exon regions from
//! a BED gene model. Algorithm ported from `classify_alignment` in
//! `oracle/upstream-src/scripts/split_bam.py` (lines 240-267).
//!
//! **Preserves the historical algorithm's specific quirks**: secondary
//! and duplicate alignments are NOT filtered (only qc-fail/unmapped route
//! to "junk"); when the mate is unmapped, only the read's own start is
//! tested (the mate start is meaningless in that case, so it's skipped
//! rather than tested against position 0); chromosome names are
//! uppercased on BOTH the exon-tree-build side and the query side
//! (`build_interval_trees`/`classify_alignment` both call `.upper()`) --
//! see the read_distribution.py case-normalization bug in project memory
//! for why both sides matter.
//!
//! Known gap: `--index-output` (BAI generation) is accepted by the CLI
//! but not implemented here -- disclosed follow-up, no BAI writer support
//! yet (same gap as divide_bam.py/RNA_fragment_size.py).

use std::io;

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};
use rseqc_formats::interval::{Bed3, MergedRegions};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    In = 0,
    Ex = 1,
    Junk = 2,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SplitCounts {
    pub total: u64,
    pub consumed: u64,
    pub excluded: u64,
    pub junk: u64,
}

/// Builds a chromosome-uppercased exon interval index from raw exon
/// `Bed3` triples (e.g. from `rseqc_formats::bed::get_exon`). Matches
/// upstream's `build_interval_trees`, which uppercases chromosome names
/// when building the tree.
pub fn build_exon_ranges(exons: &[Bed3]) -> MergedRegions {
    let uppercased: Vec<Bed3> = exons.iter().map(|(c, s, e)| (c.to_uppercase(), *s, *e)).collect();
    MergedRegions::new(&uppercased)
}

fn overlaps_exon_start(exon_ranges: &MergedRegions, chrom: &str, position: i64) -> bool {
    exon_ranges.overlap_length(chrom, position, position + 1) > 0
}

/// Classifies one alignment using the historical algorithm.
pub fn classify_alignment(
    record: &bam::Record,
    header: &sam::Header,
    exon_ranges: &MergedRegions,
) -> io::Result<Category> {
    let flags = record.flags();
    if flags.is_qc_fail() || flags.is_unmapped() {
        return Ok(Category::Junk);
    }

    let ref_id = record
        .reference_sequence_id()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "mapped record missing a reference sequence id"))?;
    let (chrom_name, _) = header
        .reference_sequences()
        .get_index(ref_id)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence id"))?;
    let chrom = String::from_utf8_lossy(chrom_name).to_uppercase();

    let read_start: i64 = record
        .alignment_start()
        .transpose()?
        .map(|p| usize::from(p) as i64 - 1)
        .unwrap_or(0);

    if flags.is_mate_unmapped() {
        return Ok(if overlaps_exon_start(exon_ranges, &chrom, read_start) {
            Category::In
        } else {
            Category::Ex
        });
    }

    let mate_start: i64 = record
        .mate_alignment_start()
        .transpose()?
        .map(|p| usize::from(p) as i64 - 1)
        .unwrap_or(0);

    if overlaps_exon_start(exon_ranges, &chrom, read_start) || overlaps_exon_start(exon_ranges, &chrom, mate_start) {
        Ok(Category::In)
    } else {
        Ok(Category::Ex)
    }
}

/// Splits `records` across the three output writers (index 0 = in, 1 =
/// ex, 2 = junk), in file order. Pure logic over generic writers so it's
/// testable with in-memory sinks; no file handling here.
pub fn split_bam<I, W>(
    records: I,
    header: &sam::Header,
    exon_ranges: &MergedRegions,
    outputs: &mut [bam::io::Writer<W>; 3],
) -> io::Result<SplitCounts>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
    W: io::Write,
{
    let mut counts = SplitCounts::default();

    for result in records {
        let record = result?;
        counts.total += 1;

        let category = classify_alignment(&record, header, exon_ranges)?;
        outputs[category as usize].write_alignment_record(header, &record)?;

        match category {
            Category::In => counts.consumed += 1,
            Category::Ex => counts.excluded += 1,
            Category::Junk => counts.junk += 1,
        }
    }

    Ok(counts)
}

/// Ports `print_report`'s literal f-string layout: each label is padded
/// (with Python's `:<55`, i.e. left-justify to a minimum width of 55,
/// never truncating a longer label) then immediately followed by the
/// count with no separator.
pub fn render_report(in_bam: &str, ex_bam: &str, junk_bam: &str, counts: &SplitCounts) -> String {
    let total_label = "Total records:".to_string();
    let in_label = format!("{in_bam} (Alignments consumed by input gene list):");
    let ex_label = format!("{ex_bam} (Alignments not consumed by input gene list):");
    let junk_label = format!("{junk_bam} (QC-failed, unmapped reads):");
    format!(
        "{total_label:<55}{}\n{in_label:<55}{}\n{ex_label:<55}{}\n{junk_label:<55}{}\n",
        counts.total, counts.consumed, counts.excluded, counts.junk,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        self as sam,
        alignment::record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind},
        alignment::record_buf::{Cigar, RecordBuf},
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()))
            .add_reference_sequence("chr2", Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()))
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

    fn exon_ranges() -> MergedRegions {
        // Exon on chr1 covering [100, 200) -- built with a lowercase
        // chrom name, verifying build_exon_ranges uppercases it.
        build_exon_ranges(&[("chr1".to_string(), 100, 200)])
    }

    #[test]
    fn junk_when_qc_fail_or_unmapped() {
        let header = test_header();
        let ranges = exon_ranges();

        let qc_fail = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::QC_FAIL)
            .build();
        let unmapped = RecordBuf::builder()
            .set_name("r2")
            .set_flags(Flags::UNMAPPED)
            .build();

        for rec in to_bam_records(&header, &[qc_fail, unmapped]) {
            assert_eq!(classify_alignment(&rec, &header, &ranges).unwrap(), Category::Junk);
        }
    }

    #[test]
    fn in_when_read_start_overlaps_exon_mate_unmapped() {
        let header = test_header();
        let ranges = exon_ranges();

        // 0-based read_start=150 (1-based alignment_start=151) is inside [100,200).
        let rec = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::MATE_UNMAPPED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(151).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let [bam_rec] = &to_bam_records(&header, &[rec])[..] else { unreachable!() };
        assert_eq!(classify_alignment(bam_rec, &header, &ranges).unwrap(), Category::In);
    }

    #[test]
    fn ex_when_neither_read_nor_mate_start_overlaps() {
        let header = test_header();
        let ranges = exon_ranges();

        // read_start (0-based 500) and mate_start (0-based 600) are both
        // outside [100, 200).
        let rec = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::SEGMENTED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(501).unwrap())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(noodles_core::Position::try_from(601).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let [bam_rec] = &to_bam_records(&header, &[rec])[..] else { unreachable!() };
        assert_eq!(classify_alignment(bam_rec, &header, &ranges).unwrap(), Category::Ex);
    }

    #[test]
    fn in_when_only_mate_start_overlaps() {
        let header = test_header();
        let ranges = exon_ranges();

        // read_start (0-based 500) is outside the exon, but mate_start
        // (0-based 150) is inside it -- the "or" branch of the check.
        let rec = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::SEGMENTED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(501).unwrap())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(noodles_core::Position::try_from(151).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let [bam_rec] = &to_bam_records(&header, &[rec])[..] else { unreachable!() };
        assert_eq!(classify_alignment(bam_rec, &header, &ranges).unwrap(), Category::In);
    }

    #[test]
    fn split_bam_counts_and_routes_records() {
        let header = test_header();
        let ranges = exon_ranges();

        let in_rec = RecordBuf::builder()
            .set_name("in_read")
            .set_flags(Flags::MATE_UNMAPPED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(151).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();
        let ex_rec = RecordBuf::builder()
            .set_name("ex_read")
            .set_flags(Flags::MATE_UNMAPPED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(501).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();
        let junk_rec = RecordBuf::builder().set_name("junk_read").set_flags(Flags::UNMAPPED).build();

        let bam_records = to_bam_records(&header, &[in_rec, ex_rec, junk_rec]);

        let mut in_buf = Vec::new();
        let mut ex_buf = Vec::new();
        let mut junk_buf = Vec::new();
        let mut in_writer = bam::io::Writer::new(&mut in_buf);
        let mut ex_writer = bam::io::Writer::new(&mut ex_buf);
        let mut junk_writer = bam::io::Writer::new(&mut junk_buf);
        in_writer.write_header(&header).unwrap();
        ex_writer.write_header(&header).unwrap();
        junk_writer.write_header(&header).unwrap();

        let mut outputs = [in_writer, ex_writer, junk_writer];
        let counts = split_bam(bam_records.into_iter().map(Ok), &header, &ranges, &mut outputs).unwrap();
        drop(outputs);

        assert_eq!(
            counts,
            SplitCounts { total: 3, consumed: 1, excluded: 1, junk: 1 }
        );

        let mut in_reader = bam::io::Reader::new(in_buf.as_slice());
        in_reader.read_header().unwrap();
        let in_records: Vec<_> = in_reader.records().map(|r| r.unwrap()).collect();
        assert_eq!(in_records.len(), 1);
        assert_eq!(in_records[0].name().unwrap().to_string(), "in_read");
    }

    #[test]
    fn render_report_exact_text() {
        let counts = SplitCounts { total: 100, consumed: 40, excluded: 55, junk: 5 };
        let report = render_report("out.in.bam", "out.ex.bam", "out.junk.bam", &counts);
        let expected = "\
Total records:                                         100
out.in.bam (Alignments consumed by input gene list):   40
out.ex.bam (Alignments not consumed by input gene list):55
out.junk.bam (QC-failed, unmapped reads):              5
";
        assert_eq!(report, expected);
    }
}
