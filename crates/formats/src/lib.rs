//! Alignment (SAM/BAM), annotation (BED12), sequence (FASTA/FASTQ), and
//! coverage (WIG/BigWig) I/O. Kept independent of command logic and CLI
//! parsing so `rseqc-commands` and `rseqc-python` share one implementation.

use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

use noodles_sam as sam;

pub mod bed;
pub mod bigwig;
pub mod cigar;
pub mod interval;

/// Opens a BAM file for sequential record reading, returning both the
/// reader (positioned at the first record) and the parsed header (needed by
/// callers that write records back out, e.g. via
/// `RecordBuf::try_from_alignment_record` or `write_alignment_record`).
pub fn open_bam(
    path: &Path,
) -> io::Result<(
    noodles_bam::io::Reader<noodles_bgzf::io::Reader<File>>,
    sam::Header,
)> {
    let mut reader = File::open(path).map(noodles_bam::io::Reader::new)?;
    let header = reader.read_header()?;
    Ok((reader, header))
}

/// Builds and writes a `<path>.bai` index for an already-written,
/// coordinate-sorted BAM file, matching upstream's `pysam.index(path)`
/// (`split_paired_bam.py`'s `index_bam`, `divide_bam.py`'s
/// `create_indexes`) -- closes DIV-0006 for any command that calls this
/// after finishing a BAM write. Requires the BAM header to declare
/// `SO:coordinate` (same requirement pysam/samtools enforce); on a
/// non-coordinate-sorted input this returns an `InvalidData` error with a
/// message deliberately worded like upstream's own
/// "could not index ...; output BAM may not be coordinate-sorted" wrapper,
/// since both sides reach the same practical conclusion (index-building
/// only works on sorted data) even though the underlying library error
/// text differs (noodles vs. htslib).
pub fn write_bai_index(path: &Path) -> io::Result<()> {
    let index = noodles_bam::fs::index(path).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("could not index {}; output BAM may not be coordinate-sorted: {err}", path.display()),
        )
    })?;
    let bai_path = format!("{}.bai", path.display());
    noodles_bam::bai::fs::write(bai_path, &index)
}

/// Opens either a BAM or plain-text SAM file for sequential record
/// reading, dispatching on the `.bam`/`.sam` extension (case-
/// insensitive) and returning a UNIFORM `(header, records)` pair
/// regardless of source format -- closing DIV-0002/0004 ("BAM or SAM"
/// input, upstream's own advertised contract via pysam) for any
/// command that switches its `open_bam` call to this function instead.
///
/// SAM-text records are converted to genuine `bam::Record`s by
/// round-tripping through an in-memory BAM byte buffer: write the
/// parsed text records out via the same `bam::io::Writer::
/// write_alignment_record` this port's own test helpers already use to
/// build BAM fixtures (it accepts anything implementing `sam::
/// alignment::Record`, including the text `sam::Record` -- no format-
/// specific encoding logic needed here), then read that buffer back
/// with a plain `bam::io::Reader`. This means every existing
/// `compute_*` function in `rseqc-commands` -- all already generic
/// over `IntoIterator<Item = io::Result<bam::Record>>` -- needs ZERO
/// changes to accept SAM-text input; only a CLI's own `open_bam` call
/// site needs to switch to this function and its module doc comment's
/// "SAM-text input is not yet supported" note needs to come off.
///
/// The whole file is decoded eagerly into memory (not streamed) for
/// both formats, unlike `open_bam`'s lazy reader -- acceptable for the
/// small/QC-scale inputs this tool targets, and it avoids needing a
/// second, owned-iterator BAM reader type alongside the existing
/// borrowed-iterator one. Extension detection (not content sniffing)
/// matches upstream's own `pysam.AlignmentFile` behavior, which also
/// dispatches BAM-vs-SAM-text by filename, not by sniffing magic bytes.
pub fn open_alignments(path: &Path) -> io::Result<(sam::Header, Vec<io::Result<noodles_bam::Record>>)> {
    let is_sam_text = path.extension().and_then(|ext| ext.to_str()).map(|ext| ext.eq_ignore_ascii_case("sam")).unwrap_or(false);

    if is_sam_text {
        use sam::alignment::io::Write as _;

        let mut text_reader = File::open(path).map(BufReader::new).map(sam::io::Reader::new)?;
        let header = text_reader.read_header()?;

        let mut buf = Vec::new();
        {
            let mut bam_writer = noodles_bam::io::Writer::new(&mut buf);
            bam_writer.write_header(&header)?;
            for result in text_reader.records() {
                let record = result?;
                bam_writer.write_alignment_record(&header, &record)?;
            }
        }

        let mut bam_reader = noodles_bam::io::Reader::new(buf.as_slice());
        bam_reader.read_header()?;
        let records: Vec<_> = bam_reader.records().collect();
        Ok((header, records))
    } else {
        let (mut reader, header) = open_bam(path)?;
        let records: Vec<_> = reader.records().collect();
        Ok((header, records))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind as CigarOpKind};
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()))
            .build()
    }

    #[test]
    fn open_alignments_decodes_sam_text_identically_to_bam() {
        use noodles_sam::alignment::io::Write as _;

        let header = test_header();
        let record = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(11).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![40; 4]))
            .build();

        let dir = std::env::temp_dir().join(format!("rseqc_formats_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        // Write the same logical record to both a real .bam and a real
        // .sam (plain text) file.
        let bam_path = dir.join("t.bam");
        {
            let mut w = noodles_bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            w.write_header(&header).unwrap();
            w.write_alignment_record(&header, &record).unwrap();
        }

        let sam_path = dir.join("t.sam");
        {
            let mut w = sam::io::Writer::new(std::fs::File::create(&sam_path).unwrap());
            w.write_header(&header).unwrap();
            w.write_alignment_record(&header, &record).unwrap();
        }

        let (bam_header, bam_records) = open_alignments(&bam_path).unwrap();
        let (sam_header, sam_records) = open_alignments(&sam_path).unwrap();

        assert_eq!(bam_header, sam_header);
        assert_eq!(bam_records.len(), 1);
        assert_eq!(sam_records.len(), 1);

        let b = bam_records.into_iter().next().unwrap().unwrap();
        let s = sam_records.into_iter().next().unwrap().unwrap();
        assert_eq!(b.name(), s.name());
        assert_eq!(b.flags(), s.flags());
        assert_eq!(b.mapping_quality(), s.mapping_quality());
        assert_eq!(b.alignment_start().unwrap().unwrap(), s.alignment_start().unwrap().unwrap());

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn write_bai_index_produces_a_real_readable_index() {
        use noodles_sam::alignment::io::Write as _;
        use noodles_sam::header::record::value::map;
        use noodles_sam::header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER};

        let header = sam::Header::builder()
            .set_header(Map::<map::Header>::builder().insert(SORT_ORDER, COORDINATE).build().unwrap())
            .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()))
            .build();

        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(11).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .build();

        let dir = std::env::temp_dir().join(format!("rseqc_formats_bai_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("sorted.bam");
        {
            let mut w = noodles_bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            w.write_header(&header).unwrap();
            w.write_alignment_record(&header, &record).unwrap();
        }

        write_bai_index(&bam_path).unwrap();

        let bai_path = dir.join("sorted.bam.bai");
        assert!(bai_path.is_file());
        let index = noodles_bam::bai::fs::read(&bai_path).unwrap();
        assert_eq!(index.reference_sequences().len(), 1);

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn write_bai_index_rejects_unsorted_header() {
        use noodles_sam::alignment::io::Write as _;

        // No SO:coordinate tag set -- matches upstream's own "may not be
        // coordinate-sorted" rejection, though via a different underlying
        // check (noodles requires the header tag; htslib/pysam inspects
        // actual record order). Both reach the same practical outcome.
        let header = test_header();
        let record = RecordBuf::builder().set_flags(Flags::UNMAPPED).build();

        let dir = std::env::temp_dir().join(format!("rseqc_formats_bai_unsorted_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("unsorted.bam");
        {
            let mut w = noodles_bam::io::Writer::new(std::fs::File::create(&bam_path).unwrap());
            w.write_header(&header).unwrap();
            w.write_alignment_record(&header, &record).unwrap();
        }

        let err = write_bai_index(&bam_path).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);

        std::fs::remove_dir_all(&dir).unwrap();
    }
}
