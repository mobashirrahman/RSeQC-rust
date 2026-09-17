//! Alignment (SAM/BAM), annotation (BED12), sequence (FASTA/FASTQ), and
//! coverage (WIG/BigWig) I/O. Kept independent of command logic and CLI
//! parsing so `rseqc-commands` and `rseqc-python` share one implementation.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use noodles_sam as sam;

pub mod bed;
pub mod bigwig;
pub mod cigar;
pub mod interval;

/// Opens a plain-text, gzip-, or bzip2-compressed file for line-buffered
/// reading, dispatching on the exact (CASE-SENSITIVE) trailing extension
/// -- matching upstream's `qcmodule.ireader.nopen`:
/// `f.endswith((".gz", ".Z", ".z"))` routes to `gzip.open` (so a literal
/// `.Z`-suffixed file, historically the `compress`(1) utility's own
/// extension and NOT actually gzip format, is still routed to a gzip
/// decoder by upstream itself -- reproduced exactly here, including that
/// it will equally fail on a genuine LZW `.Z` file on both sides, since
/// that's upstream's own established behavior, not a bug this port
/// should "fix"); `f.endswith((".bz", ".bz2", ".bzip2"))` routes to
/// `bz2.BZ2File`. Anything else is opened as plain text. Used by
/// `sc_seqQual.py`/`sc_seqLogo.py` (DIV-0002-adjacent input-format gap,
/// tracked in compatibility/divergences.yaml's compressed-input note).
///
/// `flate2`'s `MultiGzDecoder` (not plain `GzDecoder`) matches Python's
/// `gzip` module's own transparent support for concatenated
/// (multi-member) gzip streams. `bzip2-rs` is a decode-only, pure-Rust
/// bzip2 implementation (no system `libbz2` dependency, unlike the
/// `bzip2` crate's C bindings).
pub fn open_text_input(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let path_str = path.to_string_lossy();
    let file = File::open(path)?;
    if path_str.ends_with(".gz") || path_str.ends_with(".Z") || path_str.ends_with(".z") {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))))
    } else if path_str.ends_with(".bz") || path_str.ends_with(".bz2") || path_str.ends_with(".bzip2") {
        Ok(Box::new(BufReader::new(bzip2_rs::DecoderReader::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

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

/// Opens a BAM, plain-text SAM, or CRAM file for sequential record
/// reading, dispatching on the `.bam`/`.sam`/`.cram` extension (case-
/// insensitive) and returning a UNIFORM `(header, records)` pair
/// regardless of source format -- closing DIV-0002/0004 ("BAM, SAM, or
/// CRAM" input, upstream's own advertised contract via pysam) for any
/// command that switches its `open_bam` call to this function instead.
///
/// SAM-text and CRAM records are both converted to genuine
/// `bam::Record`s by round-tripping through an in-memory BAM byte
/// buffer: write the parsed records out via the same `bam::io::
/// Writer::write_alignment_record` this port's own test helpers
/// already use to build BAM fixtures (it accepts anything implementing
/// `sam::alignment::Record` -- the text `sam::Record` type AND CRAM's
/// `sam::alignment::RecordBuf`, no format-specific encoding logic
/// needed here), then read that buffer back with a plain
/// `bam::io::Reader`. This means every existing `compute_*` function
/// in `rseqc-commands` -- all already generic over
/// `IntoIterator<Item = io::Result<bam::Record>>` -- needs ZERO
/// changes to accept SAM-text or CRAM input; only a CLI's own
/// `open_bam` call site needs to switch to this function.
///
/// **CRAM's reference-sequence handling, a real scope limit, disclosed
/// rather than silently wrong**: CRAM decodes with `noodles_cram`'s
/// DEFAULT (empty) reference-sequence repository -- no external FASTA
/// is consulted. This correctly decodes CRAM written with an embedded
/// or no-reference-required encoding (confirmed via a real fixture:
/// `pysam.AlignmentFile(path, 'wc', ...)` without an explicit
/// `reference_filename` falls back to `embed_ref=2` -- htslib's own
/// term for "embed the reference in the CRAM file itself" -- when no
/// external reference is configured, which is exactly the case an
/// empty repository can decode). A CRAM file that genuinely requires
/// EXTERNAL reference resolution (encoded against a reference NOT
/// embedded and not supplied here) will surface as a decode error
/// instead of silently producing wrong sequence data. None of the 12
/// upstream commands that advertise `.cram` input expose a
/// `--reference`-style flag of their own either (checked via grep
/// across their argparse setups) -- they rely on pysam/htslib's own
/// reference resolution, which for files lacking a local/embedded
/// reference can fall back to fetching from a remote EBI/ENA reference
/// server over the network. Deliberately NOT replicated: this project
/// is offline-first by design (see README), and network-dependent,
/// non-reproducible reference fetching would be a poor fit for a QC
/// tool's I/O layer regardless of upstream's own behavior here.
///
/// The whole file is decoded eagerly into memory (not streamed) for
/// all three formats, unlike `open_bam`'s lazy reader -- acceptable
/// for the small/QC-scale inputs this tool targets, and it avoids
/// needing a second, owned-iterator BAM reader type alongside the
/// existing borrowed-iterator one. Extension detection (not content
/// sniffing) matches upstream's own `pysam.AlignmentFile` behavior,
/// which also dispatches by filename, not by sniffing magic bytes.
pub fn open_alignments(path: &Path) -> io::Result<(sam::Header, Vec<io::Result<noodles_bam::Record>>)> {
    let extension = path.extension().and_then(|ext| ext.to_str()).map(|ext| ext.to_ascii_lowercase());

    match extension.as_deref() {
        Some("sam") => {
            let mut text_reader = File::open(path).map(BufReader::new).map(sam::io::Reader::new)?;
            let header = text_reader.read_header()?;
            round_trip_through_bam(header, text_reader.records())
        }
        Some("cram") => {
            let mut cram_reader = File::open(path).map(noodles_cram::io::Reader::new)?;
            let header = cram_reader.read_header()?;
            // `records()` borrows `header`, so collect eagerly before
            // it's moved into `round_trip_through_bam`.
            let records: Vec<io::Result<sam::alignment::RecordBuf>> = cram_reader
                .records(&header)
                .map(|result| result.map(fix_unmapped_missing_mapping_quality))
                .collect();
            round_trip_through_bam(header, records.into_iter())
        }
        _ => {
            let (mut reader, header) = open_bam(path)?;
            let records: Vec<_> = reader.records().collect();
            Ok((header, records))
        }
    }
}

/// Real, live-diff-discovered CRAM/htslib interop quirk, not a bug in
/// this port's own logic: for an UNMAPPED read, `noodles_cram` decodes
/// mapping quality as genuinely MISSING (`None`) when that's what the
/// CRAM container's own MAPQ data series actually stores for that
/// record -- which is what htslib's own CRAM WRITER puts there for
/// unmapped reads (confirmed: converting `bam_stat_basic.bam`, whose
/// `unmapped1` read has an EXPLICIT MAPQ of `0`, to CRAM via
/// `pysam.AlignmentFile(..., 'wc', ...)` and back loses that `0`,
/// because htslib's writer re-encodes unmapped reads' MAPQ as missing
/// regardless of the original value). But htslib's own CRAM READER
/// (what pysam/upstream actually uses) does NOT surface that as
/// "missing" to callers -- it reports mapping quality `0` for an
/// unmapped read with a missing MAPQ data series entry, a read-time
/// convenience default. `noodles_cram` faithfully reports what's
/// actually stored (missing) instead of replicating htslib's own
/// asymmetric write/read convention. Discovered via a live diff on
/// `read_NVC.py`: this port's CRAM path counted `unmapped1` as
/// "qualifying" (missing MAPQ treated as BAM's own 255/"unavailable"
/// sentinel, which always clears any `--mapq` threshold) while real
/// upstream correctly excluded it (mapq 0 fails the default `-q 30`
/// cutoff). Reproduces htslib's own read-time default here so this
/// port's CRAM support matches what upstream ACTUALLY does, not just
/// what the CRAM container's raw bytes technically encode.
fn fix_unmapped_missing_mapping_quality(mut record: sam::alignment::RecordBuf) -> sam::alignment::RecordBuf {
    use sam::alignment::record::{Flags, MappingQuality};
    if record.flags().contains(Flags::UNMAPPED) && record.mapping_quality().is_none() {
        *record.mapping_quality_mut() = MappingQuality::new(0);
    }
    record
}

/// Shared by `open_alignments`' SAM-text and CRAM branches: writes any
/// `sam::alignment::Record`-implementing records out to an in-memory
/// BAM buffer, then reads them back as genuine `bam::Record`s.
fn round_trip_through_bam<R>(header: sam::Header, records: impl Iterator<Item = io::Result<R>>) -> io::Result<(sam::Header, Vec<io::Result<noodles_bam::Record>>)>
where
    R: sam::alignment::Record,
{
    use sam::alignment::io::Write as _;

    let mut buf = Vec::new();
    {
        let mut bam_writer = noodles_bam::io::Writer::new(&mut buf);
        bam_writer.write_header(&header)?;
        for result in records {
            let record = result?;
            bam_writer.write_alignment_record(&header, &record)?;
        }
    }

    let mut bam_reader = noodles_bam::io::Reader::new(buf.as_slice());
    bam_reader.read_header()?;
    let records: Vec<_> = bam_reader.records().collect();
    Ok((header, records))
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

    /// Real, pysam-written CRAM with no external reference configured
    /// -- htslib itself falls back to embedding the reference in this
    /// situation (confirmed via the warning it prints when writing;
    /// see the fixture's own generator script), which is exactly what
    /// `open_alignments`'s default (empty) reference-sequence
    /// repository can decode. Noodles' OWN `cram::io::Writer` was not
    /// usable to build this fixture directly: unlike pysam/htslib, it
    /// hard-requires a real `fasta::Repository` to compute `@SQ` `M5`
    /// checksums and panics without one -- a real fixture generated by
    /// the actual tool this port needs to interoperate with is more
    /// representative here anyway, matching this project's established
    /// precedent for `crates/formats/tests/fixtures/` (see
    /// `pybigwig_test.bw`'s own README entry).
    #[test]
    fn open_alignments_decodes_a_real_no_reference_cram_fixture() {
        let path = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/cram_no_reference.cram"));
        let (header, records) = open_alignments(path).unwrap();

        assert_eq!(header.reference_sequences().len(), 1);
        assert_eq!(records.len(), 1);

        let record = records.into_iter().next().unwrap().unwrap();
        assert_eq!(record.name().map(|n| n.to_vec()), Some(b"r1".to_vec()));
        assert_eq!(record.mapping_quality().map(|q| q.get()), Some(40));
        assert_eq!(record.alignment_start().unwrap().unwrap().get(), 11);
        assert_eq!(record.sequence().iter().collect::<Vec<_>>(), b"ACGT".to_vec());
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
