//! Alignment (SAM/BAM), annotation (BED12), sequence (FASTA/FASTQ), and
//! coverage (WIG/BigWig) I/O. Kept independent of command logic and CLI
//! parsing so `rseqc-commands` and `rseqc-python` share one implementation.

use std::fs::File;
use std::io;
use std::path::Path;

use noodles_sam as sam;

pub mod bed;
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
