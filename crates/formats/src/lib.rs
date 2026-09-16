//! Alignment (SAM/BAM), annotation (BED12), sequence (FASTA/FASTQ), and
//! coverage (WIG/BigWig) I/O. Kept independent of command logic and CLI
//! parsing so `rseqc-commands` and `rseqc-python` share one implementation.

use std::fs::File;
use std::io;
use std::path::Path;

/// Opens a BAM file for sequential record reading, consuming the header
/// (callers that need reference-sequence names should re-derive them from
/// the returned reader via further calls, not stored here).
pub fn open_bam(
    path: &Path,
) -> io::Result<noodles_bam::io::Reader<noodles_bgzf::io::Reader<File>>> {
    let mut reader = File::open(path).map(noodles_bam::io::Reader::new)?;
    reader.read_header()?;
    Ok(reader)
}
