//! Dispatch and flag parsing for `bam2fq.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! `-c/--compress` (DIV-0007, closed): gzips each output FASTQ file after
//! writing, matching upstream's `gzip_outputs`/`gzip_file` in
//! `oracle/upstream-src/scripts/bam2fq.py`. Uses `flate2` at compression
//! level 9 (Python's `gzip.open` default `compresslevel`). The `.gz`
//! container's own header bytes (embedded filename/mtime) are NOT
//! reproduced -- Python's `gzip.open(path, 'wb')` stamps the current wall
//! time into the header, so upstream's own output is already not
//! byte-reproducible run-to-run for the same reason; only the
//! *decompressed* content is guaranteed to match. SAM-text input
//! (DIV-0002/0004) is supported via `rseqc_formats::open_alignments`.

use std::fs::File;
use std::io::{Read as _, Write as _};
use std::path::{Path, PathBuf};

use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;
use rseqc_commands::bam2fq::{write_paired, write_single};

#[derive(Parser)]
#[command(
    name = "bam2fq.py",
    about = "Convert alignments in BAM or SAM format to FASTQ."
)]
struct Args {
    /// Input alignment file in BAM or plain-text SAM format (dispatched by
    /// the `.bam`/`.sam` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output FASTQ file(s).
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Treat the input as single-end sequencing data.
    #[arg(short = 's', long = "single-end")]
    single_end: bool,

    /// Compress output FASTQ files with gzip.
    #[arg(short = 'c', long = "compress")]
    compress: bool,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy();

    eprintln!("Convert BAM/SAM file into FASTQ format ...");

    let outputs: Vec<String> = if args.single_end {
        let path = format!("{prefix}.fastq");
        let mut out = File::create(&path)?;
        let counts = write_single(records, &mut out)?;
        eprintln!("Done");
        eprintln!("read count: {}", counts.single);
        vec![path]
    } else {
        let path1 = format!("{prefix}.R1.fastq");
        let path2 = format!("{prefix}.R2.fastq");
        let mut out1 = File::create(&path1)?;
        let mut out2 = File::create(&path2)?;
        let counts = write_paired(records, &mut out1, &mut out2)?;
        eprintln!("Done");
        eprintln!("read_1 count: {}", counts.read1);
        eprintln!("read_2 count: {}", counts.read2);
        vec![path1, path2]
    };

    if args.compress {
        for path in &outputs {
            eprintln!("Compressing {path} ...");
            gzip_file(Path::new(path))?;
        }
        eprintln!("Compression complete.");
    }

    Ok(())
}

/// Compresses `path` to `<path>.gz` and removes the original, matching
/// upstream's `gzip_file`. Compression level 9 matches Python's
/// `gzip.open` default `compresslevel`.
fn gzip_file(path: &Path) -> std::io::Result<()> {
    let mut source = File::open(path)?;
    let mut buf = Vec::new();
    source.read_to_end(&mut buf)?;

    let gz_path = format!("{}.gz", path.display());
    let mut encoder = GzEncoder::new(File::create(&gz_path)?, Compression::best());
    encoder.write_all(&buf)?;
    encoder.finish()?;

    std::fs::remove_file(path)?;
    Ok(())
}
