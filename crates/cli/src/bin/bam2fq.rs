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

// B4: link `rseqc_cli` so its `#[global_allocator]` (actionable
// out-of-memory message) applies to this binary too.
use rseqc_cli as _;
use std::fs::File;
use std::io::{BufWriter, Read as _, Write as _};
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
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
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
            eprintln!("bam2fq.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // Upstream's validate_args refuses an output prefix whose parent directory does
    // not exist, before any input is read. Omitting it here meant the whole
    // alignment was read and every metric computed, then discarded when the output
    // open failed with "No such file or directory (os error 2)" -- an error naming
    // neither the directory nor the flag, and indistinguishable from a missing
    // input. The shared helper keeps that check in one place so it cannot be
    // forgotten by the next binary.
    // DELIBERATELY ABSENT, replicating upstream. bam2wig.py and bam2fq.py are the
    // only two upstream scripts that take an output prefix and do NOT validate its
    // parent, so upstream reaches the work and then fails on the first output open
    // with `[Errno 2] No such file or directory: '.../out.wig'` and exit 1. Adding
    // the check here would be the port being better than the tool it ports: the same
    // command would then exit 2 with a different message on an input upstream
    // accepts, and every downstream comparison against upstream would diverge on a
    // case that is not a defect in the port. Recorded as DIV-0026.

    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy();

    eprintln!("Convert BAM/SAM file into FASTQ format ...");

    let outputs: Vec<String> = if args.single_end {
        let path = format!("{prefix}.fastq");
        // Buffer the output. Writing straight to a bare `File` issues one syscall per
        // `write_all`, and a FASTQ record needs 5 of them, so an unbuffered writer
        // turns 400k reads into ~5.6M write syscalls: measured at 4.6 s of system time
        // on 800k reads, which made this command 3x SLOWER than the Python reference
        // even though the reference also writes the same bytes. See
        // benchmarks/RESULTS.generated.md section 6.1.
        let mut out = BufWriter::with_capacity(1 << 20, File::create(&path)?);
        let counts = write_single(records, &mut out)?;
        out.flush()?;
        eprintln!("Done");
        eprintln!("read count: {}", counts.single);
        vec![path]
    } else {
        let path1 = format!("{prefix}.R1.fastq");
        let path2 = format!("{prefix}.R2.fastq");
        let mut out1 = BufWriter::with_capacity(1 << 20, File::create(&path1)?);
        let mut out2 = BufWriter::with_capacity(1 << 20, File::create(&path2)?);
        let counts = write_paired(records, &mut out1, &mut out2)?;
        out1.flush()?;
        out2.flush()?;
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
