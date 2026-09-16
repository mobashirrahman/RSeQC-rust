//! Dispatch and flag parsing for `bam2fq.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! `-c/--compress` (gzip output) is not implemented yet — disclosed gap, see
//! crates/commands/src/bam2fq.rs module docs.

use std::fs::File;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::bam2fq::{write_paired, write_single};

#[derive(Parser)]
#[command(
    name = "bam2fq.py",
    about = "Convert alignments in BAM or SAM format to FASTQ."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output FASTQ file(s).
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Treat the input as single-end sequencing data.
    #[arg(short = 's', long = "single-end")]
    single_end: bool,
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
    let (mut reader, _header) = rseqc_formats::open_bam(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy();

    if args.single_end {
        let mut out = File::create(format!("{prefix}.fastq"))?;
        let counts = write_single(reader.records(), &mut out)?;
        eprintln!("{counts:?}");
    } else {
        let mut out1 = File::create(format!("{prefix}.R1.fastq"))?;
        let mut out2 = File::create(format!("{prefix}.R2.fastq"))?;
        let counts = write_paired(reader.records(), &mut out1, &mut out2)?;
        eprintln!("{counts:?}");
    }

    Ok(())
}
