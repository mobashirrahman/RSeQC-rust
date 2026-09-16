//! Dispatch and flag parsing for `split_paired_bam.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10)
//! adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::path::PathBuf;

use clap::Parser;
use noodles_bam as bam;
use rseqc_commands::split_paired_bam::{SplitCounts, split_paired_bam};

#[derive(Parser)]
#[command(
    name = "split_paired_bam.py",
    about = "Split a paired-end BAM into read-1, read-2, and unmapped BAM files."
)]
struct Args {
    /// Input BAM file.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the three output BAM files
    /// (`<prefix>.R1.bam`, `<prefix>.R2.bam`, `<prefix>.unmap.bam`).
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(counts) => {
            eprintln!("{counts:?}");
            std::process::ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<SplitCounts> {
    let (reader, header) = rseqc_formats::open_bam(&args.input_file)?;

    let prefix = args.out_prefix.to_string_lossy();
    let r1_path = format!("{prefix}.R1.bam");
    let r2_path = format!("{prefix}.R2.bam");
    let unmap_path = format!("{prefix}.unmap.bam");

    let mut r1_writer = bam::io::Writer::new(File::create(&r1_path)?);
    let mut r2_writer = bam::io::Writer::new(File::create(&r2_path)?);
    let mut unmap_writer = bam::io::Writer::new(File::create(&unmap_path)?);
    r1_writer.write_header(&header)?;
    r2_writer.write_header(&header)?;
    unmap_writer.write_header(&header)?;

    let mut reader = reader;
    split_paired_bam(
        reader.records(),
        &header,
        &mut r1_writer,
        &mut r2_writer,
        &mut unmap_writer,
    )
}
