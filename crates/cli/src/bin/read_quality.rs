//! Dispatch and flag parsing for `read_quality.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! Output: .qual.r script for generating quality plots.

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_quality::{compute_quality, render_qual_r_script};

#[derive(Parser)]
#[command(
    name = "read_quality.py",
    about = "Calculate per-cycle Phred quality-score distributions for aligned reads."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Reduce data for boxplot precision (default: 1).
    #[arg(short = 'r', long = "reduce", default_value_t = 1)]
    reduce: u64,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,
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
    let hist = compute_quality(reader.records(), args.mapq)?;
    
    // Write the R script to the output file
    let r_output_path = format!("{}.qual.r", args.out_prefix.to_string_lossy());
    let mut r_output_file = File::create(&r_output_path)?;
    
    let r_script_content = render_qual_r_script(&hist, args.reduce, &args.out_prefix.to_string_lossy());
    r_output_file.write_all(r_script_content.as_bytes())?;
    
    eprintln!("R script written to: {}", r_output_path);
    
    Ok(())
}