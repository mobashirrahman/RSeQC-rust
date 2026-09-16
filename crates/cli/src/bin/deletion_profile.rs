//! Dispatch and flag parsing for `deletion_profile.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::deletion_profile::{
    compute_deletion_profile, render_deletion_r_script, render_deletion_table,
};

#[derive(Parser)]
#[command(
    name = "deletion_profile.py",
    about = "Calculate the distribution of deleted nucleotides across aligned reads."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input")]
    input_file: PathBuf,

    /// Expected aligned read length; reads whose sequence length or
    /// M/S/I-op total doesn't exactly match are skipped.
    #[arg(short = 'l', long = "read-align-length")]
    read_align_length: usize,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Maximum number of qualifying (has-deletion, right-length) reads to process.
    #[arg(short = 'n', long = "read-num", default_value_t = 1_000_000)]
    read_num: u64,

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
    let profile = compute_deletion_profile(
        reader.records(),
        args.mapq,
        args.read_align_length,
        args.read_num,
    )?;
    eprintln!("Total reads used: {}", profile.count);

    let prefix = args.out_prefix.to_string_lossy();

    let mut table = File::create(format!("{prefix}.deletion_profile.txt"))?;
    table.write_all(render_deletion_table(&profile).as_bytes())?;

    let mut r = File::create(format!("{prefix}.deletion_profile.r"))?;
    r.write_all(render_deletion_r_script(&profile, &prefix).as_bytes())?;

    Ok(())
}
