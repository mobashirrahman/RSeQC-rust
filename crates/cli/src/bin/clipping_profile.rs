//! Dispatch and flag parsing for `clipping_profile.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use rseqc_commands::clipping_profile::{
    compute_paired_end, compute_single_end, render_paired_r_script, render_paired_table,
    render_single_r_script, render_single_table,
};

#[derive(Clone, Copy, ValueEnum)]
enum Layout {
    #[value(name = "SE")]
    SingleEnd,
    #[value(name = "PE")]
    PairedEnd,
}

#[derive(Parser)]
#[command(
    name = "clipping_profile.py",
    about = "Estimate the clipping profile of RNA-seq reads from a BAM or SAM file."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Sequencing layout: SE for single-end or PE for paired-end.
    #[arg(short = 's', long = "sequencing")]
    sequencing: Layout,

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
    let prefix = args.out_prefix.to_string_lossy();

    let table_text;
    let r_script_text;

    match args.sequencing {
        Layout::SingleEnd => {
            let profile = compute_single_end(reader.records(), args.mapq, b'S')?;
            eprintln!("Total reads used: {}", profile.total_read);
            table_text = render_single_table(&profile);
            r_script_text = render_single_r_script(&profile, &prefix);
        }
        Layout::PairedEnd => {
            let profile = compute_paired_end(reader.records(), args.mapq, b'S')?;
            eprintln!(
                "Total read-1 used: {}, read-2 used: {}",
                profile.total_read1, profile.total_read2
            );
            table_text = render_paired_table(&profile);
            r_script_text = render_paired_r_script(&profile, &prefix);
        }
    }

    let mut xls = File::create(format!("{prefix}.clipping_profile.xls"))?;
    xls.write_all(table_text.as_bytes())?;

    let mut r = File::create(format!("{prefix}.clipping_profile.r"))?;
    r.write_all(r_script_text.as_bytes())?;

    Ok(())
}
