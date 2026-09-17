//! Dispatch and flag parsing for `read_duplication.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! All three outputs (.seq.DupRate.xls, .pos.DupRate.xls, and .DupRate_plot.r) are generated as upstream does.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::read_duplication::{compute_duplication, render_dup_r_script, render_pos_dup_table, render_seq_dup_table};

#[derive(Parser)]
#[command(
    name = "read_duplication.py",
    about = "Calculate sequence-based and mapping-based read duplication rates."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Upper bound of occurrence for counting.
    #[arg(short = 'u', long = "up-limit", default_value_t = 500)]
    up_limit: u64,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate duplication data but do not execute the R plotting script.
    #[arg(long = "skip-plot")]
    skip_plot: bool,

    /// Rscript executable to use.
    #[arg(long = "rscript", default_value = "Rscript")]
    rscript: String,
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
    // Upstream: `print("Load BAM file ... ", end=' ')` -- the literal's
    // own trailing space plus `end=' '` gives two spaces before "Done".
    eprint!("Load BAM file ...  ");
    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let hist = compute_duplication(reader.records(), &header, args.mapq)?;
    eprintln!("Done");

    let prefix = args.out_prefix.to_string_lossy().into_owned();

    eprintln!("report duplicte rate based on sequence ...");
    let seq_path = format!("{prefix}.seq.DupRate.xls");
    File::create(&seq_path)?.write_all(render_seq_dup_table(&hist).as_bytes())?;

    eprintln!("report duplicte rate based on mapping  ...");
    let pos_path = format!("{prefix}.pos.DupRate.xls");
    File::create(&pos_path)?.write_all(render_pos_dup_table(&hist).as_bytes())?;

    eprintln!("generate R script ...");
    let r_path = format!("{prefix}.DupRate_plot.r");
    File::create(&r_path)?.write_all(render_dup_r_script(&hist, &prefix, args.up_limit).as_bytes())?;

    if !args.skip_plot {
        let status = Command::new(&args.rscript).arg(&r_path).status()?;
        if !status.success() {
            return Err(std::io::Error::other(format!("R plotting failed for {r_path}")));
        }
    }

    Ok(())
}
