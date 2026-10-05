//! Dispatch and flag parsing for `read_quality.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! Output: .qual.r script for generating quality plots (no separate .xls
//! table -- upstream's `readsQual_boxplot` genuinely only writes the R
//! script, confirmed by reading the source: there is no other `open()`
//! call in that function).

use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_quality::run_read_quality;

#[derive(Parser)]
#[command(
    name = "read_quality.py",
    about = "Calculate per-cycle Phred quality-score distributions for aligned reads."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Ignore quality-score observations occurring fewer than this many times.
    #[arg(short = 'r', long = "reduce", default_value_t = 1)]
    reduce: u64,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate quality-profile data but do not execute the R script.
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
            eprintln!("read_quality.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // C1 multi-driver pattern: the CLI body lives in
    // `rseqc_commands::read_quality::run_read_quality` (same bytes, same
    // files, same order) so `rseqc_multi` can drive it over a record
    // broadcast with per-command stream files. The `_or_exit` parent check
    // stays here so standalone exit 2 is unchanged; the shared body
    // re-checks in `Result` form for workers.
    //
    // Upstream's validate_args refuses an output prefix whose parent directory does
    // not exist, before any input is read. Omitting it here meant the whole
    // alignment was read and every metric computed, then discarded when the output
    // open failed with "No such file or directory (os error 2)" -- an error naming
    // neither the directory nor the flag, and indistinguishable from a missing
    // input. The shared helper keeps that check in one place so it cannot be
    // forgotten by the next binary.
    rseqc_cli::require_existing_output_parent_or_exit("read_quality.py", &args.out_prefix);

    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy().into_owned();
    run_read_quality(
        records,
        args.mapq,
        &prefix,
        args.reduce,
        args.skip_plot,
        &args.rscript,
        &mut std::io::stdout(),
        &mut std::io::stderr(),
    )
}
