//! Dispatch and flag parsing for `clipping_profile.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use rseqc_commands::clipping_profile::run_clipping_profile;

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
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
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

    /// Generate the profile files but do not run the R plotting script.
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
            eprintln!("clipping_profile.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // C1 multi-driver pattern: the CLI body lives in
    // `rseqc_commands::clipping_profile::run_clipping_profile` (same bytes,
    // same files, same order) so `rseqc_multi` can drive it over a record
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
    rseqc_cli::require_existing_output_parent_or_exit("clipping_profile.py", &args.out_prefix);

    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy().into_owned();
    let sequencing = match args.sequencing {
        Layout::SingleEnd => "SE",
        Layout::PairedEnd => "PE",
    };
    run_clipping_profile(
        records,
        args.mapq,
        &prefix,
        sequencing,
        args.skip_plot,
        &args.rscript,
        &mut std::io::stdout(),
        &mut std::io::stderr(),
    )
}
