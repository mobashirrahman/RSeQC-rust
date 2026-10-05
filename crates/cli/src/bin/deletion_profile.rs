//! Dispatch and flag parsing for `deletion_profile.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::deletion_profile::run_deletion_profile;

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
            eprintln!("deletion_profile.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // C1 multi-driver pattern: the CLI body lives in
    // `rseqc_commands::deletion_profile::run_deletion_profile` (same bytes,
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
    rseqc_cli::require_existing_output_parent_or_exit("deletion_profile.py", &args.out_prefix);

    // BAM only: this command's own `-i` documents SAM-text as unsupported,
    // so it keeps the direct `open_bam` rather than the extension-dispatching
    // `open_alignments` the other profile commands use.
    let (mut reader, _header) = rseqc_formats::open_bam(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy().into_owned();
    run_deletion_profile(
        reader.records(),
        args.mapq,
        &prefix,
        args.read_align_length,
        args.read_num,
        args.skip_plot,
        &args.rscript,
        &mut std::io::stdout(),
        &mut std::io::stderr(),
    )
}
