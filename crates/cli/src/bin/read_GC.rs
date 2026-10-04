//! Dispatch and flag parsing for `read_GC.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.

use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_gc::run_read_gc;

#[derive(Parser)]
#[command(
    name = "read_GC.py",
    about = "Calculate the GC-content distribution of aligned reads."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate GC-profile data but do not execute the R plotting script.
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
            eprintln!("read_GC.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // C1: the CLI body lives in `rseqc_commands::read_gc::run_read_gc`
    // (same bytes, same files, same order) so `rseqc_multi` can drive it
    // over a record broadcast with per-command stream files. The
    // `_or_exit` parent check stays here so standalone exit 2 is
    // unchanged; the shared body re-checks in `Result` form for workers.
    rseqc_cli::require_existing_output_parent_or_exit("read_GC.py", &args.out_prefix);

    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy().into_owned();
    run_read_gc(
        records,
        args.mapq,
        &prefix,
        args.skip_plot,
        &args.rscript,
        &mut std::io::stdout(),
        &mut std::io::stderr(),
    )
}
