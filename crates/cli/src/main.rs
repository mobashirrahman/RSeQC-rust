//! Dispatch and flag parsing. Binary is named after the original `.py`
//! script (see Cargo.toml [[bin]] name) so PATH invocation matches upstream
//! once the packaging step (PORTING_PLAN Step 10) adds the `.py` alias.

// B4: link `rseqc_cli` so its `#[global_allocator]` (actionable
// out-of-memory message) applies to this binary too.
use rseqc_cli as _;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::bam_stat::run_bam_stat;

#[derive(Parser)]
#[command(
    name = "bam_stat.py",
    about = "Summarize mapping statistics for a BAM or SAM alignment file."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("bam_stat.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // C1: the CLI body lives in `rseqc_commands::bam_stat::run_bam_stat`
    // (same bytes, same order) so `rseqc_multi` can drive it over a record
    // broadcast with per-command stream files.
    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    run_bam_stat(
        records,
        args.mapq,
        &mut std::io::stdout(),
        &mut std::io::stderr(),
    )
}
