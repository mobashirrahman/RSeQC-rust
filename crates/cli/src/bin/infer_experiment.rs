//! Dispatch and flag parsing for `infer_experiment.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

// B4: link `rseqc_cli` so its `#[global_allocator]` (actionable
// out-of-memory message) applies to this binary too.
use rseqc_cli as _;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::infer_experiment::run_infer_experiment;

#[derive(Parser)]
#[command(
    name = "infer_experiment.py",
    about = "Infer RNA-seq library layout and strandedness from a SAM/BAM file."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference gene model (BED format).
    #[arg(short = 'r', long = "refgene")]
    refgene: PathBuf,

    /// Number of usable alignments to sample.
    #[arg(short = 's', long = "sample-size", default_value_t = 200_000)]
    sample_size: u64,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("infer_experiment.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // C1 multi-driver pattern: the CLI body lives in
    // `rseqc_commands::infer_experiment::run_infer_experiment` (same bytes,
    // same streams, same order) so `rseqc_multi` can drive it over a record
    // broadcast plus per-command stream files. There is no output-prefix
    // parent check to keep here: this command writes no output file, so
    // there is no parent directory to check -- upstream's shape, not a gap.
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    run_infer_experiment(
        records,
        &header,
        &args.refgene,
        args.sample_size,
        args.mapq,
        &mut std::io::stdout(),
        &mut std::io::stderr(),
    )
}
