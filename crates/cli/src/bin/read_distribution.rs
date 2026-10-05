//! Dispatch and flag parsing for `read_distribution.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

// B4: link `rseqc_cli` so its `#[global_allocator]` (actionable
// out-of-memory message) applies to this binary too.
use rseqc_cli as _;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_distribution::run_read_distribution;

#[derive(Parser)]
#[command(
    name = "read_distribution.py",
    about = "Summarize read distribution across genomic annotation categories."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference gene model (BED12 format).
    #[arg(short = 'r', long = "refgene")]
    refgene: PathBuf,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("read_distribution.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // C1 multi-driver pattern: the CLI body lives in
    // `rseqc_commands::read_distribution::run_read_distribution` (same
    // bytes, same streams, same order) so `rseqc_multi` can drive it over a
    // record broadcast plus per-command stream files. There is no
    // output-prefix parent check to keep here: this command's report IS its
    // stdout (upstream takes no `-o`), so there is no output file.
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    run_read_distribution(
        records,
        &header,
        &args.refgene,
        &args.input_file,
        &mut std::io::stdout(),
        &mut std::io::stderr(),
    )
}
