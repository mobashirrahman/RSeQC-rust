//! Dispatch for `rseqc_multi`: run several commands over one BAM read.
//!
//! Contract (card C1): `rseqc_multi -i x.bam [-r x.bed] -o <prefix>
//! --run bam_stat,read_GC` writes exactly the files each command would
//! write on its own under `-o <prefix>`, plus each command's stdout and
//! stderr to `<prefix>.<command>.stdout` / `<prefix>.<command>.stderr`.
//! See `rseqc_commands::multi` for the transport and the C2 pattern.
//!
//! Exit status: 0 when every selected command succeeds; 1 when any fails
//! (naming each failed command's error -- the same failure convention each
//! command's own binary uses); 2 for a usage error such as an unknown
//! `--run` token.

use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::multi::{MultiArgs, drive, resolve};

#[derive(Parser)]
#[command(
    name = "rseqc_multi",
    about = "Run several RSeQC commands over a single BAM read."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference BED file. Accepted and currently ignored: no pilot
    /// command reads a gene model (reserved for later commands).
    #[arg(short = 'r', long = "reference-bed")]
    reference_bed: Option<PathBuf>,

    /// Shared output prefix: data files keep the standalone names under
    /// it, and each command's streams go to `<prefix>.<command>.stdout`
    /// / `<prefix>.<command>.stderr`.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Comma-separated command tokens to run (default: every registered
    /// command). Tokens are upstream command stems (`bam_stat,read_GC`).
    #[arg(long = "run")]
    run: Option<String>,

    /// Minimum mapping quality, forwarded to every command.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Skip R plotting, forwarded to plotting commands.
    #[arg(long = "skip-plot")]
    skip_plot: bool,

    /// Rscript executable, forwarded to plotting commands.
    #[arg(long = "rscript", default_value = "Rscript")]
    rscript: String,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("rseqc_multi: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let _ = args.reference_bed.as_ref();

    let tokens: Vec<String> = match &args.run {
        Some(list) => list
            .split(',')
            .map(|token| token.trim().to_string())
            .filter(|token| !token.is_empty())
            .collect(),
        None => rseqc_commands::multi::COMMANDS
            .iter()
            .map(|entry| entry.name.to_string())
            .collect(),
    };
    if tokens.is_empty() {
        rseqc_cli::usage_exit("rseqc_multi", "--run selects no commands");
    }
    let entries = match resolve(&tokens) {
        Ok(entries) => entries,
        Err(message) => rseqc_cli::usage_exit("rseqc_multi", &message),
    };

    let prefix = args.out_prefix.to_string_lossy().into_owned();
    let multi_args = MultiArgs {
        mapq: args.mapq,
        out_prefix: prefix,
        skip_plot: args.skip_plot,
        rscript: args.rscript.clone(),
    };

    // The calling thread is the reader: open here so an unreadable input
    // fails before any stream file is created, exactly like a standalone
    // binary failing its open before any output.
    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let results = drive(entries, &multi_args, records);

    let mut failed = 0;
    let mut details = Vec::new();
    for worker in &results {
        if let Err(err) = &worker.result {
            failed += 1;
            details.push(format!("{}: {err}", worker.entry.name));
        }
    }
    if failed > 0 {
        // Name the per-command causes on multi's own stderr: when even the
        // stream files cannot be created (e.g. missing output directory),
        // there is nowhere else the reason could have gone.
        return Err(std::io::Error::other(format!(
            "{} of {} commands failed ({})",
            failed,
            results.len(),
            details.join("; ")
        )));
    }
    Ok(())
}
