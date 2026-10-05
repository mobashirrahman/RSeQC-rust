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

    /// Include ambiguous nucleotides N and X in the NVC plot (`-x`),
    /// forwarded to `read_NVC` only.
    #[arg(short = 'x', long = "nx")]
    nx: bool,

    /// Ignore quality-score observations occurring fewer than this many
    /// times, forwarded to `read_quality` only. Long-only: `-r` is
    /// already the reference BED on this driver.
    #[arg(long = "reduce", default_value_t = 1)]
    reduce: u64,

    /// Sequencing layout for `clipping_profile` and `insertion_profile`
    /// (`SE` or `PE`, the only values their standalone `--sequencing`
    /// accepts). Upstream requires the flag; the driver defaults to `SE`
    /// so unattended full runs work.
    #[arg(short = 's', long = "sequencing", default_value = "SE")]
    sequencing: String,

    /// Expected aligned read length for `deletion_profile` and
    /// `mismatch_profile` (`-l` in each), which filter reads on it.
    /// Required when either is selected: there is no sensible driver
    /// default, because a wrong length silently changes which reads
    /// qualify.
    #[arg(long = "read-align-length")]
    read_align_length: Option<usize>,

    /// Maximum qualifying reads for `deletion_profile` and
    /// `mismatch_profile` (`-n` in each). Optional; the command's own
    /// default applies when absent.
    #[arg(long = "read-num")]
    read_num: Option<u64>,
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
    // `clipping_profile` and `insertion_profile` are the only registered
    // commands whose standalone flag is restricted (clap `SE`/`PE`);
    // validate here so a bad value stays a usage error (exit 2) instead
    // of a worker failure.
    if args.sequencing != "SE"
        && args.sequencing != "PE"
        && tokens
            .iter()
            .any(|token| token == "clipping_profile" || token == "insertion_profile")
    {
        rseqc_cli::usage_exit(
            "rseqc_multi",
            "invalid --sequencing (expected SE or PE)",
        );
    }
    // `deletion_profile` and `mismatch_profile` filter on the aligned read
    // length, so there is no honest default to supply on their behalf: a
    // wrong value would silently change which reads qualify and the
    // command would still succeed.
    if tokens
        .iter()
        .any(|token| token == "deletion_profile" || token == "mismatch_profile")
        && args.read_align_length.is_none()
    {
        rseqc_cli::usage_exit(
            "rseqc_multi",
            "--read-align-length is required when deletion_profile or \
             mismatch_profile is selected (they filter reads on it)",
        );
    }

    let prefix = args.out_prefix.to_string_lossy().into_owned();
    let multi_args = MultiArgs {
        mapq: args.mapq,
        out_prefix: prefix,
        skip_plot: args.skip_plot,
        rscript: args.rscript.clone(),
        nx: args.nx,
        reduce: args.reduce,
        sequencing: args.sequencing.clone(),
        read_align_length: args.read_align_length,
        read_num: args.read_num,
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
