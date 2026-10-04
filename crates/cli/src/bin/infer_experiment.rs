//! Dispatch and flag parsing for `infer_experiment.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

// B4: link `rseqc_cli` so its `#[global_allocator]` (actionable
// out-of-memory message) applies to this binary too.
use rseqc_cli as _;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::infer_experiment::{compute_experiment, render_results, GeneRanges};

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
    // Upstream's `validate_args` prints this warning (if any) before doing
    // any real work -- ahead of even opening the refgene BED.
    if args.sample_size < 1_000 {
        eprintln!("Warning: sample size is below 1,000; the inferred protocol may be unreliable.");
    }

    // `"Reading reference gene model " + refbed + ' ...'` then `end=' '`:
    // one space from the literal's own trailing `...`+space concatenation,
    // no second space (unlike read_quality's "Read BAM file ...  Done").
    eprint!("Reading reference gene model {} ... ", args.refgene.display());
    let (gene_ranges, skipped) = GeneRanges::parse(BufReader::new(File::open(&args.refgene)?))?;
    if skipped > 0 {
        eprintln!("[NOTE: input bed must be 12-column] skipped {skipped} line(s)");
    }
    eprintln!("Done");

    // `"Loading SAM/BAM file ... "` (trailing space in the literal) plus
    // `end=' '` gives two spaces before whatever prints next.
    eprint!("Loading SAM/BAM file ...  ");
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let result = compute_experiment(records, &header, &gene_ranges, args.sample_size, args.mapq)?;
    if result.stopped_at_eof {
        eprintln!("Finished");
    }
    eprintln!("Total {} usable reads were sampled", result.sampled_count);

    println!("{}", render_results(&result));
    Ok(())
}
