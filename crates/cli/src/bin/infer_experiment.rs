//! Dispatch and flag parsing for `infer_experiment.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

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
    /// Input BAM file (SAM-text input is not yet supported).
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
            eprintln!("error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let (gene_ranges, skipped) = GeneRanges::parse(BufReader::new(File::open(&args.refgene)?))?;
    if skipped > 0 {
        eprintln!("[NOTE: input bed must be 12-column] skipped {skipped} line(s)");
    }

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let result = compute_experiment(reader.records(), &header, &gene_ranges, args.sample_size, args.mapq)?;

    println!("{}", render_results(&result));
    Ok(())
}
