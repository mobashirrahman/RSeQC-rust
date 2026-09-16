//! Dispatch and flag parsing for `read_distribution.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_distribution::{count_read_distribution, process_gene_model, render_report};

#[derive(Parser)]
#[command(
    name = "read_distribution.py",
    about = "Summarize read distribution across genomic annotation categories."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
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
            eprintln!("error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let model = process_gene_model(&args.refgene)?;
    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let counts = count_read_distribution(reader.records(), &header, &model)?;
    println!("{}", render_report(&model, &counts));
    Ok(())
}
