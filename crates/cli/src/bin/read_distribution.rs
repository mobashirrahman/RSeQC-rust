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
    // `print(f"Processing {gene_model} ...", end=" ")` then a separate
    // `print("Done")` on the same line (one space from `end=" "`, one
    // trailing newline from the "Done" print).
    eprint!("Processing {} ... ", args.refgene.display());
    let model = process_gene_model(&args.refgene)?;
    eprintln!("Done");

    // `print(f"Processing {input_file} ...", end=" ")` then a separate
    // `print("Finished\n")`: the literal "\n" plus the print's own
    // newline give a BLANK line after "Finished", not just one newline.
    eprint!("Processing {} ... ", args.input_file.display());
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let counts = count_read_distribution(records, &header, &model)?;
    eprintln!("Finished\n");

    println!("{}", render_report(&model, &counts));
    Ok(())
}
