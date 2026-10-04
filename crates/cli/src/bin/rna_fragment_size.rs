//! Dispatch and flag parsing for `RNA_fragment_size.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

// B4: link `rseqc_cli` so its `#[global_allocator]` (actionable
// out-of-memory message) applies to this binary too.
use rseqc_cli as _;
use std::fs::File;
use std::io::{BufRead, BufReader, Write as _};
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::rna_fragment_size::{compute_fragment_sizes, format_result, parse_bed12_line, IndexedReads, HEADER};

#[derive(Parser)]
#[command(
    name = "RNA_fragment_size.py",
    about = "Calculate fragment-size statistics for each transcript or gene."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input")]
    input_file: PathBuf,

    /// Reference gene model (BED12 format).
    #[arg(short = 'r', long = "refgene")]
    refgene: PathBuf,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Minimum number of fragments required to report real statistics.
    #[arg(short = 'n', long = "frag-num", default_value_t = 3)]
    frag_num: usize,

    /// Output file (defaults to stdout).
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("RNA_fragment_size.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let indexed = IndexedReads::build(reader.records(), &header)?;

    let mut out: Box<dyn std::io::Write> = match &args.output {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(std::io::stdout()),
    };
    writeln!(out, "{HEADER}")?;

    let bed_file = BufReader::new(File::open(&args.refgene)?);
    for line in bed_file.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let bed = parse_bed12_line(&line)?;
        let stats = compute_fragment_sizes(&bed, &indexed, args.mapq, args.frag_num);
        writeln!(out, "{}", format_result(&stats))?;
    }

    if let Some(path) = &args.output {
        eprintln!("Created: {}", path.display());
    }

    Ok(())
}
