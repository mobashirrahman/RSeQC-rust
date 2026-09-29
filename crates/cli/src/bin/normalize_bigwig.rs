//! Dispatch and flag parsing for `normalize_bigwig.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging
//! (PORTING_PLAN Step 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::normalize_bigwig::{calculate_wigsum, render_normalized_body};
use rseqc_commands::python_fmt::python_g12;
use rseqc_formats::bigwig::BigWigReader;

#[derive(Parser)]
#[command(name = "normalize_bigwig.py", about = "Normalize a BigWig signal to a fixed total WIG sum.")]
struct Args {
    /// Input BigWig file.
    #[arg(short = 'i', long = "bwfile")]
    bigwig_file: PathBuf,

    /// Output WIG or bedGraph file.
    #[arg(short = 'o', long = "output")]
    output_file: PathBuf,

    /// Target total WIG sum.
    #[arg(short = 't', long = "wigsum", default_value_t = 100_000_000.0)]
    total_wigsum: f64,

    /// Optional BED gene model; when supplied, the normalization factor is calculated from merged exon regions only.
    #[arg(short = 'r', long = "refgene")]
    refgene_bed: Option<PathBuf>,

    /// Chromosome chunk size in bp.
    #[arg(short = 'c', long = "chunk", default_value_t = 500_000)]
    chunk_size: i64,

    /// Output format: 'wig' for variableStep WIG or 'bgr' for bedGraph.
    #[arg(short = 'f', long = "format", default_value = "bgr")]
    out_format: String,

    /// Allow an existing output file to be replaced.
    #[arg(long = "overwrite")]
    overwrite: bool,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("normalize_bigwig.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let out_format = args.out_format.to_lowercase();
    if !matches!(out_format.as_str(), "wig" | "bgr") {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--format must be 'wig' or 'bgr'"));
    }
    if args.total_wigsum <= 0.0 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--wigsum must be greater than zero"));
    }
    if args.chunk_size <= 0 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--chunk must be greater than zero"));
    }
    if args.output_file.exists() && !args.overwrite {
        return Err(std::io::Error::new(std::io::ErrorKind::AlreadyExists, format!("output file already exists: {}; use --overwrite to replace it", args.output_file.display())));
    }

    let mut bw = BigWigReader::open(&args.bigwig_file)?;

    eprintln!("Get chromosome sizes from BigWig header ...");
    let refgene_reader = match &args.refgene_bed {
        Some(p) => Some(BufReader::new(File::open(p)?)),
        None => None,
    };
    let refgene_path = args.refgene_bed.as_ref().map(|p| p.display().to_string()).unwrap_or_default();
    if refgene_reader.is_some() {
        eprintln!("Extract exons from {refgene_path}");
    } else {
        eprintln!("Calculate WIG sum from {}", args.bigwig_file.display());
    }

    // `calculate_wigsum` prints "\nTotal WIG sum is ...\n" itself, before
    // its own zero/negative-sum check -- see its doc comment. That must
    // stay inside `calculate_wigsum` so it's emitted even on `Err`, not
    // duplicated here on `Ok`.
    let wigsum = calculate_wigsum(&mut bw, refgene_reader, args.total_wigsum, args.chunk_size, &refgene_path)?;

    eprintln!("Normalization factor: {}", python_g12(wigsum.weight));
    eprintln!("Normalizing BigWig file ...");

    let body = render_normalized_body(&mut bw, &wigsum.chrom_sizes, args.chunk_size, wigsum.weight, &out_format)?;
    File::create(&args.output_file)?.write_all(body.as_bytes())?;

    eprintln!("Created: {}", args.output_file.display());
    eprintln!("Observed WIG sum: {:.2}", wigsum.observed_wigsum);
    eprintln!("Applied normalization factor: {}", python_g12(wigsum.weight));

    Ok(())
}
