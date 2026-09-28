//! Dispatch and flag parsing for `overlay_bigwig.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging
//! (PORTING_PLAN Step 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::overlay_bigwig::{overlay_bigwigs, parse_action};
use rseqc_formats::bigwig::BigWigReader;

#[derive(Parser)]
#[command(name = "overlay_bigwig.py", about = "Apply an arithmetic operation to two BigWig signal tracks.")]
struct Args {
    /// First BigWig file.
    #[arg(short = 'i', long = "bwfile1")]
    bigwig_file1: PathBuf,

    /// Second BigWig file. Both files should use the same reference genome.
    #[arg(short = 'j', long = "bwfile2")]
    bigwig_file2: PathBuf,

    /// Arithmetic operation applied to corresponding signal values.
    #[arg(short = 'a', long = "action")]
    action: String,

    /// Output variableStep WIG file.
    #[arg(short = 'o', long = "output")]
    output_wig: PathBuf,

    /// Chromosome chunk size in bp.
    #[arg(short = 'c', long = "chunk", default_value_t = 100_000)]
    chunk_size: i64,

    /// Allow an existing output file to be replaced.
    #[arg(long = "overwrite")]
    overwrite: bool,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("overlay_bigwig.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let action = parse_action(&args.action).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("--action must be one of Add, Average, Division, Max, Min, Product, Subtract, geometricMean (got {:?})", args.action),
        )
    })?;
    if args.chunk_size <= 0 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--chunk must be greater than zero"));
    }
    if args.output_wig.exists() && !args.overwrite {
        return Err(std::io::Error::new(std::io::ErrorKind::AlreadyExists, format!("output file already exists: {}; use --overwrite to replace it", args.output_wig.display())));
    }

    let mut bw1 = BigWigReader::open(&args.bigwig_file1)?;
    let mut bw2 = BigWigReader::open(&args.bigwig_file2)?;

    eprintln!("Get chromosome sizes from BigWig headers ...");
    let body = overlay_bigwigs(&mut bw1, &mut bw2, action, args.chunk_size)?;

    File::create(&args.output_wig)?.write_all(body.as_bytes())?;
    eprintln!("Created: {}", args.output_wig.display());

    Ok(())
}
