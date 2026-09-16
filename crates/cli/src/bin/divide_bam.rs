//! Dispatch and flag parsing for `divide_bam.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! `--index` (BAM indexing) is not implemented yet — disclosed gap, see
//! crates/commands/src/divide_bam.rs module docs.

use std::fs::File;
use std::io;
use std::path::PathBuf;

use clap::Parser;
use noodles_bam as bam;
use rseqc_commands::divide_bam::{DivideCounts, divide_bam};
use rand::SeedableRng;
use rand::rngs::StdRng;

#[derive(Parser)]
#[command(
    name = "divide_bam.py",
    about = "Randomly divide a BAM file into approximately equal subsets."
)]
struct Args {
    /// Input BAM file.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Number of output BAM files to create.
    #[arg(short = 'n', long = "subset-num")]
    subset_num: usize,

    /// Prefix for the output BAM files (`<prefix>_<index>.bam`).
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Skip unmapped reads.
    #[arg(short = 's', long = "skip-unmap")]
    skip_unmap: bool,

    /// Random seed for reproducible division (optional).
    #[arg(long = "seed")]
    seed: Option<u64>,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(counts) => {
            eprintln!("{counts:?}");
            std::process::ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<DivideCounts> {
    let (reader, header) = rseqc_formats::open_bam(&args.input_file)?;

    if args.subset_num == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "subset-num must be greater than 0",
        ));
    }

    let prefix = args.out_prefix.to_string_lossy();
    let mut outputs = Vec::new();

    for i in 0..args.subset_num {
        let path = format!("{}_{}.bam", prefix, i);
        let file = File::create(&path)?;
        let writer = bam::io::Writer::new(file);
        outputs.push(writer);
    }

    // Write headers to all outputs
    for output in &mut outputs {
        output.write_header(&header)?;
    }

    let mut reader = reader;

    // Initialize RNG
    let mut rng = if let Some(seed) = args.seed {
        StdRng::seed_from_u64(seed)
    } else {
        StdRng::from_entropy()
    };

    let counts = divide_bam(
        reader.records(),
        &header,
        &mut outputs,
        args.skip_unmap,
        &mut rng,
    )?;

    // Close all writers to ensure data is flushed
    drop(outputs);
    drop(reader);

    Ok(counts)
}