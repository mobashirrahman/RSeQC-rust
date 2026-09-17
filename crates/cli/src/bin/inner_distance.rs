//! Dispatch and flag parsing for `inner_distance.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.
//!
//! `--skip-plot`/`--rscript` (running Rscript) are not implemented --
//! disclosed gap, this only generates the .r script text.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::inner_distance::{
    build_model, compute_inner_distance, histogram_buckets, render_distance_file, render_freq_table, render_r_script,
};

#[derive(Parser)]
#[command(
    name = "inner_distance.py",
    about = "Estimate the inner distance between paired-end RNA-seq reads."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Reference gene model (BED12 format).
    #[arg(short = 'r', long = "refgene")]
    refgene: PathBuf,

    /// Maximum number of qualifying read pairs to process.
    #[arg(short = 'k', long = "sample-size", default_value_t = 1_000_000)]
    sample_size: u64,

    /// Lower bound of the inner-distance histogram.
    #[arg(short = 'l', long = "lower-bound", default_value_t = -250)]
    lower_bound: i64,

    /// Upper bound of the inner-distance histogram.
    #[arg(short = 'u', long = "upper-bound", default_value_t = 250)]
    upper_bound: i64,

    /// Histogram bucket width.
    #[arg(short = 's', long = "step", default_value_t = 5)]
    step: i64,

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
    // Upstream: `print("Get exon regions from " + refbed + " ...",
    // file=sys.stderr)` -- its own line, default trailing newline.
    eprintln!("Get exon regions from {} ...", args.refgene.display());
    let model = build_model(&args.refgene)?;

    // Upstream: `print("Load BAM file ... ", end=' ', file=sys.stderr)`,
    // then (only on natural iterator exhaustion, not on hitting
    // `sample_size`) `print("Done", file=sys.stderr)`.
    eprint!("Load BAM file ...  ");
    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let result = compute_inner_distance(reader.records(), &header, &model, args.mapq, args.sample_size)?;
    if result.loop_exhausted_naturally {
        eprintln!("Done");
    }

    // Upstream: `print("Total read pairs  used " + str(pair_num),
    // file=sys.stderr)` -- literal double space before "used".
    eprintln!("Total read pairs  used {}", result.pair_num);
    if result.pair_num == 0 {
        eprintln!("Cannot find paired reads");
        return Ok(());
    }

    let prefix = args.out_prefix.to_string_lossy();

    let mut distance_file = File::create(format!("{prefix}.inner_distance.txt"))?;
    distance_file.write_all(render_distance_file(&result).as_bytes())?;

    let values: Vec<i64> = result.records.iter().filter_map(|r| r.histogram_value).collect();
    let buckets = histogram_buckets(&values, args.lower_bound, args.upper_bound, args.step);

    let mut freq_file = File::create(format!("{prefix}.inner_distance_freq.txt"))?;
    freq_file.write_all(render_freq_table(&buckets).as_bytes())?;

    let mut r_file = File::create(format!("{prefix}.inner_distance_plot.r"))?;
    r_file.write_all(render_r_script(&buckets, args.step, &prefix).as_bytes())?;

    Ok(())
}
