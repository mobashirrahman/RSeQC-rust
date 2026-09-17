//! Dispatch and flag parsing for `inner_distance.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.
//!
use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

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
    /// Input alignment file in BAM or plain-text SAM format (dispatched by
    /// the `.bam`/`.sam` extension).
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

    /// Generate inner-distance data but do not execute the R plotting script.
    #[arg(long = "skip-plot")]
    skip_plot: bool,

    /// Rscript executable to use.
    #[arg(long = "rscript", default_value = "Rscript")]
    rscript: String,
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

    // Upstream: `if self.bam_format: print("Load BAM file ... ", end=' ')
    // else: print("Load SAM file ... ", end=' ')`, then (only on natural
    // iterator exhaustion, not on hitting `sample_size`) `print("Done",
    // file=sys.stderr)`. `self.bam_format` comes from `pysam.Samfile(path,
    // 'rb')` succeeding, which it does even for genuine plain-text SAM
    // content (htslib auto-detects, ignoring the 'b' mode hint; confirmed
    // via a live diff for bam_stat.py and others, same underlying pysam
    // call here) -- the "Load SAM file" branch is practically dead code.
    eprint!("Load BAM file ...  ");
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let result = compute_inner_distance(records, &header, &model, args.mapq, args.sample_size)?;
    if result.loop_exhausted_naturally {
        eprintln!("Done");
    }

    let prefix = args.out_prefix.to_string_lossy();

    // Upstream: `mRNA_inner_distance` opens all three output files
    // UNCONDITIONALLY at the very top (`FO=open(...)`, `FQ=open(...)`,
    // `RS=open(...)`), before any read scanning happens -- so even the
    // `pair_num == 0` early exit below leaves all three as empty (but
    // existing) files. That early exit is a literal `sys.exit(0)`
    // called from INSIDE this function (not a caught exception), which
    // unwinds the whole interpreter immediately -- `main()`'s
    // `--skip-plot`/Rscript-invocation code never even runs on this
    // path, matched here by returning before reaching it.
    let distance_path = format!("{prefix}.inner_distance.txt");
    let freq_path = format!("{prefix}.inner_distance_freq.txt");
    let r_path = format!("{prefix}.inner_distance_plot.r");
    File::create(&distance_path)?;
    File::create(&freq_path)?;
    File::create(&r_path)?;

    // Upstream: `print("Total read pairs  used " + str(pair_num),
    // file=sys.stderr)` -- literal double space before "used".
    eprintln!("Total read pairs  used {}", result.pair_num);
    if result.pair_num == 0 {
        eprintln!("Cannot find paired reads");
        return Ok(());
    }

    let mut distance_file = File::create(&distance_path)?;
    distance_file.write_all(render_distance_file(&result).as_bytes())?;

    let values: Vec<i64> = result.records.iter().filter_map(|r| r.histogram_value).collect();
    let buckets = histogram_buckets(&values, args.lower_bound, args.upper_bound, args.step);

    let mut freq_file = File::create(&freq_path)?;
    freq_file.write_all(render_freq_table(&buckets).as_bytes())?;

    File::create(&r_path)?.write_all(render_r_script(&buckets, args.step, &prefix).as_bytes())?;

    if !args.skip_plot {
        let status = Command::new(&args.rscript).arg(&r_path).status()?;
        if !status.success() {
            return Err(std::io::Error::other(format!("R plotting failed for {r_path}")));
        }
    }

    Ok(())
}
