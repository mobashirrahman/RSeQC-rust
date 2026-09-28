//! Dispatch and flag parsing for `read_quality.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! Output: .qual.r script for generating quality plots (no separate .xls
//! table -- upstream's `readsQual_boxplot` genuinely only writes the R
//! script, confirmed by reading the source: there is no other `open()`
//! call in that function).

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::read_quality::{compute_quality, render_qual_r_script};

#[derive(Parser)]
#[command(
    name = "read_quality.py",
    about = "Calculate per-cycle Phred quality-score distributions for aligned reads."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Ignore quality-score observations occurring fewer than this many times.
    #[arg(short = 'r', long = "reduce", default_value_t = 1)]
    reduce: u64,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate quality-profile data but do not execute the R script.
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
    // Upstream: `if self.bam_format: print("Read BAM file ... ", end=' ')
    // else: print("Read SAM file ... ", end=' ')` -- `self.bam_format`
    // comes from `pysam.Samfile(path, 'rb')` succeeding, which it does
    // even for genuine plain-text SAM content (htslib auto-detects,
    // ignoring the 'b' mode hint; confirmed via a live diff for
    // bam_stat.py/read_NVC.py/read_GC.py, same underlying pysam.Samfile
    // call here). The "Read SAM file" branch is practically dead code
    // for any valid input. The literal's own trailing space plus
    // `end=' '` gives two spaces before "Done".
    eprint!("Read BAM file ...  ");
    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let hist = compute_quality(records, args.mapq)?;
    eprintln!("Done");

    let prefix = args.out_prefix.to_string_lossy().into_owned();
    let r_path = format!("{prefix}.qual.r");
    File::create(&r_path)?.write_all(render_qual_r_script(&hist, args.reduce, &prefix).as_bytes())?;

    if !args.skip_plot {
        let rscript_path = rseqc_commands::exec_resolve::which(&args.rscript).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("Rscript executable not found: {}", args.rscript),
            )
        })?;
        let status = Command::new(&rscript_path).arg(&r_path).status()?;
        if !status.success() {
            return Err(std::io::Error::other(format!("R plotting failed for {r_path}")));
        }
    }

    Ok(())
}
