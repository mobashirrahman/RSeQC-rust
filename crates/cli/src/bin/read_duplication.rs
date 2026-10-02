//! Dispatch and flag parsing for `read_duplication.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! All three outputs (.seq.DupRate.xls, .pos.DupRate.xls, and .DupRate_plot.r) are generated as upstream does.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::read_duplication::{compute_duplication, render_dup_r_script, render_pos_dup_table, render_seq_dup_table};

#[derive(Parser)]
#[command(
    name = "read_duplication.py",
    about = "Calculate sequence-based and mapping-based read duplication rates."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Upper bound of occurrence for counting.
    #[arg(short = 'u', long = "up-limit", default_value_t = 500)]
    up_limit: u64,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate duplication data but do not execute the R plotting script.
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
            eprintln!("read_duplication.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // Upstream's validate_args refuses an output prefix whose parent directory does
    // not exist, before any input is read. Omitting it here meant the whole
    // alignment was read and every metric computed, then discarded when the output
    // open failed with "No such file or directory (os error 2)" -- an error naming
    // neither the directory nor the flag, and indistinguishable from a missing
    // input. The shared helper keeps that check in one place so it cannot be
    // forgotten by the next binary.
    rseqc_cli::require_existing_output_parent_or_exit("read_duplication.py", &args.out_prefix);

    // Upstream: `if self.bam_format: print("Load BAM file ... ", end=' ')
    // else: print("Load SAM file ... ", end=' ')` -- `self.bam_format`
    // comes from `pysam.Samfile(path, 'rb')` succeeding, which it does
    // even for genuine plain-text SAM content (htslib auto-detects,
    // ignoring the 'b' mode hint; confirmed via a live diff for
    // bam_stat.py and others, same underlying pysam.Samfile call here).
    // The "Load SAM file" branch is practically dead code for any valid
    // input. The literal's own trailing space plus `end=' '` gives two
    // spaces before "Done".
    eprint!("Load BAM file ...  ");
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let hist = compute_duplication(records, &header, args.mapq)?;
    eprintln!("Done");

    let prefix = args.out_prefix.to_string_lossy().into_owned();

    eprintln!("report duplicte rate based on sequence ...");
    let seq_path = format!("{prefix}.seq.DupRate.xls");
    File::create(&seq_path)?.write_all(render_seq_dup_table(&hist).as_bytes())?;

    eprintln!("report duplicte rate based on mapping  ...");
    let pos_path = format!("{prefix}.pos.DupRate.xls");
    File::create(&pos_path)?.write_all(render_pos_dup_table(&hist).as_bytes())?;

    eprintln!("generate R script ...");
    let r_path = format!("{prefix}.DupRate_plot.r");
    File::create(&r_path)?.write_all(render_dup_r_script(&hist, &prefix, args.up_limit).as_bytes())?;

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
