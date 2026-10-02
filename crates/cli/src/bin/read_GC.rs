//! Dispatch and flag parsing for `read_GC.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::read_gc::{compute_gc, render_gc_r_script, render_gc_table};

#[derive(Parser)]
#[command(
    name = "read_GC.py",
    about = "Calculate the GC-content distribution of aligned reads."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate GC-profile data but do not execute the R plotting script.
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
            eprintln!("read_GC.py: error: {err}");
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
    rseqc_cli::require_existing_output_parent_or_exit("read_GC.py", &args.out_prefix);

    // Upstream: `if self.bam_format: print("Read BAM file ... ", end=' ')
    // else: print("Read SAM file ... ", end=' ')` -- `self.bam_format`
    // comes from `pysam.Samfile(path, 'rb')` succeeding, which it does
    // even for genuine plain-text SAM content (htslib auto-detects,
    // ignoring the 'b' mode hint; confirmed via a live diff for
    // bam_stat.py/read_NVC.py, same underlying pysam.Samfile call here).
    // The "Read SAM file" branch is practically dead code for any valid
    // input. The literal's own trailing space plus `end=' '` gives two
    // spaces before "Done".
    eprint!("Read BAM file ...  ");
    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let hist = compute_gc(records, args.mapq)?;
    eprintln!("Done");

    eprintln!("writing GC content ...");
    let prefix = args.out_prefix.to_string_lossy().into_owned();
    let xls_path = format!("{prefix}.GC.xls");
    File::create(&xls_path)?.write_all(render_gc_table(&hist).as_bytes())?;

    eprintln!("writing R script ...");
    let r_path = format!("{prefix}.GC_plot.r");
    File::create(&r_path)?.write_all(render_gc_r_script(&hist, &prefix).as_bytes())?;

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
