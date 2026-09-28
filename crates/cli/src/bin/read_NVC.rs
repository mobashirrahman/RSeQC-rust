//! Dispatch and flag parsing for `read_NVC.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::read_nvc::{compute_nvc, render_nvc_r_script, render_nvc_table};

#[derive(Parser)]
#[command(
    name = "read_NVC.py",
    about = "Calculate nucleotide frequency at each read cycle."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output NVC file.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Include ambiguous nucleotides N and X in the NVC plot.
    #[arg(short = 'x', long = "nx")]
    nx: bool,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate NVC data but do not execute the R plotting script.
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
    // bam_stat.py, same underlying pysam.Samfile call here). The "Read
    // SAM file" branch is practically dead code for any valid input.
    // The literal's own trailing space plus `end=' '` gives two spaces
    // before "Done".
    eprint!("Read BAM file ...  ");
    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let table = compute_nvc(records, args.mapq)?;
    eprintln!("Done");

    eprintln!("generating data matrix ...");
    let prefix = args.out_prefix.to_string_lossy().into_owned();
    let nvc_path = format!("{prefix}.NVC.xls");
    File::create(&nvc_path)?.write_all(render_nvc_table(&table).as_bytes())?;

    // Upstream: `print("generating R script  ...", ...)` -- literal has
    // two spaces between "script" and "...".
    eprintln!("generating R script  ...");
    let r_path = format!("{prefix}.NVC_plot.r");
    File::create(&r_path)?.write_all(render_nvc_r_script(&table, &prefix, args.nx).as_bytes())?;

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
