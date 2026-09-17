//! Dispatch and flag parsing for `insertion_profile.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::{Parser, ValueEnum};
use rseqc_commands::clipping_profile::{compute_paired_end, compute_single_end};
use rseqc_commands::insertion_profile::{
    render_paired_r_script, render_paired_table, render_single_r_script, render_single_table,
};

#[derive(Clone, Copy, ValueEnum)]
enum Layout {
    #[value(name = "SE")]
    SingleEnd,
    #[value(name = "PE")]
    PairedEnd,
}

#[derive(Parser)]
#[command(
    name = "insertion_profile.py",
    about = "Calculate the distribution of inserted nucleotides across RNA-seq reads."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Sequencing layout: SE for single-end or PE for paired-end.
    #[arg(short = 's', long = "sequencing")]
    sequencing: Layout,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Generate the profile files but do not run the R plotting script.
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
    let (mut reader, _header) = rseqc_formats::open_bam(&args.input_file)?;
    let prefix = args.out_prefix.to_string_lossy();

    // Upstream: `print("Load BAM file ... ", end=' ')` -- the literal's
    // own trailing space plus `end=' '` gives two spaces before "Done".
    eprint!("Load BAM file ...  ");

    let (table_text, r_script_text) = match args.sequencing {
        Layout::SingleEnd => {
            let profile = compute_single_end(reader.records(), args.mapq, b'I')?;
            eprintln!("Done");
            // Upstream: `print("Totoal reads used: %d" % ...)` -- a
            // literal upstream typo ("Totoal"), preserved exactly.
            eprintln!("Totoal reads used: {}", profile.total_read);
            (render_single_table(&profile), render_single_r_script(&profile, &prefix))
        }
        Layout::PairedEnd => {
            let profile = compute_paired_end(reader.records(), args.mapq, b'I')?;
            eprintln!("Done");
            // Upstream prints these as TWO SEPARATE lines (also with
            // the same "Totoal" typo), not one combined line.
            eprintln!("Totoal read-1 used: {}", profile.total_read1);
            eprintln!("Totoal read-2 used: {}", profile.total_read2);
            (render_paired_table(&profile), render_paired_r_script(&profile, &prefix))
        }
    };

    let xls_path = format!("{prefix}.insertion_profile.xls");
    File::create(&xls_path)?.write_all(table_text.as_bytes())?;

    let r_path = format!("{prefix}.insertion_profile.r");
    File::create(&r_path)?.write_all(r_script_text.as_bytes())?;

    if !args.skip_plot {
        let status = Command::new(&args.rscript).arg(&r_path).status()?;
        if !status.success() {
            return Err(std::io::Error::other(format!("R plotting failed for {r_path}")));
        }
    }

    Ok(())
}
