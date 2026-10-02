//! Dispatch and flag parsing for `deletion_profile.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::deletion_profile::{
    compute_deletion_profile, render_deletion_r_script, render_deletion_table,
};

#[derive(Parser)]
#[command(
    name = "deletion_profile.py",
    about = "Calculate the distribution of deleted nucleotides across aligned reads."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input")]
    input_file: PathBuf,

    /// Expected aligned read length; reads whose sequence length or
    /// M/S/I-op total doesn't exactly match are skipped.
    #[arg(short = 'l', long = "read-align-length")]
    read_align_length: usize,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Maximum number of qualifying (has-deletion, right-length) reads to process.
    #[arg(short = 'n', long = "read-num", default_value_t = 1_000_000)]
    read_num: u64,

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
            eprintln!("deletion_profile.py: error: {err}");
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
    rseqc_cli::require_existing_output_parent_or_exit("deletion_profile.py", &args.out_prefix);

    let (mut reader, _header) = rseqc_formats::open_bam(&args.input_file)?;
    // Upstream: `print("Process BAM file ... ", end=' ', file=sys.stderr)`
    // -- the string literal's own trailing space plus `end=' '` gives two
    // spaces before "Total reads used" on the same stderr line.
    eprint!("Process BAM file ...  ");
    let profile = compute_deletion_profile(
        reader.records(),
        args.mapq,
        args.read_align_length,
        args.read_num,
    )?;
    eprintln!("Total reads used: {}", profile.count);
    // Upstream's unconditional `print('\n')`: the literal "\n" plus
    // print's own trailing newline is two bytes.
    println!();
    println!();

    let prefix = args.out_prefix.to_string_lossy();

    let mut table = File::create(format!("{prefix}.deletion_profile.txt"))?;
    table.write_all(render_deletion_table(&profile).as_bytes())?;

    let r_path = format!("{prefix}.deletion_profile.r");
    File::create(&r_path)?.write_all(render_deletion_r_script(&profile, &prefix).as_bytes())?;

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
