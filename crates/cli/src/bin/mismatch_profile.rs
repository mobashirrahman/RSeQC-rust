//! Dispatch and flag parsing for `mismatch_profile.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::mismatch_profile::{
    compute_mismatch_profile, render_mismatch_r_script, render_mismatch_table,
};

#[derive(Parser)]
#[command(
    name = "mismatch_profile.py",
    about = "Calculate the distribution of mismatches across aligned reads."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported). Must contain MD tags.
    #[arg(short = 'i', long = "input")]
    input_file: PathBuf,

    /// Expected aligned read length; reads whose sequence length or
    /// matched-portion (M-op) total doesn't exactly match are skipped.
    #[arg(short = 'l', long = "read-align-length")]
    read_align_length: usize,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Maximum number of qualifying reads to process.
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
            eprintln!("mismatch_profile.py: error: {err}");
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
    rseqc_cli::require_existing_output_parent_or_exit("mismatch_profile.py", &args.out_prefix);

    let (mut reader, _header) = rseqc_formats::open_bam(&args.input_file)?;
    // Upstream: `print("Process BAM file ... ", end=' ', file=sys.stderr)`
    // -- the string literal's own trailing space plus `end=' '` gives two
    // spaces before whatever prints next on the same stderr line.
    eprint!("Process BAM file ...  ");
    let profile = compute_mismatch_profile(
        reader.records(),
        args.mapq,
        args.read_align_length,
        args.read_num,
    )?;

    if !profile.loop_exhausted_naturally {
        eprintln!("Total reads used: {}", profile.count);
    }
    // Upstream's unconditional `print('\n')`: the literal "\n" plus
    // print's own trailing newline is two bytes, always, regardless of
    // whether any mismatches were found.
    println!();
    println!();

    let prefix = args.out_prefix.to_string_lossy();

    if profile.data.is_empty() {
        // Upstream opens both output files unconditionally before this
        // check (`DOUT = open(...)`, `ROUT = open(...)`), so they exist
        // even though `sys.exit()` fires here before the table header or
        // any data rows are written. The one exception: on natural
        // iterator exhaustion, `except StopIteration: print("Total reads
        // used: " + str(count), file=DOUT)` already ran BEFORE this
        // check, so that single line is present in an otherwise-empty
        // xls file (the table header line only follows it once
        // `len(data) == 0` is known to be false, which never happens
        // here).
        let mut table = File::create(format!("{prefix}.mismatch_profile.xls"))?;
        if profile.loop_exhausted_naturally {
            writeln!(table, "Total reads used: {}", profile.count)?;
        }
        File::create(format!("{prefix}.mismatch_profile.r"))?;
        eprintln!("No mismatches found");
        return Ok(());
    }

    let mut table = File::create(format!("{prefix}.mismatch_profile.xls"))?;
    table.write_all(render_mismatch_table(&profile).as_bytes())?;

    let r_path = format!("{prefix}.mismatch_profile.r");
    File::create(&r_path)?.write_all(render_mismatch_r_script(&profile, &prefix).as_bytes())?;

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
