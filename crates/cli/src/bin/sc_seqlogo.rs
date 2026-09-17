//! Dispatch and flag parsing for `sc_seqLogo.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! **Only the `.count_matrix.csv` output is produced.** The
//! `.logo.<format>` sequence-logo image is NOT rendered -- see
//! crates/commands/src/sc_seqlogo.rs module docs and DIV-0016 in
//! compatibility/divergences.yaml. This command therefore always exits
//! non-zero (matching upstream's own behavior when logo generation
//! fails for any reason: `if not logo_path.is_file(): raise
//! RuntimeError(...)`), after writing the real, correct count matrix.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::sc_seqlogo::{compute_count_matrix, fasta_iter, fastq_seq_strings, render_count_matrix_csv};

#[derive(Parser)]
#[command(name = "sc_seqLogo.py", about = "Generate a DNA sequence logo from FASTA, FASTQ, or sequence-only input.")]
struct Args {
    /// Input FASTA or FASTQ file.
    #[arg(short = 'i', long = "infile")]
    in_file: PathBuf,

    /// Prefix for the count matrix and sequence-logo files.
    #[arg(short = 'o', long = "outfile")]
    out_file: PathBuf,

    /// Input format.
    #[arg(long = "iformat", default_value = "fq")]
    in_format: String,

    /// Sequence-logo output format.
    #[arg(long = "oformat", default_value = "pdf")]
    out_format: String,

    /// Maximum number of sequences to process.
    #[arg(short = 'n', long = "nseq-limit")]
    max_seq: Option<i64>,

    #[arg(long = "font-name", default_value = "sans")]
    font_name: String,
    #[arg(long = "stack-order", default_value = "big_on_top")]
    stack_order: String,
    #[arg(long = "flip-below")]
    flip_below: bool,
    #[arg(long = "shade-below", default_value_t = 0.0)]
    shade_below: f64,
    #[arg(long = "fade-below", default_value_t = 0.0)]
    fade_below: f64,
    #[arg(long = "exclude-N", visible_alias = "excludeN")]
    exclude_n: bool,
    #[arg(long = "highlight-start")]
    highlight_start: Option<i64>,
    #[arg(long = "highlight-end")]
    highlight_end: Option<i64>,
    #[arg(long = "step-size", default_value_t = 10_000)]
    step_size: i64,

    #[arg(long = "verbose")]
    verbose: bool,
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
    if !matches!(args.in_format.as_str(), "fq" | "fa") {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--iformat must be 'fq' or 'fa'"));
    }
    if !matches!(args.out_format.as_str(), "pdf" | "png" | "svg") {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--oformat must be 'pdf', 'png', or 'svg'"));
    }
    if let Some(n) = args.max_seq {
        if n <= 0 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--nseq-limit must be greater than zero"));
        }
    }
    if args.step_size <= 0 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--step-size must be greater than zero"));
    }
    if !(0.0..=1.0).contains(&args.shade_below) || !(0.0..=1.0).contains(&args.fade_below) {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--shade-below and --fade-below must be between 0 and 1"));
    }
    match (args.highlight_start, args.highlight_end) {
        (Some(s), Some(e)) if e < s => return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--highlight-end must be greater than or equal to --highlight-start")),
        (Some(_), None) => return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--highlight-end is required when --highlight-start is supplied")),
        (None, Some(_)) => return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--highlight-start is required when --highlight-end is supplied")),
        _ => {}
    }

    let reader = BufReader::new(File::open(&args.in_file)?);
    let seqs = if args.in_format == "fq" { fastq_seq_strings(reader)? } else { fasta_iter(reader)? };

    let matrix = compute_count_matrix(&seqs, args.max_seq, args.exclude_n)?;
    let sequence_length = matrix.rows.len() as i64;

    if let Some(s) = args.highlight_start {
        if s >= sequence_length || args.highlight_end.unwrap() >= sequence_length {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("--highlight-start/--highlight-end is outside the valid range 0..{}", sequence_length - 1)));
        }
    }

    let prefix = args.out_file.to_string_lossy().into_owned();
    let count_matrix_path = format!("{prefix}.count_matrix.csv");
    File::create(&count_matrix_path)?.write_all(render_count_matrix_csv(&matrix).as_bytes())?;
    eprintln!("Created {count_matrix_path}");

    let logo_path = format!("{prefix}.logo.{}", args.out_format);
    Err(std::io::Error::other(format!(
        "sequence logo was not created: {logo_path} (native rendering not implemented in this port -- see DIV-0016 in compatibility/divergences.yaml; {count_matrix_path} was written successfully)"
    )))
}
