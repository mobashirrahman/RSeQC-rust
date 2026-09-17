//! Dispatch and flag parsing for `sc_seqLogo.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! `--oformat svg`/`png` (DIV-0016, partially closed): renders both of
//! upstream's real logo outputs -- `<prefix>.logo.<format>` (plain,
//! frequency-based) and `<prefix>.logo.mean_centered.<format>`
//! (mean-centered, flipped-below) -- via `rseqc_render::seqlogo`/
//! `seqlogo_png`. See those modules' own doc comments for exactly what
//! is and isn't reproduced (no real font metrics/hinting, `shade_below`/
//! `fade_below` not honored, no working upstream oracle to verify
//! visual output against). `--oformat pdf` remains unimplemented (needs
//! a real PDF-writing crate, out of scope this pass) -- still fails
//! cleanly with a disclosed error, per DIV-0016.
//!
//! **Preserves a genuine upstream quirk, not "fixed"**: the
//! "Mean-centered logo saved to ..." progress line names the file as
//! `<prefix>.logo_mean_centered.<format>` (underscore) via
//! `oracle/upstream-src/src/qcmodule/fastq.py`'s own log string, but
//! the REAL file `plt.savefig(...)` actually writes is
//! `<prefix>.logo.mean_centered.<format>` (dot-separated) -- the two
//! don't match. Confirmed by reading both literal strings in
//! `make_logo()` side by side; reproduced exactly, including the
//! mismatch, not "corrected" to what was probably intended.
//!
//! Compressed (.gz/.Z/.z/.bz/.bz2/.bzip2) input IS supported via
//! `rseqc_formats::open_text_input`, matching upstream's
//! `qcmodule.ireader.nopen` extension dispatch -- independent of the
//! logo-rendering gap above.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::sc_seqlogo::{CountMatrix, compute_count_matrix, fasta_iter, fastq_seq_strings, render_count_matrix_csv};
use rseqc_render::seqlogo::{StackOrder, render_frequency_logo_svg, render_mean_centered_logo_svg};
use rseqc_render::seqlogo_png::{render_frequency_logo_png, render_mean_centered_logo_png};

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

    // Upstream: `fasta_iter`/`fastq_iter` (qcmodule/fastq.py) each open
    // with `logging.info("Reading FASTA/FASTQ file \"%s\" ...")` --
    // unconditional, before any sequence is consumed.
    if args.in_format == "fq" {
        eprintln!("Reading FASTQ file \"{}\" ...", args.in_file.display());
    } else {
        eprintln!("Reading FASTA file \"{}\" ...", args.in_file.display());
    }
    let reader = rseqc_formats::open_text_input(&args.in_file)?;
    let seqs = if args.in_format == "fq" { fastq_seq_strings(reader)? } else { fasta_iter(reader)? };

    // Upstream's `seq2countMat` also prints a `"%d sequences
    // finished\r"` progress line every `--step-size` sequences (a
    // `\r`-overwriting in-place counter, `end=' '` no newline) while
    // scanning -- not reproduced here (this port reads all sequences
    // eagerly before counting, so there's no natural mid-scan point to
    // interleave it without restructuring `compute_count_matrix`); a
    // disclosed, low-value gap given it only fires past --step-size
    // sequences (default 10,000) and is a purely cosmetic, ephemeral
    // terminal indicator, not data.
    let matrix = compute_count_matrix(&seqs, args.max_seq, args.exclude_n)?;
    let sequence_length = matrix.rows.len() as i64;

    // Upstream: `logging.info("%d sequences finished" % count)` after
    // the scan loop -- `count` there is every sequence ITERATED (even
    // ones later skipped for containing "N"), capped at --nseq-limit
    // if the loop's own `break` fired first.
    let processed_count = match args.max_seq {
        Some(limit) => seqs.len().min(limit.max(0) as usize),
        None => seqs.len(),
    };
    eprintln!("{processed_count} sequences finished");
    eprintln!("Make data frame from dict of dict ...");
    eprintln!("Filling NA as zero ...");

    if let Some(s) = args.highlight_start {
        if s >= sequence_length || args.highlight_end.unwrap() >= sequence_length {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("--highlight-start/--highlight-end is outside the valid range 0..{}", sequence_length - 1)));
        }
    }

    let prefix = args.out_file.to_string_lossy().into_owned();
    let count_matrix_path = format!("{prefix}.count_matrix.csv");
    File::create(&count_matrix_path)?.write_all(render_count_matrix_csv(&matrix).as_bytes())?;

    let logo_path = format!("{prefix}.logo.{}", args.out_format);
    if !matches!(args.out_format.as_str(), "svg" | "png") {
        return Err(std::io::Error::other(format!(
            "sequence logo was not created: {logo_path} (native rendering is only implemented for --oformat svg/png in this port -- see DIV-0016 in compatibility/divergences.yaml; {count_matrix_path} was written successfully)"
        )));
    }

    // Upstream: `logging.info("Making logo ...")` -- unconditional.
    eprintln!("Making logo ...");
    // Upstream: `if exclude_N: logging.info("'N' will be excluded.")
    // else: logging.info("'N' will be kept.")` -- unconditional.
    if args.exclude_n {
        eprintln!("'N' will be excluded.");
    } else {
        eprintln!("'N' will be kept.");
    }

    let stack_order = StackOrder::parse(&args.stack_order).ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidInput, "--stack-order must be 'big_on_top', 'small_on_top', or 'fixed'"))?;
    let highlight = match (args.highlight_start, args.highlight_end) {
        (Some(s), Some(e)) => Some((s, e)),
        _ => None,
    };

    // Upstream: `logging.info("Mean-centered logo saved to \"%s\"." %
    // (outfile + '.logo_mean_centered.' + oformat))` -- this literal
    // string does NOT match the real file `plt.savefig(...)` writes
    // (`.logo.mean_centered.<format>`, dot-separated); reproduced
    // exactly, including the mismatch (see module doc comment).
    eprintln!("Mean-centered logo saved to \"{prefix}.logo_mean_centered.{}\".", args.out_format);
    if let Some((s, e)) = highlight {
        eprintln!("Highlight logo from {s} to {e}");
    }
    let mean_centered_path = format!("{prefix}.logo.mean_centered.{}", args.out_format);
    write_logo(&mean_centered_path, &matrix, stack_order, highlight, &args.out_format, true)?;

    // Upstream: `logging.info("Logo saved to \"%s\"." % (outfile +
    // '.logo.' + oformat))` -- unconditional.
    eprintln!("Logo saved to \"{prefix}.logo.{}\".", args.out_format);
    if let Some((s, e)) = highlight {
        eprintln!("Highlight logo from {s} to {e}");
    }
    write_logo(&logo_path, &matrix, stack_order, highlight, &args.out_format, false)?;

    Ok(())
}

fn write_logo(path: &str, matrix: &CountMatrix, stack_order: StackOrder, highlight: Option<(i64, i64)>, out_format: &str, centered: bool) -> std::io::Result<()> {
    match out_format {
        "svg" => {
            let svg = if centered {
                render_mean_centered_logo_svg(&matrix.bases, &matrix.rows, stack_order, highlight)
            } else {
                render_frequency_logo_svg(&matrix.bases, &matrix.rows, stack_order, highlight)
            };
            File::create(path)?.write_all(svg.as_bytes())
        }
        "png" => {
            let png = if centered {
                render_mean_centered_logo_png(&matrix.bases, &matrix.rows, stack_order, highlight)?
            } else {
                render_frequency_logo_png(&matrix.bases, &matrix.rows, stack_order, highlight)?
            };
            File::create(path)?.write_all(&png)
        }
        _ => unreachable!("out_format already validated to be svg or png"),
    }
}
