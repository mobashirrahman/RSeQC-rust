//! Dispatch and flag parsing for `sc_bamStat.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! BAM index (`.bai`) presence is checked (matching upstream's
//! `require_index=True`) but never parsed -- see crates/commands/src/
//! sc_bamstat.rs module docs: the per-chromosome fetch loop is
//! replicated as one sequential BAM scan.

use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::sc_bamstat::{TagNames, mapping_stat, render_report};

#[derive(Parser)]
#[command(name = "sc_bamStat.py", about = "Report mapping statistics for single-cell RNA-seq BAM files.")]
struct Args {
    /// Input BAM file.
    #[arg(short = 'i', long = "infile")]
    bam_file: PathBuf,

    /// BAM tag containing the error-corrected cellular barcode.
    #[arg(long = "cb-tag", default_value = "CB")]
    cb_tag: String,

    /// BAM tag containing the alignment region type.
    #[arg(long = "re-tag", default_value = "RE")]
    re_tag: String,

    /// BAM tag marking alignments on the transcript sense strand.
    #[arg(long = "tx-tag", default_value = "TX")]
    tx_tag: String,

    /// BAM tag marking alignments on the transcript antisense strand.
    #[arg(long = "an-tag", default_value = "AN")]
    an_tag: String,

    /// BAM tag containing the error-corrected UMI.
    #[arg(long = "umi-tag", default_value = "UB")]
    umi_tag: String,

    /// BAM tag marking reads confidently assigned to a feature.
    #[arg(long = "xf-tag", default_value = "xf")]
    xf_tag: String,

    /// Mitochondrial contig name in the BAM header.
    #[arg(long = "chrM-id", default_value = "chrM")]
    mit_contig_name: String,

    /// Enable detailed progress logging.
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

fn validate_tag(option_name: &str, tag: &str) -> std::io::Result<()> {
    if tag.chars().count() != 2 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("{option_name} must be exactly two characters")));
    }
    if tag.chars().any(|c| c.is_whitespace()) {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("{option_name} cannot contain whitespace")));
    }
    Ok(())
}

fn run(args: &Args) -> std::io::Result<()> {
    for (option_name, t) in [
        ("--cb-tag", &args.cb_tag),
        ("--re-tag", &args.re_tag),
        ("--tx-tag", &args.tx_tag),
        ("--an-tag", &args.an_tag),
        ("--umi-tag", &args.umi_tag),
        ("--xf-tag", &args.xf_tag),
    ] {
        validate_tag(option_name, t)?;
    }
    if args.mit_contig_name.trim().is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--chrM-id cannot be empty"));
    }

    let bai_candidate1 = format!("{}.bai", args.bam_file.display());
    let bai_candidate2 = args.bam_file.with_extension("bai");
    if !std::path::Path::new(&bai_candidate1).is_file() && !bai_candidate2.is_file() {
        return Err(std::io::Error::new(std::io::ErrorKind::NotFound, format!("cannot find BAM index; expected either {bai_candidate1} or {}", bai_candidate2.display())));
    }

    if args.verbose {
        eprintln!("Input BAM: {}", args.bam_file.display());
    }

    let tags = TagNames { cb: args.cb_tag.clone(), umi: args.umi_tag.clone(), re: args.re_tag.clone(), tx: args.tx_tag.clone(), an: args.an_tag.clone(), xf: args.xf_tag.clone() };

    // Upstream: `logging.info("Reading BAM file \"%s\" ..." % infile)`
    // -- INFO level, unconditional (not gated by --verbose).
    eprintln!("Reading BAM file \"{}\" ...", args.bam_file.display());
    let (mut reader, header) = rseqc_formats::open_bam(&args.bam_file)?;
    let stats = mapping_stat(reader.records(), &header, &tags, &args.mit_contig_name)?;
    let report = render_report(&stats)?;

    print!("{report}");
    eprintln!("Done.");

    Ok(())
}
