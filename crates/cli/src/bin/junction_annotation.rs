//! Dispatch and flag parsing for `junction_annotation.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! `--skip-plot`/`--rscript` (running Rscript to produce the PDF) are not
//! implemented -- disclosed gap, this only generates the .r script text.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::junction_annotation::{
    build_model, compute_junction_annotation, render_bed12, render_interact, render_r_script, render_xls,
};

#[derive(Parser)]
#[command(
    name = "junction_annotation.py",
    about = "Annotate splice junctions against a reference gene model."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference gene model (BED12 format).
    #[arg(short = 'r', long = "refgene")]
    refgene: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Minimum intron size.
    #[arg(short = 'm', long = "min-intron", default_value_t = 50)]
    min_intron: i64,

    /// Minimum mapping quality for a read to be used.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Skip running Rscript to render the PDF (always skipped; disclosed gap).
    #[arg(long = "skip-plot")]
    skip_plot: bool,

    /// Path to Rscript executable (accepted, unused: no native/subprocess rendering).
    #[arg(long = "rscript", default_value = "Rscript")]
    rscript: String,

    /// Skip generating the .bed file.
    #[arg(long = "skip-bed")]
    skip_bed: bool,

    /// Skip generating the .Interact.bed file.
    #[arg(long = "skip-interact")]
    skip_interact: bool,
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
    let _ = &args.rscript;
    let _ = args.skip_plot;

    let refgene_file = File::open(&args.refgene)?;
    let model = build_model(BufReader::new(refgene_file))?;

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let result = compute_junction_annotation(reader.records(), &header, &model, args.mapq, args.min_intron)?;

    let prefix = args.out_prefix.to_string_lossy();

    let mut xls_file = File::create(format!("{prefix}.junction.xls"))?;
    xls_file.write_all(render_xls(&result).as_bytes())?;

    let mut r_file = File::create(format!("{prefix}.junction_plot.r"))?;
    r_file.write_all(render_r_script(&result, &prefix).as_bytes())?;

    if !args.skip_bed {
        let mut bed_file = File::create(format!("{prefix}.bed"))?;
        bed_file.write_all(render_bed12(&result, 1).as_bytes())?;
    }

    if !args.skip_interact {
        let bam_name = args.input_file.to_string_lossy();
        let mut interact_file = File::create(format!("{prefix}.Interact.bed"))?;
        interact_file.write_all(render_interact(&result, &bam_name, 1).as_bytes())?;
    }

    Ok(())
}
