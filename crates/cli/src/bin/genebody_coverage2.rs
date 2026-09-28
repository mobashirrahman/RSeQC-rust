//! Dispatch and flag parsing for `geneBody_coverage2.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging
//! (PORTING_PLAN Step 10) adds the `.py`-suffixed PATH alias.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::genebody_coverage2::{coverage_gene_body_bigwig, render_coverage_txt, render_r_script};
use rseqc_formats::bigwig::BigWigReader;

#[derive(Parser)]
#[command(name = "geneBody_coverage2.py", about = "Calculate RNA-seq coverage across the gene body from a BigWig file.")]
struct Args {
    /// Coverage signal file in BigWig format.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference gene model in BED12 format.
    #[arg(short = 'r', long = "refgene")]
    ref_gene_model: PathBuf,

    /// Prefix for generated output files.
    #[arg(short = 'o', long = "out-prefix")]
    output_prefix: PathBuf,

    /// Plot file type.
    #[arg(short = 't', long = "graph-type", default_value = "pdf")]
    graph_type: String,

    /// Generate data and R code but do not execute the R script.
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
            eprintln!("geneBody_coverage2.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    if !matches!(args.graph_type.as_str(), "pdf" | "png" | "bmp" | "jpeg" | "tiff") {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--graph-type must be one of pdf, png, bmp, jpeg, tiff"));
    }

    let prefix = args.output_prefix.to_string_lossy().into_owned();
    let data_path = format!("{prefix}.geneBodyCoverage.txt");
    let r_script_path = format!("{prefix}.geneBodyCoverage_plot.r");

    let mut bw = BigWigReader::open(&args.input_file)?;
    eprintln!("Calculating coverage over gene body ...");
    let (coverage, gene_count) = coverage_gene_body_bigwig(&mut bw, BufReader::new(File::open(&args.ref_gene_model)?))?;

    // Matches upstream: both output files are written unconditionally
    // (even for gene_count == 0) BEFORE the zero-gene check errors.
    File::create(&data_path)?.write_all(render_coverage_txt(&coverage).as_bytes())?;
    let script = render_r_script(&coverage, gene_count, &prefix, &args.graph_type);
    File::create(&r_script_path)?.write_all(script.as_bytes())?;

    if gene_count == 0 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "no valid BED12 records matched chromosomes in the BigWig"));
    }

    if !args.skip_plot {
        let rscript_path = rseqc_commands::exec_resolve::which(&args.rscript).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("Rscript executable not found: \"{}\"", args.rscript),
            )
        })?;
        let status = Command::new(&rscript_path).arg(&r_script_path).status()?;
        if !status.success() {
            return Err(std::io::Error::other(format!("R plotting failed for {r_script_path}")));
        }
    }

    eprintln!("Created: {data_path}");
    eprintln!("Created: {r_script_path}");

    Ok(())
}
