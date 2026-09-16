//! Dispatch and flag parsing for `junction_saturation.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! `--skip-plot`/`--rscript` (running Rscript to produce the PDF) are not
//! implemented -- disclosed gap, this only generates the .r script text.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;

use clap::Parser;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rseqc_commands::junction_saturation::{
    build_known_splice_sites, collect_splice_sites, compute_saturation, render_r_script, shuffle_sites,
};

#[derive(Parser)]
#[command(
    name = "junction_saturation.py",
    about = "Assess whether splice-junction discovery has reached sequencing saturation."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Reference gene model (BED12 format).
    #[arg(short = 'r', long = "refgene")]
    refgene: PathBuf,

    /// Starting sampling percentile.
    #[arg(short = 'l', long = "percentile-floor", default_value_t = 5)]
    percentile_low_bound: i64,

    /// Ending sampling percentile.
    #[arg(short = 'u', long = "percentile-ceiling", default_value_t = 100)]
    percentile_up_bound: i64,

    /// Sampling percentile step.
    #[arg(short = 's', long = "percentile-step", default_value_t = 5)]
    percentile_step: i64,

    /// Minimum intron length in bp.
    #[arg(short = 'm', long = "min-intron", default_value_t = 50)]
    minimum_intron_size: i64,

    /// Minimum number of supporting reads required to call a junction.
    #[arg(short = 'v', long = "min-coverage", default_value_t = 1)]
    minimum_splice_read: u32,

    /// Minimum mapping quality for an alignment to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    map_qual: u8,

    /// Generate saturation data but do not execute the R plotting script (always the case; disclosed gap).
    #[arg(long = "skip-plot")]
    skip_plot: bool,

    /// Rscript executable to use (accepted, unused: no subprocess rendering).
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
    let _ = &args.rscript;
    let _ = args.skip_plot;

    let refgene_file = File::open(&args.refgene)?;
    let (known_sites, chrom_list) = build_known_splice_sites(BufReader::new(refgene_file))?;
    eprintln!("Done! Total {} known splicing junctions.", known_sites.len());

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let mut sites = collect_splice_sites(
        reader.records(),
        &header,
        &chrom_list,
        args.map_qual,
        args.minimum_intron_size,
    )?;

    let mut rng = StdRng::from_entropy();
    shuffle_sites(&mut sites, &mut rng);

    let counts = compute_saturation(
        &sites,
        &known_sites,
        args.minimum_splice_read,
        args.percentile_low_bound,
        args.percentile_up_bound,
        args.percentile_step,
    );

    let prefix = args.out_prefix.to_string_lossy();
    let mut r_file = File::create(format!("{prefix}.junctionSaturation_plot.r"))?;
    r_file.write_all(render_r_script(&counts, &prefix).as_bytes())?;

    Ok(())
}
