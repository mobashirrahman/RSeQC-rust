//! Dispatch and flag parsing for `geneBody_coverage.py`. Binary name
//! can't contain '.'/uppercase-first-letter conventions cleanly across
//! all shells in this workspace's naming scheme (see crates/cli/
//! Cargo.toml note); packaging (PORTING_PLAN Step 10) adds the literal
//! `geneBody_coverage.py` PATH alias.
//!
//! No BAI index support (project-wide); see crates/commands/src/
//! genebody_coverage.rs module docs for the pileup-default caveats
//! shared with tin.py (DIV-0011).

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::genebody_coverage::{
    build_index, compute_coverage_for_bam, compute_coverage_windowed, genebody_percentile, get_bam_files, load_dataset,
    make_unique_sample_name, render_coverage_txt, valid_name, write_r_code, WindowedCoverage,
};

#[derive(Parser)]
#[command(name = "geneBody_coverage.py", about = "Calculate RNA-seq read coverage across the gene body.")]
struct Args {
    /// A BAM file, comma-separated BAM files, a directory containing BAM files, or a text file listing BAM paths.
    #[arg(short = 'i', long = "input")]
    input_files: String,

    /// Reference gene model in BED12 format.
    #[arg(short = 'r', long = "refgene")]
    ref_gene_model: PathBuf,

    /// Minimum transcript length in bp (must be >= 100).
    #[arg(short = 'l', long = "minimum-length", default_value_t = 100)]
    min_mrna_length: i64,

    /// Plot output format.
    #[arg(short = 'f', long = "format", default_value = "pdf")]
    output_format: String,

    /// Prefix for output files.
    #[arg(short = 'o', long = "out-prefix")]
    output_prefix: PathBuf,

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
            eprintln!("geneBody_coverage.py: error: {err}");
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
    rseqc_cli::require_existing_output_parent_or_exit("geneBody_coverage.py", &args.output_prefix);

    if args.min_mrna_length < 100 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--minimum-length cannot be smaller than 100"));
    }
    if !matches!(args.output_format.as_str(), "pdf" | "png" | "jpeg") {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--format must be one of pdf, png, jpeg"));
    }

    let prefix = args.output_prefix.to_string_lossy().into_owned();
    let coverage_path = format!("{prefix}.geneBodyCoverage.txt");
    let plot_prefix = format!("{prefix}.geneBodyCoverage");

    eprintln!("Read BED file (reference gene model) ...");
    let (transcripts, transcript_count) = genebody_percentile(BufReader::new(File::open(&args.ref_gene_model)?), args.min_mrna_length)?;
    eprintln!("Total {transcript_count} transcripts loaded");

    eprintln!("Get BAM file(s) ...");
    let bam_files = get_bam_files(&args.input_files);
    if bam_files.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "no BAM files were found"));
    }
    for f in &bam_files {
        eprintln!("\t{}", f.display());
    }

    let mut sample_names: Vec<String> = Vec::new();
    let mut samples: Vec<(String, Vec<i64>, Vec<bool>)> = Vec::new();

    for bam_file in &bam_files {
        if !bam_file.is_file() {
            eprintln!("Warning: BAM file does not exist; skipped: {}", bam_file.display());
            continue;
        }
        eprintln!("Processing {} ...", bam_file.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default());

        let (mut reader, header) = rseqc_formats::open_bam(bam_file)?;
        // Sliding-window driver: same aggregate, but only the reads that can
        // still reach an unscored transcript stay resident. Exact rather than
        // approximate because the accumulation is a sum plus an OR, both
        // commutative, so transcript order does not affect the result. Falls
        // back to the whole-file index if the input is not coordinate-sorted.
        let (coverage, float_markers) = match compute_coverage_windowed(reader.records(), &header, &transcripts)? {
            WindowedCoverage::Computed(coverage, float_markers) => (coverage, float_markers),
            WindowedCoverage::NotCoordinateSorted => {
                eprintln!("BAM is not coordinate-sorted; falling back to whole-file read index (higher memory use)");
                let (mut reader, header) = rseqc_formats::open_bam(bam_file)?;
                let reads_by_chrom = build_index(reader.records(), &header)?;
                compute_coverage_for_bam(&reads_by_chrom, &header, &transcripts)
            }
        };

        if coverage.is_empty() {
            eprintln!("\nCannot get coverage signal from {}! Skip", bam_file.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default());
            continue;
        }

        let base_name = valid_name(&bam_file.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default());
        let sample_name = make_unique_sample_name(&base_name, &sample_names);
        sample_names.push(base_name);
        samples.push((sample_name, coverage, float_markers));
    }

    File::create(&coverage_path)?.write_all(render_coverage_txt(&samples).as_bytes())?;

    // Upstream raises ZeroDivisionError here (a sample with no coverage over
    // the model) and geneBody_coverage.py reports it as
    // "geneBody_coverage.py: error: float division by zero" with exit 1. The
    // "geneBody_coverage.py: error: " prefix is added once, by main(), for
    // every error this binary reports. The coverage .txt is already written at
    // this point, matching upstream.
    let dataset = load_dataset(&samples).map_err(|e| std::io::Error::other(e.to_string()))?;
    if dataset.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, "no valid coverage profiles were generated"));
    }

    eprintln!("\n\n");
    eprintln!("\tSample\tSkewness");
    for d in &dataset {
        eprintln!("\t{}\t{}", d.name, d.skewness);
    }

    let script = write_r_code(&dataset, &plot_prefix, &args.output_format);
    let r_script_path = format!("{plot_prefix}.r");
    File::create(&r_script_path)?.write_all(script.as_bytes())?;

    if !args.skip_plot {
        eprintln!("Running R script ...");
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

    eprintln!("Created: {coverage_path}");
    eprintln!("Created: {r_script_path}");

    Ok(())
}
