//! Dispatch and flag parsing for `tin.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! Supports multiple BAM file input forms: single file, comma-separated list,
//! directory, or text file listing BAM paths (one per line), matching
//! `getBamFiles.get_bam_files` behavior. Each BAM generates separate
//! `.tin.xls` and `.summary.txt` output files.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::genebody_coverage::get_bam_files;
use rseqc_commands::tin::{build_exon_ranges, build_read_index, compute_tin, genomic_positions, render_summary, render_tin_xls};
use rseqc_formats::bed::get_exon;

#[derive(Parser)]
#[command(
    name = "tin.py",
    about = "Calculate transcript integrity number (TIN) for each transcript or gene."
)]
struct Args {
    /// Input BAM file(s): a single file, comma-separated list, directory, or text file listing BAM paths.
    #[arg(short = 'i', long = "input")]
    input_spec: PathBuf,

    /// Reference gene model in standard BED12 format.
    #[arg(short = 'r', long = "refgene")]
    ref_gene_model: PathBuf,

    /// Minimum number of distinct read starts required.
    #[arg(short = 'c', long = "minCov", default_value_t = 10)]
    minimum_coverage: i64,

    /// Number of approximately equally spaced transcript positions.
    #[arg(short = 'n', long = "sample-size", default_value_t = 100)]
    sample_size: i64,

    /// Subtract background estimated from intronic signal.
    #[arg(short = 's', long = "subtract-background")]
    subtract_bg: bool,

    /// Directory for TIN output files.
    #[arg(short = 'o', long = "output-dir", default_value = ".")]
    output_dir: PathBuf,

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

fn run(args: &Args) -> std::io::Result<()> {
    eprintln!("Get BAM file(s) ...");

    let mut bam_files = get_bam_files(&args.input_spec.to_string_lossy());
    bam_files.sort();
    if bam_files.is_empty() {
        return Err(std::io::Error::new(std::io::ErrorKind::NotFound, "No BAM files found"));
    }

    eprintln!("Total {} BAM file(s)", bam_files.len());
    for bam_file in &bam_files {
        eprintln!("  {}", bam_file.display());
    }

    let exon_ranges = if args.subtract_bg {
        let exon_file = File::open(&args.ref_gene_model)?;
        let exons = get_exon(BufReader::new(exon_file))?;
        Some(build_exon_ranges(&exons))
    } else {
        None
    };

    let refgene_file = File::open(&args.ref_gene_model)?;
    let samples = genomic_positions(BufReader::new(refgene_file), args.sample_size)?;

    for bam_file in bam_files {
        eprintln!("Processing {}", bam_file.display());

        let (mut reader, header) = rseqc_formats::open_bam(&bam_file)?;
        let reads_by_chrom = build_read_index(reader.records(), &header)?;

        let (records, summary) = compute_tin(&samples, &reads_by_chrom, args.minimum_coverage, exon_ranges.as_ref());

        let stem = bam_file
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_default();
        let stem = stem.strip_suffix(".bam").or_else(|| stem.strip_suffix(".BAM")).unwrap_or(&stem);

        let tin_path = args.output_dir.join(format!("{stem}.tin.xls"));
        let summary_path = args.output_dir.join(format!("{stem}.summary.txt"));

        let mut tin_file = File::create(&tin_path)?;
        tin_file.write_all(render_tin_xls(&records).as_bytes())?;

        let bam_name = bam_file.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
        let mut summary_file = File::create(&summary_path)?;
        summary_file.write_all(render_summary(&bam_name, &summary).as_bytes())?;

        eprintln!("Created {}", tin_path.display());
        eprintln!("Created {}", summary_path.display());
    }

    eprintln!("Done.");
    Ok(())
}
