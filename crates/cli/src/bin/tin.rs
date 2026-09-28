//! Dispatch and flag parsing for `tin.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! Supports multiple BAM file input forms: single file, comma-separated list,
//! directory, or text file listing BAM paths (one per line), matching
//! `getBamFiles.get_bam_files` behavior. Each BAM generates separate
//! `.tin.xls` and `.summary.txt` output files.

use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write as _};
use std::path::{Path, PathBuf};

use clap::Parser;
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

/// Resolve input specification to a list of BAM file paths.
/// Supports: single file, comma-separated list, directory, or text file listing.
fn resolve_bam_files(input_spec: &Path) -> std::io::Result<Vec<PathBuf>> {
    let mut bam_files = Vec::new();

    // Check if it's a directory
    if input_spec.is_dir() {
        for entry in fs::read_dir(input_spec)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_file() {
                let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
                if name.ends_with(".bam") && path.metadata()?.len() > 0 {
                    bam_files.push(path);
                }
            }
        }
        bam_files.sort();
        return Ok(bam_files);
    }

    // Check if it's a regular file
    if input_spec.is_file() {
        // Try reading as a list of BAM paths (one per line)
        if let Ok(file) = File::open(input_spec) {
            let reader = BufReader::new(file);
            let mut lines = Vec::new();
            for line_result in reader.lines() {
                if let Ok(line) = line_result {
                    let line = line.trim();
                    if !line.is_empty() && !line.starts_with('#') {
                        lines.push(line.to_string());
                    }
                }
            }

            // Check if all non-empty lines are valid BAM files
            let mut all_valid_bams = true;
            for line in &lines {
                let path = Path::new(line);
                if !path.is_file() || !line.ends_with(".bam") || path.metadata().ok().map_or(true, |m| m.len() == 0) {
                    all_valid_bams = false;
                    break;
                }
            }

            if all_valid_bams && !lines.is_empty() {
                for line in lines {
                    bam_files.push(PathBuf::from(line));
                }
                bam_files.sort();
                return Ok(bam_files);
            }
        }

        // If not a list file, treat as single BAM file
        let name = input_spec.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.ends_with(".bam") {
            bam_files.push(input_spec.to_path_buf());
            return Ok(bam_files);
        }
    }

    // Check for comma-separated list
    let input_str = input_spec.to_string_lossy();
    if input_str.contains(',') {
        for part in input_str.split(',') {
            let path = PathBuf::from(part);
            if path.is_file() && path.file_name().and_then(|n| n.to_str()).map_or(false, |n| n.ends_with(".bam")) {
                bam_files.push(path);
            }
        }
        if !bam_files.is_empty() {
            bam_files.sort();
            return Ok(bam_files);
        }
    }

    // If nothing matched, try as a single file
    if input_spec.is_file() {
        bam_files.push(input_spec.to_path_buf());
        return Ok(bam_files);
    }

    Ok(bam_files)
}

fn run(args: &Args) -> std::io::Result<()> {
    eprintln!("Get BAM file(s) ...");

    let bam_files = resolve_bam_files(&args.input_spec)?;
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
