//! Dispatch and flag parsing for `tin.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! `-i/--input` accepts every form upstream's `getBamFiles.get_bam_files`
//! does (single BAM, comma-separated list, directory, list file); the
//! resolved paths are sorted and each needs a `.bai` sidecar, as upstream
//! requires. The index is only checked for presence (region queries are
//! served from an in-memory read index -- see crates/commands/src/tin.rs).
//!
//! Upstream logs at INFO level unconditionally (`--verbose` only enables
//! DEBUG), so the progress lines below are always printed; the
//! `%(asctime)s [%(levelname)s] ` prefix is not replicated (DIV-0019).

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::{Path, PathBuf};

use clap::Parser;
use rseqc_commands::genebody_coverage::get_bam_files;
use rseqc_commands::tin::{
    build_exon_ranges, build_read_index, compute_tin, compute_tin_windowed, genomic_positions, render_summary, render_tin_xls,
    WindowedTin,
};
use rseqc_formats::interval::MergedRegions;
use rseqc_formats::bed::get_exon;

#[derive(Parser)]
#[command(
    name = "tin.py",
    about = "Calculate transcript integrity number (TIN) for each transcript or gene."
)]
struct Args {
    /// BAM input: one BAM, a comma-separated list, a directory, or a list file.
    #[arg(short = 'i', long = "input")]
    input_files: String,

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
            // Upstream: logging.error("%s", exc); return 1
            eprintln!("{err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    if args.sample_size > 1000 {
        eprintln!("--sample-size is greater than 1000; reduce it if performance is poor");
    }
    // Upstream validates the output directory BEFORE reading the alignment:
    //   if not args.output_dir.exists(): parser.error("output directory does not exist: ...")
    //   if not args.output_dir.is_dir():   parser.error("output path is not a directory: ...")
    // and exits 2 through argparse. Without this check the whole alignment is
    // scored first and the failure surfaces as a bare ENOENT when the first output
    // file is created -- after minutes of work, with no message naming the directory
    // that was wrong. Found by the command-contract checks in
    // verification/check_command_contracts.py.
    // Exit 2, as upstream's parser.error() does, not this port's usual 1 for a
    // runtime error. This port already exits 2 for every other usage error because
    // that is what clap does, so returning Err here made these two conditions the
    // only usage errors in the tree with a different status.
    if !args.output_dir.exists() {
        rseqc_cli::usage_exit(
            "tin.py",
            &format!("output directory does not exist: {}", args.output_dir.display()),
        );
    }
    if !args.output_dir.is_dir() {
        rseqc_cli::usage_exit(
            "tin.py",
            &format!("output path is not a directory: {}", args.output_dir.display()),
        );
    }
    eprintln!("Get BAM file(s) ...");
    let mut bam_files = get_bam_files(&args.input_files);
    bam_files.sort();
    if bam_files.is_empty() {
        return Err(std::io::Error::other("No BAM files found"));
    }
    for bam in &bam_files {
        if !bam.is_file() {
            return Err(std::io::Error::other(format!("BAM file does not exist: {}", bam.display())));
        }
        let mut bai = bam.as_os_str().to_owned();
        bai.push(".bai");
        if !Path::new(&bai).is_file() && !bam.with_extension("bai").is_file() {
            return Err(std::io::Error::other(format!(
                "BAM index not found for {}; expected {}.bai or {}",
                bam.display(),
                bam.display(),
                bam.with_extension("bai").display()
            )));
        }
    }
    eprintln!("Total {} BAM file(s)", bam_files.len());
    for bam in &bam_files {
        eprintln!("  {}", bam.display());
    }

    let exon_ranges = if args.subtract_bg {
        let exon_file = File::open(&args.ref_gene_model)?;
        let exons = get_exon(BufReader::new(exon_file))?;
        Some(build_exon_ranges(&exons))
    } else {
        None
    };

    for bam in &bam_files {
        eprintln!("Processing {}", bam.display());
        process_bam(args, bam, exon_ranges.as_ref())?;
    }
    eprintln!("Done.");
    Ok(())
}

fn process_bam(args: &Args, bam: &Path, exon_ranges: Option<&MergedRegions>) -> std::io::Result<()> {
    let refgene_file = File::open(&args.ref_gene_model)?;
    let samples = genomic_positions(BufReader::new(refgene_file), args.sample_size)?;

    let (mut reader, header) = rseqc_formats::open_bam(bam)?;

    // Sliding-window driver: scores the same TIN numbers as the whole-file
    // index while holding only the reads that can still reach an unscored
    // transcript, instead of every read in the BAM. Falls back to the
    // whole-file path if the input turns out not to be coordinate-sorted,
    // which is the one precondition the window relies on.
    let (records, summary) = match compute_tin_windowed(reader.records(), &header, &samples, args.minimum_coverage, exon_ranges)? {
        WindowedTin::Computed(records, summary) => (records, summary),
        WindowedTin::NotCoordinateSorted => {
            eprintln!("BAM is not coordinate-sorted; falling back to whole-file read index (higher memory use)");
            let (mut reader, header) = rseqc_formats::open_bam(bam)?;
            let reads_by_chrom = build_read_index(reader.records(), &header)?;
            compute_tin(&samples, &reads_by_chrom, args.minimum_coverage, exon_ranges)
        }
    };
    for finished in (100..=records.len()).step_by(100) {
        eprintln!("{finished} transcripts finished");
    }

    let stem = bam.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
    // Upstream strips a case-insensitive ".bam" suffix.
    let stem = if stem.to_lowercase().ends_with(".bam") { stem[..stem.len() - 4].to_string() } else { stem };

    // pathlib drops a bare "." component: Path(".") / "x" == Path("x").
    let out_path = |name: String| if args.output_dir == Path::new(".") { PathBuf::from(name) } else { args.output_dir.join(name) };
    let tin_path = out_path(format!("{stem}.tin.xls"));
    let summary_path = out_path(format!("{stem}.summary.txt"));

    let mut tin_file = File::create(&tin_path)?;
    tin_file.write_all(render_tin_xls(&records).as_bytes())?;

    let bam_name = bam.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
    let mut summary_file = File::create(&summary_path)?;
    summary_file.write_all(render_summary(&bam_name, &summary).as_bytes())?;

    eprintln!("Created {}", tin_path.display());
    eprintln!("Created {}", summary_path.display());
    Ok(())
}
