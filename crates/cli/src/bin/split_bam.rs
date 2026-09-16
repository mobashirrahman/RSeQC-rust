//! Dispatch and flag parsing for `split_bam.py`. Binary name can't contain
//! '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds
//! the `.py`-suffixed PATH alias.
//!
//! `--index-output` (BAI generation) and `--overwrite`'s pre-existing-file
//! check are not implemented -- disclosed gaps, see
//! crates/commands/src/split_bam.rs module docs.

use std::fs::File;
use std::path::PathBuf;

use clap::Parser;
use noodles_bam as bam;
use rseqc_commands::split_bam::{build_exon_ranges, render_report, split_bam};
use rseqc_formats::bed::get_exon;

#[derive(Parser)]
#[command(
    name = "split_bam.py",
    about = "Split a BAM file according to exon regions from a BED gene list."
)]
struct Args {
    /// Input BAM file.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference gene list in BED format.
    #[arg(short = 'r', long = "genelist")]
    gene_list: PathBuf,

    /// Prefix for the three output BAM files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Create BAI indexes for the output BAM files after splitting (not yet implemented).
    #[arg(long = "index-output")]
    index_output: bool,

    /// Allow existing output BAM/index files to be replaced (existence check not yet implemented).
    #[arg(long = "overwrite")]
    overwrite: bool,

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
    let _ = args.overwrite;

    let refgene_file = File::open(&args.gene_list)?;
    let exons = get_exon(std::io::BufReader::new(refgene_file))?;
    if exons.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("no exon intervals were found in {}", args.gene_list.display()),
        ));
    }
    if args.verbose {
        eprintln!("Reading gene model {}", args.gene_list.display());
    }
    let exon_ranges = build_exon_ranges(&exons);

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;

    let prefix = args.out_prefix.to_string_lossy();
    let in_path = format!("{prefix}.in.bam");
    let ex_path = format!("{prefix}.ex.bam");
    let junk_path = format!("{prefix}.junk.bam");

    let mut in_writer = bam::io::Writer::new(File::create(&in_path)?);
    let mut ex_writer = bam::io::Writer::new(File::create(&ex_path)?);
    let mut junk_writer = bam::io::Writer::new(File::create(&junk_path)?);
    in_writer.write_header(&header)?;
    ex_writer.write_header(&header)?;
    junk_writer.write_header(&header)?;

    let mut outputs = [in_writer, ex_writer, junk_writer];
    if args.verbose {
        eprintln!("Splitting {}", args.input_file.display());
    }
    let counts = split_bam(reader.records(), &header, &exon_ranges, &mut outputs)?;
    drop(outputs);

    if args.index_output {
        eprintln!("warning: --index-output is not yet implemented (no BAI writer support); skipping");
    }

    print!("{}", render_report(&in_path, &ex_path, &junk_path, &counts));
    if args.verbose {
        eprintln!("Done.");
    }

    Ok(())
}
