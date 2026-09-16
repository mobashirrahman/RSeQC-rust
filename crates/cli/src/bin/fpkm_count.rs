//! Dispatch and flag parsing for `FPKM_count.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.
//!
//! BAM index (`.bai`) presence is not checked (no BAI support at all --
//! see crates/commands/src/fpkm_count.rs module docs; region queries are
//! served from an in-memory read index instead).

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::fpkm_count::{
    build_global_exon_ranges, build_read_index, compute_fpkm_rows, count_total_fragments, parse_strand_rule,
    render_fpkm_xls,
};

#[derive(Parser)]
#[command(name = "FPKM_count.py", about = "Calculate fragment counts, FPM, and FPKM for BED12 transcript models.")]
struct Args {
    /// Alignment file in BAM format. SAM is not supported.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Reference gene model in BED12 format.
    #[arg(short = 'r', long = "refgene")]
    refgene_bed: PathBuf,

    /// Strand rule, for example '1++,1--,2+-,2-+'. Omit for unstranded RNA-seq.
    #[arg(short = 'd', long = "strand")]
    strand_rule: Option<String>,

    /// Skip alignments with mapping quality below --mapq.
    #[arg(short = 'u', long = "skip-multi-hits")]
    skip_multi: bool,

    /// Use exonic fragments rather than all fragments for normalization.
    #[arg(short = 'e', long = "only-exonic")]
    only_exon: bool,

    /// Minimum mapping quality.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    map_qual: u8,

    /// Weight for pairs with only one mapped end: 0, 0.5, or 1.
    #[arg(short = 's', long = "single-read", default_value_t = 1.0)]
    single_read: f64,
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
    let strand_map = parse_strand_rule(args.strand_rule.as_deref())
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

    eprintln!("Extract exon regions from {} ...", args.refgene_bed.display());
    let global_exon_ranges = build_global_exon_ranges(BufReader::new(File::open(&args.refgene_bed)?))?;

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let reads_by_chrom = build_read_index(reader.records(), &header, args.skip_multi, args.map_qual)?;

    eprintln!("Counting total fragment ... ");
    let (total_frags, exonic_frags) = count_total_fragments(&reads_by_chrom, &global_exon_ranges, args.single_read);
    eprintln!("Total fragment = {total_frags:<20}");
    eprintln!("Total exonic fragment = {exonic_frags:<20}");

    if total_frags <= 0.0 || exonic_frags <= 0.0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "total and exonic fragment counts must be positive",
        ));
    }
    let denominator = if args.only_exon { exonic_frags } else { total_frags };

    let rows = compute_fpkm_rows(
        BufReader::new(File::open(&args.refgene_bed)?),
        &reads_by_chrom,
        args.strand_rule.is_some(),
        &strand_map,
        args.single_read,
        denominator,
    )?;

    let output_path = format!("{}.FPKM.xls", args.out_prefix.to_string_lossy());
    let mut output_file = File::create(&output_path)?;
    output_file.write_all(render_fpkm_xls(&rows).as_bytes())?;

    eprintln!("Created {output_path}");

    Ok(())
}
