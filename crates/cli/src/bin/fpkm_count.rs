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
use rseqc_commands::python_fmt::python_str_float;
use rseqc_commands::fpkm_count::{
    build_global_exon_ranges, build_read_index, compute_fpkm_rows, compute_fpkm_rows_windowed,
    count_total_fragments_streaming, parse_strand_rule, render_fpkm_xls,
    WindowedFpkm,
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
            eprintln!("FPKM_count.py: error: {err}");
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
    rseqc_cli::require_existing_output_parent_or_exit("FPKM_count.py", &args.out_prefix);

    let strand_map = parse_strand_rule(args.strand_rule.as_deref())
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

    eprintln!("Extract exon regions from {}...", args.refgene_bed.display());
    let global_exon_ranges = build_global_exon_ranges(BufReader::new(File::open(&args.refgene_bed)?))?;

    // Streaming totals: same arithmetic on the same retained reads as
    // `count_total_fragments`, without building the whole-file index first. The
    // sums are over 0/0.5/1 -- exact in binary floating point -- so the river
    // order this makes deterministic was never load-bearing.
    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;

    // Upstream: print("Counting total fragment ... ", end=" ") ... print("Done")
    eprint!("Counting total fragment ...  ");
    let totals = count_total_fragments_streaming(
        reader.records(),
        &header,
        &global_exon_ranges,
        args.single_read,
        args.skip_multi,
        args.map_qual,
    )?;
    eprintln!("Done");
    let (total_frags, exonic_frags) = (totals.total_frags, totals.exonic_frags);
    // Both totals are Python floats (initialised to 0.0), formatted `:<20`.
    eprintln!("Total fragment = {:<20}", python_str_float(total_frags));
    eprintln!("Total exonic fragment = {:<20}", python_str_float(exonic_frags));

    if total_frags <= 0.0 || exonic_frags <= 0.0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "total and exonic fragment counts must be positive",
        ));
    }
    let denominator = if args.only_exon { exonic_frags } else { total_frags };

    // The whole-file path, used whenever the sliding-window driver cannot run.
    // This is the reference implementation: every recorded differential result
    // was produced with it, so falling back can only ever be slower, never wrong.
    let whole_file = || -> std::io::Result<Vec<rseqc_commands::fpkm_count::FpkmRow>> {
        let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
        let reads_by_chrom =
            build_read_index(reader.records(), &header, args.skip_multi, args.map_qual)?;
        // `denominator` is reused rather than recomputed: the streaming totals run
        // the same `add_total_fragment` arithmetic on the same filtered reads, and
        // the terms are exact, so the two are equal by construction.
        compute_fpkm_rows(
            BufReader::new(File::open(&args.refgene_bed)?),
            &reads_by_chrom,
            args.strand_rule.is_some(),
            &strand_map,
            args.single_read,
            denominator,
            |n| eprint!("\r{n} transcripts finished"),
        )
    };

    // Sliding-window driver: same rows, but only the reads that can still reach
    // an unscored transcript stay resident. Exact rather than approximate --
    // transcripts are visited in coordinate order, so a read ending at or before
    // the current transcript's start cannot overlap it or any later one, and
    // every read the whole-file path would select is still in the window.
    //
    // Coordinate order is established up front by the totals pass rather than
    // inferred from the window driver's own probe, which misses disorder it
    // never pulls far enough to see.
    let rows = if !totals.coordinate_sorted {
        eprintln!("BAM is not coordinate-sorted; falling back to whole-file read index (higher memory use)");
        whole_file()?
    } else {
        let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
        match compute_fpkm_rows_windowed(
            reader.records(),
            &header,
            BufReader::new(File::open(&args.refgene_bed)?),
            args.strand_rule.is_some(),
            &strand_map,
            args.single_read,
            denominator,
            args.skip_multi,
            args.map_qual,
            |n| eprint!("\r{n} transcripts finished"),
        )? {
            WindowedFpkm::Computed(rows) => rows,
            WindowedFpkm::NotCoordinateSorted => {
                eprintln!(
                    "BAM is not coordinate-sorted; falling back to whole-file read index (higher memory use)"
                );
                whole_file()?
            }
        }
    };
    eprintln!();

    let output_path = format!("{}.FPKM.xls", args.out_prefix.to_string_lossy());
    let mut output_file = File::create(&output_path)?;
    output_file.write_all(render_fpkm_xls(&rows).as_bytes())?;

    eprintln!("Created {output_path}");

    Ok(())
}
