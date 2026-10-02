//! Dispatch and flag parsing for `junction_saturation.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rseqc_commands::junction_saturation::{
    build_known_splice_sites, collect_splice_sites, compute_saturation, render_percentile_report, render_r_script, shuffle_sites,
};

#[derive(Parser)]
#[command(
    name = "junction_saturation.py",
    about = "Assess whether splice-junction discovery has reached sequencing saturation."
)]
struct Args {
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
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

    /// Generate saturation data but do not execute the R plotting script.
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
            eprintln!("junction_saturation.py: error: {err}");
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
    rseqc_cli::require_existing_output_parent_or_exit("junction_saturation.py", &args.out_prefix);

    // Upstream (SAM.py's ParseBAM.saturation_junction): `print("reading
    // reference bed file: ", refgene, " ... ", end=' ')` -- a 3-arg
    // print with the default `sep=' '` joins "reading reference bed
    // file: " + " " + refgene + " " + " ... ", then `end=' '` appends
    // one more space; byte-verified via python3 -c.
    eprint!("reading reference bed file:  {}  ...  ", args.refgene.display());
    let refgene_file = File::open(&args.refgene)?;
    let (known_sites, chrom_list) = build_known_splice_sites(BufReader::new(refgene_file))?;
    eprintln!("Done! Total {} known splicing junctions.", known_sites.len());

    // Upstream: `if self.bam_format: print("Load BAM file ... ", end=' ')
    // else: print("Load SAM file ... ", end=' ')` -- dead-code else
    // branch, same as bam_stat.py and others (pysam.Samfile(path, 'rb')
    // succeeds for genuine .sam content too).
    eprint!("Load BAM file ...  ");
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let mut sites = collect_splice_sites(
        records,
        &header,
        &chrom_list,
        args.map_qual,
        args.minimum_intron_size,
    )?;
    eprintln!("Done");

    // Upstream: `print("shuffling alignments ...", end=' ')` then a
    // separate `print("Done")` -- only ONE space before "Done" here
    // (the literal itself has no trailing space, unlike the two prints
    // above whose literals already end in a space).
    eprint!("shuffling alignments ... ");
    let mut rng = StdRng::from_entropy();
    shuffle_sites(&mut sites, &mut rng);
    eprintln!("Done");

    let counts = compute_saturation(
        &sites,
        &known_sites,
        args.minimum_splice_read,
        args.percentile_low_bound,
        args.percentile_up_bound,
        args.percentile_step,
    );
    eprint!("{}", render_percentile_report(&counts));

    let prefix = args.out_prefix.to_string_lossy();
    let r_path = format!("{prefix}.junctionSaturation_plot.r");
    File::create(&r_path)?.write_all(render_r_script(&counts, &prefix).as_bytes())?;

    if !args.skip_plot {
        let rscript_path = rseqc_commands::exec_resolve::which(&args.rscript).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("Rscript executable not found: {}", args.rscript),
            )
        })?;
        let status = Command::new(&rscript_path).arg(&r_path).status()?;
        if !status.success() {
            return Err(std::io::Error::other(format!("R plotting failed for {r_path}")));
        }
    }

    Ok(())
}
