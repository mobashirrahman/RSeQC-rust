//! Dispatch and flag parsing for `junction_annotation.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! `--skip-plot`/`--rscript` (running Rscript to produce the PDF) are not
//! implemented -- disclosed gap, this only generates the .r script text.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rseqc_commands::junction_annotation::{
    build_model, compute_junction_annotation, render_bed12, render_interact, render_r_script, render_xls,
};

#[derive(Parser)]
#[command(
    name = "junction_annotation.py",
    about = "Annotate splice junctions against a reference gene model."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference gene model (BED12 format).
    #[arg(short = 'r', long = "refgene")]
    refgene: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Minimum intron size.
    #[arg(short = 'm', long = "min-intron", default_value_t = 50)]
    min_intron: i64,

    /// Minimum mapping quality for a read to be used.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Skip running Rscript to render the PDF (always skipped; disclosed gap).
    #[arg(long = "skip-plot")]
    skip_plot: bool,

    /// Path to Rscript executable (accepted, unused: no native/subprocess rendering).
    #[arg(long = "rscript", default_value = "Rscript")]
    rscript: String,

    /// Skip generating the .bed file.
    #[arg(long = "skip-bed")]
    skip_bed: bool,

    /// Skip generating the .Interact.bed file.
    #[arg(long = "skip-interact")]
    skip_interact: bool,
}

fn main() -> ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(code) => code,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<ExitCode> {
    let _ = &args.rscript;
    let _ = args.skip_plot;

    // Upstream: `print("Reading reference bed file: ",refgene, " ... ",
    // end=' ', file=sys.stderr)` -- print's default `sep=' '` between
    // positional args, combined with each string's own leading/trailing
    // spaces, gives the exact double-space pattern reproduced literally
    // below (verified via python3 -c).
    eprint!("Reading reference bed file:  {}  ...  ", args.refgene.display());
    let refgene_file = File::open(&args.refgene)?;
    let model = build_model(BufReader::new(refgene_file))?;
    eprintln!("Done");

    eprint!("Load BAM file ...  ");
    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let result = compute_junction_annotation(reader.records(), &header, &model, args.mapq, args.min_intron)?;
    eprintln!("Done");

    // Upstream: `print("total = " + str(total_junc))` -- unconditional,
    // to stdout, regardless of whether any junctions were found.
    println!("total = {}", result.event_total);

    let prefix = args.out_prefix.to_string_lossy();

    if result.event_total == 0 {
        // Upstream opens OUT/ROUT unconditionally before this check, so
        // both output files exist but are left empty; `sys.exit()` (no
        // arg -> exit code 0) then terminates the whole process
        // immediately -- .bed/.Interact.bed are never reached.
        File::create(format!("{prefix}.junction.xls"))?;
        File::create(format!("{prefix}.junction_plot.r"))?;
        eprintln!("No splice junction found.");
        return Ok(ExitCode::SUCCESS);
    }

    let mut r_file = File::create(format!("{prefix}.junction_plot.r"))?;
    r_file.write_all(render_r_script(&result, &prefix).as_bytes())?;

    eprintln!();
    eprintln!("===================================================================");
    eprintln!("Total splicing  Events:\t{}", result.event_total);
    eprintln!("Known Splicing Events:\t{}", result.event_known);
    eprintln!("Partial Novel Splicing Events:\t{}", result.event_novel3or5);
    eprintln!("Novel Splicing Events:\t{}", result.event_novel35);
    eprintln!("Filtered Splicing Events:\t{}", result.filtered);

    let mut xls_file = File::create(format!("{prefix}.junction.xls"))?;
    xls_file.write_all(render_xls(&result).as_bytes())?;

    if result.rows.is_empty() {
        // Every event was filtered by --min-intron, leaving no unique
        // junction: upstream's second total_junc==0 check (on the
        // recomputed unique-junction count), exit code 1 this time.
        eprintln!("No splice read found");
        return Ok(ExitCode::FAILURE);
    }

    eprintln!();
    eprintln!("Total splicing  Junctions:\t{}", result.rows.len());
    eprintln!("Known Splicing Junctions:\t{}", result.junction_known);
    eprintln!("Partial Novel Splicing Junctions:\t{}", result.junction_novel3or5);
    eprintln!("Novel Splicing Junctions:\t{}", result.junction_novel35);
    eprintln!();
    eprintln!("===================================================================");

    // Upstream derives these paths from `<prefix>.junction.xls` via
    // `Path.with_suffix`/`Path.with_name(stem + ...)`, which only ever
    // strips/reuses the FINAL suffix -- since the xls path's own final
    // suffix is `.xls`, the real output names are `<prefix>.junction.bed`
    // and `<prefix>.junction.Interact.bed`, not `<prefix>.bed`/
    // `<prefix>.Interact.bed`. Verified live via `pathlib.Path(...)
    // .with_suffix('.bed')` / `.with_name(stem + '.Interact.bed')`.
    if !args.skip_bed {
        eprintln!("Create BED file ...");
        let bed_path = format!("{prefix}.junction.bed");
        let mut bed_file = File::create(&bed_path)?;
        bed_file.write_all(render_bed12(&result, 1).as_bytes())?;
        eprintln!("Created: {bed_path}");
    }

    if !args.skip_interact {
        eprintln!("Create Interact file ...");
        let bam_name = args.input_file.to_string_lossy();
        let interact_path = format!("{prefix}.junction.Interact.bed");
        let mut interact_file = File::create(&interact_path)?;
        interact_file.write_all(render_interact(&result, &bam_name, 1).as_bytes())?;
        eprintln!("Created: {interact_path}");
    }

    Ok(ExitCode::SUCCESS)
}
