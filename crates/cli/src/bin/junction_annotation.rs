//! Dispatch and flag parsing for `junction_annotation.py`. Binary name
//! can't contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::{Command, ExitCode};

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
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
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

    /// Skip running Rscript to render the PDF.
    #[arg(long = "skip-plot")]
    skip_plot: bool,

    /// Path to Rscript executable.
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
            eprintln!("junction_annotation.py: error: {err}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<ExitCode> {
    // Upstream: `print("Reading reference bed file: ",refgene, " ... ",
    // end=' ', file=sys.stderr)` -- print's default `sep=' '` between
    // positional args, combined with each string's own leading/trailing
    // spaces, gives the exact double-space pattern reproduced literally
    // below (verified via python3 -c).
    eprint!("Reading reference bed file:  {}  ...  ", args.refgene.display());
    let refgene_file = File::open(&args.refgene)?;
    let model = build_model(BufReader::new(refgene_file))?;
    eprintln!("Done");

    // Upstream: `if self.bam_format: print("Load BAM file ... ", end=' ')
    // else: print("Load SAM file ... ", end=' ')` -- dead-code else
    // branch, same as bam_stat.py and others (pysam.Samfile(path, 'rb')
    // succeeds for genuine .sam content too).
    eprint!("Load BAM file ...  ");
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let result = compute_junction_annotation(records, &header, &model, args.mapq, args.min_intron)?;
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

    let r_path = format!("{prefix}.junction_plot.r");
    File::create(&r_path)?.write_all(render_r_script(&result, &prefix).as_bytes())?;

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

    // Upstream's main() calls run_plot_script here, AFTER
    // annotate_junction() returns normally (both earlier "No splice
    // .../sys.exit()" paths bypass main() entirely via SystemExit, so
    // Rscript is never invoked on those -- matched by returning before
    // reaching this point on both of those paths above).
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
