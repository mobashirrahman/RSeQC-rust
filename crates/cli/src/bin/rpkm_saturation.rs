//! Dispatch and flag parsing for `RPKM_saturation.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step
//! 10) adds the `.py`-suffixed PATH alias.
//!
//! No BAM index (`.bai`) support -- see crates/commands/src/
//! rpkm_saturation.rs module docs; the read scan is purely sequential.
//! `--rscript`/R execution is not invoked here: this binary always
//! generates the `.saturation.r` script text and, unless `--skip-plot`,
//! shells out to the given Rscript executable exactly once (matching
//! upstream's own subprocess-wrapper behavior for this one step).

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::rpkm_saturation::{
    build_block_lists, build_quartile_plot_data, compute_saturation, parse_gene_model, render_raw_xls,
    render_rpkm_xls, render_saturation_r_script, shuffle_block_lists,
};

fn parse_strand_rule(rule: &str) -> Result<std::collections::HashMap<String, char>, String> {
    let mut map = std::collections::HashMap::new();
    let parts: Vec<&str> = rule.split(',').collect();
    if parts.len() == 4 {
        for p in &parts {
            let chars: Vec<char> = p.chars().collect();
            if chars.len() != 3 {
                return Err(format!("Unknown value of: 'strand_rule' {rule}"));
            }
            map.insert(format!("{}{}", chars[0], chars[1]), chars[2]);
        }
    } else if parts.len() == 2 {
        for p in &parts {
            let chars: Vec<char> = p.chars().collect();
            if chars.len() != 2 {
                return Err(format!("Unknown value of: 'strand_rule' {rule}"));
            }
            map.insert(chars[0].to_string(), chars[1]);
        }
    } else {
        return Err(format!("Unknown value of: 'strand_rule' {rule}"));
    }
    Ok(map)
}

#[derive(Parser)]
#[command(name = "RPKM_saturation.py", about = "Assess whether transcript RPKM estimates have reached sequencing saturation.")]
struct Args {
    /// Input alignment file in BAM or plain-text SAM format (dispatched by
    /// the `.bam`/`.sam` extension).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for generated output files.
    #[arg(short = 'o', long = "out-prefix")]
    output_prefix: PathBuf,

    /// Reference gene model in BED12 format.
    #[arg(short = 'r', long = "refgene")]
    refgene_bed: PathBuf,

    /// Strand rule, for example '1++,1--,2+-,2-+'. Omit for unstranded RNA-seq.
    #[arg(short = 'd', long = "strand")]
    strand_rule: Option<String>,

    /// Starting sampling percentile.
    #[arg(short = 'l', long = "percentile-floor", default_value_t = 5)]
    percentile_low_bound: i64,

    /// Ending sampling percentile.
    #[arg(short = 'u', long = "percentile-ceiling", default_value_t = 100)]
    percentile_up_bound: i64,

    /// Sampling percentile step.
    #[arg(short = 's', long = "percentile-step", default_value_t = 5)]
    percentile_step: i64,

    /// Omit transcripts with mean RPKM below this value from the plot.
    #[arg(short = 'c', long = "rpkm-cutoff", default_value_t = 0.01)]
    rpkm_cutoff: f64,

    /// Minimum mapping quality for an alignment to be treated as uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    map_qual: u8,

    /// Generate saturation data and R script without executing R.
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
            eprintln!("error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    let strand_map = match &args.strand_rule {
        Some(rule) => parse_strand_rule(rule).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?,
        None => std::collections::HashMap::new(),
    };
    let strand_rule_active = args.strand_rule.is_some();

    // Upstream: `if self.bam_format: print("Load BAM file ... ", end=' ')
    // else: print("Load SAM file ... ", end=' ')` -- dead-code else
    // branch, same as bam_stat.py and others (pysam.Samfile(path, 'rb')
    // succeeds for genuine .sam content too). The literal's own trailing
    // space plus `end=' '` gives TWO spaces before "Done" (verified via
    // a live python3 -c probe), then a separate `print("Done")` on the
    // same line.
    eprint!("Load BAM file ...  ");
    let (header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let mut lists = build_block_lists(records, &header, true, args.map_qual, strand_rule_active, &strand_map)?;
    eprintln!("Done");

    eprint!("shuffling alignments ... ");
    let mut rng = rand::thread_rng();
    shuffle_block_lists(&mut lists, &mut rng);
    eprintln!("Done");

    let transcripts = parse_gene_model(BufReader::new(File::open(&args.refgene_bed)?))?;

    let result = compute_saturation(
        &transcripts,
        &lists,
        strand_rule_active,
        args.percentile_low_bound,
        args.percentile_up_bound,
        args.percentile_step,
        &args.refgene_bed.display().to_string(),
    )?;

    let prefix = args.output_prefix.to_string_lossy();
    let rpkm_path = format!("{prefix}.eRPKM.xls");
    let raw_path = format!("{prefix}.rawCount.xls");
    // Upstream's main() prints no "Created" messages at all for this
    // command (unlike several other ported commands) -- confirmed by
    // grepping the actual source; do not add any here.
    File::create(&rpkm_path)?.write_all(render_rpkm_xls(&result).as_bytes())?;
    File::create(&raw_path)?.write_all(render_raw_xls(&result).as_bytes())?;

    let plot_data = build_quartile_plot_data(&result, args.rpkm_cutoff);
    let script = render_saturation_r_script(&plot_data, &prefix)?;
    let r_path = format!("{prefix}.saturation.r");
    File::create(&r_path)?.write_all(script.as_bytes())?;

    if !args.skip_plot {
        let status = Command::new(&args.rscript).arg(&r_path).status()?;
        if !status.success() {
            return Err(std::io::Error::other(format!("{} exited with status {status}", args.rscript)));
        }
    }

    Ok(())
}
