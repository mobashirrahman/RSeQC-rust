//! Dispatch and flag parsing for `bam2wig.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! BAM index presence/absence is not enforced the way upstream requires
//! it (`find_bam_index`/`.bai` check) -- see crates/commands/src/
//! bam2wig.rs module docs: no BAI support anywhere in this port, region
//! queries are served by a single sequential scan instead.
//!
//! The trailing `wigToBigWig` conversion (an external UCSC tool
//! upstream invokes best-effort, silently ignoring failure) is
//! replicated the same way: attempted once per output `.wig` file,
//! any failure (tool missing or nonzero exit) is reported to stderr
//! and otherwise ignored, matching upstream's bare `try/except: pass`.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::bam2wig::{build_wig_signal, cal_wig_sum, load_chrom_sizes, parse_strand_rule, render_stranded_wig, render_unstranded_wig};

#[derive(Parser)]
#[command(name = "bam2wig.py", about = "Convert a sorted, indexed BAM file into WIG coverage files.")]
struct Args {
    /// Coordinate-sorted BAM file.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Two-column chromosome-size file (chromosome name and length).
    #[arg(short = 's', long = "chrom-size")]
    chrom_size: PathBuf,

    /// Output prefix.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Target total WIG sum used for normalization. Omit to disable normalization.
    #[arg(short = 't', long = "wigsum")]
    total_wigsum: Option<f64>,

    /// Exclude non-unique alignments.
    #[arg(short = 'u', long = "skip-multi-hits")]
    skip_multi: bool,

    /// Strand rule, for example '1++,1--,2+-,2-+'.
    #[arg(short = 'd', long = "strand")]
    strand_rule: Option<String>,

    /// Minimum mapping quality used to identify uniquely mapped reads.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    map_qual: u8,
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
    if let Some(w) = args.total_wigsum {
        if w <= 0.0 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--wigsum must be greater than zero"));
        }
    }

    eprintln!("Skip multi-hits: {}", args.skip_multi);

    let chrom_sizes = load_chrom_sizes(BufReader::new(File::open(&args.chrom_size)?))?;
    let strand_map = parse_strand_rule(args.strand_rule.as_deref()).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
    let strand_rule_active = args.strand_rule.is_some();

    let normalization_factor = match args.total_wigsum {
        None => None,
        Some(target) => {
            eprintln!("Calcualte wigsum ... ");
            let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
            let chrom_names: HashSet<String> = header.reference_sequences().keys().map(|k| k.to_string()).collect();
            let listed: HashSet<String> = chrom_sizes.iter().map(|(c, _)| c.clone()).filter(|c| chrom_names.contains(c)).collect();
            let wig_sum = cal_wig_sum(reader.records(), &header, &listed, args.skip_multi)?;
            eprintln!("Total WIG sum: {wig_sum}");
            if wig_sum <= 0.0 {
                return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("normalization cannot be calculated because the observed WIG sum is {wig_sum:?}")));
            }
            Some(target / wig_sum)
        }
    };

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let valid_chroms: HashSet<String> = header.reference_sequences().keys().map(|k| k.to_string()).collect();
    for (chrom, _) in &chrom_sizes {
        if !valid_chroms.contains(chrom) {
            eprintln!("No alignments for {chrom}. skipped");
        } else {
            eprintln!("Processing {chrom} ...");
        }
    }

    let signal = build_wig_signal(reader.records(), &header, strand_rule_active, &strand_map, args.skip_multi, args.map_qual)?;

    let prefix = args.out_prefix.to_string_lossy();
    if !strand_rule_active {
        let wig_path = format!("{prefix}.wig");
        File::create(&wig_path)?.write_all(render_unstranded_wig(&chrom_sizes, &valid_chroms, &signal, normalization_factor).as_bytes())?;
        try_wig_to_bigwig(&wig_path, &args.chrom_size, &format!("{prefix}.bw"));
    } else {
        let (fwd, rev) = render_stranded_wig(&chrom_sizes, &valid_chroms, &signal, normalization_factor);
        let fwd_path = format!("{prefix}.Forward.wig");
        let rev_path = format!("{prefix}.Reverse.wig");
        File::create(&fwd_path)?.write_all(fwd.as_bytes())?;
        File::create(&rev_path)?.write_all(rev.as_bytes())?;
        try_wig_to_bigwig(&fwd_path, &args.chrom_size, &format!("{prefix}.Forward.bw"));
        try_wig_to_bigwig(&rev_path, &args.chrom_size, &format!("{prefix}.Reverse.bw"));
    }

    Ok(())
}

/// Best-effort `wigToBigWig -clip <wig> <chrom_size> <bw>` invocation,
/// matching upstream's bare `try/except: pass` around the same
/// subprocess call: any failure (tool not on PATH, nonzero exit) is
/// reported to stderr and otherwise ignored.
fn try_wig_to_bigwig(wig_path: &str, chrom_size: &std::path::Path, bw_path: &str) {
    let result = Command::new("wigToBigWig").arg("-clip").arg(wig_path).arg(chrom_size).arg(bw_path).status();
    match result {
        Ok(status) if status.success() => {}
        _ => eprintln!("Failed to call \"wigToBigWig\"."),
    }
}
