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
//! upstream invokes best-effort through `subprocess.call(..., shell=True)`)
//! is replicated the same way: the identical command string is run via
//! `/bin/sh -c`, preceded (unstranded mode only) by upstream's
//! "Run wigToBigWig ..." stdout line; see `try_wig_to_bigwig`.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::bam2wig::{build_wig_signal, cal_wig_sum, load_chrom_sizes, parse_strand_rule, render_stranded_wig, render_unstranded_wig};
use rseqc_commands::python_fmt::python_str_float;

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

    // Upstream: `print(f"Skip multi-hits: {args.skip_multi}", ...)` --
    // Python's bool str() is "True"/"False" (capitalized), not Rust's
    // lowercase Display impl.
    eprintln!("Skip multi-hits: {}", if args.skip_multi { "True" } else { "False" });

    let chrom_sizes = load_chrom_sizes(BufReader::new(File::open(&args.chrom_size)?))?;
    let strand_map = parse_strand_rule(args.strand_rule.as_deref()).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;
    let strand_rule_active = args.strand_rule.is_some();

    // Upstream: `calWigSum` (called during the wigsum phase) and
    // `bamTowig` (the main phase) EACH run their own independent
    // per-chromosome scan, and each prints its OWN "Processing <chrom>
    // ..." / "No alignments for <chrom>. skipped" progress line -- so
    // with `--wigsum` given, these lines appear TWICE, once per phase.
    let print_chrom_progress = |valid_chroms: &HashSet<String>| {
        for (chrom, _) in &chrom_sizes {
            if !valid_chroms.contains(chrom) {
                eprintln!("No alignments for {chrom}. skipped");
            } else {
                eprintln!("Processing {chrom} ...");
            }
        }
    };

    let normalization_factor = match args.total_wigsum {
        None => None,
        Some(target) => {
            eprintln!("Calcualte wigsum ... ");
            let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
            let chrom_names: HashSet<String> = header.reference_sequences().keys().map(|k| k.to_string()).collect();
            print_chrom_progress(&chrom_names);
            let listed: HashSet<String> = chrom_sizes.iter().map(|(c, _)| c.clone()).filter(|c| chrom_names.contains(c)).collect();
            let wig_sum = cal_wig_sum(reader.records(), &header, &listed, args.skip_multi)?;
            // Upstream's `wig_sum` is a Python float, so an f-string
            // always shows it as e.g. "100.0", never bare "100".
            eprintln!("Total WIG sum: {}", python_str_float(wig_sum));
            if wig_sum <= 0.0 {
                return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, format!("normalization cannot be calculated because the observed WIG sum is {wig_sum:?}")));
            }
            Some(target / wig_sum)
        }
    };

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let valid_chroms: HashSet<String> = header.reference_sequences().keys().map(|k| k.to_string()).collect();
    print_chrom_progress(&valid_chroms);

    let signal = build_wig_signal(reader.records(), &header, strand_rule_active, &strand_map, args.skip_multi, args.map_qual)?;

    let prefix = args.out_prefix.to_string_lossy();
    if !strand_rule_active {
        let wig_path = format!("{prefix}.wig");
        File::create(&wig_path)?.write_all(render_unstranded_wig(&chrom_sizes, &valid_chroms, &signal, normalization_factor).as_bytes())?;
        // Upstream prints the (flag-less) command line to stdout first.
        println!("Run wigToBigWig {prefix}.wig {} {prefix}.bw ", args.chrom_size.display());
        std::io::stdout().flush()?;
        try_wig_to_bigwig(&[format!("wigToBigWig -clip {prefix}.wig {} {prefix}.bw ", args.chrom_size.display())]);
    } else {
        let (fwd, rev) = render_stranded_wig(&chrom_sizes, &valid_chroms, &signal, normalization_factor);
        let fwd_path = format!("{prefix}.Forward.wig");
        let rev_path = format!("{prefix}.Reverse.wig");
        File::create(&fwd_path)?.write_all(fwd.as_bytes())?;
        File::create(&rev_path)?.write_all(rev.as_bytes())?;
        let cs = args.chrom_size.display();
        try_wig_to_bigwig(&[
            format!("wigToBigWig -clip {fwd_path} {cs} {prefix}.Forward.bw "),
            format!("wigToBigWig -clip {rev_path} {cs} {prefix}.Reverse.bw "),
        ]);
    }

    Ok(())
}

/// Best-effort `wigToBigWig -clip <wig> <chrom_size> <bw>` invocation(s),
/// reproducing upstream's `subprocess.call(<string>, shell=True)` inside a
/// bare `try/except: pass`: each command string is handed to `/bin/sh -c`
/// verbatim, so a missing tool yields the shell's own "not found" message
/// (and a nonzero status, which upstream ignores). "Failed to call" is only
/// printed when the shell itself cannot be spawned -- the only case in which
/// `subprocess.call` raises -- and then aborts the remaining calls, as the
/// exception would.
fn try_wig_to_bigwig(commands: &[String]) {
    for cmd in commands {
        if Command::new("/bin/sh").arg("-c").arg(cmd).status().is_err() {
            eprintln!("Failed to call \"wigToBigWig\".");
            return;
        }
    }
}
