//! Dispatch and flag parsing. Binary is named after the original `.py`
//! script (see Cargo.toml [[bin]] name) so PATH invocation matches upstream
//! once the packaging step (PORTING_PLAN Step 10) adds the `.py` alias.

use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::bam_stat::{BamStatCounts, compute_stats};

#[derive(Parser)]
#[command(
    name = "bam_stat.py",
    about = "Summarize mapping statistics for a BAM or SAM alignment file."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported; see
    /// compatibility/divergences.yaml DIV-0002).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,
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
    let mut reader = rseqc_formats::open_bam(&args.input_file)?;
    let counts = compute_stats(reader.records(), args.mapq)?;
    print_report(&counts);
    Ok(())
}

fn print_report(c: &BamStatCounts) {
    println!();
    println!("#==================================================");
    println!("#All numbers are READ count");
    println!("#==================================================");
    println!();
    println!("{:<40}{}", "Total records:", c.total);
    println!();
    println!("{:<40}{}", "QC failed:", c.qc_fail);
    println!("{:<40}{}", "Optical/PCR duplicate:", c.duplicate);
    println!("{:<40}{}", "Non primary hits", c.non_primary);
    println!("{:<40}{}", "Unmapped reads:", c.unmapped);
    println!("{:<40}{}", "mapq < mapq_cut (non-unique):", c.multi_hit);
    println!();
    println!("{:<40}{}", "mapq >= mapq_cut (unique):", c.uniq_hit);
    println!("{:<40}{}", "Read-1:", c.read1);
    println!("{:<40}{}", "Read-2:", c.read2);
    println!("{:<40}{}", "Reads map to '+':", c.forward);
    println!("{:<40}{}", "Reads map to '-':", c.reverse);
    println!("{:<40}{}", "Non-splice reads:", c.non_splice);
    println!("{:<40}{}", "Splice reads:", c.splice);
    println!("{:<40}{}", "Reads mapped in proper pairs:", c.proper_pair);
    println!(
        "{:<40}{}",
        "Proper-paired reads map to different chrom:", c.proper_pair_diff_chrom
    );
}
