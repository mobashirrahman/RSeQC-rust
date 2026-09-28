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
    /// Input alignment file in BAM, plain-text SAM, or CRAM format
    /// (dispatched by the `.bam`/`.sam`/`.cram` extension).
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
            eprintln!("bam_stat.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> std::io::Result<()> {
    // Upstream: `if self.bam_format: print("Load BAM file ... ",
    // end=' ') else: print("Load SAM file ... ", end=' ')` --
    // `self.bam_format` comes from trying `pysam.Samfile(path, 'rb')`
    // FIRST and only falling back to `'r'` (bam_format=False) if that
    // raises. Confirmed via a live diff against real upstream: htslib's
    // `'rb'` open is lenient about actual content and succeeds for a
    // genuine plain-text SAM file too (it auto-detects format,
    // effectively ignoring the 'b' mode hint) -- so `bam_format` is
    // `True`, and "Load BAM file" prints, EVEN for `.sam` input. The
    // "Load SAM file" branch is practically dead code for any valid
    // input, not something this port needs a format check to trigger.
    // The literal's own trailing space plus `end=' '` gives two spaces
    // before "Done".
    eprint!("Load BAM file ...  ");
    let (_header, records) = rseqc_formats::open_alignments(&args.input_file)?;
    let counts = compute_stats(records, args.mapq)?;
    eprintln!("Done");
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
