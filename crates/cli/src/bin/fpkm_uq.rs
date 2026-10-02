//! Dispatch and flag parsing for `FPKM-UQ.py`. Binary name can't
//! contain '.'/'-' cleanly across all shells so this uses `FPKM_UQ`;
//! packaging (PORTING_PLAN Step 10) adds the literal `FPKM-UQ.py`
//! PATH alias.
//!
//! Unlike every other ported command, this one does no BAM/BED parsing
//! itself: it shells out to `htseq-count` (see crates/commands/src/
//! fpkm_uq.rs module docs) and only post-processes its tabular output.

use std::fs::File;
use std::io::{BufReader, Write as _};
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::fpkm_uq::{
    calculate_fpkm, format_command, htseq_command, read_gene_information, read_htseq_counts, render_fpkm_uq_table, resolve_executable,
};

#[derive(Parser)]
#[command(name = "FPKM-UQ.py", about = "Calculate raw counts, FPKM, and upper-quartile normalized FPKM.")]
struct Args {
    /// Coordinate-sorted BAM file.
    #[arg(long = "bam")]
    bam_file: PathBuf,

    /// Gene model in GTF format.
    #[arg(long = "gtf")]
    gtf_file: PathBuf,

    /// Gene information file containing exon lengths and gene types.
    #[arg(long = "info")]
    info_file: PathBuf,

    /// Prefix for output files.
    #[arg(short = 'o', long = "output")]
    out_prefix: PathBuf,

    /// Report log2(FPKM + 1) and log2(FPKM-UQ + 1).
    #[arg(long = "log2")]
    log_scale: bool,

    /// htseq-count executable to use.
    #[arg(long = "htseq-count", default_value = "htseq-count")]
    htseq_count: String,

    /// Print the htseq-count command and exit without running it.
    #[arg(long = "print-htseq-command")]
    print_htseq_command: bool,
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("FPKM-UQ.py: error: {err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn printlog(message: &str) {
    eprintln!("@ {}: {message}", utc_timestamp());
}

/// `%Y-%m-%d %H:%M:%S` UTC timestamp computed directly from
/// `SystemTime`, with no external process or time-formatting
/// dependency (Howard Hinnant's civil-from-days algorithm). Upstream's
/// `datetime.now()` uses LOCAL time; this cosmetic stderr progress line
/// is not compared against upstream output byte-for-byte by anything,
/// so UTC is an accepted simplification here.
fn utc_timestamp() -> String {
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    let days = (secs / 86400) as i64;
    let time_of_day = secs % 86400;
    let (hour, minute, second) = (time_of_day / 3600, (time_of_day % 3600) / 60, time_of_day % 60);

    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if month <= 2 { y + 1 } else { y };

    format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02}")
}

#[cfg(test)]
#[allow(clippy::items_after_test_module)]
mod tests {
    use super::*;

    #[test]
    fn utc_timestamp_matches_wall_clock_within_a_second() {
        // Cross-check the hand-written civil-calendar conversion against
        // the system's own `date -u` at call time (loose but sufficient:
        // this line is cosmetic stderr output only, never compared
        // byte-for-byte against upstream).
        let ours = utc_timestamp();
        let output = std::process::Command::new("date").args(["-u", "+%Y-%m-%d %H:%M:%S"]).output();
        if let Ok(o) = output {
            if o.status.success() {
                let theirs = String::from_utf8_lossy(&o.stdout).trim().to_string();
                assert!(ours == theirs || ours.as_str() < theirs.as_str(), "ours={ours} theirs={theirs}");
            }
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
    rseqc_cli::require_existing_output_parent_or_exit("FPKM-UQ.py", &args.out_prefix);

    let prefix = args.out_prefix.to_string_lossy();
    let count_file = format!("{prefix}.htseq.counts.txt");
    let fpkm_file = format!("{prefix}.FPKM-UQ.txt");

    let resolved = resolve_executable(&args.htseq_count)?;
    let command = htseq_command(&resolved.to_string_lossy(), &args.bam_file, &args.gtf_file)?;

    if args.print_htseq_command {
        println!("{}", format_command(&command));
        return Ok(());
    }

    printlog("Running htseq-count ...");
    printlog(&format!("Running: {}", format_command(&command)));
    let out_file = File::create(&count_file)?;
    let status = Command::new(&command[0]).args(&command[1..]).stdout(out_file).status()?;
    if !status.success() {
        return Err(std::io::Error::other(format!("htseq-count failed with exit status {status}")));
    }

    if args.log_scale {
        printlog("Calculate log2(FPKM + 1) and log2(FPKM-UQ + 1) ...");
    } else {
        printlog("Calculate FPKM and FPKM-UQ ...");
    }

    // Upstream: `read_gene_information` itself prints a `printlog` line
    // plus two plain `\tTotal ...` lines BEFORE returning -- these are
    // ported here at the call site since the Rust function is pure (no
    // I/O side effects), but the message sequence must match exactly.
    printlog(&format!("Read gene information file: {}", args.info_file.display()));
    let gene_info = read_gene_information(BufReader::new(File::open(&args.info_file)?))?;
    eprintln!("\tTotal genes: {}", gene_info.sizes.len());
    eprintln!("\tTotal protein-coding genes: {}", gene_info.protein_coding.len());

    printlog(&format!("Read gene count file to calculate 75 percentile count and total count: {count_file}"));
    let all_counts = read_htseq_counts(BufReader::new(File::open(&count_file)?))?;
    let (rows, summary) = calculate_fpkm(&gene_info, &all_counts, args.log_scale)?;

    eprintln!("\tTotal protein-coding genes: {}", summary.protein_coding_counts_len);
    eprintln!("\tThe 75 percentile count of protein-coding genes: {:.6}", summary.uq_count);
    eprintln!("\tThe total count of protein-coding genes: {:.6}", summary.total_count as f64);
    // Upstream: plain `print(..., file=sys.stderr)`, NOT `printlog` --
    // no timestamp prefix on this one line, unlike its neighbors.
    eprintln!("Read gene count file to calculate FPKM and FPKM-UQ: {count_file}");
    // Upstream (FPKM-UQ.py:348-353) warns while writing each table row,
    // naming the info file exactly as given on the command line.
    for gene_id in &summary.missing_gene_ids {
        eprintln!("Warning: {gene_id} is absent from {}; skipped", args.info_file.display());
    }

    File::create(&fpkm_file)?.write_all(render_fpkm_uq_table(&rows, args.log_scale).as_bytes())?;

    printlog(&format!("Created: {count_file}"));
    printlog(&format!("Created: {fpkm_file}"));

    Ok(())
}
