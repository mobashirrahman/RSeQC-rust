//! Dispatch and flag parsing for `read_duplication.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! All three outputs (.seq.DupRate.xls, .pos.DupRate.xls, and .DupRate_plot.r) are generated as upstream does.

use std::fs::File;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_duplication::{compute_duplication, render_seq_dup_table, render_pos_dup_table, render_dup_r_script};

#[derive(Parser)]
#[command(
    name = "read_duplication.py",
    about = "Calculate sequence-based and mapping-based read duplication rates."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Upper bound of occurrence for counting.
    #[arg(short = 'u', long = "up-limit", default_value_t = 500)]
    up_limit: u64,

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
    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;
    let hist = compute_duplication(reader.records(), &header, args.mapq)?;
    
    // Write the sequence duplication table to the output file
    let seq_xls_output_path = format!("{}.seq.DupRate.xls", args.out_prefix.to_string_lossy());
    let mut seq_xls_output_file = File::create(&seq_xls_output_path)?;
    
    let seq_table_content = render_seq_dup_table(&hist);
    use std::io::Write;
    seq_xls_output_file.write_all(seq_table_content.as_bytes())?;
    
    eprintln!("Sequence duplication table written to: {}", seq_xls_output_path);
    
    // Write the position duplication table to the output file
    let pos_xls_output_path = format!("{}.pos.DupRate.xls", args.out_prefix.to_string_lossy());
    let mut pos_xls_output_file = File::create(&pos_xls_output_path)?;
    
    let pos_table_content = render_pos_dup_table(&hist);
    pos_xls_output_file.write_all(pos_table_content.as_bytes())?;
    
    eprintln!("Position duplication table written to: {}", pos_xls_output_path);
    
    // Write the R script to the output file
    let r_output_path = format!("{}.DupRate_plot.r", args.out_prefix.to_string_lossy());
    let mut r_output_file = File::create(&r_output_path)?;
    
    let r_script_content = render_dup_r_script(&hist, &args.out_prefix.to_string_lossy(), args.up_limit);
    r_output_file.write_all(r_script_content.as_bytes())?;
    
    eprintln!("R script written to: {}", r_output_path);
    
    Ok(())
}