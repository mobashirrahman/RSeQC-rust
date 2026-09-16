//! Dispatch and flag parsing for `read_GC.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! Both outputs (.GC.xls and .GC_plot.r) are generated as upstream does.

use std::fs::File;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_gc::{compute_gc, render_gc_table, render_gc_r_script};

#[derive(Parser)]
#[command(
    name = "read_GC.py",
    about = "Calculate the GC-content distribution of aligned reads."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

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
    let (mut reader, _header) = rseqc_formats::open_bam(&args.input_file)?;
    let hist = compute_gc(reader.records(), args.mapq)?;
    
    // Write the GC table to the output file
    let xls_output_path = format!("{}.GC.xls", args.out_prefix.to_string_lossy());
    let mut xls_output_file = File::create(&xls_output_path)?;
    
    let table_content = render_gc_table(&hist);
    use std::io::Write;
    xls_output_file.write_all(table_content.as_bytes())?;
    
    eprintln!("GC table written to: {}", xls_output_path);
    
    // Write the R script to the output file
    let r_output_path = format!("{}.GC_plot.r", args.out_prefix.to_string_lossy());
    let mut r_output_file = File::create(&r_output_path)?;
    
    let r_script_content = render_gc_r_script(&hist, &args.out_prefix.to_string_lossy());
    r_output_file.write_all(r_script_content.as_bytes())?;
    
    eprintln!("R script written to: {}", r_output_path);
    
    Ok(())
}