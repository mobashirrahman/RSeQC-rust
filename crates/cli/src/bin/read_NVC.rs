//! Dispatch and flag parsing for `read_NVC.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! `-x/--nx` and plot generation are not implemented yet — disclosed gap,
//! see crates/commands/src/read_nvc.rs module docs.

use std::fs::File;
use std::path::PathBuf;

use clap::Parser;
use rseqc_commands::read_nvc::{compute_nvc, render_nvc_table};

#[derive(Parser)]
#[command(
    name = "read_NVC.py",
    about = "Calculate nucleotide frequency at each read cycle."
)]
struct Args {
    /// Input BAM file (SAM-text input is not yet supported).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the output NVC file.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Minimum mapping quality for a read to be considered uniquely mapped.
    #[arg(short = 'q', long = "mapq", default_value_t = 30)]
    mapq: u8,

    /// Include ambiguous nucleotides in the calculation (affects plot generation only).
    #[arg(short = 'x', long = "nx")]
    nx: bool,
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
    let table = compute_nvc(reader.records(), args.mapq)?;
    
    // Write the table to the output file
    let output_path = format!("{}.NVC.xls", args.out_prefix.to_string_lossy());
    let mut output_file = File::create(&output_path)?;
    
    let table_content = render_nvc_table(&table);
    use std::io::Write;
    output_file.write_all(table_content.as_bytes())?;
    
    eprintln!("NVC table written to: {}", output_path);
    
    // If nx flag is provided, note that it's not implemented for plotting yet
    if args.nx {
        eprintln!("Note: -x/--nx flag affects plot generation only (not implemented yet)");
    }
    
    Ok(())
}