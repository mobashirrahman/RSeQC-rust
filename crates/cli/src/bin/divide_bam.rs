//! Dispatch and flag parsing for `divide_bam.py`. Binary name can't contain '.'
//! (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds the
//! `.py`-suffixed PATH alias.
//!
//! `--index` (DIV-0006, closed): writes a `.bai` index for each output
//! BAM via `rseqc_formats::write_bai_index`, matching upstream's
//! `pysam.index(path)` call in `create_indexes`. Unlike
//! `split_paired_bam.py`, upstream's `create_indexes` has no logging
//! calls at all, so no extra stderr output accompanies this.
//!
//! Upstream has no `--overwrite` flag for this command: it always refuses
//! to run if any `<prefix>_<index>.bam` output already exists
//! (`main()`'s `existing = [... path.exists() ...]` check before any
//! writer is opened).

use std::fs::File;
use std::io;
use std::path::PathBuf;

use clap::Parser;
use noodles_bam as bam;
use rseqc_commands::divide_bam::{divide_bam, render_report};
use rand::SeedableRng;
use rand::rngs::StdRng;

#[derive(Parser)]
#[command(
    name = "divide_bam.py",
    about = "Randomly divide a BAM file into approximately equal subsets."
)]
struct Args {
    /// Input BAM file. Upstream (`divide_bam.py`'s own `validate_args`)
    /// explicitly rejects any non-`.bam` extension via `parser.error` before
    /// ever reading the file -- unlike `bam_stat.py`/`bam2fq.py`, this
    /// command genuinely is BAM-only by upstream's own design, not a
    /// disclosed gap (see DIV-0004).
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Number of output BAM files to create.
    #[arg(short = 'n', long = "subset-num")]
    subset_num: usize,

    /// Prefix for the output BAM files (`<prefix>_<index>.bam`).
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Skip unmapped reads.
    #[arg(short = 's', long = "skip-unmap")]
    skip_unmap: bool,

    /// Random seed for reproducible division (optional).
    #[arg(long = "seed")]
    seed: Option<u64>,

    /// Create a BAM index for each output file after writing.
    #[arg(long = "index")]
    index: bool,
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

fn output_paths(prefix: &str, subset_num: usize) -> Vec<String> {
    (0..subset_num).map(|i| format!("{prefix}_{i}.bam")).collect()
}

fn run(args: &Args) -> std::io::Result<()> {
    if args.subset_num == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "subset-num must be greater than 0",
        ));
    }

    let prefix = args.out_prefix.to_string_lossy();
    let paths = output_paths(&prefix, args.subset_num);

    // Upstream's `main()` refuses to run (before opening any writer) if
    // any output path already exists; there is no `--overwrite` escape
    // hatch for this command.
    let existing: Vec<&String> = paths.iter().filter(|p| std::path::Path::new(p).exists()).collect();
    if !existing.is_empty() {
        let listed = existing.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ");
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("refusing to overwrite existing output file(s): {listed}"),
        ));
    }

    let (reader, header) = rseqc_formats::open_bam(&args.input_file)?;

    let mut outputs = Vec::new();
    for path in &paths {
        let file = File::create(path)?;
        let writer = bam::io::Writer::new(file);
        outputs.push(writer);
    }

    for output in &mut outputs {
        output.write_header(&header)?;
    }

    let mut reader = reader;

    let mut rng = if let Some(seed) = args.seed {
        StdRng::seed_from_u64(seed)
    } else {
        StdRng::from_entropy()
    };

    // Upstream: `print(f"Dividing {input_file} ...", file=sys.stderr)`
    // before the read loop, `print("Done", file=sys.stderr)` after it.
    eprintln!("Dividing {} ...", args.input_file.display());
    let counts = divide_bam(
        reader.records(),
        &header,
        &mut outputs,
        args.skip_unmap,
        &mut rng,
    )?;
    eprintln!("Done");

    drop(outputs);
    drop(reader);

    if args.index {
        for path in &paths {
            rseqc_formats::write_bai_index(std::path::Path::new(path))?;
        }
    }

    let (stdout, stderr) = render_report(&paths, &counts, args.skip_unmap);
    print!("{stdout}");
    eprint!("{stderr}");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_paths_use_prefix_index_suffix() {
        let paths = output_paths("out/sample", 3);
        assert_eq!(paths, vec!["out/sample_0.bam", "out/sample_1.bam", "out/sample_2.bam"]);
    }
}