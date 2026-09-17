//! Dispatch and flag parsing for `split_paired_bam.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10)
//! adds the `.py`-suffixed PATH alias.
//!
//! `--index-output` (BAI generation) remains rejected: no BAI writer exists
//! yet, same disclosed gap as `split_bam.py`. `--verbose` (upstream's extra
//! logging) is accepted but not yet wired to any extra output.

use std::fs::File;
use std::path::{Path, PathBuf};

use clap::Parser;
use noodles_bam as bam;
use rseqc_commands::split_paired_bam::{SplitCounts, render_report, split_paired_bam};

#[derive(Parser)]
#[command(
    name = "split_paired_bam.py",
    about = "Split a paired-end BAM into read-1, read-2, and unmapped BAM files."
)]
struct Args {
    /// Input BAM file.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Prefix for the three output BAM files
    /// (`<prefix>.R1.bam`, `<prefix>.R2.bam`, `<prefix>.unmap.bam`).
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Create BAI indexes for the output BAM files after splitting (currently unsupported).
    #[arg(long = "index-output")]
    index_output: bool,

    /// Allow existing output BAM/index files to be replaced.
    #[arg(long = "overwrite")]
    overwrite: bool,

    /// Enable detailed progress logging (currently a no-op).
    #[arg(long = "verbose")]
    verbose: bool,
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

fn output_paths(prefix: &Path) -> [PathBuf; 3] {
    [
        PathBuf::from(format!("{}.R1.bam", prefix.display())),
        PathBuf::from(format!("{}.R2.bam", prefix.display())),
        PathBuf::from(format!("{}.unmap.bam", prefix.display())),
    ]
}

fn check_output_paths(paths: &[PathBuf], overwrite: bool) -> std::io::Result<()> {
    if overwrite {
        return Ok(());
    }

    let mut existing = Vec::new();
    for path in paths {
        if path.exists() {
            existing.push(path.clone());
        }
        let sidecar = PathBuf::from(format!("{}.bai", path.display()));
        if sidecar.exists() {
            existing.push(sidecar);
        }
        if let Some(stem) = path.file_stem() {
            let sibling = path.with_file_name(format!("{}.bai", stem.to_string_lossy()));
            if sibling.exists() {
                existing.push(sibling);
            }
        }
    }
    if existing.is_empty() {
        return Ok(());
    }

    let listed = existing.iter().map(|path| path.display().to_string()).collect::<Vec<_>>().join(", ");
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        format!("output file already exists; use --overwrite to replace: {listed}"),
    ))
}

fn run(args: &Args) -> std::io::Result<()> {
    if args.index_output {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "--index-output is not supported by this build (no BAI writer)",
        ));
    }

    let (reader, header) = rseqc_formats::open_bam(&args.input_file)?;

    let [r1_path, r2_path, unmap_path] = output_paths(&args.out_prefix);
    check_output_paths(&[r1_path.clone(), r2_path.clone(), unmap_path.clone()], args.overwrite)?;

    let mut r1_writer = bam::io::Writer::new(File::create(&r1_path)?);
    let mut r2_writer = bam::io::Writer::new(File::create(&r2_path)?);
    let mut unmap_writer = bam::io::Writer::new(File::create(&unmap_path)?);
    r1_writer.write_header(&header)?;
    r2_writer.write_header(&header)?;
    unmap_writer.write_header(&header)?;

    let mut reader = reader;
    let counts: SplitCounts = split_paired_bam(
        reader.records(),
        &header,
        &mut r1_writer,
        &mut r2_writer,
        &mut unmap_writer,
    )?;
    drop(r1_writer);
    drop(r2_writer);
    drop(unmap_writer);

    print!(
        "{}",
        render_report(&r1_path.to_string_lossy(), &r2_path.to_string_lossy(), &unmap_path.to_string_lossy(), &counts)
    );
    if args.verbose {
        eprintln!("Done.");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_paths_use_upstream_suffixes() {
        let paths = output_paths(&PathBuf::from("out/sample"));
        assert_eq!(paths[0], PathBuf::from("out/sample.R1.bam"));
        assert_eq!(paths[1], PathBuf::from("out/sample.R2.bam"));
        assert_eq!(paths[2], PathBuf::from("out/sample.unmap.bam"));
    }

    #[test]
    fn existing_output_requires_overwrite() {
        let root = std::env::temp_dir().join(format!("rseqc_split_paired_bam_test_{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let paths = output_paths(&root.join("sample"));
        std::fs::write(&paths[0], b"existing").unwrap();
        let err = check_output_paths(&paths, false).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::AlreadyExists);
        assert!(check_output_paths(&paths, true).is_ok());
        std::fs::remove_file(&paths[0]).unwrap();
        std::fs::remove_dir(&root).unwrap();
    }
}
