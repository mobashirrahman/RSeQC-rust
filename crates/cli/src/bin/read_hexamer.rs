//! Dispatch and flag parsing for `read_hexamer.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufReader, Write as _};
use std::path::{Path, PathBuf};

use clap::Parser;
use rseqc_commands::read_hexamer::{kmer_freq_file, render_report, unique_display_name};

#[derive(Parser)]
#[command(name = "read_hexamer.py", about = "Calculate normalized hexamer frequencies from FASTA/FASTQ files.")]
struct Args {
    /// Comma-separated FASTA/FASTQ read files, for example 'reads_1.fq,reads_2.fa'.
    #[arg(short = 'i', long = "input")]
    input_reads: String,

    /// Optional reference-genome FASTA file.
    #[arg(short = 'r', long = "refgenome")]
    ref_genome: Option<PathBuf>,

    /// Optional reference mRNA/transcript FASTA file.
    #[arg(short = 'g', long = "refgene")]
    ref_gene: Option<PathBuf>,

    /// Write the frequency table to this file instead of standard output.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Skip missing read files with a warning. By default, any missing
    /// input file is treated as an error.
    #[arg(long = "skip-missing")]
    skip_missing: bool,
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

/// Parses and validates the comma-separated read-file list. Ports
/// `parse_read_files`: a missing file is a hard error unless
/// `skip_missing`, in which case it's a stderr warning -- but if
/// `skip_missing` leaves NO valid files at all, that is still a hard
/// error (`"none of the requested input read files exists"`).
fn parse_read_files(value: &str, skip_missing: bool) -> io::Result<Vec<PathBuf>> {
    let files: Vec<PathBuf> = value.split(',').map(str::trim).filter(|s| !s.is_empty()).map(PathBuf::from).collect();
    if files.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "--input must contain at least one file"));
    }

    let mut valid_files = Vec::new();
    for path in files {
        if path.is_file() {
            valid_files.push(path);
            continue;
        }
        if skip_missing {
            eprintln!("Warning: input file does not exist and will be skipped: {}", path.display());
        } else {
            return Err(io::Error::new(io::ErrorKind::NotFound, format!("input file does not exist: {}", path.display())));
        }
    }

    if valid_files.is_empty() {
        return Err(io::Error::new(io::ErrorKind::NotFound, "none of the requested input read files exists"));
    }

    Ok(valid_files)
}

fn validate_optional_file(path: &Option<PathBuf>, label: &str) -> io::Result<()> {
    if let Some(p) = path {
        if !p.is_file() {
            return Err(io::Error::new(io::ErrorKind::NotFound, format!("{label} file does not exist: {}", p.display())));
        }
    }
    Ok(())
}

#[allow(clippy::type_complexity)]
fn collect_inputs(
    read_files: &[PathBuf],
    ref_genome: &Option<PathBuf>,
    ref_gene: &Option<PathBuf>,
) -> io::Result<(Vec<String>, HashMap<String, HashMap<String, i64>>, HashMap<String, f64>)> {
    let mut names = Vec::new();
    let mut tables = HashMap::new();
    let mut totals = HashMap::new();
    let mut existing_names: HashSet<String> = HashSet::new();

    let mut ordered_paths = read_files.to_vec();
    if let Some(rg) = ref_genome {
        ordered_paths.push(rg.clone());
    }
    if let Some(rg) = ref_gene {
        ordered_paths.push(rg.clone());
    }

    for path in ordered_paths {
        let name = unique_display_name(&path, &mut existing_names);
        eprint!("Calculate hexamer frequencies for {} ... ", path.display());
        let counts = kmer_freq_file(BufReader::new(File::open(&path)?))?;
        let total = counts.values().sum::<i64>() as f64;
        eprintln!("Done");

        names.push(name.clone());
        tables.insert(name.clone(), counts);
        totals.insert(name, total);
    }

    Ok((names, tables, totals))
}

fn run(args: &Args) -> io::Result<()> {
    let read_files = parse_read_files(&args.input_reads, args.skip_missing)?;
    validate_optional_file(&args.ref_genome, "reference genome")?;
    validate_optional_file(&args.ref_gene, "reference transcript")?;

    let (names, tables, totals) = collect_inputs(&read_files, &args.ref_genome, &args.ref_gene)?;
    let report = render_report(&names, &tables, &totals);

    if let Some(output_path) = &args.output {
        // Upstream: `output_parent = args.output.parent`. For a bare
        // relative filename like "result.tsv" (no directory component),
        // Python's `Path.parent` is `Path('.')` -- the current directory,
        // which exists and is a directory -- NOT an error. Confirmed via
        // a live `python3 -c` probe. Map an empty `Path::parent()` the
        // same way instead of treating it as "no parent directory".
        let output_parent = match output_path.parent() {
            Some(p) if p.as_os_str().is_empty() => Path::new("."),
            Some(p) => p,
            None => Path::new("."),
        };
        if !output_parent.is_dir() {
            return Err(io::Error::new(io::ErrorKind::NotFound, format!("output directory does not exist: {}", output_parent.display())));
        }
        File::create(output_path)?.write_all(report.as_bytes())?;
        eprintln!("Created: {}", output_path.display());
    } else {
        io::stdout().write_all(report.as_bytes())?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_read_files_valid() {
        let dir = std::env::temp_dir().join(format!("read_hexamer_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f1 = dir.join("a.fa");
        let f2 = dir.join("b.fa");
        std::fs::write(&f1, ">x\nACGT\n").unwrap();
        std::fs::write(&f2, ">y\nACGT\n").unwrap();

        let input = format!("{},{}", f1.display(), f2.display());
        let result = parse_read_files(&input, false).unwrap();
        assert_eq!(result, vec![f1.clone(), f2.clone()]);

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn parse_read_files_empty_input_errors() {
        assert!(parse_read_files("", false).is_err());
    }

    #[test]
    fn parse_read_files_skip_missing_keeps_existing_and_skips_missing() {
        let dir = std::env::temp_dir().join(format!("read_hexamer_test2_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let exists = dir.join("exists.fa");
        std::fs::write(&exists, ">x\nACGT\n").unwrap();
        let missing = dir.join("missing.fa");

        let input = format!("{},{}", exists.display(), missing.display());
        let result = parse_read_files(&input, true).unwrap();
        assert_eq!(result, vec![exists.clone()]);

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn parse_read_files_skip_missing_all_missing_still_errors() {
        // Upstream: skip_missing only skips INDIVIDUAL missing files with a
        // warning; if NONE of them exist, it's still a hard error.
        let result = parse_read_files("nonexistent1.fa,nonexistent2.fa", true);
        assert!(result.is_err());
    }

    #[test]
    fn parse_read_files_without_skip_missing_errors_on_first_missing() {
        assert!(parse_read_files("nonexistent.fa", false).is_err());
    }
}
