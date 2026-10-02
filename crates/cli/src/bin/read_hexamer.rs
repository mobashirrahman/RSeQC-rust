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

/// The compression format of `path`, from its magic bytes, or None if it is not a
/// compressed file. Extension is not enough: a `.fastq` symlinked to a `.fastq.gz`
/// is gzipped, which is exactly how the benchmark workload supplied this command's
/// input the first time, and the extension said otherwise.
fn compressed_format(path: &Path) -> Option<&'static str> {
    use std::io::Read as _;
    let mut head = [0u8; 4];
    let mut file = File::open(path).ok()?;
    file.read_exact(&mut head).ok()?;
    if head[0] == 0x1f && head[1] == 0x8b {
        Some("gzip")
    } else if &head[..4] == [0x28, 0xb5, 0x2f, 0xfd] {
        Some("zstd")
    } else if &head[..3] == [0x42, 0x5a, 0x68] {
        Some("bzip2")
    } else {
        None
    }
}

fn main() -> std::process::ExitCode {
    let args = Args::parse();
    match run(&args) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("read_hexamer.py: error: {err}");
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
        // Name the cause before opening, because the symptom is unreadable. Upstream
        // reads its inputs as text, so a gzipped FASTQ -- which is what a sequencing
        // facility actually hands you, and what this repository's own held-out data
        // is -- reaches the reader as binary. Upstream then reports a UnicodeDecode
        // error naming a byte offset; this would have reported "stream did not
        // contain valid UTF-8", which names neither the file format nor the fix. Both
        // refuse the input and exit 1, so this changes no outcome; it changes only
        // whether the diagnostic is actionable.
        if let Some(format) = compressed_format(&path) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{} is {} compressed, and this command reads its inputs as plain \
                     text (as upstream does). Decompress it first, e.g. \
                     `gunzip -c {} > {}.plain`",
                    path.display(),
                    format,
                    path.display(),
                    path.display(),
                ),
            ));
        }
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

    #[test]
    fn compressed_format_detects_magic_bytes_not_extensions() {
        // The failure this exists for: a `.fastq` symlink pointing at a `.fastq.gz`
        // is gzipped, and the extension says otherwise. Reading it as text then
        // produces a Unicode error naming neither the format nor the fix.
        let dir = std::env::temp_dir();
        let gz_path = dir.join("rseqc-hexamer-sniff.fastq");
        std::fs::write(&gz_path, [0x1f, 0x8b, 0x08, 0x00]).unwrap();
        assert_eq!(compressed_format(&gz_path), Some("gzip"));

        let bz2_path = dir.join("rseqc-hexamer-sniff-bz2.fastq");
        std::fs::write(&bz2_path, [0x42, 0x5a, 0x68, 0x39]).unwrap();
        assert_eq!(compressed_format(&bz2_path), Some("bzip2"));

        let zs_path = dir.join("rseqc-hexamer-sniff-zstd.fastq");
        std::fs::write(&zs_path, [0x28, 0xb5, 0x2f, 0xfd]).unwrap();
        assert_eq!(compressed_format(&zs_path), Some("zstd"));

        let plain_path = dir.join("rseqc-hexamer-sniff-plain.fastq");
        std::fs::write(&plain_path, b"@r1\nACGT\n+\nIIII\n").unwrap();
        assert_eq!(compressed_format(&plain_path), None);

        for p in [&gz_path, &bz2_path, &zs_path, &plain_path] {
            let _ = std::fs::remove_file(p);
        }
    }

    #[test]
    fn a_compressed_input_is_refused_before_any_work() {
        // Both implementations refuse a gzipped input and exit 1, so no outcome
        // depends on this; the test is that the refusal happens immediately and
        // names the format and the fix, instead of after a Unicode error mid-file.
        let dir = std::env::temp_dir();
        let path = dir.join("rseqc-hexamer-refuse.fastq");
        std::fs::write(&path, [0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00]).unwrap();
        let format = compressed_format(&path);
        assert_eq!(format, Some("gzip"));
        let _ = std::fs::remove_file(&path);
    }

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
