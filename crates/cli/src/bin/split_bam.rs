//! Dispatch and flag parsing for `split_bam.py`. Binary name can't contain
//! '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN Step 10) adds
//! the `.py`-suffixed PATH alias.
//!
//! `--index-output` (same DIV-0006 pattern as `split_paired_bam.py`,
//! `divide_bam.py`): writes a `.bai` per output BAM via
//! `rseqc_formats::write_bai_index`. Upstream's own `index_bam()` here has
//! the identical `pysam.index("-f", path)` bug (`-f` is not a valid
//! `samtools index` option in any version) that makes `--index-output`
//! crash on every invocation of the real upstream CLI -- this port
//! implements the working, intended behavior instead; see DIV-0006 in
//! compatibility/divergences.yaml for the full account.
//!
//! Upstream's logging is unconditional at INFO level (`--verbose` only
//! raises the threshold to DEBUG, which INFO messages already clear) --
//! "Reading gene model <path>", "Loaded exon intervals for N
//! chromosome(s)", "Splitting <input>", one "Indexing <path>" per
//! `--index-output` file, and "Done." all print regardless of
//! `--verbose`. No timestamp prefix is emitted (same approach as
//! `sc_bamStat.py`/DIV-0019 -- upstream's own prefix can never be
//! byte-reproduced regardless of content).

use std::collections::HashSet;
use std::fs::File;
use std::path::{Path, PathBuf};

use clap::Parser;
use noodles_bam as bam;
use rseqc_commands::split_bam::{build_exon_ranges, render_report, split_bam};
use rseqc_formats::bed::get_exon;

#[derive(Parser)]
#[command(
    name = "split_bam.py",
    about = "Split a BAM file according to exon regions from a BED gene list."
)]
struct Args {
    /// Input BAM file.
    #[arg(short = 'i', long = "input-file")]
    input_file: PathBuf,

    /// Reference gene list in BED format.
    #[arg(short = 'r', long = "genelist")]
    gene_list: PathBuf,

    /// Prefix for the three output BAM files.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: PathBuf,

    /// Create BAI indexes for the output BAM files after splitting.
    #[arg(long = "index-output")]
    index_output: bool,

    /// Allow existing output BAM/index files to be replaced.
    #[arg(long = "overwrite")]
    overwrite: bool,

    /// Enable detailed progress logging (accepted; upstream's messages are
    /// always emitted at INFO level regardless of this flag -- verbose
    /// only lowers the threshold, it does not gate these specific lines).
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

fn run(args: &Args) -> std::io::Result<()> {
    // Upstream: `logging.info("Reading gene model %s", gene_list)` --
    // unconditional.
    eprintln!("Reading gene model {}", args.gene_list.display());
    let refgene_file = File::open(&args.gene_list)?;
    let exons = get_exon(std::io::BufReader::new(refgene_file))?;
    if exons.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("no exon intervals were found in {}", args.gene_list.display()),
        ));
    }
    let chrom_count = exons.iter().map(|(chrom, _, _)| chrom.to_uppercase()).collect::<HashSet<_>>().len();
    // Upstream: `logging.info("Loaded exon intervals for %d
    // chromosome(s)", len(exon_ranges))` -- unconditional.
    eprintln!("Loaded exon intervals for {chrom_count} chromosome(s)");
    let exon_ranges = build_exon_ranges(&exons);

    let (mut reader, header) = rseqc_formats::open_bam(&args.input_file)?;

    let [in_path, ex_path, junk_path] = output_paths(&args.out_prefix);
    check_output_paths(&[in_path.clone(), ex_path.clone(), junk_path.clone()], args.overwrite)?;

    let mut in_writer = bam::io::Writer::new(File::create(&in_path)?);
    let mut ex_writer = bam::io::Writer::new(File::create(&ex_path)?);
    let mut junk_writer = bam::io::Writer::new(File::create(&junk_path)?);
    in_writer.write_header(&header)?;
    ex_writer.write_header(&header)?;
    junk_writer.write_header(&header)?;

    let mut outputs = [in_writer, ex_writer, junk_writer];
    // Upstream: `logging.info("Splitting %s", input_file)` -- unconditional.
    eprintln!("Splitting {}", args.input_file.display());
    let counts = split_bam(reader.records(), &header, &exon_ranges, &mut outputs)?;
    drop(outputs);

    print!("{}", render_report(&in_path.to_string_lossy(), &ex_path.to_string_lossy(), &junk_path.to_string_lossy(), &counts));

    if args.index_output {
        for path in [&in_path, &ex_path, &junk_path] {
            // Upstream: `logging.info("Indexing %s", path)` -- unconditional.
            eprintln!("Indexing {}", path.display());
            rseqc_formats::write_bai_index(path)?;
        }
    }

    // Upstream: `logging.info("Done.")` -- unconditional, not gated by
    // --verbose (see the module doc comment).
    eprintln!("Done.");

    Ok(())
}

fn output_paths(prefix: &Path) -> [PathBuf; 3] {
    [
        PathBuf::from(format!("{}.in.bam", prefix.display())),
        PathBuf::from(format!("{}.ex.bam", prefix.display())),
        PathBuf::from(format!("{}.junk.bam", prefix.display())),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_paths_use_upstream_suffixes() {
        let paths = output_paths(&PathBuf::from("out/sample"));
        assert_eq!(paths[0], PathBuf::from("out/sample.in.bam"));
        assert_eq!(paths[1], PathBuf::from("out/sample.ex.bam"));
        assert_eq!(paths[2], PathBuf::from("out/sample.junk.bam"));
    }

    #[test]
    fn existing_output_requires_overwrite() {
        let root = std::env::temp_dir().join(format!("rseqc_split_bam_test_{}", std::process::id()));
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
