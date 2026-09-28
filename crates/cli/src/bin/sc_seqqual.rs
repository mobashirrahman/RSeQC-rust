//! Dispatch and flag parsing for `sc_seqQual.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! Compressed (.gz/.Z/.z/.bz/.bz2/.bzip2) FASTQ input is supported via
//! `rseqc_formats::open_text_input`, matching upstream's
//! `qcmodule.ireader.nopen` extension dispatch.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::sc_editmatrix::render_heatmap_r_script;
use rseqc_commands::sc_seqqual::{fastq_qual_strings, qual2count_mat, render_quality_matrices};

#[derive(Parser)]
#[command(name = "sc_seqQual.py", about = "Generate sequencing-quality matrices and a heatmap from a FASTQ file.")]
struct Args {
    /// Input FASTQ file.
    #[arg(short = 'i', long = "infile")]
    in_file: PathBuf,

    /// Prefix for generated matrix and heatmap files.
    #[arg(short = 'o', long = "outfile")]
    out_file: PathBuf,

    /// Maximum number of sequences to process.
    #[arg(short = 'n', long = "nseq-limit")]
    max_seq: Option<i64>,

    #[arg(long = "cell-width", default_value_t = 12)]
    cell_width: i64,
    #[arg(long = "cell-height", default_value_t = 10)]
    cell_height: i64,
    #[arg(long = "font-size", default_value_t = 6)]
    font_size: i64,
    #[arg(long = "angle", default_value_t = 45)]
    col_angle: i64,
    #[arg(long = "text-color", default_value = "black")]
    text_color: String,
    #[arg(long = "file-type", default_value = "pdf")]
    file_type: String,
    #[arg(long = "no-num")]
    no_num: bool,
    #[arg(long = "skip-heatmap")]
    skip_heatmap: bool,

    #[arg(long = "rscript", default_value = "Rscript")]
    rscript: String,
    #[arg(long = "install-r-deps")]
    install_r_deps: bool,
    #[arg(long = "cran-mirror", default_value = "https://cloud.r-project.org")]
    cran_mirror: String,

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

fn r_package_available(rscript: &str, package: &str) -> bool {
    Command::new(rscript)
        .arg("-e")
        .arg(format!("quit(status=ifelse(requireNamespace('{package}', quietly=TRUE), 0, 1))"))
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn ensure_r_dependencies(rscript: &str, install_missing: bool, cran_mirror: &str) -> std::io::Result<()> {
    let rscript_path = rseqc_commands::exec_resolve::which(rscript).ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("Rscript executable not found: {rscript}"),
        )
    })?;
    let rscript = rscript_path.to_string_lossy();
    if r_package_available(&rscript, "pheatmap") {
        return Ok(());
    }
    if !install_missing {
        return Err(std::io::Error::other(
            "R package 'pheatmap' is not installed. Install it with \"Rscript -e \\\"install.packages('pheatmap', repos='https://cloud.r-project.org')\\\"\", rerun with --install-r-deps, or use --skip-heatmap.",
        ));
    }
    eprintln!("Installing R package pheatmap from {cran_mirror}");
    let status = Command::new(rscript.as_ref()).arg("-e").arg(format!("install.packages('pheatmap', repos='{cran_mirror}')")).status()?;
    if !status.success() {
        return Err(std::io::Error::other("R dependency installation failed"));
    }
    if !r_package_available(&rscript, "pheatmap") {
        return Err(std::io::Error::other("R package 'pheatmap' is still unavailable after installation"));
    }
    Ok(())
}

fn run(args: &Args) -> std::io::Result<()> {
    if let Some(n) = args.max_seq {
        if n <= 0 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--nseq-limit must be greater than zero"));
        }
    }
    if args.cell_width <= 0 || args.cell_height <= 0 || args.font_size <= 0 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--cell-width, --cell-height, and --font-size must be greater than zero"));
    }

    let quals = fastq_qual_strings(rseqc_formats::open_text_input(&args.in_file)?)?;
    let dat = qual2count_mat(&quals, args.max_seq);
    let (count_csv, percent_csv) = render_quality_matrices(&dat)?;

    let prefix = args.out_file.to_string_lossy().into_owned();
    let count_path = format!("{prefix}.qual_count.csv");
    let percent_path = format!("{prefix}.qual_percent.csv");
    File::create(&count_path)?.write_all(count_csv.as_bytes())?;
    File::create(&percent_path)?.write_all(percent_csv.as_bytes())?;

    eprintln!("Created {count_path}");
    eprintln!("Created {percent_path}");

    if !args.skip_heatmap {
        ensure_r_dependencies(&args.rscript, args.install_r_deps, &args.cran_mirror)?;

        let heatmap_prefix = format!("{prefix}.qual_heatmap");
        let heatmap_path = format!("{heatmap_prefix}.{}", args.file_type);
        if std::path::Path::new(&heatmap_path).exists() {
            std::fs::remove_file(&heatmap_path)?;
        }

        // Matches upstream: log2_scale is NOT passed here (defaults to
        // False), unlike sc_editMatrix.py's heatmaps.
        let script = render_heatmap_r_script(&percent_path, &heatmap_prefix, &args.file_type, args.cell_width, args.cell_height, args.col_angle, args.font_size, &args.text_color, args.no_num, false);
        let r_path = format!("{heatmap_prefix}.r");
        File::create(&r_path)?.write_all(script.as_bytes())?;

        // Matches upstream's heatmap.make_heatmap: hardcoded "Rscript",
        // not args.rscript.
        let status = Command::new("Rscript").arg(&r_path).status();
        if !matches!(status, Ok(s) if s.success()) {
            eprintln!("Failed to run Rscript file \"{r_path}\"");
        }

        if !std::path::Path::new(&heatmap_path).is_file() {
            return Err(std::io::Error::other(format!("heatmap generation failed; expected output was not created: {heatmap_path}")));
        }
        eprintln!("Created {heatmap_path}");
    }

    eprintln!("Done.");
    Ok(())
}
