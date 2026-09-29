//! Dispatch and flag parsing for `sc_editMatrix.py`. Binary name can't
//! contain '.' (see crates/cli/Cargo.toml); packaging (PORTING_PLAN
//! Step 10) adds the `.py`-suffixed PATH alias.
//!
//! No BAI parsing (project-wide policy); the BAM-index staleness
//! warning upstream prints is a courtesy check, replicated as a
//! best-effort mtime comparison.

use std::fs::File;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::Command;

use clap::Parser;
use rseqc_commands::pylog::log_line;
use rseqc_commands::sc_editmatrix::{BarcodeTagNames, barcode_edits, render_edit_matrix_csv, render_freq_tsv, render_heatmap_r_script};

#[derive(Parser)]
#[command(name = "sc_editMatrix.py", about = "Visualize error-correction edits in cellular barcodes and UMIs.")]
struct Args {
    /// Input BAM file.
    #[arg(short = 'i', long = "infile")]
    in_file: PathBuf,

    /// Prefix for generated matrix and heatmap files.
    #[arg(short = 'o', long = "outfile")]
    out_file: PathBuf,

    /// Maximum number of alignments to process.
    #[arg(long = "limit")]
    reads_num: Option<i64>,

    #[arg(long = "cr-tag", default_value = "CR")]
    cr_tag: String,
    #[arg(long = "cb-tag", default_value = "CB")]
    cb_tag: String,
    #[arg(long = "ur-tag", default_value = "UR")]
    ur_tag: String,
    #[arg(long = "ub-tag", default_value = "UB")]
    ub_tag: String,

    #[arg(long = "cell-width", default_value_t = 15)]
    cell_width: i64,
    #[arg(long = "cell-height", default_value_t = 10)]
    cell_height: i64,
    #[arg(long = "font-size", default_value_t = 8)]
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
            // Upstream: logging.error("%s", exc) -- the message alone, after
            // the (unreplicated, DIV-0019) timestamp/level prefix.
            eprintln!("{err}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn validate_tag(option_name: &str, tag: &str) -> std::io::Result<()> {
    if tag.chars().count() != 2 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("{option_name} must be exactly two characters")));
    }
    if tag.chars().any(|c| c.is_whitespace()) {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, format!("{option_name} cannot contain whitespace")));
    }
    Ok(())
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
    eprintln!("{}", log_line("INFO", &format!("Installing R package pheatmap from {cran_mirror}")));
    let status = Command::new(rscript.as_ref()).arg("-e").arg(format!("install.packages('pheatmap', repos='{cran_mirror}')")).status()?;
    if !status.success() {
        return Err(std::io::Error::other("R dependency installation failed"));
    }
    if !r_package_available(&rscript, "pheatmap") {
        return Err(std::io::Error::other("R package 'pheatmap' is still unavailable after installation"));
    }
    Ok(())
}

/// Ports `generate_heatmap` + `heatmap.make_heatmap`'s R-execution half.
/// The actual Rscript invocation is hardcoded to the literal `"Rscript"`
/// (not `args.rscript`), matching upstream's own inconsistency -- see
/// crates/commands/src/sc_editmatrix.rs module docs.
#[allow(clippy::too_many_arguments)]
fn generate_heatmap(matrix_file: &str, out_prefix: &str, file_type: &str, cell_width: i64, cell_height: i64, col_angle: i64, font_size: i64, text_color: &str, no_numbers: bool) -> std::io::Result<()> {
    let script = render_heatmap_r_script(matrix_file, out_prefix, file_type, cell_width, cell_height, col_angle, font_size, text_color, no_numbers, true);
    let r_path = format!("{out_prefix}.r");
    File::create(&r_path)?.write_all(script.as_bytes())?;

    let plot_path = format!("{out_prefix}.{file_type}");
    if std::path::Path::new(&plot_path).exists() {
        std::fs::remove_file(&plot_path)?;
    }

    let status = Command::new("Rscript").arg(&r_path).status();
    match status {
        Ok(s) if s.success() => {}
        _ => eprintln!("Failed to run Rscript file \"{r_path}\""),
    }

    if !std::path::Path::new(&plot_path).is_file() {
        return Err(std::io::Error::other(format!("heatmap generation failed; expected output was not created: {plot_path}")));
    }
    Ok(())
}

fn run(args: &Args) -> std::io::Result<()> {
    if let Some(n) = args.reads_num {
        if n <= 0 {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--limit must be greater than zero"));
        }
    }
    for (option_name, t) in [("--cr-tag", &args.cr_tag), ("--cb-tag", &args.cb_tag), ("--ur-tag", &args.ur_tag), ("--ub-tag", &args.ub_tag)] {
        validate_tag(option_name, t)?;
    }
    if args.cell_width <= 0 || args.cell_height <= 0 || args.font_size <= 0 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "--cell-width, --cell-height, and --font-size must be greater than zero"));
    }

    let tags = BarcodeTagNames { cr: args.cr_tag.clone(), cb: args.cb_tag.clone(), ur: args.ur_tag.clone(), ub: args.ub_tag.clone() };

    // Ports `barcode_edits`' own `logging.info` calls (scbam.py), which
    // this port's pure `barcode_edits` compute function doesn't print
    // itself -- reproduced here, in the same order, around/after the
    // equivalent call.
    eprintln!("{}", log_line("INFO", &format!("Reading BAM file \"{}\" ...", args.in_file.display())));
    let (mut reader, _header) = rseqc_formats::open_bam(&args.in_file)?;
    let stats = barcode_edits(reader.records(), &tags, args.reads_num)?;

    eprintln!("{}", log_line("INFO", &format!("Total alignments processed: {}", stats.total_alignments)));
    eprintln!("{}", log_line("INFO", &format!("Number of alignmenets with <cell barcode> kept AS IS: {}", stats.cb.same)));
    eprintln!("{}", log_line("INFO", &format!("Number of alignmenets with <cell barcode> edited: {}", stats.cb.diff)));
    eprintln!("{}", log_line("INFO", &format!("Number of alignmenets with <cell barcode> missing: {}", stats.cb.miss)));
    eprintln!("{}", log_line("INFO", &format!("Number of alignmenets with UMI kept AS IS: {}", stats.umi.same)));
    eprintln!("{}", log_line("INFO", &format!("Number of alignmenets with UMI edited: {}", stats.umi.diff)));
    eprintln!("{}", log_line("INFO", &format!("Number of alignmenets with UMI missing: {}", stats.umi.miss)));

    let prefix = args.out_file.to_string_lossy().into_owned();
    let cb_freq_path = format!("{prefix}.CB_freq.tsv");
    let umi_freq_path = format!("{prefix}.UMI_freq.tsv");
    eprintln!("{}", log_line("INFO", &format!("Writing cell barcode frequencies to \"{cb_freq_path}\"")));
    File::create(&cb_freq_path)?.write_all(render_freq_tsv(&stats.cb).as_bytes())?;
    eprintln!("{}", log_line("INFO", &format!("Writing UMI frequencies to \"{umi_freq_path}\"")));
    File::create(&umi_freq_path)?.write_all(render_freq_tsv(&stats.umi).as_bytes())?;

    let cb_matrix_path = format!("{prefix}.CB_edits_count.csv");
    let umi_matrix_path = format!("{prefix}.UMI_edits_count.csv");
    eprintln!("{}", log_line("INFO", &format!("Writing the nucleotide editing matrix (count) of cell barcode to \"{cb_matrix_path}\"")));
    File::create(&cb_matrix_path)?.write_all(render_edit_matrix_csv(&stats.cb.corrected_bases).as_bytes())?;
    eprintln!("{}", log_line("INFO", &format!("Writing the nucleotide editing matrix of molecular barcode (UMI) to \"{umi_matrix_path}\"")));
    File::create(&umi_matrix_path)?.write_all(render_edit_matrix_csv(&stats.umi.corrected_bases).as_bytes())?;

    eprintln!("{}", log_line("INFO", &format!("Created {cb_matrix_path}")));
    eprintln!("{}", log_line("INFO", &format!("Created {umi_matrix_path}")));

    if !args.skip_heatmap {
        ensure_r_dependencies(&args.rscript, args.install_r_deps, &args.cran_mirror)?;

        generate_heatmap(&cb_matrix_path, &format!("{prefix}.CB_edits_heatmap"), &args.file_type, args.cell_width, args.cell_height, args.col_angle, args.font_size, &args.text_color, args.no_num)?;
        generate_heatmap(&umi_matrix_path, &format!("{prefix}.UMI_edits_heatmap"), &args.file_type, args.cell_width, args.cell_height, args.col_angle, args.font_size, &args.text_color, args.no_num)?;
    }

    eprintln!("{}", log_line("INFO", "Done."));
    Ok(())
}
