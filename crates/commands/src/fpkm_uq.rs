//! Port of `FPKM-UQ.py`: run `htseq-count` as a subprocess, then compute
//! FPKM and upper-quartile-normalized FPKM (FPKM-UQ) from its output.
//! Ported from `oracle/upstream-src/scripts/FPKM-UQ.py` -- unlike every
//! other command in this port, the actual per-alignment counting is
//! delegated entirely to the external `htseq-count` tool; this module
//! only builds that command line and post-processes its tabular output.
//! No BAM/CIGAR/BED parsing lives here at all.
//!
//! **Preserves several deliberate upstream quirks, do not "fix"**:
//! - `gene_sizes` is populated for EVERY gene in the info file, but
//!   `uq_count`/`total_count` (the two normalization denominators) are
//!   computed from ONLY the protein-coding subset of counts -- yet the
//!   per-gene FPKM/FPKM-UQ loop below iterates over ALL genes in the
//!   htseq-count output (not just protein-coding ones), applying the
//!   protein-coding-derived normalizers to every gene's raw count.
//! - `fpkm_uq` is computed BEFORE `fpkm` in the same try block; if
//!   `gene_size == 0` (or, independently, `total_count == 0`), the
//!   SECOND statement (`fpkm`'s pure-integer division) is what actually
//!   raises Python's `ZeroDivisionError` -- `gene_size * uq_count`
//!   dividing count*1e9 involves a numpy float64 and so silently
//!   produces `inf`/`nan` rather than raising, even when `uq_count == 0`
//!   (a real, reachable low-coverage scenario, not just a defensive
//!   edge case). Both fields are overwritten to `"NA"` only when the
//!   SECOND statement's plain-Python-int division raises -- i.e. only
//!   when `gene_size == 0 || total_count == 0`. Replicated here as: if
//!   that condition holds, both cells are `"NA"`; otherwise `fpkm_uq` is
//!   computed via ordinary IEEE-754 float division (naturally yielding
//!   `inf`/`nan` on a zero `uq_count`, exactly like the numpy path).
//! - A gene present in the htseq-count output but absent from the info
//!   file is skipped with a stderr warning, not a hard failure.
#![allow(non_snake_case)]

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};

use crate::python_fmt::python_str_float;

/// Builds the original TCGA-compatible `htseq-count` command line.
/// `executable` should already be a resolved absolute path (see
/// `resolve_executable`); `bam_file`/`gtf_file` are canonicalized here,
/// matching upstream's `Path.resolve()` calls.
pub fn htseq_command(executable: &str, bam_file: &Path, gtf_file: &Path) -> io::Result<Vec<String>> {
    let bam_resolved = std::fs::canonicalize(bam_file)?;
    let gtf_resolved = std::fs::canonicalize(gtf_file)?;
    Ok(vec![
        executable.to_string(),
        "-f".to_string(),
        "bam".to_string(),
        "-r".to_string(),
        "pos".to_string(),
        "-s".to_string(),
        "no".to_string(),
        "-a".to_string(),
        "10".to_string(),
        "-t".to_string(),
        "exon".to_string(),
        "-i".to_string(),
        "gene_id".to_string(),
        "-m".to_string(),
        "intersection-nonempty".to_string(),
        bam_resolved.to_string_lossy().into_owned(),
        gtf_resolved.to_string_lossy().into_owned(),
    ])
}

/// Shell-quotes a single argument, replicating Python's `shlex.quote`
/// (POSIX-safe-character allowlist `[A-Za-z0-9_@%+=:,./-]`).
fn shlex_quote(s: &str) -> String {
    if !s.is_empty() && s.bytes().all(|b| b.is_ascii_alphanumeric() || b"_@%+=:,./-".contains(&b)) {
        return s.to_string();
    }
    format!("'{}'", s.replace('\'', r#"'"'"'"#))
}

/// Renders a command as a shell-readable string, for display only.
pub fn format_command(command: &[String]) -> String {
    command.iter().map(|s| shlex_quote(s)).collect::<Vec<_>>().join(" ")
}

/// Resolves an executable name to an absolute path by searching `PATH`
/// (or checking directly if `name` already contains a `/`), replicating
/// `shutil.which`'s behavior on Unix (executable-bit check).
pub fn resolve_executable(name: &str) -> io::Result<PathBuf> {
    use std::os::unix::fs::PermissionsExt;

    let is_executable = |p: &Path| -> bool { p.metadata().map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0).unwrap_or(false) };

    if name.contains('/') {
        let p = PathBuf::from(name);
        return if is_executable(&p) {
            Ok(p)
        } else {
            Err(io::Error::new(io::ErrorKind::NotFound, format!("cannot find htseq-count executable \"{name}\"")))
        };
    }

    if let Ok(path_var) = std::env::var("PATH") {
        for dir in std::env::split_paths(&path_var) {
            let candidate = dir.join(name);
            if is_executable(&candidate) {
                return Ok(candidate);
            }
        }
    }
    Err(io::Error::new(io::ErrorKind::NotFound, format!("cannot find htseq-count executable \"{name}\"")))
}

#[derive(Debug, Default, Clone)]
pub struct GeneInfo {
    pub sizes: HashMap<String, i64>,
    pub info: HashMap<String, String>,
    pub protein_coding: HashSet<String>,
}

/// Parses the gene information file: `gene_id` plus at least 10 more
/// columns, of which columns 1-5 (0-based) are display metadata, column
/// 6 is the gene type (`"protein_coding"` membership test), and column
/// 10 is the exon length. Hard-fails (unlike most BED12 parsers in this
/// port) on any line with fewer than 11 columns or a non-integer exon
/// length -- upstream's uncaught `ValueError`.
pub fn read_gene_information(reader: impl BufRead) -> io::Result<GeneInfo> {
    let mut out = GeneInfo::default();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let stripped = line.trim();
        if stripped.is_empty() || stripped.starts_with("gene_id") {
            continue;
        }

        let fields: Vec<&str> = stripped.split_whitespace().collect();
        if fields.len() < 11 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("line {}: expected at least 11 columns", line_number + 1)));
        }

        let gene_id = fields[0].to_string();
        let exon_length: i64 = fields[10]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("line {}: exon_length must be an integer", line_number + 1)))?;

        out.sizes.insert(gene_id.clone(), exon_length);
        out.info.insert(gene_id.clone(), fields[1..6].join("\t"));
        if fields[6] == "protein_coding" {
            out.protein_coding.insert(gene_id);
        }
    }

    Ok(out)
}

/// Parses non-summary records (lines not starting with `__`) from an
/// `htseq-count` output file, preserving row order.
pub fn read_htseq_counts(reader: impl BufRead) -> io::Result<Vec<(String, i64)>> {
    let mut out = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let stripped = line.trim();
        if stripped.is_empty() || stripped.starts_with("__") {
            continue;
        }

        let fields: Vec<&str> = stripped.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("line {}: expected two columns", line_number + 1)));
        }
        let count: i64 = fields[1]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("line {}: count must be an integer", line_number + 1)))?;
        out.push((fields[0].to_string(), count));
    }

    Ok(out)
}

/// Computes the q-th percentile of an already-sorted slice using numpy's
/// default `'linear'` interpolation method.
pub fn percentile_linear(sorted_values: &[f64], q: f64) -> f64 {
    let n = sorted_values.len();
    if n == 1 {
        return sorted_values[0];
    }
    let idx = (q / 100.0) * (n - 1) as f64;
    let lower = idx.floor() as usize;
    let upper = idx.ceil() as usize;
    let frac = idx - lower as f64;
    sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * frac
}

#[derive(Debug, Clone)]
pub struct FpkmUqRow {
    pub gene_id: String,
    pub info: String,
    pub count: i64,
    pub fpkm: String,
    pub fpkm_uq: String,
}

#[derive(Debug, Clone)]
pub struct FpkmUqSummary {
    pub protein_coding_counts_len: usize,
    pub uq_count: f64,
    pub total_count: i64,
}

/// Computes FPKM/FPKM-UQ for every gene present in `all_counts`. Returns
/// the rendered rows plus the summary values needed for the upstream
/// stderr progress messages. Errors (propagated `io::Error`) match
/// upstream's uncaught `ValueError("no protein-coding gene counts were
/// found")`.
pub fn calculate_fpkm(gene_info: &GeneInfo, all_counts: &[(String, i64)], log2_flag: bool) -> io::Result<(Vec<FpkmUqRow>, FpkmUqSummary)> {
    let protein_coding_counts: Vec<i64> = all_counts.iter().filter(|(id, _)| gene_info.protein_coding.contains(id)).map(|&(_, c)| c).collect();

    if protein_coding_counts.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "no protein-coding gene counts were found"));
    }

    let mut sorted_counts: Vec<f64> = protein_coding_counts.iter().map(|&c| c as f64).collect();
    sorted_counts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let uq_count = percentile_linear(&sorted_counts, 75.0);
    let total_count: i64 = protein_coding_counts.iter().sum();

    let mut rows = Vec::new();
    for (gene_id, count) in all_counts {
        let Some(&gene_size) = gene_info.sizes.get(gene_id) else {
            eprintln!("Warning: {gene_id} is absent from info file; skipped");
            continue;
        };
        let Some(info) = gene_info.info.get(gene_id) else { continue };

        let (fpkm_str, fpkm_uq_str) = if gene_size == 0 || total_count == 0 {
            ("NA".to_string(), "NA".to_string())
        } else {
            let numerator = *count as f64 * 1_000_000_000.0;
            let fpkm_uq_ratio = numerator / (gene_size as f64 * uq_count);
            let fpkm_ratio = numerator / (gene_size as f64 * total_count as f64);
            if log2_flag {
                (python_str_float((fpkm_ratio + 1.0).log2()), python_str_float((fpkm_uq_ratio + 1.0).log2()))
            } else {
                (python_str_float(fpkm_ratio), python_str_float(fpkm_uq_ratio))
            }
        };

        rows.push(FpkmUqRow { gene_id: gene_id.clone(), info: info.clone(), count: *count, fpkm: fpkm_str, fpkm_uq: fpkm_uq_str });
    }

    Ok((rows, FpkmUqSummary { protein_coding_counts_len: protein_coding_counts.len(), uq_count, total_count }))
}

/// Renders the `.FPKM-UQ.txt` table, including the log2-labeled header
/// variant.
pub fn render_fpkm_uq_table(rows: &[FpkmUqRow], log2_flag: bool) -> String {
    let mut out = String::new();
    let header = if log2_flag {
        ["gene_ID", "symbol", "chrom", "start", "end", "strand", "raw_count", "FPKM(log2(x+1))", "FPKM-UQ(log2(x+1))"]
    } else {
        ["gene_ID", "symbol", "chrom", "start", "end", "strand", "raw_count", "FPKM", "FPKM-UQ"]
    };
    out.push_str(&header.join("\t"));
    out.push('\n');
    for row in rows {
        out.push_str(&row.gene_id);
        out.push('\t');
        out.push_str(&row.info);
        out.push('\t');
        out.push_str(&row.count.to_string());
        out.push('\t');
        out.push_str(&row.fpkm);
        out.push('\t');
        out.push_str(&row.fpkm_uq);
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn percentile_linear_matches_numpy_default_method() {
        // Cross-checked against `numpy.percentile([0,0,0,3,5,10,20], 75)` == 7.5.
        let data = [0.0, 0.0, 0.0, 3.0, 5.0, 10.0, 20.0];
        assert_eq!(percentile_linear(&data, 75.0), 7.5);
        assert_eq!(percentile_linear(&[42.0], 75.0), 42.0);
    }

    #[test]
    fn shlex_quote_matches_python_shlex_quote() {
        // Cross-checked against `shlex.quote(...)` in python3.
        assert_eq!(shlex_quote("/tmp/plain.bam"), "/tmp/plain.bam");
        assert_eq!(shlex_quote("/tmp/some file.bam"), "'/tmp/some file.bam'");
        assert_eq!(shlex_quote("it's"), r#"'it'"'"'s'"#);
    }

    #[test]
    fn read_gene_information_parses_and_skips_header() {
        let text = "\
gene_id\tsymbol\tchrom\tstart\tend\tstrand\tgene_type\tx\ty\tz\texon_length
ENSG1\tGENE1\tchr1\t0\t100\t+\tprotein_coding\ta\tb\tc\t150
ENSG2\tGENE2\tchr1\t200\t300\t-\tlincRNA\ta\tb\tc\t80
";
        let info = read_gene_information(Cursor::new(text)).unwrap();
        assert_eq!(info.sizes["ENSG1"], 150);
        assert_eq!(info.sizes["ENSG2"], 80);
        assert_eq!(info.info["ENSG1"], "GENE1\tchr1\t0\t100\t+");
        assert!(info.protein_coding.contains("ENSG1"));
        assert!(!info.protein_coding.contains("ENSG2"));
    }

    #[test]
    fn read_htseq_counts_skips_blank_and_summary_lines() {
        let text = "\
ENSG1\t10
ENSG2\t0

__no_feature\t5
__ambiguous\t2
";
        let counts = read_htseq_counts(Cursor::new(text)).unwrap();
        assert_eq!(counts, vec![("ENSG1".to_string(), 10), ("ENSG2".to_string(), 0)]);
    }

    fn gene_info_fixture() -> GeneInfo {
        let mut info = GeneInfo::default();
        info.sizes.insert("g1".to_string(), 1000);
        info.sizes.insert("g2".to_string(), 0); // zero exon length -> forces NA
        info.info.insert("g1".to_string(), "SYM1\tchr1\t0\t1000\t+".to_string());
        info.info.insert("g2".to_string(), "SYM2\tchr1\t2000\t2000\t-".to_string());
        info.protein_coding.insert("g1".to_string());
        info.protein_coding.insert("g2".to_string());
        info
    }

    #[test]
    fn calculate_fpkm_basic_values() {
        let gene_info = gene_info_fixture();
        let counts = vec![("g1".to_string(), 100), ("g2".to_string(), 50)];
        let (rows, summary) = calculate_fpkm(&gene_info, &counts, false).unwrap();
        assert_eq!(summary.total_count, 150);
        // protein_coding_counts sorted = [50,100], n=2, idx=75/100*1=0.75,
        // result = 50 + (100-50)*0.75 = 87.5
        assert_eq!(summary.uq_count, 87.5);

        let g1 = rows.iter().find(|r| r.gene_id == "g1").unwrap();
        // fpkm = 100*1e9/(1000*150) = 1e11/150000 = 666666.6666666666
        // fpkm_uq = 100*1e9/(1000*87.5) = 1e11/87500 = 1142857.142857143
        assert_eq!(g1.fpkm, python_str_float(100.0 * 1e9 / (1000.0 * 150.0)));
        assert_eq!(g1.fpkm_uq, python_str_float(100.0 * 1e9 / (1000.0 * 87.5)));

        // g2 has gene_size==0 -> both NA regardless of counts.
        let g2 = rows.iter().find(|r| r.gene_id == "g2").unwrap();
        assert_eq!(g2.fpkm, "NA");
        assert_eq!(g2.fpkm_uq, "NA");
    }

    #[test]
    fn calculate_fpkm_zero_uq_count_yields_inf_not_na() {
        // 6 zero-count protein-coding genes + 1 nonzero -- 75th percentile
        // of [0,0,0,0,0,0,9] (numpy 'linear' method, n=7) is exactly 0,
        // while total_count is still nonzero. This is the real reachable
        // low-coverage scenario described in the module docs: fpkm_uq
        // divides by zero silently via numpy semantics and does NOT
        // become "NA" (only fpkm's pure-int division triggers that).
        let mut gene_info = GeneInfo::default();
        let genes = ["g1", "g2", "g3", "g4", "g5", "g6", "g7"];
        for g in genes {
            gene_info.sizes.insert(g.to_string(), 100);
            gene_info.info.insert(g.to_string(), "S\tchr1\t0\t100\t+".to_string());
            gene_info.protein_coding.insert(g.to_string());
        }
        let counts: Vec<(String, i64)> = genes.iter().map(|&g| (g.to_string(), 0)).collect();
        let mut counts = counts;
        counts.last_mut().unwrap().1 = 9; // g7 = 9, rest = 0
        let (rows, summary) = calculate_fpkm(&gene_info, &counts, false).unwrap();
        assert_eq!(summary.uq_count, 0.0);
        assert_eq!(summary.total_count, 9);
        let g7 = rows.iter().find(|r| r.gene_id == "g7").unwrap();
        assert_eq!(g7.fpkm_uq, "inf");
        assert_ne!(g7.fpkm, "NA");
    }

    #[test]
    fn render_fpkm_uq_table_header_and_row_shape() {
        let rows = vec![FpkmUqRow {
            gene_id: "g1".to_string(),
            info: "SYM1\tchr1\t0\t100\t+".to_string(),
            count: 10,
            fpkm: "1.5".to_string(),
            fpkm_uq: "2.5".to_string(),
        }];
        let table = render_fpkm_uq_table(&rows, false);
        assert_eq!(table, "gene_ID\tsymbol\tchrom\tstart\tend\tstrand\traw_count\tFPKM\tFPKM-UQ\ng1\tSYM1\tchr1\t0\t100\t+\t10\t1.5\t2.5\n");

        let log2_table = render_fpkm_uq_table(&rows, true);
        assert!(log2_table.starts_with("gene_ID\tsymbol\tchrom\tstart\tend\tstrand\traw_count\tFPKM(log2(x+1))\tFPKM-UQ(log2(x+1))\n"));
    }
}
