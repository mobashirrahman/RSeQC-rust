//! Port of `sc_seqLogo.py`'s count-matrix computation. Ports
//! `qcmodule.fastq.fasta_iter`/`seq2countMat`.
//!
//! **DIV-0005 applies, and applies differently here than anywhere else
//! in this port**: every other "plot" command deferred under DIV-0005
//! still ships a complete, real artifact -- literal R-script TEXT that
//! would produce the plot if `Rscript` were run against it. This
//! command has no such escape hatch: upstream's `qcmodule.fastq.
//! make_logo` calls the Python `logomaker` plotting library directly
//! (matplotlib-based), not R, so there is no intermediate text artifact
//! to generate in its place. The `.count_matrix.csv` output (this
//! module's whole scope) is fully computed and correct; the
//! `.logo.<format>` image is NOT produced by this port and requires
//! real native rendering (`crates/render`, not started) to implement.
//! Disclosed as DIV-0016 rather than silently skipped.
//!
//! Compressed (`.gz`/`.Z`/`.z`/`.bz`/`.bz2`/`.bzip2`) input is supported
//! at the CLI layer via `rseqc_formats::open_text_input`, same as
//! `sc_seqQual.py` -- independent of the logo-rendering gap above.
//!
//! **Preserves the `pandas.DataFrame.from_dict` column/row ORDER
//! quirk**, verified against a real `pandas` run in `oracle/venv`: with
//! no explicit `sort_index()` call anywhere in this code path (unlike
//! `sc_seqQual.py`'s quality matrices, which ARE explicitly sorted),
//! the count matrix's row order (positions) is natural ascending
//! (0,1,2,... -- the outer dict's insertion order), but its COLUMN
//! order (nucleotide bases) is neither alphabetical nor a fixed
//! ACGT(N) order -- it's the order each base is FIRST ENCOUNTERED while
//! scanning positions 0,1,2,... in turn (and, within a position, in
//! that position's own first-seen order across all input sequences).
//! Replicated with an explicit per-position insertion-order list, not
//! a `HashMap`/`BTreeMap` (whose iteration order would not match).
use std::io::{self, BufRead};

use crate::python_fmt::python_str_float;

/// Ports `fastq.fasta_iter`: yields every non-blank, non-`>`-prefixed
/// line as its OWN sequence (does NOT concatenate multi-line FASTA
/// records -- upstream assumes one line per record, matching this
/// tool's short-barcode use case).
pub fn fasta_iter(reader: impl BufRead) -> io::Result<Vec<String>> {
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let l = line.trim();
        if l.is_empty() || l.starts_with('>') {
            continue;
        }
        out.push(l.to_string());
    }
    Ok(out)
}

/// Ports `fastq.fastq_iter(mode='seq')`: yields the 2nd line of each
/// 4-line FASTQ record (the sequence line).
pub fn fastq_seq_strings(reader: impl BufRead) -> io::Result<Vec<String>> {
    let mut out = Vec::new();
    let mut count = 0u8;
    let mut seq = String::new();
    for line in reader.lines() {
        let line = line?;
        let l = line.trim().replace('\r', "");
        if l.is_empty() {
            continue;
        }
        count += 1;
        if count == 2 {
            seq = l;
        }
        if count == 4 {
            out.push(std::mem::take(&mut seq));
            count = 0;
        }
    }
    Ok(out)
}

/// One position's per-base counts, preserving FIRST-SEEN base order
/// (see module docs -- required to replicate pandas' column ordering).
#[derive(Debug, Default, Clone)]
struct PositionCounts {
    order: Vec<char>,
    counts: std::collections::HashMap<char, i64>,
}

impl PositionCounts {
    fn add(&mut self, base: char) {
        if !self.counts.contains_key(&base) {
            self.order.push(base);
        }
        *self.counts.entry(base).or_insert(0) += 1;
    }
}

/// Ports `seq2countMat`: per-position nucleotide histograms. `limit`
/// stops after that many sequences; `exclude_n` skips any sequence
/// containing the literal character `'N'` entirely (not just at that
/// position).
fn seq2count_mat(seqs: &[String], limit: Option<i64>, exclude_n: bool) -> Vec<PositionCounts> {
    let mut mat: Vec<PositionCounts> = Vec::new();
    let mut count = 0i64;

    for s in seqs {
        count += 1;
        if exclude_n && s.contains('N') {
            continue;
        }
        for (i, base) in s.chars().enumerate() {
            if mat.len() <= i {
                mat.resize_with(i + 1, PositionCounts::default);
            }
            mat[i].add(base);
        }
        if let Some(l) = limit {
            if count >= l {
                break;
            }
        }
    }

    mat
}

pub struct CountMatrix {
    /// Column (base) order, first-seen while scanning positions ascending.
    pub bases: Vec<char>,
    /// One row per position (already in position order 0..len).
    pub rows: Vec<Vec<i64>>,
}

/// Computes the full count matrix and its column order from raw input
/// sequences. This is the whole `.count_matrix.csv` deliverable (see
/// module docs for why the `.logo.<format>` image is out of scope).
pub fn compute_count_matrix(seqs: &[String], limit: Option<i64>, exclude_n: bool) -> io::Result<CountMatrix> {
    let mat = seq2count_mat(seqs, limit, exclude_n);
    if mat.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "no usable sequences were found; check the input format, sequence content, and --exclude-N setting"));
    }

    let mut bases: Vec<char> = Vec::new();
    for pos in &mat {
        for &b in &pos.order {
            if !bases.contains(&b) {
                bases.push(b);
            }
        }
    }

    let rows: Vec<Vec<i64>> = mat.iter().map(|pos| bases.iter().map(|b| pos.counts.get(b).copied().unwrap_or(0)).collect()).collect();

    Ok(CountMatrix { bases, rows })
}

/// Renders the `.count_matrix.csv` text: `Index,<base>,<base>,...`
/// header, then one `<position>,<count>,...` row per position.
///
/// Cell dtype is decided GLOBALLY (one flag for the whole matrix), the
/// same rule as `sc_seqqual::render_quality_matrices` and for the same
/// reason: `seq2countMat` (`fastq.py`) calls `pandas.DataFrame.T` before
/// returning, and transposing a DataFrame with heterogeneous per-column
/// dtypes forces pandas to upcast EVERY column to a common dtype
/// (float64) -- confirmed via a live pandas probe. Only when every
/// position's observed-base set is a subset of the union (no position
/// is sparse relative to some other position) does the whole matrix
/// stay integer; this is DIFFERENT from `sc_editmatrix::
/// render_edit_matrix_csv`, whose pipeline never transposes and so
/// decides dtype per COLUMN instead.
pub fn render_count_matrix_csv(matrix: &CountMatrix) -> String {
    let fully_dense = matrix.rows.iter().all(|row| row.iter().all(|&c| c > 0));

    let mut out = String::from("Index");
    for b in &matrix.bases {
        out.push(',');
        out.push(*b);
    }
    out.push('\n');
    for (pos, row) in matrix.rows.iter().enumerate() {
        out.push_str(&pos.to_string());
        for &c in row {
            out.push(',');
            if fully_dense {
                out.push_str(&c.to_string());
            } else {
                out.push_str(&python_str_float(c as f64));
            }
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn fasta_iter_yields_each_non_header_line_separately() {
        let text = ">seq1\nACGT\n>seq2\n\nACGA\n";
        let seqs = fasta_iter(Cursor::new(text)).unwrap();
        assert_eq!(seqs, vec!["ACGT".to_string(), "ACGA".to_string()]);
    }

    #[test]
    fn fastq_seq_strings_extracts_second_line() {
        let text = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n";
        let seqs = fastq_seq_strings(Cursor::new(text)).unwrap();
        assert_eq!(seqs, vec!["ACGT".to_string(), "TTTT".to_string()]);
    }

    #[test]
    fn compute_count_matrix_matches_real_pandas_column_order() {
        // Cross-checked byte-for-byte against a real pandas
        // DataFrame.from_dict + fillna(0) + .T pipeline run in
        // oracle/venv with this exact input (seqs: ACGT,ACGA,TCGT,ACGT).
        let seqs = vec!["ACGT".to_string(), "ACGA".to_string(), "TCGT".to_string(), "ACGT".to_string()];
        let matrix = compute_count_matrix(&seqs, None, false).unwrap();
        // Column order: A (pos0 first seq), T (pos0 third seq),
        // C (pos1), G (pos2) -- NOT alphabetical, NOT fixed ACGT.
        assert_eq!(matrix.bases, vec!['A', 'T', 'C', 'G']);
        assert_eq!(matrix.rows, vec![vec![3, 1, 0, 0], vec![0, 0, 4, 0], vec![0, 0, 0, 4], vec![1, 3, 0, 0],]);
    }

    #[test]
    fn render_count_matrix_csv_matches_real_pandas_output() {
        let seqs = vec!["ACGT".to_string(), "ACGA".to_string(), "TCGT".to_string(), "ACGT".to_string()];
        let matrix = compute_count_matrix(&seqs, None, false).unwrap();
        let csv = render_count_matrix_csv(&matrix);
        let expected = "Index,A,T,C,G\n0,3.0,1.0,0.0,0.0\n1,0.0,0.0,4.0,0.0\n2,0.0,0.0,0.0,4.0\n3,1.0,3.0,0.0,0.0\n";
        assert_eq!(csv, expected);
    }

    #[test]
    fn render_count_matrix_csv_stays_integer_when_fully_dense() {
        // Cross-checked against a real pandas run: when every position
        // observes every base in the union (no position is sparse
        // relative to another), the transpose inside seq2countMat
        // doesn't need to upcast any column, so the whole matrix stays
        // int64 -- "3", not "3.0". This is the shape that was
        // previously mis-rendered (the bug this test guards against).
        let matrix = CountMatrix { bases: vec!['A', 'T'], rows: vec![vec![3, 1], vec![2, 2], vec![1, 3]] };
        let csv = render_count_matrix_csv(&matrix);
        assert_eq!(csv, "Index,A,T\n0,3,1\n1,2,2\n2,1,3\n");
    }

    #[test]
    fn compute_count_matrix_exclude_n_skips_whole_sequence() {
        let seqs = vec!["ACGT".to_string(), "ACGN".to_string()];
        let with_n = compute_count_matrix(&seqs, None, false).unwrap();
        // bases order: A(pos0),C(pos1),G(pos2),T(pos3 from ACGT),N(pos3 from ACGN)
        assert_eq!(with_n.bases, vec!['A', 'C', 'G', 'T', 'N']);
        assert_eq!(with_n.rows[3], vec![0, 0, 0, 1, 1]); // T then N both counted at pos3

        let without_n = compute_count_matrix(&seqs, None, true).unwrap();
        // Only "ACGT" contributes -- pos3 has just T, no N column at all.
        assert_eq!(without_n.bases, vec!['A', 'C', 'G', 'T']);
        assert_eq!(without_n.rows[3], vec![0, 0, 0, 1]);
    }

    #[test]
    fn compute_count_matrix_errors_on_no_usable_sequences() {
        assert!(compute_count_matrix(&[], None, false).is_err());
        let all_excluded = vec!["ACGN".to_string()];
        assert!(compute_count_matrix(&all_excluded, None, true).is_err());
    }
}
