//! Port of `sc_seqQual.py`: per-FASTQ-cycle Phred quality count/
//! percentage matrices plus the `pheatmap` heatmap R script (shared
//! generator with `sc_editMatrix.py`, see `crate::sc_editmatrix::
//! render_heatmap_r_script`). Ports `qcmodule.fastq.fastq_iter`/
//! `qual2countMat`.
//!
//! Compressed (`.gz`/`.Z`/`.z`/`.bz`/`.bz2`/`.bzip2`) input is supported
//! at the CLI layer via `rseqc_formats::open_text_input`, matching
//! upstream's `qcmodule.ireader.nopen` extension dispatch; this module's
//! own functions are unaffected (already generic over `impl BufRead`).
//!
//! **Preserves a `pandas.DataFrame.from_dict(...).fillna(0)` dtype
//! quirk, but NOT the same one as `sc_editmatrix.rs`**: unlike
//! `sc_editMatrix.py` (no transpose), this command's pipeline calls
//! `.T` on the matrix before `to_csv`, which forces pandas to upcast
//! EVERY column to a common dtype when the pre-transpose columns had
//! heterogeneous dtypes -- so the count matrix's float-vs-int cell
//! rendering is decided GLOBALLY (one flag for the whole matrix), not
//! per column. See `render_quality_matrices`'s doc comment for the full
//! reasoning, confirmed via live pandas probes. The percentage matrix's
//! cells are genuine fractions regardless, rendered with full
//! `str(float)` precision (`python_str_float`), verified against a real
//! `pandas` run in `oracle/venv`.
use std::collections::{BTreeSet, HashMap};
use std::io::{self, BufRead};

use crate::pandas_repr::{DebugColumn, DebugMatrix};
use crate::python_fmt::python_str_float;

/// Ports `fastq_iter(infile, mode='qual')`: yields the 4th line of each
/// 4-line FASTQ record (the quality string). Blank lines (after
/// trimming) are skipped and do NOT count toward the 4-line cycle,
/// matching upstream's `if len(l)==0: continue` before the line-count
/// increment.
pub fn fastq_qual_strings(reader: impl BufRead) -> io::Result<Vec<String>> {
    let mut out = Vec::new();
    let mut count = 0u8;
    for line in reader.lines() {
        let line = line?;
        let l = line.trim().replace('\r', "");
        if l.is_empty() {
            continue;
        }
        count += 1;
        if count == 4 {
            out.push(l);
            count = 0;
        }
    }
    Ok(out)
}

/// Ports `qual2countMat`'s counting half: per-position histogram of
/// Phred quality scores (`ord(char) - 33`). `limit` stops after that
/// many quality strings, matching upstream's `--nseq-limit`.
pub fn qual2count_mat(quals: &[String], limit: Option<i64>) -> Vec<HashMap<i64, i64>> {
    let mut dat: Vec<HashMap<i64, i64>> = Vec::new();
    let mut count = 0i64;

    for q in quals {
        count += 1;
        for (i, ch) in q.chars().enumerate() {
            let score = ch as i64 - 33;
            if dat.len() <= i {
                dat.resize_with(i + 1, HashMap::new);
            }
            *dat[i].entry(score).or_insert(0) += 1;
        }
        if let Some(l) = limit {
            if count >= l {
                break;
            }
        }
    }

    dat
}

/// Renders `(count_csv, percent_csv)` from the per-position histogram.
/// Rows are quality scores sorted DESCENDING (matching the CLI's
/// `sort_index(ascending=False)` after transposing); columns are
/// positions in natural ascending order. Ports the CLI script's own
/// matrix-shaping and division logic, including the hard error on any
/// zero-total column.
///
/// `count_csv`'s cell dtype is decided GLOBALLY, not per column, unlike
/// `sc_editmatrix::render_edit_matrix_csv`: `qual2countMat` (`fastq.py`)
/// and the CLI both call `pandas.DataFrame.T` on the matrix before
/// `to_csv`, and transposing a DataFrame with heterogeneous per-column
/// dtypes forces pandas to upcast EVERY column to a common dtype (here,
/// float64) -- confirmed via a live pandas probe comparing pre- and
/// post-transpose `.dtypes`. So: if every read-cycle's observed-score
/// set is a subset of some OTHER cycle's observed-score set (i.e. no
/// cycle is sparse relative to the union), the whole count matrix stays
/// integer; if even one cycle is missing a score present elsewhere, the
/// ENTIRE matrix (every cell, not just that column) renders with a
/// trailing `.0`. `percent_csv`'s cells are always genuine floats
/// regardless (a division result), so this doesn't apply there.
pub fn render_quality_matrices(dat: &[HashMap<i64, i64>]) -> io::Result<(String, String)> {
    if dat.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "no usable quality records were found in the FASTQ file"));
    }

    let mut scores: BTreeSet<i64> = BTreeSet::new();
    for pos_map in dat {
        scores.extend(pos_map.keys().copied());
    }
    let scores: Vec<i64> = scores.into_iter().rev().collect(); // descending

    let totals: Vec<i64> = dat.iter().map(|m| m.values().sum()).collect();
    let zero_cols: Vec<String> = totals.iter().enumerate().filter(|(_, &t)| t == 0).map(|(i, _)| i.to_string()).collect();
    if !zero_cols.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, format!("one or more read cycles contain zero total observations: {}", zero_cols.join(", "))));
    }

    let fully_dense = dat.iter().all(|pos_map| scores.iter().all(|s| pos_map.contains_key(s)));

    let header: String = (0..dat.len()).map(|i| format!(",{i}")).collect();

    let mut count_csv = format!("Index{header}\n");
    let mut percent_csv = format!("Index{header}\n");

    for &score in &scores {
        count_csv.push_str(&score.to_string());
        percent_csv.push_str(&score.to_string());
        for (i, pos_map) in dat.iter().enumerate() {
            let c = pos_map.get(&score).copied().unwrap_or(0);
            if fully_dense {
                count_csv.push_str(&format!(",{c}"));
            } else {
                count_csv.push_str(&format!(",{}", python_str_float(c as f64)));
            }
            percent_csv.push_str(&format!(",{}", python_str_float(c as f64 / totals[i] as f64)));
        }
        count_csv.push('\n');
        percent_csv.push('\n');
    }

    Ok((count_csv, percent_csv))
}

/// Builds the `(raw_counts, fraction)` `DebugMatrix` pair for the CLI's
/// `--verbose` `logging.debug("Sequence quality score matrix (...):\n%s",
/// quality_matrix)` lines -- same row/column shape and same
/// `fully_dense`-decided int-vs-float dtype as [`render_quality_matrices`]'s
/// CSV output (both come from the SAME upstream `quality_matrix`/
/// `quality_percent` objects), just rendered as a pandas `repr()` instead
/// of a CSV. Returns `None` if `dat` is empty (mirrors
/// `render_quality_matrices`'s own empty check; the CLI only calls this
/// after that check has already passed, so `None` shouldn't occur in
/// practice, but this avoids a second copy of the "empty" error).
pub fn build_debug_matrices(dat: &[HashMap<i64, i64>]) -> Option<(DebugMatrix, DebugMatrix)> {
    if dat.is_empty() {
        return None;
    }

    let mut scores: BTreeSet<i64> = BTreeSet::new();
    for pos_map in dat {
        scores.extend(pos_map.keys().copied());
    }
    let scores: Vec<i64> = scores.into_iter().rev().collect(); // descending

    let totals: Vec<i64> = dat.iter().map(|m| m.values().sum()).collect();
    let fully_dense = dat.iter().all(|pos_map| scores.iter().all(|s| pos_map.contains_key(s)));

    let row_labels: Vec<String> = scores.iter().map(|s| s.to_string()).collect();

    let mut count_columns = Vec::with_capacity(dat.len());
    let mut percent_columns = Vec::with_capacity(dat.len());
    for (i, pos_map) in dat.iter().enumerate() {
        let counts: Vec<f64> = scores.iter().map(|s| pos_map.get(s).copied().unwrap_or(0) as f64).collect();
        let fractions: Vec<f64> = counts.iter().map(|&c| c / totals[i] as f64).collect();
        count_columns.push(DebugColumn { label: i.to_string(), values: counts });
        percent_columns.push(DebugColumn { label: i.to_string(), values: fractions });
    }

    let count_matrix = DebugMatrix { index_name: "pos".to_string(), row_labels: row_labels.clone(), columns: count_columns, all_integer: fully_dense };
    let percent_matrix = DebugMatrix { index_name: "pos".to_string(), row_labels, columns: percent_columns, all_integer: false };

    Some((count_matrix, percent_matrix))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn fastq_qual_strings_extracts_fourth_line_and_skips_blanks() {
        let text = "@r1\nACGT\n+\nIIII\n\n@r2\nGGGG\n+\n##II\n";
        let quals = fastq_qual_strings(Cursor::new(text)).unwrap();
        assert_eq!(quals, vec!["IIII".to_string(), "##II".to_string()]);
    }

    #[test]
    fn qual2count_mat_builds_per_position_histogram() {
        let quals = vec!["III".to_string(), "II#".to_string(), "#II".to_string()];
        let dat = qual2count_mat(&quals, None);
        assert_eq!(dat.len(), 3);
        // position 0: 'I'(40) x2, '#'(2) x1
        assert_eq!(dat[0][&40], 2);
        assert_eq!(dat[0][&2], 1);
        // position 1: 'I'(40) x3
        assert_eq!(dat[1][&40], 3);
        // position 2: 'I'(40) x2, '#'(2) x1
        assert_eq!(dat[2][&40], 2);
        assert_eq!(dat[2][&2], 1);
    }

    #[test]
    fn render_quality_matrices_matches_real_pandas_output() {
        // Cross-checked byte-for-byte against a real pandas
        // qual2countMat + CLI division/to_csv pipeline run in
        // oracle/venv with this exact input (reads "III","II#","#II").
        let quals = vec!["III".to_string(), "II#".to_string(), "#II".to_string()];
        let dat = qual2count_mat(&quals, None);
        let (count_csv, percent_csv) = render_quality_matrices(&dat).unwrap();
        assert_eq!(count_csv, "Index,0,1,2\n40,2.0,3.0,2.0\n2,1.0,0.0,1.0\n");
        assert_eq!(percent_csv, "Index,0,1,2\n40,0.6666666666666666,1.0,0.6666666666666666\n2,0.3333333333333333,0.0,0.3333333333333333\n");
    }

    #[test]
    fn render_quality_matrices_stays_integer_when_fully_dense() {
        // Cross-checked against a real pandas run: when every position's
        // observed-score set is the SAME (no position is sparse relative
        // to the union), the transpose doesn't need to upcast any
        // column, so the whole count matrix stays int64 -- "2", not
        // "2.0". This is the shape that was previously mis-rendered
        // (the bug this test guards against).
        let dat = vec![
            HashMap::from([(40, 2), (2, 1)]),
            HashMap::from([(40, 5), (2, 3)]),
            HashMap::from([(40, 9), (2, 7)]),
        ];
        let (count_csv, _) = render_quality_matrices(&dat).unwrap();
        assert_eq!(count_csv, "Index,0,1,2\n40,2,5,9\n2,1,3,7\n");
    }

    #[test]
    fn render_quality_matrices_errors_on_empty_input() {
        assert!(render_quality_matrices(&[]).is_err());
    }
}
