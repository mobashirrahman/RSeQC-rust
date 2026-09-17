//! Port of `sc_seqQual.py`: per-FASTQ-cycle Phred quality count/
//! percentage matrices plus the `pheatmap` heatmap R script (shared
//! generator with `sc_editMatrix.py`, see `crate::sc_editmatrix::
//! render_heatmap_r_script`). Ports `qcmodule.fastq.fastq_iter`/
//! `qual2countMat`.
//!
//! **Scope gap, disclosed rather than silently ignored**: upstream's
//! `qcmodule.ireader.reader` transparently decompresses `.gz`/`.bz2`
//! input by file extension. This port only reads plain-text FASTQ --
//! compressed input is not yet supported (no compression crate has been
//! added; every other "format" in this port turned out to need no new
//! dependency, so this is a genuine, disclosed gap rather than an
//! assumed non-issue).
//!
//! **Preserves the same `pandas.DataFrame.from_dict(...).fillna(0)`
//! all-float-cells quirk documented in `sc_editmatrix.rs`**: every count
//! matrix cell prints with a trailing `.0`. The percentage matrix's
//! cells are genuine fractions, rendered with full `str(float)`
//! precision (`python_str_float`), verified against a real `pandas` run
//! in `oracle/venv`.
use std::collections::{BTreeSet, HashMap};
use std::io::{self, BufRead};

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

    let header: String = (0..dat.len()).map(|i| format!(",{i}")).collect();

    let mut count_csv = format!("Index{header}\n");
    let mut percent_csv = format!("Index{header}\n");

    for &score in &scores {
        count_csv.push_str(&score.to_string());
        percent_csv.push_str(&score.to_string());
        for (i, pos_map) in dat.iter().enumerate() {
            let c = pos_map.get(&score).copied().unwrap_or(0);
            count_csv.push_str(&format!(",{}", python_str_float(c as f64)));
            percent_csv.push_str(&format!(",{}", python_str_float(c as f64 / totals[i] as f64)));
        }
        count_csv.push('\n');
        percent_csv.push('\n');
    }

    Ok((count_csv, percent_csv))
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
    fn render_quality_matrices_errors_on_empty_input() {
        assert!(render_quality_matrices(&[]).is_err());
    }
}
