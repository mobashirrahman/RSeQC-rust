//! A narrow, purpose-built reimplementation of pandas' `DataFrame`
//! `repr()`/`to_string()` layout engine -- specifically, and ONLY, the
//! shape `sc_seqQual.py`'s two `logging.debug("...\n%s", quality_matrix)`
//! calls need under `--verbose`: a small numeric matrix with a named,
//! left-justified integer row index; sequential integer column labels;
//! right-justified numeric cells (either all-`int64` or all-`float64`,
//! decided by the caller); and column truncation with a `...` separator
//! when the rendered width exceeds the (piped-non-tty) terminal width.
//!
//! This is NOT a general DataFrame formatter. It was reverse-engineered
//! empirically by reading `pandas/io/formats/{format,string,printing}.py`
//! (this repo's `oracle/venv` pandas 3.0.5) and cross-checking every step
//! against real `pandas.DataFrame.__repr__()` output on the ACTUAL
//! `quality_matrix` object `sc_seqQual.py` builds (reconstructed via
//! `oracle/venv`'s `qcmodule.fastq` against a kept synthetic-sweep
//! fixture) -- this text is compared byte-for-byte by
//! `verification/synthetic_sweep.py`'s `sc_seqqual_n` case.
//!
//! Key, non-obvious facts confirmed via live pandas probes (see git log
//! for this file for the investigation):
//! - `repr(df)` uses `display.width=80`/`display.max_rows=60` as
//!   OPTION VALUES (not a real terminal query) when running as a
//!   non-interactive script (`__main__` has `__file__`); but the
//!   ACTUAL column-fit truncation is driven by a SEPARATE mechanism
//!   (`StringFormatter._fit_strcols_to_terminal_width`) that DOES query
//!   `shutil.get_terminal_size()` (COLUMNS/LINES env vars, else the
//!   (80, 24) fallback used whenever stdout isn't a real tty -- always
//!   true here, since both python and this port's binaries run under
//!   `subprocess.PIPE` in the sweep).
//! - Every numeric cell AND every column-label header gets a single
//!   leading "sign" space (`f"{value: .6f}"`/`f"{x: d}"` in Python --
//!   reserves a column for a minus sign; always a space here since
//!   quality counts/fractions are never negative).
//! - Float columns are formatted at `display.precision=6` decimals,
//!   then the COMMON number of trailing zero digits shared by every
//!   value in that column is trimmed (never trimming past 1 remaining
//!   digit) -- confirmed against a dozen synthetic columns, e.g.
//!   `[0.5, 0.25]` -> 2 decimals (`"0.50"`/`"0.25"`), `[1.0,2.0,0.0]` ->
//!   1 decimal (`"1.0"`/`"2.0"`/`"0.0"`), `[1/3, 2/3, 1.0]` -> full 6
//!   decimals (no common trailing zero).
//! - Column labels are first LEFT-justified to the widest label's digit
//!   width across ALL columns (e.g. `"0 "`.."49"` for 50 columns), THEN
//!   get the single leading sign-style space prepended -- this two-step
//!   order matters and was only found by tracing real `adj.justify`
//!   calls, not by reading the source alone.
//!
//! Vertical (row-count) truncation is implemented analogously to the
//! horizontal case (same source-derived logic, `display.max_rows=60`/
//! `display.min_rows=10`) but is NOT independently verified against real
//! pandas output on live data: every real workload this port's sweep
//! exercises stays at 38 distinct quality-score rows (well under 60), so
//! this path is a from-source best-effort safety net, not a
//! byte-verified one.

/// One numeric column: `label` is the plain column label (e.g. `"7"`,
/// NOT yet padded/spaced), `values` are the raw numeric cells for every
/// row, in the SAME order as [`DebugMatrix::row_labels`].
#[derive(Debug, Clone)]
pub struct DebugColumn {
    pub label: String,
    pub values: Vec<f64>,
}

/// The full matrix to render, already shaped exactly as it will appear
/// (row order, column order) -- this module does no sorting itself.
#[derive(Debug, Clone)]
pub struct DebugMatrix {
    pub index_name: String,
    /// Row labels in DISPLAY order (already sorted by the caller).
    pub row_labels: Vec<String>,
    pub columns: Vec<DebugColumn>,
    /// True if pandas would have kept this matrix as `int64` (no
    /// `NaN`-fill ever happened -- see `sc_seqqual::render_quality_matrices`'s
    /// `fully_dense` flag, which decides the SAME thing for `to_csv`).
    /// False means every cell renders as a float (the `display.precision`
    /// fixed-then-trimmed rule).
    pub all_integer: bool,
}

const MAX_ROWS_OPTION: usize = 60;
const MIN_ROWS_OPTION: usize = 10;

/// `shutil.get_terminal_size()`'s fallback logic: COLUMNS env var if a
/// valid positive integer, else the (80, 24) fallback pandas always
/// lands on when stdout isn't a real terminal -- which is always true
/// for a process launched with piped stdout/stderr, as both sides are
/// in `verification/synthetic_sweep.py`.
fn terminal_width() -> usize {
    std::env::var("COLUMNS").ok().and_then(|s| s.parse::<i64>().ok()).filter(|&n| n > 0).map(|n| n as usize).unwrap_or(80)
}

/// Python's `round()`: round half to even ("banker's rounding"). Only
/// ever called here on `n / 2` for integer `n`, so the input is always
/// either an integer or an exact `.5`.
fn round_half_to_even_div2(n: usize) -> usize {
    if n.is_multiple_of(2) {
        n / 2
    } else {
        let k = n / 2; // floor((n-1)/2 implied by integer division)
        if k.is_multiple_of(2) { k } else { k + 1 }
    }
}

/// Formats one float column exactly as pandas' `FloatArrayFormatter`
/// does when `fixed_width` is set (the default for `repr()`): fixed
/// `precision` decimals with a leading sign space, then trim the
/// longest COMMON run of trailing zero digits shared by every value in
/// the column (capped so at least one digit after the decimal point
/// always remains).
fn format_float_column(values: &[f64], precision: usize) -> Vec<String> {
    let full: Vec<String> = values.iter().map(|v| format_signed_fixed(*v, precision)).collect();

    let trailing_zeros = |s: &str| -> usize {
        let frac = s.rsplit('.').next().unwrap_or("");
        frac.chars().rev().take_while(|&c| c == '0').count()
    };
    let common = full.iter().map(|s| trailing_zeros(s)).min().unwrap_or(0);
    let trim = common.min(precision.saturating_sub(1));
    let new_precision = precision - trim;

    values.iter().map(|v| format_signed_fixed(*v, new_precision)).collect()
}

/// `f"{value: .{precision}f}"`: a leading space for non-negative values
/// (reserving the sign column), the literal `-` for negative ones.
fn format_signed_fixed(value: f64, precision: usize) -> String {
    if value.is_sign_negative() && value != 0.0 {
        format!("{value:.precision$}")
    } else {
        format!(" {value:.precision$}")
    }
}

/// `f"{x: d}"` for the all-`int64` case.
fn format_signed_int(value: i64) -> String {
    if value < 0 { format!("{value}") } else { format!(" {value}") }
}

/// Pads every string in `strings` to the widest one's length (LEFT- or
/// RIGHT-justified per `left`), matching `_make_fixed_width` /
/// `adj.justify` for our always-ASCII, no-`col_space`, no-`max_colwidth`
/// case.
fn justify_uniform(strings: &[String], left: bool) -> (Vec<String>, usize) {
    let width = strings.iter().map(|s| s.chars().count()).max().unwrap_or(0);
    let padded = strings
        .iter()
        .map(|s| {
            if left {
                format!("{s:<width$}")
            } else {
                format!("{s:>width$}")
            }
        })
        .collect();
    (padded, width)
}

/// pandas' `adjoin(space=1, *strcols)`: left-justifies (pads with
/// trailing spaces) every column except the last to `own_width + space`,
/// the last to just `own_width`, then joins row-wise.
fn adjoin(strcols: &[Vec<String>]) -> String {
    let n = strcols.len();
    let widths: Vec<usize> = strcols.iter().map(|c| c.iter().map(|s| s.chars().count()).max().unwrap_or(0)).collect();
    let lengths: Vec<usize> = widths.iter().enumerate().map(|(i, &w)| if i + 1 < n { w + 1 } else { w }).collect();
    let nrows = strcols.first().map(|c| c.len()).unwrap_or(0);
    let mut lines = Vec::with_capacity(nrows);
    for r in 0..nrows {
        let mut line = String::new();
        for (i, col) in strcols.iter().enumerate() {
            line.push_str(&format!("{:<width$}", col[r], width = lengths[i]));
        }
        lines.push(line);
    }
    lines.join("\n")
}

/// Builds the index column's `strcol`: `["", index_name, row_label_0, ...]`,
/// left-justified to the widest of `index_name`/any row label (matching
/// `_get_formatted_index`'s `justify="left"`).
fn build_index_strcol(index_name: &str, row_labels: &[String]) -> Vec<String> {
    let mut all: Vec<String> = Vec::with_capacity(row_labels.len() + 1);
    all.push(index_name.to_string());
    all.extend(row_labels.iter().cloned());
    let (padded, width) = justify_uniform(&all, true);
    let mut out = Vec::with_capacity(padded.len() + 1);
    out.push(" ".repeat(width)); // blank column-header row (adjoin pads "" itself, but this keeps widths honest for our own width table)
    out.extend(padded);
    out
}

/// Builds one data column's `strcol`: `[header_row, blank_row, value_0, ...]`.
/// `label_width` is the shared digit-width across ALL columns (numeric
/// column labels are left-padded to it before the leading sign space).
fn build_data_strcol(col: &DebugColumn, label_width: usize, all_integer: bool) -> Vec<String> {
    let header = format!(" {:<label_width$}", col.label);
    let value_strs = if all_integer {
        col.values.iter().map(|v| format_signed_int(*v as i64)).collect::<Vec<_>>()
    } else {
        format_float_column(&col.values, 6)
    };
    let header_colwidth = header.chars().count();
    let value_width = value_strs.iter().map(|s| s.chars().count()).max().unwrap_or(0);
    let width = header_colwidth.max(value_width);
    let padded_values: Vec<String> = value_strs.iter().map(|s| format!("{s:>width$}")).collect();
    let padded_header = format!("{header:>width$}");
    let blank = " ".repeat(width);
    let mut out = Vec::with_capacity(padded_values.len() + 2);
    out.push(padded_header);
    out.push(blank);
    out.extend(padded_values);
    out
}

/// Renders `matrix` exactly as `repr(quality_matrix)` would, including
/// the `\n\n[R rows x C columns]` suffix when truncated.
pub fn render_debug_matrix(matrix: &DebugMatrix) -> String {
    let n_rows = matrix.row_labels.len();
    let n_cols = matrix.columns.len();

    let label_width = matrix.columns.iter().map(|c| c.label.chars().count()).max().unwrap_or(0);
    let index_strcol = build_index_strcol(&matrix.index_name, &matrix.row_labels);
    let data_strcols: Vec<Vec<String>> = matrix.columns.iter().map(|c| build_data_strcol(c, label_width, matrix.all_integer)).collect();

    let mut strcols: Vec<Vec<String>> = Vec::with_capacity(n_cols + 1);
    strcols.push(index_strcol.clone());
    strcols.extend(data_strcols.iter().cloned());

    let full_text = adjoin(&strcols);
    let max_len = full_text.lines().map(|l| l.chars().count()).max().unwrap_or(0);
    let width = terminal_width();

    let mut truncated_h = false;
    let mut truncated_v = false;
    let mut final_strcols = strcols;

    if max_len > width {
        truncated_h = true;
        let mut lens: Vec<usize> = final_strcols.iter().map(|c| c.iter().map(|s| s.chars().count()).max().unwrap_or(0)).collect();
        let mut adj_dif = (max_len as i64 - width as i64) + 1;
        while adj_dif > 0 && lens.len() > 1 {
            let mid = round_half_to_even_div2(lens.len());
            let removed = lens.remove(mid);
            adj_dif -= removed as i64 + 1;
        }
        let n_cols_after = lens.len();
        let max_cols_fitted = (n_cols_after.saturating_sub(1)).max(2); // -1 for the (always-shown) index column
        let col_num = max_cols_fitted / 2;

        let dot_col: Vec<String> = std::iter::repeat_n(" ...".to_string(), index_strcol.len()).collect();
        let mut rebuilt: Vec<Vec<String>> = Vec::with_capacity(2 * col_num + 2);
        rebuilt.push(index_strcol);
        rebuilt.extend(data_strcols[..col_num].iter().cloned());
        rebuilt.push(dot_col);
        rebuilt.extend(data_strcols[data_strcols.len() - col_num..].iter().cloned());
        final_strcols = rebuilt;
    }

    // Vertical (row) truncation: same source-derived logic as above, but
    // not independently byte-verified against real pandas (see module
    // doc comment) since no real workload this port targets ever
    // produces more than 38 distinct quality-score rows.
    let max_rows_fitted = if n_rows > MAX_ROWS_OPTION { MIN_ROWS_OPTION.min(MAX_ROWS_OPTION) } else { MAX_ROWS_OPTION };
    if n_rows > max_rows_fitted {
        truncated_v = true;
        let row_num = max_rows_fitted / 2;
        let header_rows = 2; // column-header row + index-name row
        let mut rebuilt: Vec<Vec<String>> = Vec::with_capacity(final_strcols.len());
        for col in &final_strcols {
            let width = col.iter().map(|s| s.chars().count()).max().unwrap_or(0);
            let mut new_col = Vec::with_capacity(header_rows + 2 * row_num + 1);
            new_col.extend(col[..header_rows].iter().cloned());
            new_col.extend(col[header_rows..header_rows + row_num].iter().cloned());
            new_col.push(format!("{:>width$}", ".."));
            new_col.extend(col[col.len() - row_num..].iter().cloned());
            rebuilt.push(new_col);
        }
        final_strcols = rebuilt;
    }

    let mut text = adjoin(&final_strcols);
    if truncated_h || truncated_v {
        text.push_str(&format!("\n\n[{n_rows} rows x {n_cols} columns]"));
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix_from_grid(row_labels: &[i64], values: &[Vec<f64>], all_integer: bool) -> DebugMatrix {
        let n_cols = values[0].len();
        let columns = (0..n_cols)
            .map(|j| DebugColumn {
                label: j.to_string(),
                values: values.iter().map(|row| row[j]).collect(),
            })
            .collect();
        DebugMatrix {
            index_name: "pos".to_string(),
            row_labels: row_labels.iter().map(|v| v.to_string()).collect(),
            columns,
            all_integer,
        }
    }

    #[test]
    fn small_untruncated_matrix_matches_real_pandas() {
        // Cross-checked byte-for-byte against a real
        // `repr(pd.DataFrame({0:[1.5,2.25,3.0], 1:..., 2:...}))` run in
        // oracle/venv (3 columns, 3 rows, non-integral values -- no
        // truncation needed): `repr(repr(df))` ==
        // `'        0     1     2\\npos                  \\n2    1.50  1.50  1.50\\n1    2.25  2.25  2.25\\n0    3.00  3.00  3.00'`.
        let m = matrix_from_grid(&[2, 1, 0], &[vec![1.5, 1.5, 1.5], vec![2.25, 2.25, 2.25], vec![3.0, 3.0, 3.0]], false);
        let out = render_debug_matrix(&m);
        let expected = "        0     1     2\npos                  \n2    1.50  1.50  1.50\n1    2.25  2.25  2.25\n0    3.00  3.00  3.00";
        assert_eq!(out, expected);
    }

    #[test]
    fn all_integer_column_has_no_decimal_point() {
        // Cross-checked against a real `repr(pd.DataFrame({0:[1,3],
        // 1:[2,4]}))` (plain Python `int`s -> genuine `int64` dtype) run
        // in oracle/venv: `repr(repr(df))` ==
        // `'     0  1\\npos      \\n1    1  2\\n0    3  4'`.
        let m = matrix_from_grid(&[1, 0], &[vec![1.0, 2.0], vec![3.0, 4.0]], true);
        let out = render_debug_matrix(&m);
        let expected = "     0  1\npos      \n1    1  2\n0    3  4";
        assert_eq!(out, expected);
    }
}
