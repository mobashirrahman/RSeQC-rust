//! Port of `overlay_bigwig.py`: combine two BigWig signal tracks
//! position-by-position and write the result as variableStep WIG.
//! Ported from `oracle/upstream-src/scripts/overlay_bigwig.py` and
//! `qcmodule/twoList.py`.
//!
//! **Preserves a real, confirmed-broken upstream action, do not "fix"**:
//! `twoList.Division` calls `(v1+1).__div__(v2+1)` -- `__div__` is a
//! Python 2 dunder method that numpy arrays do NOT implement under
//! Python 3 (only `__truediv__`/`__floordiv__` exist), so `-a Division`
//! ALWAYS raises an uncaught `AttributeError` in upstream, crashing the
//! whole program with a traceback (not caught by the CLI's own
//! `except (OSError, ValueError, RuntimeError)` clause, since
//! `AttributeError` isn't one of those). Confirmed by actually calling
//! `twoList.Division` in `oracle/venv`. This is dead, uncallable code in
//! upstream today -- NOT a working division a user could rely on.
//! Documented as DIV-0015 and replicated as an immediate hard error
//! (not a silent divide, and not a literal AttributeError message,
//! since that would misleadingly suggest a Rust-side bug).
use std::collections::{HashMap, HashSet};
use std::io;

use rseqc_formats::bigwig::BigWigReader;

use crate::normalize_bigwig::chromosome_chunks;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    Add,
    Subtract,
    Product,
    Division,
    Average,
    GeometricMean,
    Max,
    Min,
}

/// Ports `resolve_action`'s name lookup against `qcmodule.twoList`.
pub fn parse_action(name: &str) -> Option<Action> {
    match name {
        "Add" => Some(Action::Add),
        "Subtract" => Some(Action::Subtract),
        "Product" => Some(Action::Product),
        "Division" => Some(Action::Division),
        "Average" => Some(Action::Average),
        "geometricMean" => Some(Action::GeometricMean),
        "Max" => Some(Action::Max),
        "Min" => Some(Action::Min),
        _ => None,
    }
}

/// Ports the `qcmodule.twoList` operation for `action`. Inputs are
/// already NaN-replaced-with-0 by the caller (matching upstream's own
/// `np.nan_to_num` before calling into `twoList`), so no NaN handling
/// is needed here.
pub fn apply_action(action: Action, v1: &[f64], v2: &[f64]) -> io::Result<Vec<f64>> {
    if v1.len() != v2.len() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "the length of both arrays must be the same"));
    }
    let result = match action {
        Action::Add => v1.iter().zip(v2).map(|(&a, &b)| a + b).collect(),
        Action::Subtract => v1.iter().zip(v2).map(|(&a, &b)| a - b).collect(),
        Action::Product => v1.iter().zip(v2).map(|(&a, &b)| a * b).collect(),
        Action::Average => v1.iter().zip(v2).map(|(&a, &b)| (a + b) / 2.0).collect(),
        Action::GeometricMean => v1.iter().zip(v2).map(|(&a, &b)| (a * b).sqrt()).collect(),
        Action::Max => v1.iter().zip(v2).map(|(&a, &b)| a.max(b)).collect(),
        Action::Min => v1.iter().zip(v2).map(|(&a, &b)| a.min(b)).collect(),
        Action::Division => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "the 'Division' action is unusable in upstream RSeQC (qcmodule.twoList.Division calls a Python-2-only ndarray method that raises AttributeError under Python 3) -- see DIV-0015",
            ));
        }
    };
    Ok(result)
}

/// Merges two BigWigs' chromosome headers (bigwig1's chroms first, in
/// its own order, then any bigwig2-only chroms appended after, in its
/// own order), erroring if a shared chromosome has mismatched sizes.
/// Ports `combined_chromosome_sizes`.
pub fn combined_chromosome_sizes(first: &[(String, i64)], second: &[(String, i64)]) -> io::Result<Vec<(String, i64)>> {
    let mut combined: Vec<(String, i64)> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();

    for (chrom, size) in first {
        index.insert(chrom.clone(), combined.len());
        combined.push((chrom.clone(), *size));
    }
    for (chrom, size) in second {
        if let Some(&i) = index.get(chrom) {
            if combined[i].1 != *size {
                return Err(io::Error::new(io::ErrorKind::InvalidData, format!("chromosome size mismatch for {chrom:?}: {} versus {}", combined[i].1, size)));
            }
        } else {
            index.insert(chrom.clone(), combined.len());
            combined.push((chrom.clone(), *size));
        }
    }
    Ok(combined)
}

fn interval_has_signal(bw: &mut BigWigReader, chrom_set: &HashSet<String>, chrom: &str, start: i64, end: i64) -> io::Result<bool> {
    if !chrom_set.contains(chrom) {
        return Ok(false);
    }
    Ok(!bw.intervals(chrom, start as u32, end as u32)?.is_empty())
}

fn interval_values(bw: &mut BigWigReader, chrom_set: &HashSet<String>, chrom: &str, start: i64, end: i64) -> io::Result<Vec<f64>> {
    if !chrom_set.contains(chrom) {
        return Ok(vec![0.0; (end - start) as usize]);
    }
    let values = bw.values(chrom, start as u32, end as u32)?;
    Ok(values.iter().map(|&v| if v.is_nan() { 0.0 } else { v as f64 }).collect())
}

/// Ports Python's `f"{value:.2f}"` for the one case where it diverges
/// from Rust's default `{:.2}`: NaN. Python renders lowercase "nan";
/// Rust's Display renders "NaN". (Infinity is NOT special-cased: both
/// languages render "inf"/"-inf" lowercase already, confirmed via a
/// live `python3 -c` probe.) A NaN geometricMean result (sqrt of a
/// negative product) is genuinely reachable and printed by upstream,
/// not filtered out by the `value != 0.0` check (NaN never equals
/// anything, including 0.0).
fn format_value_2dp(value: f64) -> String {
    if value.is_nan() { "nan".to_string() } else { format!("{value:.2}") }
}

/// Runs the full chunked overlay pass, returning the rendered
/// variableStep WIG text. Ports `overlay_bigwigs`.
pub fn overlay_bigwigs(bw1: &mut BigWigReader, bw2: &mut BigWigReader, action: Action, chunk_size: i64) -> io::Result<String> {
    let sizes1: Vec<(String, i64)> = bw1.chroms().into_iter().map(|(n, l)| (n, l as i64)).collect();
    let sizes2: Vec<(String, i64)> = bw2.chroms().into_iter().map(|(n, l)| (n, l as i64)).collect();
    let chrom_sizes = combined_chromosome_sizes(&sizes1, &sizes2)?;

    let chrom_set1: HashSet<String> = sizes1.iter().map(|(n, _)| n.clone()).collect();
    let chrom_set2: HashSet<String> = sizes2.iter().map(|(n, _)| n.clone()).collect();

    let mut out = String::new();
    for (chrom, size) in &chrom_sizes {
        eprintln!("Processing {chrom} ...");
        out.push_str(&format!("variableStep chrom={chrom}\n"));

        for (start, end) in chromosome_chunks(*size, chunk_size) {
            let has1 = interval_has_signal(bw1, &chrom_set1, chrom, start, end)?;
            let has2 = interval_has_signal(bw2, &chrom_set2, chrom, start, end)?;
            if !has1 && !has2 {
                continue;
            }

            let v1 = interval_values(bw1, &chrom_set1, chrom, start, end)?;
            let v2 = interval_values(bw2, &chrom_set2, chrom, start, end)?;
            let result = apply_action(action, &v1, &v2)?;

            let mut coordinate = start;
            for value in result {
                coordinate += 1;
                // `value != 0.0` is true for NaN too (NaN never equals
                // anything, including 0.0), so a NaN geometricMean
                // result (sqrt of a negative product) DOES get printed
                // -- matching upstream.
                if value != 0.0 {
                    out.push_str(&format!("{coordinate}\t{}\n", format_value_2dp(value)));
                }
            }
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn apply_action_arithmetic_ops() {
        let v1 = vec![4.0, 9.0, 6.0];
        let v2 = vec![2.0, 3.0, 8.0];
        assert_eq!(apply_action(Action::Add, &v1, &v2).unwrap(), vec![6.0, 12.0, 14.0]);
        assert_eq!(apply_action(Action::Subtract, &v1, &v2).unwrap(), vec![2.0, 6.0, -2.0]);
        assert_eq!(apply_action(Action::Product, &v1, &v2).unwrap(), vec![8.0, 27.0, 48.0]);
        assert_eq!(apply_action(Action::Average, &v1, &v2).unwrap(), vec![3.0, 6.0, 7.0]);
        assert_eq!(apply_action(Action::GeometricMean, &v1, &v2).unwrap(), vec![8.0_f64.sqrt(), 27.0_f64.sqrt(), 48.0_f64.sqrt()]);
        assert_eq!(apply_action(Action::Max, &v1, &v2).unwrap(), vec![4.0, 9.0, 8.0]);
        assert_eq!(apply_action(Action::Min, &v1, &v2).unwrap(), vec![2.0, 3.0, 6.0]);
    }

    #[test]
    fn apply_action_division_is_a_hard_error() {
        // DIV-0015: Division is broken upstream, replicated as an error.
        assert!(apply_action(Action::Division, &[1.0], &[2.0]).is_err());
    }

    #[test]
    fn apply_action_rejects_mismatched_lengths() {
        assert!(apply_action(Action::Add, &[1.0, 2.0], &[1.0]).is_err());
    }

    #[test]
    fn parse_action_recognizes_all_choices_and_rejects_unknown() {
        for (name, expected) in [
            ("Add", Action::Add),
            ("Subtract", Action::Subtract),
            ("Product", Action::Product),
            ("Division", Action::Division),
            ("Average", Action::Average),
            ("geometricMean", Action::GeometricMean),
            ("Max", Action::Max),
            ("Min", Action::Min),
        ] {
            assert_eq!(parse_action(name), Some(expected));
        }
        assert_eq!(parse_action("Nope"), None);
    }

    #[test]
    fn combined_chromosome_sizes_orders_first_then_second_only_and_checks_mismatch() {
        let first = vec![("chr1".to_string(), 1000), ("chr2".to_string(), 2000)];
        let second = vec![("chr2".to_string(), 2000), ("chr3".to_string(), 3000)];
        let combined = combined_chromosome_sizes(&first, &second).unwrap();
        assert_eq!(combined, vec![("chr1".to_string(), 1000), ("chr2".to_string(), 2000), ("chr3".to_string(), 3000)]);

        let mismatched = vec![("chr2".to_string(), 9999), ("chr3".to_string(), 3000)];
        assert!(combined_chromosome_sizes(&first, &mismatched).is_err());
    }

    fn fixture_reader() -> BigWigReader {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../formats/tests/fixtures/pybigwig_test.bw");
        BigWigReader::open(&path).unwrap()
    }

    #[test]
    fn format_value_2dp_renders_nan_lowercase_matching_python() {
        // Cross-checked against a real `python3 -c "print(f'{float(\"nan\"):.2f}')"`.
        assert_eq!(format_value_2dp(f64::NAN), "nan");
        assert_eq!(format_value_2dp(1.5), "1.50");
        assert_eq!(format_value_2dp(f64::INFINITY), "inf");
        assert_eq!(format_value_2dp(f64::NEG_INFINITY), "-inf");
    }

    #[test]
    fn geometric_mean_of_a_negative_product_is_nan() {
        // (-1)*(4) = -4, sqrt(-4) = NaN -- a genuinely reachable
        // geometricMean result upstream prints (not filtered by the
        // `value != 0.0` check, since NaN never equals anything).
        let v1 = vec![-1.0];
        let v2 = vec![4.0];
        let result = apply_action(Action::GeometricMean, &v1, &v2).unwrap();
        assert!(result[0].is_nan());
    }

    #[test]
    fn overlay_bigwigs_add_self_doubles_signal() {
        let mut bw1 = fixture_reader();
        let mut bw2 = fixture_reader();
        let wig = overlay_bigwigs(&mut bw1, &mut bw2, Action::Add, 500_000).unwrap();
        assert!(wig.starts_with("variableStep chrom=1\n"));
        // fixture position 1 (1-based) = value[0] = 0.1 -> doubled = 0.20
        assert!(wig.contains("1\t0.20\n"));
        assert!(wig.contains("3\t0.60\n")); // 0.3 doubled
        assert!(wig.contains("variableStep chrom=10\n"));
    }
}
