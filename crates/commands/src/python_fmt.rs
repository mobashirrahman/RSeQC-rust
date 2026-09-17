//! Small shared text-formatting helpers that replicate specific Python
//! `str()`/`print()` behaviors not matched by Rust's defaults. Kept
//! separate from any one command since multiple ports need the same
//! fix -- see each command's module docs for where this matters.

/// Renders a Python `str(float)`-equivalent string: Rust's default `f64`
/// `Display` already produces the shortest round-tripping decimal (like
/// Python's `repr`/`str` since 3.1) for non-integral values, but drops
/// the trailing `.0` for whole numbers that Python always keeps (e.g.
/// `5.0` prints as `"5"` in Rust, `"5.0"` in Python).
///
/// Known gap: Python's `str(float)` switches to scientific notation
/// (`"1e+16"`) for magnitudes >= 1e16 or < 1e-4; Rust's `{}` never does,
/// so this always renders the full decimal expansion instead. Not
/// replicated -- values in that range are not expected from any current
/// caller (TIN/FPKM/FPM scores stay well within ordinary magnitudes).
/// Python 3's `round()`: round-half-to-even, NOT round-half-away-from-zero.
pub fn python_round(x: f64) -> i64 {
    let floor = x.floor();
    let diff = x - floor;
    if diff < 0.5 {
        floor as i64
    } else if diff > 0.5 {
        floor as i64 + 1
    } else if (floor as i64) % 2 == 0 {
        floor as i64
    } else {
        floor as i64 + 1
    }
}

pub fn python_str_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }
    let s = format!("{x}");
    if s.contains('.') || s.contains('e') || s.contains('E') { s } else { format!("{s}.0") }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn whole_numbers_get_a_trailing_dot_zero() {
        assert_eq!(python_str_float(5.0), "5.0");
        assert_eq!(python_str_float(0.0), "0.0");
        assert_eq!(python_str_float(-0.0), "-0.0");
        assert_eq!(python_str_float(100.0), "100.0");
    }

    #[test]
    fn non_integral_values_pass_through() {
        assert_eq!(python_str_float(42.75), "42.75");
        assert_eq!(python_str_float(58.651499544403826), "58.651499544403826");
    }
}
