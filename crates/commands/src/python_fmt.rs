//! Small shared text-formatting helpers that replicate specific Python
//! `str()`/`print()` behaviors not matched by Rust's defaults. Kept
//! separate from any one command since multiple ports need the same
//! fix -- see each command's module docs for where this matters.

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

/// Formats `value` the way Python's `f"{value:.12g}"` does (C's
/// `%.12g`): 12 significant digits, fixed notation when the decimal
/// exponent `X` satisfies `-4 <= X < 12`, otherwise scientific notation
/// with a lowercase `e`, a sign, and a minimum-2-digit (not
/// zero-truncated beyond that) exponent. Trailing zeros (and a bare
/// trailing decimal point) are stripped in both cases.
pub fn python_g12(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }

    const PRECISION: i32 = 12;
    let neg = value < 0.0;
    let abs = value.abs();

    // Use Rust's own scientific-notation rounding (via LowerExp) to
    // determine the correctly-rounded decimal exponent and mantissa,
    // rather than hand-rolling log10-based exponent math (which is
    // prone to off-by-one errors right at power-of-ten boundaries).
    let sci = format!("{:.*e}", (PRECISION - 1) as usize, abs);
    let e_pos = sci.find('e').expect("LowerExp always emits 'e'");
    let mantissa_str = &sci[..e_pos];
    let exp: i32 = sci[e_pos + 1..].parse().expect("LowerExp exponent is always a valid integer");

    let sign_str = if neg { "-" } else { "" };

    if exp >= -4 && exp < PRECISION {
        let decimals = (PRECISION - 1 - exp).max(0) as usize;
        let fixed = format!("{abs:.decimals$}");
        let trimmed = if fixed.contains('.') { fixed.trim_end_matches('0').trim_end_matches('.') } else { fixed.as_str() };
        format!("{sign_str}{trimmed}")
    } else {
        let mantissa = if mantissa_str.contains('.') { mantissa_str.trim_end_matches('0').trim_end_matches('.') } else { mantissa_str };
        let exp_sign = if exp < 0 { '-' } else { '+' };
        format!("{sign_str}{mantissa}e{exp_sign}{:02}", exp.abs())
    }
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

    #[test]
    fn python_g12_matches_cpython_format_spec() {
        // Cross-checked against real `python3 -c "print(f'{value:.12g}')"`.
        assert_eq!(python_g12(0.0), "0");
        assert_eq!(python_g12(1.0), "1");
        assert_eq!(python_g12(0.030029296875), "0.030029296875");
        assert_eq!(python_g12(0.3333333333333333), "0.333333333333");
        assert_eq!(python_g12(123456789.123456), "123456789.123");
        assert_eq!(python_g12(0.0001234567891234), "0.000123456789123");
        assert_eq!(python_g12(100.0), "100");
        assert_eq!(python_g12(0.5), "0.5");
        assert_eq!(python_g12(1.0 / 30000.0), "3.33333333333e-05");
        assert_eq!(python_g12(1e15), "1e+15");
        assert_eq!(python_g12(1e20), "1e+20");
        assert_eq!(python_g12(1e-15), "1e-15");
        assert_eq!(python_g12(1.23456 * 1e-9), "1.23456e-09");
        assert_eq!(python_g12(1.23456 * 1e11), "123456000000");
        assert_eq!(python_g12(1.23456 * 1e12), "1.23456e+12");
    }
}
