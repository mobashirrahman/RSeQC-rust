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

fn fixed_from_scientific(value: &str) -> String {
    let (sign, unsigned) = match value.strip_prefix('-') {
        Some(v) => ("-", v),
        None => ("", value.strip_prefix('+').unwrap_or(value)),
    };
    let Some((mantissa, exponent)) = unsigned.split_once(['e', 'E']) else {
        return value.to_string();
    };
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = format!("{whole}{fraction}");
    let decimal_position = whole.len() as i32 + exponent;
    let body = if decimal_position <= 0 {
        format!("0.{}{}", "0".repeat((-decimal_position) as usize), digits)
    } else if decimal_position >= digits.len() as i32 {
        format!("{}{}", digits, "0".repeat((decimal_position as usize).saturating_sub(digits.len())))
    } else {
        let position = decimal_position as usize;
        format!("{}.{}", &digits[..position], &digits[position..])
    };
    format!("{sign}{body}")
}

fn scientific_from_fixed(value: &str) -> String {
    let (sign, unsigned) = match value.strip_prefix('-') {
        Some(v) => ("-", v),
        None => ("", value.strip_prefix('+').unwrap_or(value)),
    };
    let (whole, fraction) = unsigned.split_once('.').unwrap_or((unsigned, ""));
    let digits = format!("{whole}{fraction}");
    let Some(first_nonzero) = digits.bytes().position(|b| b != b'0') else {
        return "0".to_string();
    };
    let mut significant = digits[first_nonzero..].to_string();
    while significant.ends_with('0') {
        significant.pop();
    }
    let exponent = whole.len() as i32 - first_nonzero as i32 - 1;
    let mut chars = significant.chars();
    let first = chars.next().expect("non-zero significant digits");
    let rest: String = chars.collect();
    let mantissa = if rest.is_empty() { first.to_string() } else { format!("{first}.{rest}") };
    let exponent_sign = if exponent < 0 { '-' } else { '+' };
    let exponent_abs = exponent.unsigned_abs();
    let exponent_digits = if exponent_abs < 10 { format!("0{exponent_abs}") } else { exponent_abs.to_string() };
    format!("{sign}{mantissa}e{exponent_sign}{exponent_digits}")
}

/// Adjusts the last digit of a plain decimal digit string (optional
/// leading `-`, digits, optional `.`) by `delta` (`+1`/`-1`),
/// propagating any carry/borrow leftward through the digits (skipping
/// over `.`). Returns `None` if the carry/borrow would need to grow the
/// digit count (e.g. `"9.99"` incremented needs a new leading digit) --
/// that case is astronomically rare here and callers just skip the
/// adjustment.
fn adjust_last_digit(raw: &str, delta: i32) -> Option<String> {
    let (sign, unsigned) = match raw.strip_prefix('-') {
        Some(u) => ("-", u),
        None => ("", raw),
    };
    let mut chars: Vec<char> = unsigned.chars().collect();
    let mut carry = delta;
    let mut idx = chars.len();
    while carry != 0 {
        if idx == 0 {
            return None;
        }
        idx -= 1;
        if chars[idx] == '.' {
            continue;
        }
        let d = (chars[idx] as u8 - b'0') as i32 + carry;
        if d >= 10 {
            chars[idx] = '0';
            carry = 1;
        } else if d < 0 {
            chars[idx] = '9';
            carry = -1;
        } else {
            chars[idx] = (b'0' + d as u8) as char;
            carry = 0;
        }
    }
    Some(format!("{sign}{}", chars.into_iter().collect::<String>()))
}

fn last_digit_is_even(s: &str) -> bool {
    s.chars().rev().find(|c| c.is_ascii_digit()).map(|c| (c as u8 - b'0').is_multiple_of(2)).unwrap_or(false)
}

/// Number of digits after the decimal point in a plain (non-scientific)
/// decimal string, or `0` if there is no `.`.
fn frac_digit_count(s: &str) -> usize {
    s.split_once('.').map(|(_, frac)| frac.len()).unwrap_or(0)
}

/// True only if `x`'s EXACT binary value (not merely "round-trips to a
/// decimal near here") is precisely the arithmetic midpoint between
/// `lo` and `lo`'s next representable decimal at its own digit count
/// (i.e. `lo` with a `5` appended one place further right) -- a genuine
/// tie, not just "this neighboring decimal also happens to round-trip".
/// `lo` must be the unsigned (no leading `-`), numerically smaller of
/// the two tied candidates.
///
/// Every `f64` has a finite exact decimal expansion (its value is
/// mantissa * 2^exponent, and 2^k always terminates in decimal), so
/// formatting `x` to enough fractional digits reproduces that exact
/// expansion verbatim, padded with trailing zeros -- no approximation,
/// this is exact by construction of Rust's fixed-precision float
/// formatting. `MARGIN` digits beyond the midpoint's own precision is
/// far more than any realistic normal-magnitude double from this
/// codebase's commands needs (bounded by mantissa bits plus a handful
/// of doublings/halvings from summation, well under 100); it would not
/// suffice for subnormals near the denormal boundary, which this
/// codebase never produces.
fn is_exact_midpoint(lo: &str, x: f64) -> bool {
    const MARGIN: usize = 400;
    let frac = frac_digit_count(lo);
    let precision = frac + 1 + MARGIN;

    let mut expected = lo.to_string();
    if frac == 0 {
        expected.push('.');
    }
    expected.push('5');
    expected.push_str(&"0".repeat(MARGIN));

    format!("{:.*}", precision, x.abs()) == expected
}

/// Rust's `Display` for `f64` (used as `raw` below) gives a
/// correctly-rounded, shortest round-tripping decimal, but Rust and
/// Python can disagree on which decimal to pick when a double's exact
/// binary value sits precisely halfway between the two nearest
/// candidates at that digit count (a genuine tie): Python's
/// `repr`/`str` (its own dtoa) rounds such ties to the EVEN last digit;
/// Rust's does not consistently. Concrete example that surfaced this:
/// `f64::from_bits(0x4056ae1480000000)` is exactly `90.720001220703125`
/// (equidistant between `...312` and `...313` at 16 significant
/// digits) -- Rust's `Display` prints `"90.72000122070313"`, Python's
/// `str()` prints `"90.72000122070312"` (the even one).
///
/// Merely checking "does the neighboring decimal ALSO round-trip to the
/// same bits" is NOT sufficient to detect a tie (tried first, caused
/// regressions): a double's round-trip acceptance interval is often
/// wider than one decimal step at the last digit, so a neighbor can
/// round-trip while still being strictly farther from `x`'s exact value
/// than `raw` -- e.g. `85.63000106811523` (exact value
/// `85.630001068115234375`) and `85.63000106811524` both round-trip to
/// the same bits, but the first is genuinely closer (Rust's `Display`
/// already gets this one right); swapping to the neighbor there would
/// have been a bug, not a fix. Genuine ties are instead confirmed via
/// `is_exact_midpoint`.
fn round_ties_to_even(raw: &str, x: f64) -> String {
    let bits = x.to_bits();
    for delta in [1i32, -1i32] {
        let Some(neighbor) = adjust_last_digit(raw, delta) else { continue };
        if neighbor.parse::<f64>().map(|v| v.to_bits()) != Ok(bits) {
            continue;
        }
        let lo = if delta > 0 { raw } else { neighbor.as_str() };
        let lo_unsigned = lo.strip_prefix('-').unwrap_or(lo);
        if !is_exact_midpoint(lo_unsigned, x) {
            continue;
        }
        return if last_digit_is_even(&neighbor) && !last_digit_is_even(raw) { neighbor } else { raw.to_string() };
    }
    raw.to_string()
}

/// Renders Python 3's `str(float)` representation, including its notation
/// switch at `1e-4`/`1e16` and two-digit signed exponents. Rust's default
/// Display supplies the shortest round-tripping decimal; the helpers above
/// only normalize its notation and Python's required trailing `.0`, and
/// `round_ties_to_even` corrects the rare exact-tie last digit.
pub fn python_str_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0".to_string() } else { "0.0".to_string() };
    }
    let raw = round_ties_to_even(&format!("{x}"), x);
    let scientific = x.abs() >= 1e16 || x.abs() < 1e-4;
    let normalized = if scientific {
        scientific_from_fixed(&fixed_from_scientific(&raw))
    } else if raw.contains('e') || raw.contains('E') {
        fixed_from_scientific(&raw)
    } else {
        raw
    };
    if scientific || normalized.contains('.') { normalized } else { format!("{normalized}.0") }
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

    if (-4..PRECISION).contains(&exp) {
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

/// numpy's `pairwise_sum_DOUBLE` (the float64 `np.add.reduce` kernel used
/// by `np.sum`/`np.mean`/`np.std`): 8-way unrolled blocks of up to 128
/// elements, recursively halved above that. Plain left-to-right summation
/// differs from it in the last bits for long inputs.
pub fn numpy_sum(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 8 {
        let mut res = -0.0;
        for &v in values {
            res += v;
        }
        res
    } else if n <= 128 {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&values[..8]);
        let mut i = 8;
        while i < n - (n % 8) {
            for j in 0..8 {
                r[j] += values[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += values[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        numpy_sum(&values[..n2]) + numpy_sum(&values[n2..])
    }
}

/// `np.mean` of a float64 array.
pub fn numpy_mean(values: &[f64]) -> f64 {
    numpy_sum(values) / values.len() as f64
}

/// `np.std(values, ddof=ddof)`: pairwise mean, squared deviations
/// (`x * x`), pairwise sum, divided by `n - ddof`, square root.
pub fn numpy_std(values: &[f64], ddof: usize) -> f64 {
    let mean = numpy_mean(values);
    let squares: Vec<f64> = values.iter().map(|&v| (v - mean) * (v - mean)).collect();
    (numpy_sum(&squares) / (values.len() - ddof) as f64).sqrt()
}

#[cfg(test)]
mod tests {
    #[test]
    fn python_str_float_breaks_genuine_ties_to_even() {
        // Regression for `genebody_coverage2_png`/`_sig1`/`_sig2`:
        // f64::from_bits(0x4056ae1480000000) is exactly
        // 90.720001220703125 -- precisely equidistant between
        // "90.72000122070312" and "90.72000122070313" at 16 significant
        // digits. Python's str() picks the even one ("...312"); Rust's
        // native Display picked "...313".
        let x = f64::from_bits(0x4056ae1480000000);
        assert_eq!(super::python_str_float(x), "90.72000122070312");

        // Second real case from the same sweep (`geneBodyCoverage.txt`
        // line 71, size1500_seed3): f64::from_bits(0x3ffb5c2800000000)
        // is exactly 1.70999908447265625 -- equidistant between
        // "1.7099990844726562" and "1.7099990844726563".
        let y = f64::from_bits(0x3ffb5c2800000000);
        assert_eq!(super::python_str_float(y), "1.7099990844726562");
    }

    #[test]
    fn python_str_float_does_not_flip_genuinely_closer_neighbor() {
        // Regression for a bug in an earlier version of this fix: merely
        // checking "does the decimal neighbor also round-trip" is not
        // enough to detect a tie. 85.63000106811523's exact binary value
        // is 85.630001068115234375 -- strictly closer to "...23" than to
        // "...24" even though BOTH decimal strings round-trip to the
        // same bits. Rust's Display already gets this one right and
        // must not be "corrected" to the farther neighbor.
        let x: f64 = "85.63000106811523".parse().unwrap();
        assert_eq!(
            x.to_bits(),
            "85.63000106811524".parse::<f64>().unwrap().to_bits(),
            "test premise: both decimal strings must round-trip to the same bits"
        );
        assert_eq!(super::python_str_float(x), "85.63000106811523");
    }

    #[test]
    fn python_str_float_ties_still_respect_notation_switch_and_sign() {
        // The tie fix runs before the sci/fixed notation branch and the
        // trailing-".0" branch in python_str_float; make sure ordinary
        // non-tied values (including negative and whole-number ones)
        // are completely unaffected.
        assert_eq!(super::python_str_float(-f64::from_bits(0x4056ae1480000000)), "-90.72000122070312");
        assert_eq!(super::python_str_float(5.0), "5.0");
        assert_eq!(super::python_str_float(1e-5), "1e-05");
    }

    #[test]
    fn numpy_sum_uses_pairwise_blocks() {
        // 0.1 summed 1000 times: numpy's pairwise kernel gives
        // 100.00000000000001 (np.sum(np.full(1000, 0.1))), naive left-to-right
        // summation 99.9999999999986.
        let v = vec![0.1f64; 1000];
        assert_eq!(numpy_sum(&v), 100.00000000000001);
        assert_ne!(v.iter().sum::<f64>(), 100.00000000000001);
        assert_eq!(numpy_std(&[1.0, 2.0, 3.0, 4.0], 0), 1.118033988749895);
    }

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
    fn python_float_notation_switch_matches_cpython() {
        assert_eq!(python_str_float(1e-5), "1e-05");
        assert_eq!(python_str_float(1e-4), "0.0001");
        assert_eq!(python_str_float(1e16), "1e+16");
        assert_eq!(python_str_float(1e15), "1000000000000000.0");
        assert_eq!(python_str_float(1.23456e-9), "1.23456e-09");
        assert_eq!(python_str_float(1.23456e12), "1234560000000.0");
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
