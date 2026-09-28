// Fuzz tests for Python numeric formatting functions
use proptest::prelude::*;
use rseqc_commands::python_fmt;

proptest! {
    #[test]
    fn fuzz_python_round(x in any::<f64>()) {
        // Just ensure python_round doesn't panic on any f64
        let _result = python_fmt::python_round(x);
    }

    #[test]
    fn fuzz_python_str_float(x in any::<f64>()) {
        // Just ensure python_str_float doesn't panic and returns a valid string
        let result = python_fmt::python_str_float(x);
        prop_assert!(!result.is_empty());
    }

    #[test]
    fn fuzz_python_g12(x in any::<f64>()) {
        // Just ensure python_g12 doesn't panic and returns a valid string
        let result = python_fmt::python_g12(x);
        prop_assert!(!result.is_empty());
    }

    // Test specific edge cases for python_round
    #[test]
    fn test_python_round_half_to_even(x in 0.0f64..=10.0) {
        let rounded = python_fmt::python_round(x);
        // Test that half-way cases round to even
        if x.fract() == 0.5 {
            // Rounded result should be even when x is halfway
            prop_assert_eq!(rounded % 2, 0, "half-way case should round to even");
        }
    }

    // Test python_str_float with various magnitude extremes
    #[test]
    fn test_python_str_float_magnitudes(
        exponent in -308i32..=308i32,
    ) {
        let value = 10.0f64.powi(exponent);
        let result = python_fmt::python_str_float(value);
        // Check that it doesn't crash and produces valid output
        prop_assert!(!result.is_empty());
        // For very small numbers, should use scientific notation (unless value is 0 or denormal)
        if exponent < -4 && value != 0.0 && value.is_normal() {
            prop_assert!(result.contains('e'), "small normal numbers should use scientific notation");
        }
    }

    // Test python_str_float preserves sign
    #[test]
    fn test_python_str_float_sign_preservation(
        x in 0.0001f64..=1000.0f64,
    ) {
        let result_pos = python_fmt::python_str_float(x);
        let result_neg = python_fmt::python_str_float(-x);

        if x != 0.0 {
            prop_assert!(!result_pos.starts_with('-'));
            prop_assert!(result_neg.starts_with('-'));
        }
    }

    // Test that NaN, Inf are handled correctly
    #[test]
    fn test_python_str_float_special_values(_dummy in 0..1) {
        let nan_str = python_fmt::python_str_float(f64::NAN);
        prop_assert_eq!(nan_str, "nan");

        let inf_str = python_fmt::python_str_float(f64::INFINITY);
        prop_assert_eq!(inf_str, "inf");

        let neg_inf_str = python_fmt::python_str_float(f64::NEG_INFINITY);
        prop_assert_eq!(neg_inf_str, "-inf");

        let pos_zero_str = python_fmt::python_str_float(0.0);
        prop_assert_eq!(pos_zero_str, "0.0");

        let neg_zero_str = python_fmt::python_str_float(-0.0);
        prop_assert_eq!(neg_zero_str, "-0.0");
    }

    // Test g12 format boundary
    #[test]
    fn test_python_g12_boundary(
        exponent in -10i32..=20i32,
    ) {
        let value = 10.0f64.powi(exponent);
        let result = python_fmt::python_g12(value);

        // g12 format uses fixed notation for -4 <= exp < 12
        if (-4..12).contains(&exponent) {
            prop_assert!(!result.contains('e') && !result.contains('E'),
                "fixed notation expected for exponent {}", exponent);
        } else if !(-4..12).contains(&exponent) {
            // Should use scientific notation outside the range
            // (for most cases; zero is special)
            if value != 0.0 {
                prop_assert!(result.contains('e') || result.contains('E'),
                    "scientific notation expected for exponent {}", exponent);
            }
        }
    }

    // Test very large and very small numbers don't overflow
    #[test]
    fn test_python_fmt_no_overflow(_dummy in 0..1) {
        let very_large = 1e308;
        let very_small = 1e-308;
        let subnormal = 1e-320;

        let _ = python_fmt::python_str_float(very_large);
        let _ = python_fmt::python_str_float(very_small);
        let _ = python_fmt::python_str_float(subnormal);

        let _ = python_fmt::python_g12(very_large);
        let _ = python_fmt::python_g12(very_small);
        let _ = python_fmt::python_g12(subnormal);
    }

    // Test common FPKM-like values
    #[test]
    fn test_python_fmt_fpkm_values(
        fpkm in 0.0f64..=1_000_000.0,
    ) {
        let str_result = python_fmt::python_str_float(fpkm);
        let g12_result = python_fmt::python_g12(fpkm);

        prop_assert!(!str_result.is_empty());
        prop_assert!(!g12_result.is_empty());
    }

    // Round-trip test for common values
    #[test]
    fn test_python_round_consistency(x in -1e10f64..=1e10f64) {
        let rounded = python_fmt::python_round(x);
        // Just ensure it doesn't panic and returns a reasonable value
        prop_assert!(rounded >= (x.floor() as i64) - 1);
        prop_assert!(rounded <= (x.ceil() as i64) + 1);
    }

    // Test with integers (should format as X.0)
    #[test]
    fn test_python_str_float_integers(x in 0i32..=1000i32) {
        let result = python_fmt::python_str_float(x as f64);
        prop_assert!(result.contains(".0"), "integers should have .0");
    }

    // Test zero-variance case for percentile calculations
    #[test]
    fn test_python_fmt_zeros(_dummy in 0..1) {
        let result = python_fmt::python_str_float(0.0);
        prop_assert_eq!(result, "0.0");

        let g12_result = python_fmt::python_g12(0.0);
        prop_assert_eq!(g12_result, "0");
    }
}
