#!/usr/bin/env python3
"""Unit tests for the differential runner's failure-closed semantics.

These tests deliberately exercise the harness itself, rather than any RSeQC
algorithm.  They run without the oracle environment or a compiled Rust
workspace, so a broken comparator cannot hide behind missing dependencies.
"""
from __future__ import annotations

import contextlib
import io
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("run_diff.py")
SPEC = importlib.util.spec_from_file_location("run_diff", MODULE_PATH)
assert SPEC and SPEC.loader
run_diff = importlib.util.module_from_spec(SPEC)
sys.modules["run_diff"] = run_diff
SPEC.loader.exec_module(run_diff)


def _case(**kwargs):
    defaults = dict(
        name="unit",
        ensure_fixture=lambda: None,
        py_script="unused.py",
        rust_bin="unused",
        py_args=lambda _scratch: [],
        rust_args=lambda _scratch: [],
    )
    defaults.update(kwargs)
    return run_diff.Case(**defaults)


class RunnerComparatorTests(unittest.TestCase):
    def compare(self, case, py, rust):
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            with contextlib.redirect_stdout(io.StringIO()):
                return run_diff.compare_results(case, py, rust, Path(py_tmp), Path(rust_tmp))

    def test_scientific_notation_is_not_truncated(self):
        parsed = run_diff.extract_labeled_counts("small: 1e-3\nlarge: 1e+3\n")
        self.assertEqual(parsed, {"small": "1e-3", "large": "1e+3"})

    def test_duplicate_labels_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate report label"):
            run_diff.extract_labeled_counts("value: 1\nvalue: 2\n")

    def test_equal_nonzero_exit_codes_do_not_pass_positive_case(self):
        case = _case(stream_format="exact")
        failed_py = run_diff.RunResult(1, "same", "traceback")
        failed_rust = run_diff.RunResult(1, "same", "error")
        self.assertFalse(self.compare(case, failed_py, failed_rust))

    def test_empty_label_stream_does_not_pass(self):
        case = _case()
        success = run_diff.RunResult(0, "", "")
        self.assertFalse(self.compare(case, success, success))

    def test_required_labels_are_checked_on_both_sides(self):
        case = _case(required_labels=("Total records",))
        py = run_diff.RunResult(0, "Total records: 1\n", "")
        rust = run_diff.RunResult(0, "Other: 1\n", "")
        self.assertFalse(self.compare(case, py, rust))

    def test_decimal_and_scientific_values_compare_semantically(self):
        case = _case()
        py = run_diff.RunResult(0, "value: 0.001\n", "")
        rust = run_diff.RunResult(0, "value: 1e-3\n", "")
        self.assertTrue(self.compare(case, py, rust))

    def test_timeout_is_a_failed_run(self):
        case = _case(stream_format="exact")
        timed_out = run_diff.RunResult(None, "", "timeout", timed_out=True)
        self.assertFalse(self.compare(case, timed_out, timed_out))

    def test_run_enforces_timeout(self):
        result = run_diff.run(["python3", "-c", "import time; time.sleep(1)"], timeout_s=0.05)
        self.assertIsNone(result.exit_code)
        self.assertTrue(result.timed_out)


if __name__ == "__main__":
    unittest.main()
