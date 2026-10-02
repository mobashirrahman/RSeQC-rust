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

    def test_per_side_exit_overrides_pass_when_each_side_meets_its_own(self):
        # sc_seqlogo_basic contract: upstream exits 1 on its own pandas
        # bug, candidate succeeds. Diminished-exit must be explicit and
        # both sides must still meet their declared expectation.
        case = _case(compare_stream="none", py_expected_exit=1, rust_expected_exit=0)
        py = run_diff.RunResult(1, "", "upstream pandas crash")
        rust = run_diff.RunResult(0, "", "")
        self.assertTrue(self.compare(case, py, rust))

    def test_per_side_exit_override_rejects_side_failing_its_own_expectation(self):
        case = _case(compare_stream="none", py_expected_exit=1, rust_expected_exit=0)
        py = run_diff.RunResult(1, "", "ok")
        rust = run_diff.RunResult(1, "", "candidate crashed")
        self.assertFalse(self.compare(case, py, rust))

    def test_no_overrides_still_requires_identical_exits(self):
        case = _case(compare_stream="none")
        py = run_diff.RunResult(0, "", "")
        rust = run_diff.RunResult(1, "", "candidate crashed")
        self.assertFalse(self.compare(case, py, rust))

    def test_rust_expect_files_passes_for_nonempty_candidate_artifact(self):
        case = _case(compare_stream="none", rust_expect_files=("out.logo.pdf",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            (Path(rust_tmp) / "out.logo.pdf").write_bytes(b"%PDF-1.4")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertTrue(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))

    def test_rust_expect_files_fails_when_artifact_missing(self):
        case = _case(compare_stream="none", rust_expect_files=("out.logo.pdf",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))

    def test_nan_versus_nan_does_not_pass(self):
        """Two arms both collapsing to NaN is a defect, not agreement."""
        case = _case(required_labels=("Rate",))
        py = run_diff.RunResult(0, "Rate: nan\n", "")
        rust = run_diff.RunResult(0, "Rate: nan\n", "")
        self.assertFalse(self.compare(case, py, rust))

    def test_infinite_values_do_not_pass(self):
        case = _case(required_labels=("Rate",))
        for token in ("inf", "-inf", "Infinity"):
            with self.subTest(token=token):
                py = run_diff.RunResult(0, f"Rate: {token}\n", "")
                self.assertFalse(self.compare(case, py, py))

    def test_finite_replaced_by_nan_in_a_file_does_not_pass(self):
        """The file comparator had the same NaN-is-agreement path."""
        case = _case(compare_stream="none", compare_files=("out.txt",),
                     numeric_files=("out.txt",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            (Path(py_tmp) / "out.txt").write_text("value 1.0\n")
            (Path(rust_tmp) / "out.txt").write_text("value nan\n")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))

    def test_zero_sized_compared_file_does_not_pass(self):
        """Two empty files are byte-identical, which is not a result."""
        case = _case(compare_stream="none", compare_files=("out.txt",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            (Path(py_tmp) / "out.txt").write_text("")
            (Path(rust_tmp) / "out.txt").write_text("")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))

    def test_declared_empty_file_passes(self):
        """A case may declare that empty is the correct artifact for a named file.

        Three real cases need this: commands with nothing to report, where upstream
        and the port both write a zero-byte file. The exemption names the file, so it
        cannot silently extend to an artifact nobody looked at.
        """
        case = _case(compare_stream="none", compare_files=("out.txt", "other.txt"),
                     allow_empty_files=("out.txt",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            (Path(py_tmp) / "out.txt").write_text("")
            (Path(rust_tmp) / "out.txt").write_text("")
            for d in (py_tmp, rust_tmp):
                (Path(d) / "other.txt").write_text("")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))

    def test_exemption_does_not_extend_to_undeclared_files(self):
        """The exemption is per-file, not per-case."""
        case = _case(compare_stream="none", compare_files=("a.txt", "b.txt"),
                     allow_empty_files=("a.txt",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            for d in (py_tmp, rust_tmp):
                (Path(d) / "a.txt").write_text("")
                (Path(d) / "b.txt").write_text("")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))

    def test_matching_nonempty_file_still_passes(self):
        """Positive control: the file comparator must remain capable of passing."""
        case = _case(compare_stream="none", compare_files=("out.txt",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            (Path(py_tmp) / "out.txt").write_text("value 1.0\n")
            (Path(rust_tmp) / "out.txt").write_text("value 1.0\n")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertTrue(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))

    def test_expected_file_present_but_empty_does_not_pass(self):
        case = _case(compare_stream="none", expect_files=("out.txt",))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            for d in (py_tmp, rust_tmp):
                (Path(d) / "out.txt").write_text("")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))


class PinnedDivergenceTests(unittest.TestCase):
    """A recorded divergence must be PINNED, not merely tolerated.

    `divergent_files` exists so that a known difference cannot be forgotten. It is
    only worth having if it fails in both directions: when the divergence is fixed
    (files become identical) and when it changes shape (difference falls below the
    recorded floor). A check that only ever passes is decoration, and DIV-0024 is
    exactly the kind of entry that goes stale unnoticed -- which is how DIV-0005 was
    found to be out of date in the first place.
    """

    HEADER = "Percentile\t1\t2\t3\n"
    ROW = "{}\t{}\t{}\t{}\n"

    def _run(self, py_row, rust_row, floor):
        case = _case(compare_stream="none",
                     divergent_files=(("out.geneBodyCoverage.txt", floor),))
        ok = run_diff.RunResult(0, "", "")
        buffer = io.StringIO()
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            for d, row in ((py_tmp, py_row), (rust_tmp, rust_row)):
                (Path(d) / "out.geneBodyCoverage.txt").write_text(self.HEADER + self.ROW.format(*row))
            with contextlib.redirect_stdout(buffer):
                verdict = run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp))
        return verdict, buffer.getvalue()

    def test_identical_files_fail_because_the_divergence_is_gone(self):
        row = ("sample", 8000.0, 11.0, 20.0)
        verdict, text = self._run(row, row, 10.0)
        self.assertFalse(verdict, "identical files must fail a pinned divergence")
        self.assertIn("now IDENTICAL", text)
        self.assertIn("Retire the ledger", text)

    def test_difference_at_or_above_the_floor_passes(self):
        verdict, text = self._run(("s", 8000.0, 11.0, 20.0), ("s", 8000.0, 0.0, 0.0), 10.0)
        self.assertTrue(verdict)
        self.assertIn("largest difference 20", text)

    def test_difference_below_the_floor_fails_as_changed_shape(self):
        verdict, text = self._run(("s", 8000.0, 11.0, 20.0), ("s", 8000.0, 9.0, 20.0), 10.0)
        self.assertFalse(verdict, "a shrunken divergence must fail, not pass quietly")
        self.assertIn("below the recorded floor", text)

    def test_leading_sample_name_cell_is_not_mistaken_for_a_numeric_row(self):
        # The first cell is a text label. If it caused the row to be skipped, no
        # numeric cell would be compared and the divergence could not be pinned --
        # so the case would report success without having looked at anything.
        verdict, text = self._run(("sample_name", 1.0, 2.0, 3.0), ("sample_name", 1.0, 2.0, 9.0), 5.0)
        self.assertTrue(verdict)
        self.assertIn("largest difference 6", text)

    def test_unparseable_table_is_not_silently_accepted(self):
        case = _case(compare_stream="none",
                     divergent_files=(("out.geneBodyCoverage.txt", 10.0),))
        ok = run_diff.RunResult(0, "", "")
        with tempfile.TemporaryDirectory() as py_tmp, tempfile.TemporaryDirectory() as rust_tmp:
            for d in (py_tmp, rust_tmp):
                (Path(d) / "out.geneBodyCoverage.txt").write_text("no numbers at all\n")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(run_diff.compare_results(case, ok, ok, Path(py_tmp), Path(rust_tmp)))


if __name__ == "__main__":
    unittest.main()
