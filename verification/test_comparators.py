#!/usr/bin/env python3
"""Shared comparators used by the differential runner and the benchmark harness.

Both harnesses compare artifacts produced by two independent implementations, so
a weakness here is shared by both: the audit found the benchmark's comparators
accepting changed BAM quality scores, duplicate flags and NM tags, a FASTQ with a
trailing partial record, and finite text replaced by NaN. Those are the same
class of defect this module exists to prevent in one place.

Every function here fails closed. A comparator that cannot decide must return a
failure with a reason, never a pass: a false pass converts a silent regression
into a reported speedup, while a false failure costs one debugging session.

Needs only pysam for the BAM path; the text and FASTQ paths are standard library.

Usage:
    oracle/venv/bin/python3 verification/test_comparators.py
"""
from __future__ import annotations

import importlib.util
import re
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "rseqc_comparators", REPO / "verification" / "comparators.py")
assert SPEC and SPEC.loader
cmp_mod = importlib.util.module_from_spec(SPEC)
sys.modules["rseqc_comparators"] = cmp_mod
SPEC.loader.exec_module(cmp_mod)

try:
    import pysam
    HAVE_PYSAM = True
except ImportError:
    HAVE_PYSAM = False

write_bam = cmp_mod.write_bam


class TextComparatorTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def compare(self, left, right, **kw):
        (self.dir / "a.txt").write_text(left)
        (self.dir / "b.txt").write_text(right)
        return cmp_mod.compare_text(self.dir / "a.txt", self.dir / "b.txt", **kw)

    def test_identical_text_passes(self):
        self.assertEqual(self.compare("a\t1.0\n", "a\t1.0\n"), (True, ""))

    def test_one_changed_character_fails(self):
        ok, msg = self.compare("a\t1.0\n", "a\t1.1\n")
        self.assertFalse(ok)
        self.assertTrue(msg, "a failure must name the differing line")

    def test_line_count_difference_fails(self):
        ok, msg = self.compare("a\nb\n", "a\n")
        self.assertFalse(ok)
        self.assertIn("line count", msg)

    def test_nonfinite_is_rejected(self):
        for token in ("nan", "inf", "-inf", "Infinity", "NaN"):
            for pair in (("1.0", token), (token, "1.0"), (token, token)):
                with self.subTest(token=token, pair=pair):
                    ok, msg = self.compare(f"a\t{pair[0]}\n", f"a\t{pair[1]}\n")
                    self.assertFalse(ok, f"{pair} must not pass")
                    self.assertIn("nonfinite", msg)

    def test_nonfinite_in_a_later_cell_is_rejected(self):
        ok, msg = self.compare("a\t1.0\t2.0\n", "a\t1.0\tnan\n")
        self.assertFalse(ok)
        self.assertIn("nonfinite", msg)

    def test_label_named_nan_is_not_a_value(self):
        """A word containing 'nan' as a substring is a label, not a number."""
        ok, _ = self.compare("nanotube\t1\n", "nanotube\t1\n")
        self.assertTrue(ok)
        ok, _ = self.compare("nanotube\t1\n", "nanotube\t2\n")
        self.assertFalse(ok)


class FastqComparatorTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.good = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n"

    def compare(self, left, right):
        (self.dir / "a.fq").write_text(left)
        (self.dir / "b.fq").write_text(right)
        return cmp_mod.compare_fastq(self.dir / "a.fq", self.dir / "b.fq")

    def test_identical_fastq_passes(self):
        self.assertEqual(self.compare(self.good, self.good), (True, ""))

    def test_trailing_partial_record_fails(self):
        ok, msg = self.compare(self.good, self.good + "@trunc\nAC")
        self.assertFalse(ok)
        self.assertIn("multiple of 4", msg)

    def test_partial_record_fails_in_either_arm(self):
        for left, right in ((self.good, self.good + "@x\nAC"),
                            (self.good + "@x\nAC", self.good)):
            with self.subTest(side="both"):
                ok, _ = self.compare(left, right)
                self.assertFalse(ok)

    def test_length_mismatch_fails(self):
        ok, msg = self.compare(self.good, "@r1\nACGTA\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")
        self.assertFalse(ok)
        self.assertIn("length", msg)

    def test_bad_header_fails(self):
        ok, msg = self.compare(self.good, self.good.replace("@r1", "r1", 1))
        self.assertFalse(ok)
        self.assertIn("@", msg)

    def test_bad_separator_fails(self):
        ok, msg = self.compare(self.good, self.good.replace("+", "-", 1))
        self.assertFalse(ok)
        self.assertIn("+", msg)

    def test_changed_sequence_fails(self):
        ok, msg = self.compare(self.good, self.good.replace("ACGT", "ACGA"))
        self.assertFalse(ok)
        self.assertIn("sequence", msg)

    def test_changed_quality_fails(self):
        ok, msg = self.compare(self.good, self.good.replace("IIII", "!!!!"))
        self.assertFalse(ok)
        self.assertIn("quality", msg)

    def test_empty_file_fails(self):
        ok, msg = self.compare("", self.good)
        self.assertFalse(ok)
        self.assertTrue(msg)

    def test_record_count_difference_fails(self):
        ok, msg = self.compare(self.good, self.good + "@r3\nAC\n+\nII\n")
        self.assertFalse(ok)
        self.assertIn("record count", msg)


@unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
class BamComparatorTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.base = [
            {"name": "r1", "start": 100, "qual": "IIIIIIIIII", "tags": {"NM": 0}},
            {"name": "r2", "start": 200, "qual": "JJJJJJJJJJ", "tags": {"NM": 1}},
        ]

    def pair(self, mutated, require_sorted=False, header=None):
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, self.base)
        write_bam(b, mutated, header=header)
        return cmp_mod.compare_bam(a, b, require_sorted=require_sorted)

    def test_identical_passes(self):
        self.assertEqual(self.pair([dict(r) for r in self.base]), (True, ""))

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_quality_mutation_fails(self):
        mutated = [dict(self.base[0], qual="!!!!!!!!!!"), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_duplicate_flag_mutation_fails(self):
        mutated = [dict(self.base[0], flag=1024), self.base[1]]
        ok, msg = self.pair(mutated)
        self.assertFalse(ok)
        self.assertTrue(msg)

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_secondary_flag_mutation_fails(self):
        mutated = [dict(self.base[0], flag=256), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_tag_mutation_fails(self):
        mutated = [dict(self.base[0], tags={"NM": 99}), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_added_tag_fails(self):
        mutated = [dict(self.base[0], tags={"NM": 0, "MD": 10}), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_coordinate_shift_fails(self):
        mutated = [dict(self.base[0], start=101), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_cigar_change_fails(self):
        mutated = [dict(self.base[0], cigar=[(0, 9), (1, 1)]), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_sequence_change_fails(self):
        # The default written sequence is ACGTACGTAC, so this must actually differ.
        mutated = [dict(self.base[0], seq="ACGTACGTAG"), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_mapq_change_fails(self):
        mutated = [dict(self.base[0], mapq=1), self.base[1]]
        self.assertFalse(self.pair(mutated)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_mate_position_change_fails(self):
        paired = [
            {"name": "p1", "flag": 99, "next_ref": 0, "next_start": 300, "tlen": 200},
            {"name": "p1", "flag": 147, "next_ref": 0, "next_start": 100, "tlen": -200},
        ]
        mutated = [paired[0], dict(paired[1], next_start=999)]
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, paired)
        write_bam(b, mutated)
        self.assertFalse(cmp_mod.compare_bam(a, b)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_insert_size_change_fails(self):
        paired = [
            {"name": "p1", "flag": 99, "next_ref": 0, "next_start": 300, "tlen": 200},
            {"name": "p1", "flag": 147, "next_ref": 0, "next_start": 100, "tlen": -200},
        ]
        mutated = [paired[0], dict(paired[1], tlen=-201)]
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, paired)
        write_bam(b, mutated)
        self.assertFalse(cmp_mod.compare_bam(a, b)[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_record_count_difference_fails(self):
        self.assertFalse(self.pair([dict(self.base[0])])[0])

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_reference_length_change_fails(self):
        other = {"HD": {"VN": "1.6", "SO": "coordinate"},
                 "SQ": [{"SN": "chr1", "LN": 200000}]}
        ok, msg = self.pair([dict(r) for r in self.base], header=other)
        self.assertFalse(ok)
        self.assertIn("reference dictionary", msg)

    @unittest.skipUnless(HAVE_PYSAM, "pysam")
    def test_sort_order_change_is_ignored_by_default_and_caught_when_required(self):
        """Order-sensitivity is a per-command contract, so it must be a parameter.

        Several benchmarked commands legitimately emit unsorted subsets, so the
        default must stay order-insensitive; a command whose contract is
        coordinate-sorted output passes require_sorted=True and must then catch it.
        """
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, self.base)
        write_bam(b, [self.base[1], self.base[0]])
        self.assertTrue(cmp_mod.compare_bam(a, b)[0],
                        "default must tolerate order differences")
        ok, msg = cmp_mod.compare_bam(a, b, require_sorted=True)
        self.assertFalse(ok)
        self.assertIn("sorted order", msg)


class NumericTableComparatorTests(unittest.TestCase):
    """The numeric table comparator the differential runner now shares.

    These exist because the shared copy replaced two divergent private ones, and
    "we deduplicated" is not evidence that the merged version still rejects
    corruption. Every test feeds a defect and requires a FAILING result.
    """

    def cmp(self, a, b):
        return cmp_mod.compare_numeric_table(a, b)

    def test_equal_tables_pass(self):
        self.assertEqual(self.cmp(b"1 2\n3 4\n", b"1 2\n3 4\n"), (True, ""))

    def test_changed_value_fails(self):
        ok, msg = self.cmp(b"1 2\n", b"1 3\n")
        self.assertFalse(ok)
        self.assertIn("'2' vs '3'", msg)

    def test_changed_row_count_fails(self):
        ok, _ = self.cmp(b"1 2\n3 4\n", b"1 2\n")
        self.assertFalse(ok)

    def test_changed_cell_count_fails(self):
        ok, msg = self.cmp(b"1 2\n", b"1 2 3\n")
        self.assertFalse(ok)
        self.assertIn("cell count", msg)

    def test_symmetric_nan_fails(self):
        """The defect the old `Decimal` path got right by accident; keep it asserted.

        A table whose metric collapsed to NaN on both arms is a broken command, not
        two implementations agreeing on a broken answer.
        """
        for spelling in (b"nan", b"NaN", b"-nan", b"inf", b"-inf", b"Infinity"):
            with self.subTest(spelling=spelling):
                ok, msg = self.cmp(b"1 " + spelling + b"\n", b"1 " + spelling + b"\n")
                self.assertFalse(ok, f"{spelling!r} must not compare equal to itself")
                self.assertIn("nonfinite", msg)

    def test_asymmetric_nan_fails(self):
        ok, _ = self.cmp(b"nan 1\n", b"5 1\n")
        self.assertFalse(ok)

    def test_nonfinite_detected_in_either_table(self):
        for side in (0, 1):
            with self.subTest(side=side):
                a, b = (b"nan 1\n", b"5 1\n") if side == 0 else (b"5 1\n", b"nan 1\n")
                self.assertFalse(self.cmp(a, b)[0])

    def test_label_columns_fall_back_to_string_equality(self):
        self.assertTrue(self.cmp(b"total 10\n", b"total 10\n")[0])
        ok, _ = self.cmp(b"total 10\n", b"other 10\n")
        self.assertFalse(ok)

    def test_float_spelling_difference_is_a_difference(self):
        """`1.0` and `1.00` are the same number; these tables must still agree."""
        self.assertTrue(self.cmp(b"1.0 2.50\n", b"1.00 2.5\n")[0])

    def test_precise_difference_is_not_absorbed(self):
        """Decimal, not float: a real numeric difference must not round away."""
        ok, _ = self.cmp(b"1.0000000001\n", b"1.0\n")
        self.assertFalse(ok)

    def test_accepts_str_as_well_as_bytes(self):
        self.assertTrue(self.cmp("1 2\n", "1 2\n")[0])
        self.assertTrue(self.cmp("1 2\n", b"1 2\n")[0])
        self.assertFalse(self.cmp("1 2\n", "1 3\n")[0])

    def test_numeric_cell_equal_directly(self):
        self.assertTrue(cmp_mod.numeric_cell_equal("1.0", "1.00"))
        self.assertFalse(cmp_mod.numeric_cell_equal("1", "2"))
        self.assertFalse(cmp_mod.numeric_cell_equal("nan", "nan"))
        self.assertFalse(cmp_mod.numeric_cell_equal("inf", "inf"))
        self.assertTrue(cmp_mod.numeric_cell_equal("abc", "abc"))
        self.assertFalse(cmp_mod.numeric_cell_equal("abc", "abd"))


class SharedModuleUseTests(unittest.TestCase):
    """The differential runner must actually import the shared comparators.

    The claim "the comparators are shared" is only true if the import exists. A
    renamed function or a deleted `sys.path` entry would leave the claim standing in
    a docstring while each harness went back to its own private copy -- which is the
    state the audit found.
    """

    def test_run_diff_imports_comparators(self):
        source = (REPO / "verification" / "run_diff.py").read_text()
        self.assertIn("from comparators import", source,
                      "run_diff.py must import the shared comparators")
        self.assertNotIn("def _numeric_table_equal", source.split("from comparators import")[0],
                         "the private copy must be gone, not merely supplemented")

    def test_bench_imports_comparators(self):
        source = (REPO / "benchmarks" / "bench.py").read_text()
        self.assertIn("from comparators import", source)

    def test_numeric_helper_is_not_redefined_locally(self):
        """Only the thin `_numeric_equal` alias may remain in run_diff.py."""
        source = (REPO / "verification" / "run_diff.py").read_text()
        for name in ("numeric_cell_equal", "compare_numeric_table"):
            redefinitions = [
                line for line in source.splitlines()
                if line.startswith(f"def {name}")
            ]
            self.assertEqual(redefinitions, [],
                             f"{name} must come from comparators, not be redefined")


if __name__ == "__main__":
    unittest.main(verbosity=2)
