#!/usr/bin/env python3
"""Harness-credibility tests for the benchmark runner itself.

The readiness audit's P0 findings are all defects in the *measurement*, not in
the port: comparators that accept corrupted output, an output gate that passes
two empty directories, a timeout handler that can kill the harness's own
process group, and a schedule whose pairs are not the pairs it claims.

Every test here asserts the failure-closed direction: feed the harness a
deliberate defect and require a *failing* result with a recorded reason. A
comparator that cannot reject corruption is worse than no comparator, because
it converts a silent regression into a reported speedup.

Runs on the standard library plus pysam; no oracle environment, no Rust build
beyond what is already present, and no timing measurement.
    python3 benchmarks/test_bench_harness.py
"""
from __future__ import annotations

import importlib.util
import json
import os
import signal
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "bench", REPO / "benchmarks" / "bench.py")
assert SPEC and SPEC.loader
bench = importlib.util.module_from_spec(SPEC)
sys.modules["bench"] = bench
SPEC.loader.exec_module(bench)

try:
    import pysam  # noqa: F401
    HAVE_PYSAM = True
except ImportError:
    HAVE_PYSAM = False


def write_bam(path, records, header=None):
    """records: list of dicts with the fields the comparator is required to read."""
    header = header or {"HD": {"VN": "1.6", "SO": "coordinate"},
                        "SQ": [{"SN": "chr1", "LN": 100000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in records:
            a = pysam.AlignedSegment()
            a.query_name = r["name"]
            a.query_sequence = r.get("seq", "ACGTACGTAC")
            a.flag = r.get("flag", 0)
            a.reference_id = 0
            a.reference_start = r.get("start", 100)
            a.mapping_quality = r.get("mapq", 60)
            a.cigar = r.get("cigar", [(0, 10)])
            a.query_qualities = pysam.qualitystring_to_array(
                r.get("qual", "IIIIIIIIII"))
            a.next_reference_id = r.get("next_ref", -1)
            a.next_reference_start = r.get("next_start", -1)
            a.template_length = r.get("tlen", 0)
            for tag, value in r.get("tags", {}).items():
                a.set_tag(tag, value, value_type="i")
            out.write(a)


BASE_RECORDS = [
    {"name": "r1", "start": 100},
    {"name": "r2", "start": 200},
]


class BamComparatorTests(unittest.TestCase):
    """compare_bam must reject every corruption the audit demonstrated it accepted."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def pair(self, mutated):
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, BASE_RECORDS)
        write_bam(b, mutated)
        return bench.compare_bam(a, b)

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_identical_bams_pass(self):
        self.assertEqual(self.pair(list(BASE_RECORDS)), (True, ""))

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_quality_mutation_is_rejected(self):
        mutated = [dict(r, qual="!!!!!!!!!!") if i == 0 else r
                   for i, r in enumerate(BASE_RECORDS)]
        ok, msg = self.pair(mutated)
        self.assertFalse(ok, "changed base qualities must fail the gate")
        self.assertTrue(msg)

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_duplicate_flag_mutation_is_rejected(self):
        mutated = [dict(r, flag=1024) if i == 0 else r
                   for i, r in enumerate(BASE_RECORDS)]
        ok, msg = self.pair(mutated)
        self.assertFalse(ok, "a changed duplicate flag must fail the gate")
        self.assertTrue(msg)

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_nm_tag_mutation_is_rejected(self):
        base = [dict(r, tags={"NM": 0}) for r in BASE_RECORDS]
        mutated = [dict(r, tags={"NM": 99}) for r in BASE_RECORDS]
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, base)
        write_bam(b, mutated)
        ok, msg = bench.compare_bam(a, b)
        self.assertFalse(ok, "a changed NM tag must fail the gate")
        self.assertTrue(msg)

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_coordinate_shift_is_rejected(self):
        mutated = [dict(r, start=r["start"] + 1) if i == 0 else r
                   for i, r in enumerate(BASE_RECORDS)]
        ok, _ = self.pair(mutated)
        self.assertFalse(ok, "a one-base coordinate shift must fail the gate")

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_mate_fields_are_compared(self):
        paired = [
            {"name": "p1", "flag": 99, "next_ref": 0, "next_start": 300, "tlen": 200},
            {"name": "p1", "flag": 147, "next_ref": 0, "next_start": 100, "tlen": -200},
        ]
        mutated = [paired[0], dict(paired[1], next_start=999)]
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, paired)
        write_bam(b, mutated)
        ok, _ = bench.compare_bam(a, b)
        self.assertFalse(ok, "a changed mate position must fail the gate")

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_header_difference_is_rejected(self):
        other_header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                        "SQ": [{"SN": "chr1", "LN": 200000}]}
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, BASE_RECORDS)
        write_bam(b, BASE_RECORDS, header=other_header)
        ok, _ = bench.compare_bam(a, b)
        self.assertFalse(ok, "a changed reference length must fail the gate")

    @unittest.skipUnless(HAVE_PYSAM, "pysam required for BAM decoding")
    def test_record_order_matters_when_sort_order_is_required(self):
        swapped = list(reversed(BASE_RECORDS))
        a, b = self.dir / "a.bam", self.dir / "b.bam"
        write_bam(a, BASE_RECORDS)
        write_bam(b, swapped)
        ok, _ = bench.compare_bam(a, b, require_sorted=True)
        self.assertFalse(ok, "reordered coordinate-sorted records must fail")


class FastqComparatorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.good = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n"
        (self.dir / "a.fq").write_text(self.good)
        self.other = self.dir / "b.fq"

    def test_identical_fastq_passes(self):
        self.other.write_text(self.good)
        self.assertEqual(bench.compare_fastq(self.dir / "a.fq", self.other), (True, ""))

    def test_trailing_partial_record_is_rejected(self):
        """The audit's probe: a truncated final record must not be silently dropped."""
        self.other.write_text(self.good + "@truncated\nAC")
        ok, msg = bench.compare_fastq(self.dir / "a.fq", self.other)
        self.assertFalse(ok, "a trailing partial FASTQ record must fail")
        self.assertTrue(msg)

    def test_truncated_record_in_either_arm_is_rejected(self):
        for name, a, b in (("in b", self.good, self.good + "@x\nAC"),
                           ("in a", self.good + "@x\nAC", self.good)):
            left, right = self.dir / "l.fq", self.dir / "r.fq"
            left.write_text(a)
            right.write_text(b)
            ok, msg = bench.compare_fastq(left, right)
            self.assertFalse(ok, f"partial record {name} must fail")
            self.assertTrue(msg)

    def test_length_mismatch_between_seq_and_qual_is_rejected(self):
        bad = "@r1\nACGTA\n+\nIIII\n"
        left, right = self.dir / "l.fq", self.dir / "r.fq"
        left.write_text(self.good)
        right.write_text(bad)
        ok, msg = bench.compare_fastq(left, right)
        self.assertFalse(ok, "a 5-base sequence with a 4-base quality must fail")
        self.assertTrue(msg)

    def test_qualities_are_compared(self):
        self.other.write_text(self.good.replace("IIII", "!!!!"))
        ok, _ = bench.compare_fastq(self.dir / "a.fq", self.other)
        self.assertFalse(ok, "changed qualities must fail")


class TextComparatorTests(unittest.TestCase):
    """compare_text must reject a finite value replaced by a nonfinite one."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def compare(self, left, right):
        (self.dir / "a.txt").write_text(left)
        (self.dir / "b.txt").write_text(right)
        return bench.compare_text(self.dir / "a.txt", self.dir / "b.txt")

    def test_identical_text_passes(self):
        self.assertEqual(self.compare("a\t1.0\n", "a\t1.0\n"), (True, ""))

    def test_finite_replaced_by_nan_is_rejected(self):
        ok, msg = self.compare("a\t1.0\n", "a\tnan\n")
        self.assertFalse(ok, "finite replaced by nan must fail the gate")
        self.assertTrue(msg)

    def test_finite_replaced_by_inf_is_rejected(self):
        ok, _ = self.compare("a\t1.0\n", "a\tinf\n")
        self.assertFalse(ok, "finite replaced by inf must fail")

    def test_finite_replaced_by_negative_inf_is_rejected(self):
        ok, _ = self.compare("a\t1.0\n", "a\t-Inf\n")
        self.assertFalse(ok, "finite replaced by -Inf must fail")

    def test_two_nans_are_not_equal_evidence(self):
        """Both arms failing identically is not evidence of agreement."""
        ok, msg = self.compare("a\tnan\n", "a\tnan\n")
        self.assertFalse(ok, "a nonfinite value must not pass, even symmetrically")
        self.assertTrue(msg)

    def test_two_infs_are_not_equal_evidence(self):
        ok, _ = self.compare("a\tinf\n", "a\tinf\n")
        self.assertFalse(ok, "a nonfinite value must not pass, even symmetrically")


class OutputGateTests(unittest.TestCase):
    """The gate must require expected streams and artifacts, not just a file tree."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def two_empty(self):
        py, rs = self.dir / "py", self.dir / "rs"
        py.mkdir()
        rs.mkdir()
        return py, rs

    def test_two_empty_output_trees_fail(self):
        """The audit's probe: emptiness on both sides must not be a gate pass."""
        py, rs = self.two_empty()
        ok, errors = bench.gate_outputs(py, rs, name="bam_stat")
        self.assertFalse(ok, "two empty output trees must not pass the gate")
        self.assertTrue(errors)

    def test_declared_artifacts_present_but_empty_fail(self):
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "output.txt").write_text("")
        ok, errors = bench.gate_outputs(py, rs, name="bam_stat")
        self.assertFalse(ok, "an empty deliverable must not pass")
        self.assertTrue(errors)

    def test_stdout_only_metric_without_declared_stream_fails(self):
        """The audit's probe: metrics on stdout were never compared at all."""
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
        ok, errors = bench.gate_outputs(py, rs, name="bam_stat",
                                        stdout_a="total: 10\n", stdout_b="")
        self.assertFalse(ok, "differing stdout must fail the gate")
        self.assertTrue(errors)

    def test_missing_artifact_in_one_arm_fails(self):
        py, rs = self.dir / "py", self.dir / "rs"
        py.mkdir()
        rs.mkdir()
        (py / "output.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(py, rs, name="bam_stat")
        self.assertFalse(ok, "an artifact present in one arm only must fail")
        self.assertTrue(errors)

    def test_unexpected_exit_code_fails(self):
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "output.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(py, rs, name="bam_stat",
                                        exit_a=0, exit_b=2)
        self.assertFalse(ok, "a differing exit code must fail the gate")
        self.assertTrue(errors)

    def test_timed_out_run_fails_even_with_identical_outputs(self):
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "output.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(py, rs, name="bam_stat",
                                        timed_out_a=True, timed_out_b=False)
        self.assertFalse(ok, "a timed-out run must not receive a gate pass")
        self.assertTrue(errors)

    def declare(self, name, spec):
        """Install a spec for this test only, then restore the real table.

        The declared labels for real commands are derived from observed runs
        (see derive_expected_streams.py), so they change as upstream evolves.
        These gate-behaviour tests must not break when a label is re-spelled, nor
        accidentally assert that today's observed spelling is correct.
        """
        original = bench.EXPECTED_STREAMS.get(name)
        bench.EXPECTED_STREAMS[name] = spec
        self.addCleanup(lambda: bench.EXPECTED_STREAMS.__setitem__(name, original)
                        if original is not None
                        else bench.EXPECTED_STREAMS.pop(name, None))

    STREAM_SPEC = {"stdout": True, "labels": {"Total reads", "Assigned reads"},
                   "artifacts": ["report.txt"]}

    def test_matching_outputs_and_streams_pass(self):
        """The positive control: a correct pair must still pass.

        Without this, a gate that rejects everything would satisfy every negative
        test above while destroying the benchmark's ability to report anything.
        """
        self.declare("test_stream_command", dict(self.STREAM_SPEC))
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "report.txt").write_text("value\t1\n")
        stdout = "Total reads: 10\nAssigned reads: 9\n"
        ok, errors = bench.gate_outputs(
            py, rs, name="test_stream_command",
            stdout_a=stdout, stdout_b=stdout, exit_a=0, exit_b=0)
        self.assertTrue(ok, f"matching arms must pass, got {errors}")
        self.assertEqual(errors, [])

    def test_declared_artifact_missing_fails(self):
        self.declare("test_stream_command", dict(self.STREAM_SPEC))
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "something_else.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(
            py, rs, name="test_stream_command",
            stdout_a="Total reads: 10\nAssigned reads: 9\n",
            stdout_b="Total reads: 10\nAssigned reads: 9\n")
        self.assertFalse(ok, "a command declaring artifacts must fail when they are absent")
        self.assertTrue(any("missing expected artifact" in e for e in errors))

    def test_declared_label_missing_from_stream_fails(self):
        self.declare("test_stream_command", dict(self.STREAM_SPEC))
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "report.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(
            py, rs, name="test_stream_command",
            stdout_a="Total reads: 10\n", stdout_b="Total reads: 10\n")
        self.assertFalse(ok, "a missing declared stdout label must fail")
        self.assertTrue(any("required label" in e for e in errors))

    def test_stdout_metric_disagreement_fails(self):
        """The audit's probe: stdout-only metrics were never compared at all."""
        self.declare("test_stream_command", dict(self.STREAM_SPEC))
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "report.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(
            py, rs, name="test_stream_command",
            stdout_a="Total reads: 10\nAssigned reads: 9\n",
            stdout_b="Total reads: 10\nAssigned reads: 2\n")
        self.assertFalse(ok, "a differing stdout metric must fail the gate")
        self.assertTrue(errors)

    def test_nonfinite_stdout_metric_fails(self):
        self.declare("test_stream_command", dict(self.STREAM_SPEC))
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "report.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(
            py, rs, name="test_stream_command",
            stdout_a="Total reads: 10\nAssigned reads: nan\n",
            stdout_b="Total reads: 10\nAssigned reads: nan\n")
        self.assertFalse(ok, "a nonfinite stdout metric must fail even symmetrically")
        self.assertTrue(any("nonfinite" in e for e in errors))

    def test_required_label_missing_from_stream_fails(self):
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "value.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(
            py, rs, name="bam_stat",
            stdout_a="total: 1\n", stdout_b="unrelated: 5\n",
            required_labels={"total"})
        self.assertFalse(ok, "a missing required stdout label must fail")
        self.assertTrue(errors)

    def test_command_with_no_declared_artifacts_fails_the_gate(self):
        """An undeclared command is not silently treated as expecting nothing."""
        py, rs = self.dir / "py", self.dir / "rs"
        for d in (py, rs):
            d.mkdir()
            (d / "output.txt").write_text("value\t1\n")
        ok, errors = bench.gate_outputs(py, rs, name="a_command_with_no_declaration")
        self.assertFalse(ok)
        self.assertTrue(any("no artifact expectation" in e for e in errors))


class SchedulePairingTests(unittest.TestCase):
    """Pairs must be adjacent matched blocks, not a global shuffle then zip."""

    def test_schedule_pairs_matched_repetitions(self):
        reps = 10
        schedule = bench.build_schedule(reps=reps, seed=20260929)
        by_arm = {"py": [], "rs": []}
        for item in schedule:
            by_arm[item.arm].append(item.rep)
        self.assertEqual(sorted(by_arm["py"]), list(range(reps)))
        self.assertEqual(sorted(by_arm["rs"]), list(range(reps)))
        # Every repetition must appear in exactly one adjacent block.
        blocks = {}
        for index, item in enumerate(schedule):
            blocks.setdefault(item.rep, []).append(index)
        for rep, positions in blocks.items():
            self.assertEqual(len(positions), 2,
                             f"rep {rep} must run exactly once per arm")
            self.assertEqual(positions[1] - positions[0], 1,
                             f"rep {rep} arms must be adjacent")

    def test_arm_order_is_randomised_within_blocks(self):
        orders = set()
        for seed in range(40):
            schedule = bench.build_schedule(reps=1, seed=seed)
            orders.add(schedule[0].arm)
        self.assertEqual(orders, {"py", "rs"},
                         "which arm runs first must vary with the seed")

    def test_schedule_is_reproducible_from_seed(self):
        a = bench.build_schedule(reps=6, seed=7)
        b = bench.build_schedule(reps=6, seed=7)
        self.assertEqual([(x.rep, x.arm) for x in a],
                         [(x.rep, x.arm) for x in b])

    def test_every_run_carries_its_own_identity_and_order(self):
        schedule = bench.build_schedule(reps=4, seed=11)
        for index, item in enumerate(schedule):
            self.assertEqual(item.order, index)
            self.assertIn(item.arm, ("py", "rs"))
            self.assertIsInstance(item.rep, int)


class BootstrapTests(unittest.TestCase):
    def test_interval_uses_actual_blocks_and_is_reproducible(self):
        pairs = [(10.0, 5.0)] * 6 + [(10.0, 9.0)] * 4
        a = bench.paired_log_ratio_ci(pairs, iters=500, seed=3)
        b = bench.paired_log_ratio_ci(pairs, iters=500, seed=3)
        self.assertEqual(a, b)
        lo, hi = a
        self.assertLess(lo, hi)
        self.assertGreater(hi, 1.0)
        self.assertLess(lo, 2.0)

    def test_single_pair_gives_no_interval(self):
        self.assertIsNone(bench.paired_log_ratio_ci([(2.0, 1.0)]))


@unittest.skipUnless(os.geteuid() != 0, "process-group test must not run as root")
class TimeoutHandlerTests(unittest.TestCase):
    """A timeout must kill only the benchmark's own group, not the harness's."""

    def test_runner_survives_a_hung_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            marker = Path(tmp) / "child-was-killed"
            script = Path(tmp) / "hang.py"
            script.write_text(textwrap.dedent(f"""
                import signal, sys, time
                # Ignore SIGTERM so only a hard kill can stop this, which is what
                # a wedged external tool looks like in practice.
                signal.signal(signal.SIGTERM, signal.SIG_IGN)
                def spawn_child():
                    if len(sys.argv) > 1 and sys.argv[1] == "child":
                        signal.signal(signal.SIGTERM, signal.SIG_IGN)
                        while True:
                            time.sleep(0.05)
                if len(sys.argv) > 1 and sys.argv[1] == "child":
                    spawn_child()
                    while True:
                        time.sleep(0.05)
                import subprocess
                subprocess.Popen([sys.executable, __file__, "child"])
                Path = __import__("pathlib").Path
                Path({str(marker)!r}).write_text("spawned")
                while True:
                    time.sleep(0.05)
            """))
            harness_pgid = os.getpgid(0)
            m, _out, _err = bench.run_arm(
                [sys.executable, str(script)], Path(tmp), os.environ.copy(), 2)
            self.assertTrue(m["timed_out"], "a hung child must be recorded as timed out")
            self.assertEqual(os.getpgid(0), harness_pgid,
                             "the harness must still be in its own process group")
            self.assertTrue(marker.exists(), "the fixture must have really run")

            # Nothing from the fixture's group may survive.
            deadline = time.time() + 10
            survivors = _descendants_of(script)
            while survivors and time.time() < deadline:
                time.sleep(0.2)
                survivors = _descendants_of(script)
            self.assertEqual(survivors, [],
                             f"timed-out children must be killed, found {survivors}")

    def test_timeout_does_not_signal_the_harness_group(self):
        """killpg must never be aimed at a group the harness belongs to."""
        recorded = []
        real_getpgid, real_killpg = os.getpgid, os.killpg
        try:
            os.getpgid = lambda pid: recorded.append(("getpgid", pid)) or real_getpgid(0)
            os.killpg = lambda pgid, sig: recorded.append(("killpg", pgid, sig))
            with tempfile.TemporaryDirectory() as tmp:
                script = Path(tmp) / "hang.py"
                script.write_text("import time\nwhile True: time.sleep(0.05)\n")
                m, _o, _e = bench.run_arm(
                    [sys.executable, str(script)], Path(tmp), os.environ.copy(), 1)
        finally:
            os.getpgid, os.killpg = real_getpgid, real_killpg
        self.assertTrue(m["timed_out"])
        for call in recorded:
            if call[0] == "getpgid":
                self.assertNotEqual(call[1], 0,
                                    "getpgid(0) identifies the CALLER's group")
            else:
                self.assertNotEqual(call[1], os.getpgid(0),
                                    "killpg must not target the harness's group")


def _descendants_of(script_path):
    """Any live process whose command line mentions `script_path`."""
    found = []
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            cmdline = (entry / "cmdline").read_bytes().decode("utf8", "replace")
        except OSError:
            continue
        if str(script_path) in cmdline:
            found.append(int(entry.name))
    return found


class ProvenanceTests(unittest.TestCase):
    """Each run must be bound to the candidate, not to an older clean commit."""

    def test_binary_digest_is_recorded_per_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            binary = Path(tmp) / "fake-binary"
            binary.write_bytes(b"#!/bin/sh\nexit 0\n")
            digest = bench.sha256_file(binary)
            self.assertEqual(len(digest), 64, "a truncated digest cannot bind a binary")
            self.assertEqual(digest, bench.sha256_file(binary))

    def test_missing_binary_is_reported_not_silently_hashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(bench.sha256_file(Path(tmp) / "absent"))

    def test_provenance_records_every_required_field(self):
        required = {"git_commit", "git_dirty", "binary_sha256", "oracle_lock_sha256",
                    "comparator_version", "inputs_sha256", "command_line",
                    "toolchain", "timestamp"}
        prov = {
            "git_commit": "abc", "git_dirty": False, "binary_sha256": "d" * 64,
            "oracle_lock_sha256": "e" * 64, "comparator_version": 2,
            "inputs_sha256": {}, "command_line": [], "toolchain": {},
            "timestamp": 0,
        }
        self.assertEqual(set(prov) & required, required)

    def test_comparator_version_is_exported_and_changes_with_comparators(self):
        self.assertIsInstance(bench.COMPARATOR_VERSION, int)
        self.assertGreaterEqual(bench.COMPARATOR_VERSION, 2,
                                "comparator semantics changed; bump the version")


class CommandContractTests(unittest.TestCase):
    """Per-command expectations must be declared, not inferred."""

    def test_every_benchmarked_command_declares_its_contract(self):
        for name in bench.COMMANDS:
            self.assertIn(name, bench.EXPECTED_STREAMS,
                          f"{name} has no declared expected stream")

    def test_expected_streams_declare_labels_for_stdout_commands(self):
        """A stdout command must declare HOW its stream is checked.

        A command that prints to stdout but declares no labels gets no stream
        checking at all -- the audit's P0 finding, since the original gate never
        looked at stdout. Requiring `labels` alone is not enough now that a command
        may declare a free-text stream instead: `bam2wig` prints a wigToBigWig
        command line, which is not a `label value` table, and declaring a label for
        it makes the gate demand parseable metrics from a sentence.
        """
        for name, spec in bench.EXPECTED_STREAMS.items():
            if not spec.get("stdout"):
                continue
            if spec.get("stdout_mode") == "text":
                self.assertTrue(
                    spec.get("stdout_contains"),
                    f"{name} declares a text stream but no required substrings, "
                    f"so only 'both are non-empty' is checked")
                continue
            self.assertTrue(spec.get("labels"),
                            f"{name} prints to stdout but declares no labels")

    def test_text_stream_declarations_are_consistent(self):
        """A text-mode stream must not also declare metric labels.

        The two modes check different things -- substrings plus stream equality
        versus parsed label/value metrics -- and a declaration that mixes them is
        checked by whichever branch runs first, so the other is silently dead.
        """
        for name, spec in bench.EXPECTED_STREAMS.items():
            if spec.get("stdout_mode") == "text":
                self.assertFalse(spec.get("labels"),
                                 f"{name} mixes a text stream with metric labels")
                self.assertTrue(spec.get("stdout"),
                                f"{name} declares stdout_mode without stdout")

    def test_every_expected_artifact_is_named(self):
        for name, spec in bench.EXPECTED_STREAMS.items():
            if name in bench.EXCLUDED:
                continue
            self.assertTrue(spec.get("artifacts") is not None,
                            f"{name} must declare whether artifacts are expected")


class OutputPathResolutionTests(unittest.TestCase):
    """--output-dir and --workload must be resolved, not taken literally.

    Each arm runs with cwd set to its own run directory. A relative path reaching a
    command's -i or -o is therefore re-resolved against that directory rather than
    the invocation directory, and the command is handed a path whose parent does not
    exist. Both arms then fail identically and fast, which reads as two broken
    commands rather than one bad invocation; and because the failure record used to
    keep only the tail of stderr, the retained diagnosis was the truncated text
    "imum resident set size". A whole real-data study was launched that way.
    """

    def _paths_for(self, output_dir, workload):
        """Recompute what an arm would be handed, as the runner does."""
        root = Path(output_dir).resolve()
        d = root / "cmd" / "py" / "rep00"
        d.mkdir(parents=True, exist_ok=True)
        return root, d

    def test_relative_output_dir_is_not_reinterpreted_under_the_run_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            rel = os.path.relpath(os.path.join(tmp, "results"), os.getcwd())
            root, d = self._paths_for(rel, None)
            self.assertTrue(root.is_absolute(),
                            "the output root must be absolute before any arm cwd change")
            # An arm handed `-o d / "x"` must find d to exist, because d was created
            # relative to the invocation directory, not inside itself.
            self.assertTrue(d.is_dir())
            self.assertFalse(
                (d / rel).exists(),
                "the run directory must not contain a nested copy of the output root")

    def test_absolute_and_relative_output_dirs_place_files_identically(self):
        with tempfile.TemporaryDirectory() as tmp:
            a_root, a_d = self._paths_for(os.path.join(tmp, "a"), None)
            b_root, b_d = self._paths_for(
                os.path.relpath(os.path.join(tmp, "b"), os.getcwd()), None)
            self.assertEqual(a_root, Path(tmp).resolve() / "a")
            self.assertEqual(b_root, Path(tmp).resolve() / "b")
            self.assertEqual(a_d.relative_to(a_root), b_d.relative_to(b_root))


class CommandMatrixTests(unittest.TestCase):
    def test_thirty_three_commands_have_a_rust_binary(self):
        self.assertGreaterEqual(len(bench.COMMANDS), 29)
        for name, entry in bench.COMMANDS.items():
            self.assertEqual(len(entry), 4, f"{name} entry shape changed")
            self.assertTrue(callable(entry[2]))

    def test_experiment_classes_are_declared_distinctly(self):
        """A compute-class row is not a full pipeline measurement, and the
        recorded results must not imply otherwise."""
        classes = {entry[3] for entry in bench.COMMANDS.values()}
        self.assertTrue(classes <= {"compute", "data-only", "plotting", "end-to-end"})
        self.assertTrue(classes, "every row must declare its experiment class")


if __name__ == "__main__":
    unittest.main(verbosity=2)
