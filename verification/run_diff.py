#!/usr/bin/env python3
"""Differential verification: runs the real upstream Python CLI and this
port's Rust binary against the same fixture and diffs their output.

This is a deliberately minimal first slice of PORTING_PLAN.md's Steps 3-4
(a genuine runner that executes reference and candidate against the same
case and reports actionable differences), not the full spec (no public
dataset downloads, no BAM/BigWig-decoded structural comparators, no
metamorphic invariants yet). Extend CASES below as more commands get
fixtures.

Usage:
    oracle/venv/bin/python3 verification/run_diff.py [case_name ...]
    (with no arguments, runs every case)

Requires: the Rust workspace already built in release mode
(`cargo build --workspace --release`) and `oracle/venv`/`oracle/upstream-src`
present (both gitignored, set up separately -- see docs/PORTING_PLAN.md).
Uses only the Python standard library at runtime.
"""
from __future__ import annotations

import dataclasses
import os
import re
import signal
import subprocess
import sys
from pathlib import Path
from typing import Callable

REPO_ROOT = Path(__file__).resolve().parent.parent
ORACLE_SCRIPTS = REPO_ROOT / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHONPATH = str(REPO_ROOT / "oracle" / "upstream-src" / "src")
ORACLE_PYTHON = REPO_ROOT / "oracle" / "venv" / "bin" / "python3"
RUST_BIN_DIR = REPO_ROOT / "target" / "release"


@dataclasses.dataclass
class Case:
    name: str
    # Ensures the fixture file(s) this case needs exist; returns nothing.
    ensure_fixture: Callable[[], None]
    py_script: str  # filename under oracle/upstream-src/scripts/
    rust_bin: str  # filename under target/release/
    # CLI args as a function of the run's own scratch directory (so
    # `-o`-style output-prefix args can point INTO it). Separate
    # callables for each side since the two CLIs' flag surfaces are not
    # always identical (e.g. a flag genuinely not implemented yet on
    # one side) -- default both to the same function when they match.
    py_args: Callable[[Path], list[str]]
    rust_args: Callable[[Path], list[str]]
    # Which stream(s) carry the numeric report to compare: "stdout",
    # "stderr", "both" (concatenated), or "none" (skip stream
    # comparison, e.g. when the command's real output is in files).
    compare_stream: str = "stdout"
    # Relative file paths (relative to the run's own scratch directory)
    # to byte-compare between the two sides' output, e.g. ["out.NVC.xls"].
    compare_files: tuple[str, ...] = ()
    # When True, each side's own scratch-directory absolute path is
    # stripped from compare_files' content AND from the compared stream
    # text before comparing -- needed for R scripts that embed their own
    # output path (e.g. `pdf('...')`) or console reports that echo their
    # own output file paths (e.g. split_bam-style "<path> (Read 1): N"
    # lines), which legitimately differ between the two sides' separate
    # scratch directories even when the actual data is identical. This is
    # the "explicitly named normalization for paths" PORTING_PLAN.md's
    # Step 4 table allows, not a way to hide a real difference.
    normalize_paths: bool = False
    # Files whose whitespace-delimited table cells should be compared as
    # numbers where possible.  This is explicit per case: byte comparison
    # remains the default, while numeric tables do not fail on Python's
    # harmless `0` versus `0.0` rendering difference.
    numeric_files: tuple[str, ...] = ()
    # Exit status expected from both implementations.  A positive
    # compatibility case must therefore not pass merely because both sides
    # failed in the same way.
    expected_exit_code: int = 0
    # Stream comparison is either semantic labelled-number comparison or
    # exact text comparison.  Cases with file outputs normally use "none".
    stream_format: str = "labels"
    # Optional labels that must be present on both sides.  This closes the
    # empty-stream false-pass path while allowing commands whose reports are
    # intentionally file-only.
    required_labels: tuple[str, ...] = ()
    allow_empty_stream: bool = False
    # Wall-clock limit for each implementation.  A timeout is a failed run,
    # never an equivalent result.
    timeout_s: float = 120.0


@dataclasses.dataclass
class RunResult:
    exit_code: int | None
    stdout: str
    stderr: str
    timed_out: bool = False


def run(
    argv: list[str],
    pythonpath: str | None = None,
    *,
    cwd: Path = REPO_ROOT,
    timeout_s: float = 120.0,
) -> RunResult:
    env = dict(os.environ)
    if pythonpath is not None:
        env["PYTHONPATH"] = pythonpath
    try:
        proc = subprocess.Popen(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=env,
            start_new_session=True,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
            return RunResult(exit_code=proc.returncode, stdout=stdout, stderr=stderr)
        except subprocess.TimeoutExpired as exc:
            # subprocess.run does not kill grandchildren.  Since each
            # command is a process group, terminate the group before
            # returning a failed result.
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = proc.communicate()
            stdout = stdout or exc.stdout or ""
            stderr = stderr or exc.stderr or ""
            return RunResult(
                exit_code=None,
                stdout=stdout if isinstance(stdout, str) else stdout.decode(errors="replace"),
                stderr=(stderr if isinstance(stderr, str) else stderr.decode(errors="replace"))
                + f"\n[verification timeout after {timeout_s:g}s]",
                timed_out=True,
            )
    except OSError as exc:
        return RunResult(exit_code=None, stdout="", stderr=f"[verification could not execute command: {exc}]")


# Keep the token boundary explicit: the old expression parsed ``1e-3`` as
# ``1`` and ``1e+3`` as ``1``.  Decimal and scientific notation are accepted,
# but arbitrary trailing prose is still tolerated because upstream reports
# often append percentages or units.
NUMBER_TOKEN = r"[-+]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|(?:inf|nan))"
LABEL_COUNT_RE = re.compile(rf"^([A-Za-z][^:]*?):\s*({NUMBER_TOKEN})(?=\s|$)", re.IGNORECASE)


def extract_labeled_counts(text: str) -> dict[str, str]:
    """Pulls `<label>: <number>` pairs out of free-form report text,
    tolerant of trailing percentages/whitespace/decoration -- enough to
    compare the STATISTIC VALUES two differently-formatted reports
    print, without requiring byte-identical text (timestamps, spacing,
    and decoration legitimately differ between the two CLIs)."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        m = LABEL_COUNT_RE.match(line.strip())
        if m:
            label = m.group(1).strip()
            if label in out:
                raise ValueError(f"duplicate report label: {label!r}")
            out[label] = m.group(2)
    return out


def stream_for(result: RunResult, which: str) -> str:
    if which == "stdout":
        return result.stdout
    if which == "stderr":
        return result.stderr
    return result.stdout + result.stderr


def _numeric_equal(left: str, right: str) -> bool:
    from decimal import Decimal, InvalidOperation

    try:
        lval = Decimal(left)
        rval = Decimal(right)
    except InvalidOperation:
        return left == right
    if lval.is_nan() or rval.is_nan():
        return lval.is_nan() and rval.is_nan()
    return lval == rval


def _numeric_table_equal(left: bytes, right: bytes) -> bool:
    """Compare a whitespace-delimited text table with numeric cell semantics."""
    left_rows = left.decode("utf-8", errors="replace").splitlines()
    right_rows = right.decode("utf-8", errors="replace").splitlines()
    if len(left_rows) != len(right_rows):
        return False
    from decimal import Decimal, InvalidOperation

    for left_row, right_row in zip(left_rows, right_rows):
        left_cells = left_row.split()
        right_cells = right_row.split()
        if len(left_cells) != len(right_cells):
            return False
        for left_cell, right_cell in zip(left_cells, right_cells):
            try:
                Decimal(left_cell)
                Decimal(right_cell)
            except InvalidOperation:
                if left_cell != right_cell:
                    return False
            else:
                if not _numeric_equal(left_cell, right_cell):
                    return False
    return True


def compare_results(case: Case, py_result: RunResult, rust_result: RunResult, py_dir: Path, rust_dir: Path) -> bool:
    """Compare two completed runs and print actionable diagnostics."""
    ok = True
    expected = case.expected_exit_code
    for side, result in (("python", py_result), ("rust", rust_result)):
        if result.exit_code != expected:
            suffix = " (timed out)" if result.timed_out else ""
            print(f"  FAIL {side} exit code: expected={expected} actual={result.exit_code}{suffix}")
            if result.stderr:
                print(f"  --- {side} stderr ---")
                print(result.stderr)
            ok = False
    if py_result.exit_code != rust_result.exit_code:
        print(f"  FAIL exit code differs: python={py_result.exit_code} rust={rust_result.exit_code}")
        ok = False
    if not ok:
        # Never treat output produced by a failed process as a valid
        # compatibility result.
        return False

    if case.compare_stream != "none":
        py_text = stream_for(py_result, case.compare_stream)
        rust_text = stream_for(rust_result, case.compare_stream)
        if case.normalize_paths:
            py_text = py_text.replace(str(py_dir), "<SCRATCH_DIR>")
            rust_text = rust_text.replace(str(rust_dir), "<SCRATCH_DIR>")
        if case.stream_format == "exact":
            if not py_text and not case.allow_empty_stream:
                print("  FAIL stream comparison: both streams are empty")
                ok = False
            elif py_text != rust_text:
                print("  FAIL exact stream content differs")
                print("  --- python output ---")
                print(py_text)
                print("  --- rust output ---")
                print(rust_text)
                ok = False
            else:
                print(f"  stream comparison PASS (exact, {len(py_text)} characters)")
        elif case.stream_format == "labels":
            try:
                py_counts = extract_labeled_counts(py_text)
                rust_counts = extract_labeled_counts(rust_text)
            except ValueError as exc:
                print(f"  FAIL stream parser: {exc}")
                return False
            if not case.allow_empty_stream and (not py_counts or not rust_counts):
                print(
                    "  FAIL stream comparison: expected labelled values, "
                    f"got python={len(py_counts)} rust={len(rust_counts)}"
                )
                ok = False
            missing = sorted(set(case.required_labels) - set(py_counts))
            missing += sorted(set(case.required_labels) - set(rust_counts))
            if missing:
                print(f"  FAIL required report labels missing: {sorted(set(missing))}")
                ok = False
            all_labels = sorted(set(py_counts) | set(rust_counts))
            for label in all_labels:
                pv = py_counts.get(label)
                rv = rust_counts.get(label)
                if pv is None or rv is None or not _numeric_equal(pv, rv):
                    print(f"  FAIL '{label}': python={pv!r} rust={rv!r}")
                    ok = False
            if ok:
                print(f"  stream comparison PASS ({len(all_labels)} labelled values matched)")
            else:
                print("  --- python output ---")
                print(py_text)
                print("  --- rust output ---")
                print(rust_text)
        else:
            print(f"  FAIL unsupported stream format: {case.stream_format!r}")
            ok = False

    for rel_path in case.compare_files:
        py_file = py_dir / rel_path
        rust_file = rust_dir / rel_path
        if not py_file.is_file() or not rust_file.is_file():
            print(f"  FAIL file '{rel_path}': python_exists={py_file.is_file()} rust_exists={rust_file.is_file()}")
            ok = False
            continue
        py_bytes = py_file.read_bytes()
        rust_bytes = rust_file.read_bytes()
        if case.normalize_paths:
            py_bytes = py_bytes.replace(str(py_dir).encode(), b"<SCRATCH_DIR>")
            rust_bytes = rust_bytes.replace(str(rust_dir).encode(), b"<SCRATCH_DIR>")
        numeric = rel_path in case.numeric_files
        file_equal = _numeric_table_equal(py_bytes, rust_bytes) if numeric else py_bytes == rust_bytes
        if file_equal:
            qualifier = "numeric cells" if numeric else f"byte-identical, {len(py_bytes)} bytes"
            print(f"  file '{rel_path}' PASS ({qualifier})")
        else:
            print(f"  FAIL file '{rel_path}': byte content differs ({len(py_bytes)} vs {len(rust_bytes)} bytes)")
            print(f"  --- python {rel_path} ---")
            print(py_bytes.decode("utf-8", errors="replace"))
            print(f"  --- rust {rel_path} ---")
            print(rust_bytes.decode("utf-8", errors="replace"))
            ok = False
    return ok


def run_case(case: Case) -> bool:
    import shutil
    import tempfile

    print(f"=== {case.name} ===")
    case.ensure_fixture()

    ok = True
    keep_failures = os.environ.get("RSEQC_KEEP_FAILURES") == "1"
    if keep_failures:
        failure_root = REPO_ROOT / "verification" / "failures"
        failure_root.mkdir(parents=True, exist_ok=True)
        py_dir = Path(tempfile.mkdtemp(prefix=f"{case.name}_py_", dir=failure_root))
        rust_dir = Path(tempfile.mkdtemp(prefix=f"{case.name}_rust_", dir=failure_root))
    else:
        py_dir = Path(tempfile.mkdtemp(prefix="rseqc_verify_py_"))
        rust_dir = Path(tempfile.mkdtemp(prefix="rseqc_verify_rust_"))
    try:
        py_args = case.py_args(py_dir)
        rust_args = case.rust_args(rust_dir)

        py_result = run(
            [str(ORACLE_PYTHON), str(ORACLE_SCRIPTS / case.py_script)] + py_args,
            pythonpath=ORACLE_PYTHONPATH,
            cwd=py_dir,
            timeout_s=case.timeout_s,
        )
        rust_result = run(
            [str(RUST_BIN_DIR / case.rust_bin)] + rust_args,
            cwd=rust_dir,
            timeout_s=case.timeout_s,
        )
        ok = compare_results(case, py_result, rust_result, py_dir, rust_dir)
    finally:
        if not keep_failures or ok:
            shutil.rmtree(py_dir, ignore_errors=True)
            shutil.rmtree(rust_dir, ignore_errors=True)
        else:
            print(f"  kept failure evidence under {py_dir} and {rust_dir}")

    return ok


def ensure_bam_stat_fixture() -> None:
    fixture = REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam"
    if fixture.is_file():
        return
    generator = REPO_ROOT / "verification" / "fixtures" / "make_bam_stat_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture)], cwd=REPO_ROOT, check=True)


def ensure_regression_fixtures() -> None:
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    required = (
        fixture_dir / "regression_single_exon.bed12",
        fixture_dir / "regression_fpkm_mate_overlap.bam",
        fixture_dir / "regression_fpkm_mate_overlap.bam.bai",
        fixture_dir / "regression_fpkm_fetch_span.bam",
        fixture_dir / "regression_fpkm_fetch_span.bam.bai",
        fixture_dir / "regression_splice_fetch.bed12",
        fixture_dir / "regression_rna_equals.bam",
        fixture_dir / "regression_rna_equals.bam.bai",
        fixture_dir / "regression_overlap_pair.bam",
        fixture_dir / "regression_overlap_pair.bam.bai",
        fixture_dir / "regression_genebody_depth.bam",
        fixture_dir / "regression_genebody_depth.bam.bai",
    )
    if all(path.is_file() for path in required):
        return
    generator = fixture_dir / "make_regression_fixtures.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture_dir)], cwd=REPO_ROOT, check=True)


def _bam_stat_args(_scratch_dir: Path) -> list[str]:
    return ["-i", str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam")]


def _nvc_fixture_path() -> str:
    return str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam")


def _regression_fixture(name: str) -> str:
    return str(REPO_ROOT / "verification" / "fixtures" / name)


def ensure_hexamer_fixtures() -> None:
    # Plain-text FASTA fixtures, committed directly (no generator needed,
    # same pattern as the .bed12 fixtures).
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    required = (
        fixture_dir / "regression_hexamer_reads.fa",
        fixture_dir / "regression_hexamer_ref.fa",
    )
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"missing committed fixture(s): {missing}")


CASES: list[Case] = [
    Case(
        name="bam_stat_basic",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="bam_stat.py",
        rust_bin="bam_stat",
        py_args=_bam_stat_args,
        rust_args=_bam_stat_args,
        compare_stream="stdout",
        required_labels=("Total records", "Unmapped reads", "Read-1"),
    ),
    Case(
        name="read_NVC_basic",
        # Reuses bam_stat_basic.bam -- read_NVC.py only needs a BAM with
        # some real sequence content, no need for a dedicated fixture.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_NVC.py",
        rust_bin="read_NVC",
        # --skip-plot is needed on the python side to avoid it trying to
        # invoke Rscript (not installed here) -- the Rust CLI does not
        # accept this flag at all, because it never implements plot
        # generation in the first place (a pre-existing, disclosed gap;
        # see crates/cli/src/bin/read_NVC.rs's own module docs), so it's
        # simply omitted there rather than needing to be "skipped".
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out")],
        compare_stream="none",
        compare_files=("out.NVC.xls",),
    ),
    Case(
        name="read_GC_basic",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_GC.py",
        rust_bin="read_GC",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out")],
        compare_stream="none",
        compare_files=("out.GC.xls",),
    ),
    Case(
        name="read_duplication_basic",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_duplication.py",
        rust_bin="read_duplication",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out")],
        compare_stream="none",
        compare_files=("out.pos.DupRate.xls", "out.seq.DupRate.xls"),
    ),
    Case(
        name="read_quality_basic",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_quality.py",
        rust_bin="read_quality",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out")],
        compare_stream="none",
        # The .qual.r script embeds each side's own absolute scratch-
        # directory path in its pdf('...') line -- normalize it away,
        # see normalize_paths' docstring above.
        compare_files=("out.qual.r",),
        normalize_paths=True,
    ),
    Case(
        name="FPKM_count_exonic_mate_overlap",
        ensure_fixture=ensure_regression_fixtures,
        py_script="FPKM_count.py",
        rust_bin="FPKM_count",
        py_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_fpkm_mate_overlap.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "-e",
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_fpkm_mate_overlap.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "-e",
        ],
        compare_stream="none",
        compare_files=("out.FPKM.xls",),
    ),
    Case(
        name="RNA_fragment_size_equals_cigar",
        ensure_fixture=ensure_regression_fixtures,
        py_script="RNA_fragment_size.py",
        rust_bin="RNA_fragment_size",
        py_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_rna_equals.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-n",
            "1",
            "-o",
            str(scratch_dir / "out.tsv"),
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_rna_equals.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-n",
            "1",
            "-o",
            str(scratch_dir / "out.tsv"),
        ],
        compare_stream="none",
        compare_files=("out.tsv",),
    ),
    Case(
        name="FPKM_count_fetch_reference_span",
        ensure_fixture=ensure_regression_fixtures,
        py_script="FPKM_count.py",
        rust_bin="FPKM_count",
        py_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_fpkm_fetch_span.bam"),
            "-r",
            _regression_fixture("regression_splice_fetch.bed12"),
            "-o",
            str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_fpkm_fetch_span.bam"),
            "-r",
            _regression_fixture("regression_splice_fetch.bed12"),
            "-o",
            str(scratch_dir / "out"),
        ],
        compare_stream="none",
        compare_files=("out.FPKM.xls",),
    ),
    Case(
        name="geneBody_coverage_max_depth",
        ensure_fixture=ensure_regression_fixtures,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_genebody_depth.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_genebody_depth.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "--skip-plot",
        ],
        compare_stream="none",
        compare_files=("out.geneBodyCoverage.txt",),
    ),
    Case(
        name="geneBody_coverage_pair_overlap",
        ensure_fixture=ensure_regression_fixtures,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_overlap_pair.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_overlap_pair.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "--skip-plot",
        ],
        compare_stream="none",
        compare_files=("out.geneBodyCoverage.txt",),
        numeric_files=("out.geneBodyCoverage.txt",),
    ),
    Case(
        name="tin_pair_overlap",
        ensure_fixture=ensure_regression_fixtures,
        py_script="tin.py",
        rust_bin="tin",
        py_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_overlap_pair.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-c",
            "0",
            "-n",
            "100",
            "-o",
            str(scratch_dir),
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_overlap_pair.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-c",
            "0",
            "-n",
            "100",
            "-o",
            str(scratch_dir),
        ],
        compare_stream="none",
        compare_files=("regression_overlap_pair.tin.xls", "regression_overlap_pair.summary.txt"),
    ),
    Case(
        name="clipping_profile_basic",
        # Regression case for the defaultdict(int)-vs-float duck-typing
        # bug: an untouched position's Clipped_nt must render as bare "0",
        # not "0.0" -- see crates/commands/src/clipping_profile.rs.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="clipping_profile.py",
        rust_bin="clipping_profile",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE"],
        compare_stream="none",
        compare_files=("out.clipping_profile.xls", "out.clipping_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="insertion_profile_basic",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="insertion_profile.py",
        rust_bin="insertion_profile",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE"],
        compare_stream="none",
        compare_files=("out.insertion_profile.xls", "out.insertion_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="bam2fq_single_end",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="bam2fq.py",
        rust_bin="bam2fq",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s"],
        compare_stream="none",
        compare_files=("out.fastq",),
    ),
    Case(
        name="split_bam_basic",
        ensure_fixture=ensure_regression_fixtures,
        py_script="split_bam.py",
        rust_bin="split_bam",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        compare_stream="stdout",
        required_labels=("Total records",),
        # BAM byte output is not a fair comparison here: noodles' BGZF
        # writer and pysam/htslib's produce different compressed bytes
        # (block boundaries, @PG header line) for logically identical
        # content -- the stdout record-count report is the real
        # compatibility signal for this case.
    ),
    Case(
        name="infer_experiment_basic",
        ensure_fixture=ensure_regression_fixtures,
        py_script="infer_experiment.py",
        rust_bin="infer_experiment",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
        ],
        compare_stream="both",
        allow_empty_stream=True,
    ),
    Case(
        name="split_paired_bam_basic",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="split_paired_bam.py",
        rust_bin="split_paired_bam",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out")],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out")],
        compare_stream="stdout",
        stream_format="exact",
        normalize_paths=True,
    ),
    Case(
        name="deletion_profile_no_deletions",
        # bam_stat_basic.bam's reads are all plain 20M (no deletions), so
        # this exercises the "0 qualifying reads" path -- still a real,
        # oracle-verified case, not a dummy. stdout is intentionally not
        # compared: deletion_profile.py has no --skip-plot flag, so the
        # real upstream CLI always spawns Rscript and its console noise
        # ("null device"/"1") leaks into stdout, an already-accepted
        # architecture difference (this port never invokes Rscript).
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="deletion_profile.py",
        rust_bin="deletion_profile",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out")],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out")],
        compare_stream="none",
        compare_files=("out.deletion_profile.txt", "out.deletion_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="mismatch_profile_no_mismatches",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="mismatch_profile.py",
        rust_bin="mismatch_profile",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out")],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out")],
        compare_stream="stdout",
        stream_format="exact",
        compare_files=("out.mismatch_profile.xls", "out.mismatch_profile.r"),
    ),
    Case(
        name="junction_annotation_no_junctions",
        ensure_fixture=ensure_regression_fixtures,
        py_script="junction_annotation.py",
        rust_bin="junction_annotation",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        compare_stream="both",
        stream_format="exact",
        compare_files=("out.junction.xls", "out.junction_plot.r"),
    ),
    Case(
        name="junction_annotation_with_junction",
        # Real splice read (20M100N20M) against a single-exon model with
        # no annotated introns -- exercises the "complete_novel"
        # classification path plus the .bed/.Interact.bed outputs.
        ensure_fixture=ensure_regression_fixtures,
        py_script="junction_annotation.py",
        rust_bin="junction_annotation",
        py_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_fpkm_fetch_span.bam"),
            "-r", _regression_fixture("regression_splice_fetch.bed12"),
            "-m", "10", "-o", str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_fpkm_fetch_span.bam"),
            "-r", _regression_fixture("regression_splice_fetch.bed12"),
            "-m", "10", "-o", str(scratch_dir / "out"),
        ],
        compare_stream="none",
        compare_files=("out.junction.xls", "out.junction.bed", "out.junction.Interact.bed"),
    ),
    Case(
        name="inner_distance_basic",
        ensure_fixture=ensure_regression_fixtures,
        py_script="inner_distance.py",
        rust_bin="inner_distance",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.inner_distance.txt", "out.inner_distance_freq.txt", "out.inner_distance_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_distribution_basic",
        # stdout only (not "both"): upstream additionally writes
        # "Processing ... Done"/"Finished" progress lines to stderr that
        # this port's CLI does not -- a disclosed, accepted difference,
        # not part of the actual report this case checks. The report
        # itself has no "label: value" colons (fixed-width columns), so
        # it needs exact text comparison rather than the labels parser.
        ensure_fixture=ensure_regression_fixtures,
        py_script="read_distribution.py",
        rust_bin="read_distribution",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
        ],
        compare_stream="stdout",
        stream_format="exact",
    ),
    Case(
        name="divide_bam_basic",
        # Both sides seed their own RNG deterministically, but Python's
        # random.Random and Rust's rand::StdRng are different algorithms
        # (DIV-0017) -- the per-file query-name assignment legitimately
        # differs, so stdout (the "<path>\t<count>" per-file table) is NOT
        # compared, and the output BAMs are not byte-compared for the same
        # reason split_bam_basic's aren't (different bgzf writers, and
        # here also different record distribution). What IS a fair,
        # bug-catching check: the stderr progress/summary lines, which
        # don't depend on the RNG stream at all -- this is exactly the
        # text the CLI used to omit entirely (a raw `{:?}` Debug dump of
        # the counts struct instead of upstream's report).
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="divide_bam.py",
        rust_bin="divide_bam",
        py_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam"),
            "-n", "2", "-o", str(scratch_dir / "out"), "--seed", "42", "-s",
        ],
        rust_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam"),
            "-n", "2", "-o", str(scratch_dir / "out"), "--seed", "42", "-s",
        ],
        compare_stream="stderr",
        stream_format="exact",
    ),
    Case(
        name="read_hexamer_basic",
        # Found by this case (before it was formalized as a harness
        # entry): a bare relative `-o` filename (no directory component,
        # e.g. "result.tsv") errored out in the Rust CLI ("output path
        # has no parent directory") but succeeds in upstream, whose
        # `Path.parent` for such a path is `Path('.')` -- the current
        # directory, which exists. This case uses stdout output (no -o)
        # so it doesn't exercise that specific path, but does cover the
        # header-name-uppercasing fix in seq_generator (dead code for
        # this command's own output, but shared/documented-as-exact) and
        # the whole report/progress-line pipeline end to end.
        ensure_fixture=ensure_hexamer_fixtures,
        py_script="read_hexamer.py",
        rust_bin="read_hexamer",
        py_args=lambda scratch_dir: [
            "-i", f"{_regression_fixture('regression_hexamer_reads.fa')},{_regression_fixture('regression_hexamer_ref.fa')}",
        ],
        rust_args=lambda scratch_dir: [
            "-i", f"{_regression_fixture('regression_hexamer_reads.fa')},{_regression_fixture('regression_hexamer_ref.fa')}",
        ],
        compare_stream="both",
        stream_format="exact",
    ),
    Case(
        name="junction_saturation_basic",
        # Reuses the same BAM/BED12 pair as junction_annotation_with_
        # junction: exactly one qualifying spliced read (20M100N20M
        # against a single-exon model, no annotated introns), so the
        # random-shuffle-order dependence of the middle percentile steps
        # is moot -- there is only one splice site total, so every
        # percentile step's counts are deterministic regardless of RNG
        # stream. This exercises the previously entirely-missing
        # progress-message pipeline (reading bed file / Load BAM file /
        # shuffling / per-percentile summary lines) end to end.
        ensure_fixture=ensure_regression_fixtures,
        py_script="junction_saturation.py",
        rust_bin="junction_saturation",
        py_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_fpkm_fetch_span.bam"),
            "-r", _regression_fixture("regression_splice_fetch.bed12"),
            "-m", "10", "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_fpkm_fetch_span.bam"),
            "-r", _regression_fixture("regression_splice_fetch.bed12"),
            "-m", "10", "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        compare_stream="both",
        stream_format="exact",
        compare_files=("out.junctionSaturation_plot.r",),
        normalize_paths=True,
    ),
]


def main() -> int:
    requested = sys.argv[1:]
    case_by_name = {case.name: case for case in CASES}
    unknown = sorted(set(requested) - set(case_by_name))
    if unknown:
        print(f"Unknown case(s): {', '.join(unknown)}", file=sys.stderr)
        return 2
    cases = [case_by_name[name] for name in requested] if requested else CASES

    if not cases:
        print(f"No matching cases for: {requested}", file=sys.stderr)
        return 2

    if not RUST_BIN_DIR.is_dir():
        print(f"error: {RUST_BIN_DIR} not found -- run `cargo build --workspace --release` first", file=sys.stderr)
        return 2

    all_ok = True
    for case in cases:
        all_ok &= run_case(case)
        print()

    if all_ok:
        print(f"All {len(cases)} case(s) PASSED")
        return 0
    print("One or more cases FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
