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
import re
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
    # stripped from compare_files' content before comparing -- needed
    # for R scripts that embed their own output path (e.g. `pdf('...')`),
    # which legitimately differs between the two sides' separate scratch
    # directories even when the actual data is identical. This is the
    # "explicitly named normalization for paths" PORTING_PLAN.md's Step 4
    # table allows, not a way to hide a real difference.
    normalize_paths: bool = False


@dataclasses.dataclass
class RunResult:
    exit_code: int
    stdout: str
    stderr: str


def run(argv: list[str], pythonpath: str | None = None) -> RunResult:
    env = None
    if pythonpath is not None:
        import os

        env = dict(os.environ)
        env["PYTHONPATH"] = pythonpath
    proc = subprocess.run(argv, capture_output=True, text=True, cwd=REPO_ROOT, env=env)
    return RunResult(exit_code=proc.returncode, stdout=proc.stdout, stderr=proc.stderr)


LABEL_COUNT_RE = re.compile(r"^([A-Za-z][^:]*?):\s*(-?\d+(?:\.\d+)?)")


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
            out[m.group(1).strip()] = m.group(2)
    return out


def stream_for(result: RunResult, which: str) -> str:
    if which == "stdout":
        return result.stdout
    if which == "stderr":
        return result.stderr
    return result.stdout + result.stderr


def run_case(case: Case) -> bool:
    import shutil
    import tempfile

    print(f"=== {case.name} ===")
    case.ensure_fixture()

    ok = True
    py_dir = Path(tempfile.mkdtemp(prefix="rseqc_verify_py_"))
    rust_dir = Path(tempfile.mkdtemp(prefix="rseqc_verify_rust_"))
    try:
        py_args = case.py_args(py_dir)
        rust_args = case.rust_args(rust_dir)

        py_result = run([str(ORACLE_PYTHON), str(ORACLE_SCRIPTS / case.py_script)] + py_args, pythonpath=ORACLE_PYTHONPATH)
        rust_result = run([str(RUST_BIN_DIR / case.rust_bin)] + rust_args)

        if py_result.exit_code != rust_result.exit_code:
            print(f"  FAIL exit code: python={py_result.exit_code} rust={rust_result.exit_code}")
            ok = False
        else:
            print(f"  exit code matches ({py_result.exit_code})")

        if case.compare_stream != "none":
            py_text = stream_for(py_result, case.compare_stream)
            rust_text = stream_for(rust_result, case.compare_stream)
            py_counts = extract_labeled_counts(py_text)
            rust_counts = extract_labeled_counts(rust_text)
            all_labels = sorted(set(py_counts) | set(rust_counts))
            for label in all_labels:
                pv = py_counts.get(label)
                rv = rust_counts.get(label)
                if pv != rv:
                    print(f"  FAIL '{label}': python={pv!r} rust={rv!r}")
                    ok = False
            if ok:
                print(f"  stream comparison PASS ({len(all_labels)} labeled values matched)")
            else:
                print("  --- python output ---")
                print(py_text)
                print("  --- rust output ---")
                print(rust_text)

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
            if py_bytes == rust_bytes:
                print(f"  file '{rel_path}' PASS (byte-identical, {len(py_bytes)} bytes)")
            else:
                print(f"  FAIL file '{rel_path}': byte content differs ({len(py_bytes)} vs {len(rust_bytes)} bytes)")
                print(f"  --- python {rel_path} ---")
                print(py_bytes.decode("utf-8", errors="replace"))
                print(f"  --- rust {rel_path} ---")
                print(rust_bytes.decode("utf-8", errors="replace"))
                ok = False
    finally:
        shutil.rmtree(py_dir, ignore_errors=True)
        shutil.rmtree(rust_dir, ignore_errors=True)

    return ok


def ensure_bam_stat_fixture() -> None:
    fixture = REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam"
    if fixture.is_file():
        return
    generator = REPO_ROOT / "verification" / "fixtures" / "make_bam_stat_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture)], cwd=REPO_ROOT, check=True)


def _bam_stat_args(_scratch_dir: Path) -> list[str]:
    return ["-i", str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam")]


def _nvc_fixture_path() -> str:
    return str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam")


CASES: list[Case] = [
    Case(
        name="bam_stat_basic",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="bam_stat.py",
        rust_bin="bam_stat",
        py_args=_bam_stat_args,
        rust_args=_bam_stat_args,
        compare_stream="stdout",
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
]


def main() -> int:
    requested = sys.argv[1:]
    cases = [c for c in CASES if not requested or c.name in requested]

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
