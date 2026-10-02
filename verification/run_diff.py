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
# comparators.py lives beside this file; make it importable regardless of the caller's
# working directory, since the benchmark harness is loaded by path from elsewhere.
sys.path.insert(0, str(REPO_ROOT / "verification"))
ORACLE_SCRIPTS = REPO_ROOT / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHONPATH = str(REPO_ROOT / "oracle" / "upstream-src" / "src")
ORACLE_PYTHON = REPO_ROOT / "oracle" / "venv" / "bin" / "python3"
RUST_BIN_DIR = Path(os.environ.get("RSEQC_RUST_BIN_DIR", REPO_ROOT / "target" / "release"))


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
    # Files recorded as DIVERGENT rather than equal, each with the smallest
    # largest-absolute-cell-difference the divergence is known to have.
    # Asserting a divergence is stronger than tolerating one: the case passes only
    # while the files still differ by at least the recorded amount, so a fix that
    # makes them identical fails here and points at the ledger entry to retire, and
    # a divergence that changes shape fails here too. `file: floor` entries.
    divergent_files: tuple[tuple[str, float], ...] = ()
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
    # Files whose content should be gunzip-decompressed before comparison
    # (e.g. bam2fq.py's -c/--compress output). The gzip CONTAINER itself
    # (compressed bytes, header) is never byte-comparable between the two
    # implementations regardless of correctness: Python's `gzip.open`
    # stamps the current wall-clock time into the header on every run, so
    # even upstream's own output isn't reproducible run-to-run for that
    # reason -- only the decompressed content is a fair comparison.
    gzip_files: tuple[str, ...] = ()
    # Files that must exist (non-empty) on BOTH sides but whose content is
    # never byte-comparable in principle -- e.g. a `.bai` index, whose
    # binning-index encoder detail differs between noodles and htslib even
    # though both produce a functionally correct, real-htslib-readable
    # index (verified separately via a live pysam.fetch() probe during
    # development, not by this harness). Existence is still a genuine,
    # bug-catching check: a command that silently failed to write the file
    # at all would otherwise pass unnoticed.
    expect_files: tuple[str, ...] = ()
    # Exit status expected from both implementations.  A positive
    # compatibility case must therefore not pass merely because both sides
    # failed in the same way.
    expected_exit_code: int = 0
    # Per-side exit overrides for capability-scoped cases where one side is
    # a known-degraded upstream environment (e.g. sc_seqlogo_basic: the
    # pinned upstream crashes on its own pandas bug after writing the CSV,
    # while this port is expected to succeed). When an override is set on a
    # side, that side's exit must equal ITS override and the identical-exit
    # assertion between the two sides is skipped (the exits legitimately
    # differ). Overrides must always be declared explicitly; both-side
    # equal-failure is never a positive pass either way.
    py_expected_exit: int | None = None
    rust_expected_exit: int | None = None
    # Candidate-only artifacts that must exist and be non-empty on the Rust
    # side, without content comparison -- used when the pinned upstream
    # cannot produce the file at all in its environment (DIV-0016's rendered
    # logos) yet the candidate is required to. Content is checked separately
    # (e.g. structural unit tests in the render crate); here the existence
    # check catches a command that silently skips the artifact.
    rust_expect_files: tuple[str, ...] = ()
    # Stream comparison is either semantic labelled-number comparison or
    # exact text comparison.  Cases with file outputs normally use "none".
    stream_format: str = "labels"
    # Optional labels that must be present on both sides.  This closes the
    # empty-stream false-pass path while allowing commands whose reports are
    # intentionally file-only.
    required_labels: tuple[str, ...] = ()
    allow_empty_stream: bool = False
    # When True, timestamp-bearing log prefixes (Python logging's
    # "YYYY-MM-DD HH:MM:SS [LEVEL] " and printlog()'s "@ YYYY-MM-DD HH:MM:SS: ")
    # are removed from the start of each line of the compared stream --
    # timestamps can never be byte-reproduced (DIV-0019/DIV-0022 context).
    strip_log_prefixes: bool = False
    # Compared files that are permitted to be zero bytes ON BOTH SIDES.
    #
    # Two empty files are byte-identical, so a blanket "empty must fail" rule would
    # reject cases where an empty artifact is the correct answer -- and the rule
    # exists because the opposite failure is real: a command that writes nothing at
    # all, in either arm, used to compare equal. So the default stays closed and a
    # case must name the specific file and say why empty is right.
    #
    # Every current use is a command that legitimately has nothing to report, or
    # takes an early-exit path before writing:
    #   * mismatch_profile on an input with no mismatches -- upstream's own
    #     `sys.exit()` inside mismatchProfile bypasses the R-script writer on both
    #     sides, leaving a zero-byte file.
    #   * junction_annotation on an alignment with no junctions -- there are no rows
    #     to tabulate, and the stream comparison already asserts both arms reported
    #     the same "no junctions" text.
    allow_empty_files: tuple[str, ...] = ()
    # Wall-clock limit for each implementation.  A timeout is a failed run,
    # never an equivalent result. None means "use the default for the current
    # input class" -- 120s for the synthetic fixtures, much longer for the real
    # panel, which is orders of magnitude bigger. It must be None here and not
    # a number: this field is always passed to run() explicitly, so a numeric
    # default would silently win over any real-data-aware default.
    timeout_s: float | None = None


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
    timeout_s: float | None = None,
) -> RunResult:
    if timeout_s is None:
        timeout_s = 3600.0 if REAL_DATA_DIR else 120.0
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


LOG_PREFIX_RE = re.compile(
    r"^(?:\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:,\d+)? \[[A-Z]+\] |@ \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}: )",
    re.M,
)


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


# Numeric cell and table comparison are shared with the benchmark harness
# (verification/comparators.py) rather than reimplemented here. Two private copies is
# how the two harnesses drifted apart in the first place: the benchmark's copy
# accepted changed BAM qualities, flags and tags and a truncated FASTQ record, and the
# differential suite's own weaker copy meant nothing in this suite noticed.
#
# The rule both copies obey: NONFINITE IS NOT A VALUE. `NaN == NaN` is False, so a
# comparator written with plain equality rejects it for free; a comparator written
# with a tolerance or a "both look the same" shortcut accepts it. Two arms both
# dividing by zero is a defect that reached the output, not agreement, so it fails
# even when both sides produce the identical token.
from comparators import compare_numeric_table, numeric_cell_equal


def _numeric_equal(left: str, right: str) -> bool:
    """Whether one reported value agrees; nonfinite on either side never agrees."""
    return numeric_cell_equal(left, right)


def _numeric_table_equal(left: bytes, right: bytes) -> tuple[bool, str]:
    """Whether two whitespace-delimited numeric tables agree, with a reason."""
    return compare_numeric_table(left, right)


def _numeric_rows(data: bytes) -> list[list[float | None]] | None:
    """Parse a whitespace-delimited table whose cells are mostly numbers.

    A cell that does not parse becomes ``None`` rather than discarding its row:
    these reports begin with a text cell (a sample name, a ``Percentile`` header), and
    dropping the whole row instead is how a real numeric row silently stops being
    compared. Placeholders keep columns aligned, and only a table with no numeric
    cell anywhere is reported as unparseable.
    """
    rows = []
    for raw in data.decode("utf-8", "replace").splitlines():
        cells = raw.split()
        if not cells:
            continue
        values: list[float | None] = []
        for cell in cells:
            try:
                values.append(float(cell))
            except ValueError:
                values.append(None)
        if any(v is not None for v in values):
            rows.append(values)
    return rows or None


def compare_results(case: Case, py_result: RunResult, rust_result: RunResult, py_dir: Path, rust_dir: Path) -> bool:
    """Compare two completed runs and print actionable diagnostics."""
    ok = True
    expected = case.expected_exit_code
    for side, result in (("python", py_result), ("rust", rust_result)):
        side_expected = case.py_expected_exit if side == "python" else case.rust_expected_exit
        if side_expected is not None:
            expected = side_expected
        if result.exit_code != expected:
            suffix = " (timed out)" if result.timed_out else ""
            expected_desc = f"expected={expected}" + (f" [override for {side}]" if side_expected is not None else "")
            print(f"  FAIL {side} exit code: {expected_desc} actual={result.exit_code}{suffix}")
            if result.stderr:
                print(f"  --- {side} stderr ---")
                print(result.stderr)
            ok = False
    if case.py_expected_exit is None and case.rust_expected_exit is None and py_result.exit_code != rust_result.exit_code:
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
        if case.strip_log_prefixes:
            py_text = LOG_PREFIX_RE.sub("", py_text)
            rust_text = LOG_PREFIX_RE.sub("", rust_text)
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
        elif case.stream_format == "error_line":
            # Compare only the final non-empty line -- the operative
            # `prog: error: <message>` -- plus the exit status, which the caller
            # checks separately. The usage block printed above it is deliberately
            # excluded: upstream emits argparse's `usage: ...` summary and the port
            # emits clap's, and hand-replicating 33 argparse usage strings is a way
            # to introduce 33 chances to be subtly wrong (DIV-0025). What must match
            # is WHICH condition was reported and that both sides refused.
            def _last_error_line(text):
                for line in reversed([ln for ln in text.splitlines() if ln.strip()]):
                    return line.strip()
                return ""
            py_last = _last_error_line(py_text)
            rs_last = _last_error_line(rust_text)
            if not py_last or not rs_last:
                print("  FAIL error_line comparison: an arm printed no error line "
                      f"(python={py_last!r} rust={rs_last!r})")
                ok = False
            elif py_last != rs_last:
                print("  FAIL error line differs")
                print(f"  --- python: {py_last}")
                print(f"  --- rust  : {rs_last}")
                ok = False
            else:
                print(f"  error line PASS ({py_last!r})")
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

    for rel_path, floor in case.divergent_files:
        py_file = py_dir / rel_path
        rust_file = rust_dir / rel_path
        if not py_file.is_file() or not rust_file.is_file():
            print(f"  FAIL divergent file '{rel_path}': python_exists={py_file.is_file()} rust_exists={rust_file.is_file()}")
            ok = False
            continue
        py_rows = _numeric_rows(py_file.read_bytes())
        rust_rows = _numeric_rows(rust_file.read_bytes())
        if py_rows is None or rust_rows is None:
            print(f"  FAIL divergent file '{rel_path}': not parseable as numeric rows "
                  f"(python={py_rows is not None} rust={rust_rows is not None})")
            ok = False
            continue
        widest, where = 0.0, None
        compared = 0
        for r, (py_row, rust_row) in enumerate(zip(py_rows, rust_rows)):
            for c, (pv, rv) in enumerate(zip(py_row, rust_row)):
                if pv is None or rv is None:
                    # A label cell, not a numeric cell: comparing it as text is what
                    # compare_files is for, and there are none here.
                    continue
                compared += 1
                delta = abs(pv - rv)
                if delta > widest:
                    widest, where = delta, (r, c)
        if compared == 0:
            print(f"  FAIL divergent file '{rel_path}': no numeric cell was "
                  "compared, so the divergence cannot be pinned")
            ok = False
            continue
        if widest == 0.0:
            print(f"  FAIL divergent file '{rel_path}': the files are now IDENTICAL.")
            print("       The recorded divergence no longer holds. Retire the ledger")
            print("       entry, re-run the command contract tests, and remove this")
            print("       case from divergent_files in the SAME change.")
            ok = False
        elif widest < floor:
            print(f"  FAIL divergent file '{rel_path}': largest difference is {widest:g}, "
                  f"below the recorded floor {floor:g}.")
            print("       The divergence changed shape. Re-derive it and update the")
            print("       ledger entry and this floor together.")
            ok = False
        else:
            print(f"  divergent file '{rel_path}' PASS (largest difference {widest:g} "
                  f"at row {where[0]} col {where[1]}, floor {floor:g})")

    for rel_path in case.compare_files:
        py_file = py_dir / rel_path
        rust_file = rust_dir / rel_path
        if not py_file.is_file() or not rust_file.is_file():
            print(f"  FAIL file '{rel_path}': python_exists={py_file.is_file()} rust_exists={rust_file.is_file()}")
            ok = False
            continue
        py_bytes = py_file.read_bytes()
        rust_bytes = rust_file.read_bytes()
        if rel_path in case.gzip_files:
            import gzip as _gzip

            py_bytes = _gzip.decompress(py_bytes)
            rust_bytes = _gzip.decompress(rust_bytes)
        if case.normalize_paths:
            py_bytes = py_bytes.replace(str(py_dir).encode(), b"<SCRATCH_DIR>")
            rust_bytes = rust_bytes.replace(str(rust_dir).encode(), b"<SCRATCH_DIR>")
        # Two zero-byte files are byte-identical, so without this a command that
        # truncated its output identically in both arms would pass. A compared file
        # must contain something unless the case declares empty correct for that
        # specific file -- naming the file rather than the whole case, so the
        # exemption cannot quietly extend to an artifact nobody examined.
        if not py_bytes and rel_path not in case.allow_empty_files:
            print(f"  FAIL file '{rel_path}': both sides are empty (0 bytes)")
            print("       If an empty artifact is the correct result here, name the "
                  "file in allow_empty_files with the reason.")
            ok = False
            continue
        numeric = rel_path in case.numeric_files
        reason = ""
        if numeric:
            file_equal, reason = _numeric_table_equal(py_bytes, rust_bytes)
        else:
            file_equal = py_bytes == rust_bytes
        if file_equal:
            qualifier = "numeric cells" if numeric else f"byte-identical, {len(py_bytes)} bytes"
            print(f"  file '{rel_path}' PASS ({qualifier})")
        else:
            print(f"  FAIL file '{rel_path}': {reason or f'content differs ({len(py_bytes)} vs {len(rust_bytes)} bytes)'}")
            print(f"       {len(py_bytes)} vs {len(rust_bytes)} bytes")
            print(f"  --- python {rel_path} ---")
            print(py_bytes.decode("utf-8", errors="replace"))
            print(f"  --- rust {rel_path} ---")
            print(rust_bytes.decode("utf-8", errors="replace"))
            ok = False

    for rel_path in case.expect_files:
        py_file = py_dir / rel_path
        rust_file = rust_dir / rel_path
        py_size = py_file.stat().st_size if py_file.is_file() else None
        rust_size = rust_file.stat().st_size if rust_file.is_file() else None
        if py_size and rust_size:
            print(f"  file '{rel_path}' PASS (exists on both sides, {py_size}/{rust_size} bytes)")
        else:
            # `if py_size` treats a zero-byte file as absent, which is the intent:
            # a truncated deliverable is not a deliverable.
            print(f"  FAIL file '{rel_path}': python_size={py_size!r} rust_size={rust_size!r}")
            ok = False

    for rel_path in case.rust_expect_files:
        rust_file = rust_dir / rel_path
        rust_size = rust_file.stat().st_size if rust_file.is_file() else None
        if rust_size:
            print(f"  rust-only file '{rel_path}' PASS (exists, {rust_size} bytes)")
        else:
            print(f"  FAIL rust-only file '{rel_path}': rust_size={rust_size!r}")
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


def ensure_malformed_bam_fixtures() -> None:
    """Builds the degenerate-input fixtures the negative branch cases need.

    Three shapes, all in `verification/fixtures/`:

    - `empty.bam` (+`.bai`): a real BAM header copied from the synthetic
      panel with ZERO alignment records. This is not a corrupt file -- it is
      what a filtering or post-processing pipeline legitimately produces,
      and it is the input class that exposed DIV-0023 (upstream's
      `readsNVC` crashes with an unhandled `UnboundLocalError` on it).
    - `junk.bam`: a syntactically valid but empty gzip member. The
      container decompresses fine and then the reader hits EOF where a BAM
      header is required, so this exercises header handling specifically
      rather than gzip parsing.
    - `noidx.bam`: a normal BAM with its `.bai` removed, for the commands
      whose CLI requires an index sidecar.

    Generated rather than committed because `empty.bam` must share the
    synthetic panel's header exactly; both are deterministic, so a
    regeneration is byte-identical.
    """
    import pysam

    out_dir = REPO_ROOT / "verification" / "fixtures"
    empty = out_dir / "empty.bam"
    junk = out_dir / "junk.bam"
    noidx = out_dir / "noidx.bam"
    if empty.is_file() and junk.is_file() and noidx.is_file():
        return
    ensure_synthetic_fixtures()
    with pysam.AlignmentFile(SYNTHETIC_DIR / "pe.bam") as src:
        header = src.header
    if not empty.is_file():
        with pysam.AlignmentFile(empty, "wb", header=header):
            pass
        pysam.index(str(empty))
    if not junk.is_file():
        # A genuine BGZF member holding no alignment block at all: the
        # container parses and then the reader hits EOF where a BAM header
        # is required. Chosen over random bytes because it fails as a BAM
        # rather than as a gzip.
        import gzip as _gzip
        junk.write_bytes(_gzip.compress(b"", mtime=0))
    if not noidx.is_file():
        import shutil as _shutil
        _shutil.copy(SYNTHETIC_DIR / "pe.bam", noidx)
    # A BED12 whose every line is unparseable (too few columns, then a
    # non-integer coordinate). Upstream skips such lines individually with
    # a stderr note rather than aborting, so the command still succeeds
    # with an empty transcript set; the port must match that, and in
    # particular must not abort on the first bad line.
    malformed = out_dir / "malformed.bed12"
    if not malformed.is_file():
        malformed.write_text(
            "chr1\t0\t100\n"
            "chr1\t0\t100\tg\t0\t+\n"
            "chr1\tX\t100\tg\t0\t+\t0\t100\t0\t1\t100,\t0,\n"
        )


def _generator_cwd():
    """A scratch working directory for fixture generators.

    Upstream RSeQC writes `log.txt` into the process's current directory rather than
    next to its output (`geneBody_coverage.py:57` passes `log_file=Path("log.txt")`).
    Running the generators with `cwd=REPO_ROOT` therefore dropped a `log.txt` in the
    repository root -- an untracked file that looks like a stray artefact and has to be
    cleaned up by hand after every suite run. Confining generator output to a
    temporary directory keeps the checkout clean without hiding the behaviour, which
    is separately accounted for by NON_DELIVERABLE_UPSTREAM on the comparison side.
    """
    import atexit
    import shutil
    import tempfile

    global _GENERATOR_CWD
    if _GENERATOR_CWD is None or not Path(_GENERATOR_CWD).is_dir():
        _GENERATOR_CWD = tempfile.mkdtemp(prefix="rseqc-fixtures-")
        atexit.register(shutil.rmtree, _GENERATOR_CWD, True)
    return _GENERATOR_CWD


_GENERATOR_CWD = None


def ensure_bam_stat_fixture() -> None:
    fixture = REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam"
    if fixture.is_file():
        return
    generator = REPO_ROOT / "verification" / "fixtures" / "make_bam_stat_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture)], cwd=_generator_cwd(), check=True)


def ensure_mismatch_profile_fixture() -> None:
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    required = (fixture_dir / "mismatch_profile_basic.bam", fixture_dir / "mismatch_profile_basic.bam.bai")
    if all(path.is_file() for path in required):
        return
    generator = fixture_dir / "make_mismatch_profile_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture_dir / "mismatch_profile_basic.bam")], cwd=_generator_cwd(), check=True)


def ensure_deletion_profile_fixture() -> None:
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    required = (fixture_dir / "deletion_profile_basic.bam", fixture_dir / "deletion_profile_basic.bam.bai")
    if all(path.is_file() for path in required):
        return
    generator = fixture_dir / "make_deletion_profile_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture_dir / "deletion_profile_basic.bam")], cwd=_generator_cwd(), check=True)


def ensure_bam_stat_sam_fixture() -> None:
    # Same alignments as bam_stat_basic.bam, re-encoded as plain-text SAM
    # via pysam -- exercises open_alignments()'s SAM-text round-trip path
    # (DIV-0002/0004) against the exact same data as the .bam case.
    ensure_bam_stat_fixture()
    fixture = REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.sam"
    if fixture.is_file():
        return
    generator = REPO_ROOT / "verification" / "fixtures" / "make_bam_stat_sam_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture)], cwd=_generator_cwd(), check=True)


def ensure_bam_stat_cram_fixture() -> None:
    # Same alignments as bam_stat_basic.bam, re-encoded as CRAM (no
    # external reference -- htslib embeds it) via pysam -- exercises
    # open_alignments()'s CRAM round-trip path (DIV-0002/0004) against
    # the exact same data as the .bam case. This specific fixture
    # (including unmapped1, whose explicit MAPQ 0 doesn't survive
    # htslib's own CRAM write/read round-trip faithfully) caught a real
    # htslib/noodles-cram interop discrepancy -- see
    # rseqc_formats::fix_unmapped_missing_mapping_quality.
    ensure_bam_stat_fixture()
    fixture = REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.cram"
    if fixture.is_file():
        return
    generator = REPO_ROOT / "verification" / "fixtures" / "make_bam_stat_cram_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture)], cwd=_generator_cwd(), check=True)


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
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture_dir)], cwd=_generator_cwd(), check=True)


def _bam_stat_args(_scratch_dir: Path) -> list[str]:
    return ["-i", str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.bam")]


def _bam_stat_sam_fixture_path() -> str:
    return str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.sam")


def _bam_stat_sam_args(_scratch_dir: Path) -> list[str]:
    return ["-i", _bam_stat_sam_fixture_path()]


def _bam_stat_cram_fixture_path() -> str:
    return str(REPO_ROOT / "verification" / "fixtures" / "bam_stat_basic.cram")


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


def ensure_sc_bamstat_fixture() -> None:
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    required = (fixture_dir / "sc_bamstat_basic.bam", fixture_dir / "sc_bamstat_basic.bam.bai")
    if all(path.is_file() for path in required):
        return
    generator = fixture_dir / "make_sc_bamstat_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture_dir / "sc_bamstat_basic.bam")], cwd=_generator_cwd(), check=True)


def ensure_sc_editmatrix_fixture() -> None:
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    fixture = fixture_dir / "sc_editmatrix_basic.bam"
    if fixture.is_file():
        return
    generator = fixture_dir / "make_sc_editmatrix_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture)], cwd=_generator_cwd(), check=True)


def ensure_sc_seqqual_fixture() -> None:
    fixture = REPO_ROOT / "verification" / "fixtures" / "regression_sc_seqqual.fq"
    if not fixture.is_file():
        raise FileNotFoundError(f"missing committed fixture: {fixture}")


def ensure_sc_seqqual_gz_fixture() -> None:
    # Same content as regression_sc_seqqual.fq, gzip-compressed, for
    # exercising open_text_input's .gz dispatch. Regenerated with
    # Python's own gzip module (deterministic content; the gzip
    # CONTAINER's own mtime/filename header is not reproducible run-to-
    # run regardless, same as DIV-0007's bam2fq.py -c output -- callers
    # only need this file to decompress to the right content, not to be
    # byte-identical to any previous run of this generator).
    ensure_sc_seqqual_fixture()
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    src = fixture_dir / "regression_sc_seqqual.fq"
    dst = fixture_dir / "regression_sc_seqqual.fq.gz"
    if dst.is_file():
        return
    import gzip as _gzip

    with open(src, "rb") as f_in, open(dst, "wb") as raw, _gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f_out:
        f_out.write(f_in.read())


def ensure_sc_seqlogo_fixture() -> None:
    fixture = REPO_ROOT / "verification" / "fixtures" / "regression_sc_seqlogo.fa"
    if not fixture.is_file():
        raise FileNotFoundError(f"missing committed fixture: {fixture}")


def ensure_fpkm_uq_fixtures() -> None:
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    required = (
        fixture_dir / "mock_htseq_count.sh",
        fixture_dir / "regression_fpkm_uq_genes.info.txt",
        fixture_dir / "regression_fpkm_uq_dummy.bam",
        fixture_dir / "regression_fpkm_uq_dummy.gtf",
    )
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"missing committed fixture(s): {missing}")
    (fixture_dir / "mock_htseq_count.sh").chmod(0o755)


def _track_fixture(name: str) -> str:
    return str(REPO_ROOT / "verification" / "fixtures" / "track" / name)


def ensure_rpkm_saturation_fixture() -> None:
    fixture_dir = REPO_ROOT / "verification" / "fixtures"
    required = (
        fixture_dir / "rpkm_saturation_basic.bam",
        fixture_dir / "rpkm_saturation_basic.bam.bai",
        fixture_dir / "rpkm_saturation_model.bed12",
    )
    if all(path.is_file() for path in required):
        return
    generator = fixture_dir / "make_rpkm_saturation_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture_dir / "rpkm_saturation_basic.bam")], cwd=_generator_cwd(), check=True)


def ensure_rpkm_saturation_sam_fixture() -> None:
    # Same single alignment as rpkm_saturation_basic.bam, re-encoded as
    # plain-text SAM -- preserves the "exactly one qualifying alignment"
    # shuffle-order-independence property (see make_rpkm_saturation_
    # fixture.py's docstring).
    ensure_rpkm_saturation_fixture()
    fixture = REPO_ROOT / "verification" / "fixtures" / "rpkm_saturation_basic.sam"
    if fixture.is_file():
        return
    generator = REPO_ROOT / "verification" / "fixtures" / "make_rpkm_saturation_sam_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture)], cwd=_generator_cwd(), check=True)


def ensure_track_fixtures() -> None:
    """BigWig/WIG-family fixtures (bam2wig.py, geneBody_coverage2.py,
    normalize_bigwig.py, overlay_bigwig.py). Regenerated via pyBigWig +
    pysam if missing; verified reproducible/deterministic byte-for-byte
    before being relied on here."""
    track_dir = REPO_ROOT / "verification" / "fixtures" / "track"
    required = (
        track_dir / "track_chrom.sizes",
        track_dir / "track_signal.bw",
        track_dir / "track_signal2.bw",
        track_dir / "track_reads.bam",
        track_dir / "track_reads.bam.bai",
        track_dir / "track_model.bed12",
    )
    if all(path.is_file() for path in required):
        return
    generator = REPO_ROOT / "verification" / "fixtures" / "make_track_fixtures.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(track_dir)], cwd=_generator_cwd(), check=True)


def ensure_genebody_coverage_float_fixture() -> None:
    """Fixture for testing geneBody_coverage.py's float vs int duck-typing
    behavior. Generates a deterministic workload using the benchmarks/
    generate_workload.py script with a fixed seed."""
    fixture_dir = REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_float"
    required = (
        fixture_dir / "reads.bam",
        fixture_dir / "reads.bam.bai",
        fixture_dir / "model.bed12",
    )
    if all(path.is_file() for path in required):
        return
    fixture_dir.mkdir(parents=True, exist_ok=True)
    generator = REPO_ROOT / "benchmarks" / "generate_workload.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), "--size", "300", "--output-dir", str(fixture_dir)], cwd=_generator_cwd(), check=True)


def ensure_genebody_coverage_edge_fixture() -> None:
    """Edge-case fixture for geneBody_coverage.py CIGAR D/N operations.
    Verifies that deletion-only and skip-only spans appear in pileup and
    render as int 0 (visited), not float 0.0 (unvisited)."""
    fixture_dir = REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_edge"
    required = (
        fixture_dir / "edge.bam",
        fixture_dir / "edge.bam.bai",
        fixture_dir / "edge.bed12",
    )
    if all(path.is_file() for path in required):
        return
    fixture_dir.mkdir(parents=True, exist_ok=True)
    generator = REPO_ROOT / "verification" / "fixtures" / "make_genebody_coverage_edge_fixture.py"
    subprocess.run([str(ORACLE_PYTHON), str(generator), str(fixture_dir)], cwd=_generator_cwd(), check=True)


SYNTHETIC_DIR = REPO_ROOT / "verification" / "fixtures" / "synthetic"


# Real-data panel redirection (testing.md section 11 / T4).
#
# Setting RSEQC_REAL_DATA=<dir> re-points the *bulk* fixtures at the real
# aligned panel, so the same 90 cases -- same comparators, same tolerances,
# same expected divergences -- run against real sequenced reads instead of
# generated ones. The point is not a parallel second harness that could drift
# from the first; it is the existing matrix, unchanged, fed different input.
#
# Only the bulk fixtures have real equivalents. The single-cell and BigWig
# fixtures do not, so cases needing them are SKIPPED with a printed reason
# rather than silently falling back to synthetic -- silently falling back is
# how a "real-data validation" ends up validating nothing real.
REAL_DATA_DIR_ENV = "RSEQC_REAL_DATA"
REAL_DATA_DIR = os.environ.get(REAL_DATA_DIR_ENV) or None
# Every synthetic fixture a case may ask for. Used to detect, per case, whether
# it needs something the real panel does not provide.
REAL_FIXTURE_NAMES = (
    "pe.bam", "se.bam", "sc.bam", "pe_placed.bam", "se_placed.bam",
    "model.bed12", "model_hdr.bed12", "sig1.bw", "sig2.bw",
    "barcodes.fq", "barcodes.fa", "barcodes.fq.gz",
    "reads.fq", "reads.fq.gz", "reads.fa", "mrna.fa", "genome.fa",
    "chrom.sizes", "genes.info.txt", "htseq_counts.txt",
)
# Maps a synthetic fixture name to the real panel file that replaces it.
#
# The panel directory IS a fixture directory: it uses the same names, so a
# case's expected outputs (which are derived from the BAM's stem, e.g.
# "pe.tin.xls") keep matching. An earlier revision renamed them to
# real_pe.bam etc. and then, because the lookup was a plain `is_file()` that
# returned None on a miss, silently fell back to the SYNTHETIC fixtures -- so a
# run that printed "real-data panel:" in its banner was in fact validating
# nothing real. `require_real_panel` below exists to make that impossible.
REAL_EQUIVALENTS = {
    "pe.bam": "pe.bam",
    "se.bam": "se.bam",
    "model.bed12": "model.bed12",
    "chrom.sizes": "chrom.sizes",
}


def require_real_panel() -> None:
    """Refuses to run under RSEQC_REAL_DATA unless every required file is present.

    Without this, a typo'd panel path or a misnamed file produces a run whose
    header says "real-data panel: ..." while every case silently reads the
    synthetic fixtures. That is the worst possible failure for this harness:
    it looks like a pass.
    """
    if not REAL_DATA_DIR:
        return
    root = Path(REAL_DATA_DIR)
    if not root.is_dir():
        sys.exit(f"error: RSEQC_REAL_DATA={root} is not a directory -- refusing to fall "
                 f"back to the synthetic fixtures, which would make this run validate "
                 f"nothing real")
    missing = [f"{k} (expected {v})" for k, v in sorted(REAL_EQUIVALENTS.items())
               if not (root / v).is_file()]
    if missing:
        joined = "; ".join(missing)
        sys.exit(
            f"error: real-data panel {root} is missing required fixture(s): {joined}\n"
            f"       Refusing to fall back to the synthetic fixtures: a run that "
            f"silently used them would report 'real-data' results it never produced."
        )


def _real_path(name: str) -> str | None:
    """Real-panel path for a synthetic fixture name, or None if there isn't one."""
    if not REAL_DATA_DIR:
        return None
    replacement = REAL_EQUIVALENTS.get(name)
    if replacement is None:
        return None
    candidate = Path(REAL_DATA_DIR) / replacement
    return str(candidate) if candidate.is_file() else None


def _synthetic(name: str) -> str:
    real = _real_path(name)
    if real is not None:
        return real
    return str(SYNTHETIC_DIR / name)


def _real_data_gaps(case: Case) -> list[str]:
    """Which of a case's fixtures have no real-panel equivalent.

    Non-empty means the case cannot be run under RSEQC_REAL_DATA without
    quietly falling back to synthetic input, so it is skipped instead.

    The arg builders are CALLED, not stringified. An earlier revision did
    `str(case.py_args)`, which for a lambda is `"<function <lambda> at 0x...>"`
    and contains no fixture names at all -- so this function never detected a
    single gap and every case ran, the ones needing absent fixtures included.
    They then failed for a data-shape reason (a synthetic 8 KB BAM against a
    real 3,000-transcript model: zero exonic fragments) which read like a port
    defect and was not one.
    """
    if not REAL_DATA_DIR:
        return []
    probe = Path(os.environ.get("TMPDIR", "/tmp")) / "rseqc_realdata_probe"
    try:
        argv = list(case.py_args(probe)) + list(case.rust_args(probe))
    except Exception:  # a builder that cannot run without real inputs
        return ["<args unavailable>"]
    # Match on the resolved fixture PATH, not on a substring of the argument.
    # A substring test flagged read_hexamer's cases as needing "reads.fa" because
    # their argument contains "regression_hexamer_reads.fa" -- a different file,
    # from a different fixture family, that the real panel is not expected to
    # provide. Two fixtures skipped for a bogus reason is worse than one missed.
    # Multi-input commands take a COMMA-JOINED list in a single argument
    # (`"-i", pe + "," + se + "," + sc`), so each argument has to be split
    # before comparing. Without the split, no multi-BAM case was ever
    # recognised as needing a fixture, and genebody_coverage's three-sample
    # case ran the real pe/se against the SYNTHETIC sc.bam -- which has no
    # coverage in the real model's window, so both arms raised the same
    # ZeroDivisionError and the row was reported as a failure of both.
    tokens = [tok for a in argv for tok in a.split(",") if tok]
    missing = []
    for name in sorted(REAL_FIXTURE_NAMES):
        if _real_path(name) is not None:
            continue
        synthetic_path = str(SYNTHETIC_DIR / name)
        if any(tok == synthetic_path for tok in tokens):
            missing.append(name)
    return missing


def ensure_synthetic_fixtures() -> None:
    """Committed output of `verification/synthetic_data.py --seed 1 --size 40`
    (minus the machine-specific bams.txt): multi-chromosome, spliced/clipped/
    indel reads with MD tags, flags of every kind, PE+SE+single-cell BAMs,
    BigWigs, FASTA/FASTQ. Regenerated byte-identically if missing."""
    if (SYNTHETIC_DIR / "pe.bam").is_file():
        return
    subprocess.run(
        [str(ORACLE_PYTHON), str(REPO_ROOT / "verification" / "synthetic_data.py"),
         "--seed", "1", "--size", "40", "--output-dir", str(SYNTHETIC_DIR)],
        check=True,
    )
    (SYNTHETIC_DIR / "bams.txt").unlink(missing_ok=True)


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
        name="bam_stat_sam_text",
        # DIV-0002/0004: SAM-text input support, closed for bam_stat.py via
        # rseqc_formats::open_alignments (round-trips SAM-text records
        # through an in-memory BAM byte buffer). Confirmed via live diff
        # against real upstream that pysam's Samfile(path, 'rb') succeeds
        # even for genuine plain-text SAM content (htslib auto-detects,
        # ignoring the 'b' mode hint) -- so upstream ALSO prints "Load BAM
        # file" for .sam input, not "Load SAM file"; comparing full stderr
        # here (not just stdout) pins that behavior down.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="bam_stat.py",
        rust_bin="bam_stat",
        py_args=_bam_stat_sam_args,
        rust_args=_bam_stat_sam_args,
        compare_stream="both",
        stream_format="exact",
        required_labels=("Total records", "Unmapped reads", "Read-1"),
    ),
    Case(
        name="read_NVC_basic",
        # Reuses bam_stat_basic.bam -- read_NVC.py only needs a BAM with
        # some real sequence content, no need for a dedicated fixture.
        # Found by this case (before it was formalized): read_NVC.py's
        # port never generated an R script or invoked Rscript at all
        # (missing --skip-plot/--rscript flags entirely, plus invented
        # stderr messages instead of upstream's real ones) -- a genuine
        # DIV-0005 instance, now fixed. --skip-plot on both sides here
        # since real Rscript execution is separately spot-checked
        # manually (both sides produce a real, valid, byte-identical-
        # sized single-page PDF via the actual Rscript binary), not
        # something this harness needs to assert every run.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_NVC.py",
        rust_bin="read_NVC",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.NVC.xls", "out.NVC_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_NVC_sam_text",
        # DIV-0002/0004: SAM-text input, closed for read_NVC.py via
        # open_alignments.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="read_NVC.py",
        rust_bin="read_NVC",
        py_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.NVC.xls", "out.NVC_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_NVC_cram",
        # DIV-0002/0004 CRAM remnant, closed for read_NVC.py via
        # open_alignments. Found and fixed a real bug via this exact
        # case: unmapped1's explicit MAPQ 0 doesn't survive htslib's
        # own CRAM write/read round-trip faithfully -- see
        # rseqc_formats::fix_unmapped_missing_mapping_quality. stderr is
        # NOT compared: htslib always attempts to load a .crai index
        # for CRAM input, even for pure sequential access, and prints
        # "[E::cram_index_load] Could not retrieve index file ..." to
        # stderr when none exists (this fixture has none, by design --
        # this port's CRAM support never needs an index). Environment/
        # library noise, not a correctness signal; the .xls/.r file
        # outputs are the real check here.
        ensure_fixture=ensure_bam_stat_cram_fixture,
        py_script="read_NVC.py",
        rust_bin="read_NVC",
        py_args=lambda scratch_dir: ["-i", _bam_stat_cram_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _bam_stat_cram_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="none",
        compare_files=("out.NVC.xls", "out.NVC_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_NVC_with_nx",
        # Same fixture, -x/--nx branch: the R script's total/ym/yn
        # expressions and legend include the N/X count series too.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_NVC.py",
        rust_bin="read_NVC",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-x", "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-x", "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.NVC.xls", "out.NVC_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_GC_basic",
        # Found by this case (before it was formalized): read_GC.py's
        # port never generated an R script or invoked Rscript at all
        # (missing --skip-plot/--rscript flags entirely, plus invented
        # stderr messages "GC table written to: ..."/"R script written
        # to: ..." instead of upstream's real "Read BAM file ...  Done"/
        # "writing GC content ..."/"writing R script ..." -- the exact
        # same DIV-0005 pattern read_NVC.py had). --skip-plot on both
        # sides here; real Rscript execution (both sides produce a
        # real, valid PDF via the actual Rscript binary) spot-checked
        # manually, not asserted every run.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_GC.py",
        rust_bin="read_GC",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.GC.xls", "out.GC_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_GC_sam_text",
        # DIV-0002/0004: SAM-text input, closed for read_GC.py via
        # open_alignments.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="read_GC.py",
        rust_bin="read_GC",
        py_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.GC.xls", "out.GC_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_duplication_basic",
        # Found by this case (before it was formalized): (1) same
        # DIV-0005 pattern as read_NVC.py/read_GC.py -- no --skip-plot/
        # --rscript, no Rscript invocation, invented stderr messages;
        # (2) a genuine, separate bug in render_dup_r_script: the
        # pdf(...) line used double quotes, but upstream's
        # `print("pdf(\'%s\')" % ...)` -- the `\'` is an ESCAPED SINGLE
        # QUOTE inside the double-quoted Python literal -- renders
        # single quotes. An earlier "cross-checked" unit test asserted
        # the wrong (double-quote) text; only a live diff against the
        # real upstream CLI caught it. --skip-plot on both sides here;
        # real Rscript execution (both sides produce a real, valid PDF)
        # spot-checked manually.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_duplication.py",
        rust_bin="read_duplication",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.pos.DupRate.xls", "out.seq.DupRate.xls", "out.DupRate_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_duplication_sam_text",
        # DIV-0002/0004: SAM-text input, closed for read_duplication.py via
        # open_alignments.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="read_duplication.py",
        rust_bin="read_duplication",
        py_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.pos.DupRate.xls", "out.seq.DupRate.xls", "out.DupRate_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="read_quality_basic",
        # Found by this case (before it was formalized): read_quality.py's
        # port never accepted --skip-plot/--rscript or invoked Rscript
        # (a DIV-0005 instance) and was missing upstream's "Read BAM
        # file ...  Done" progress line. --skip-plot on both sides here;
        # real Rscript execution (both sides exit 1 and produce the same
        # 2 PDFs -- an upstream-inherent R error with real data, not a
        # divergence) spot-checked manually, not asserted every run.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="read_quality.py",
        rust_bin="read_quality",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        # The .qual.r script embeds each side's own absolute scratch-
        # directory path in its pdf('...') line -- normalize it away,
        # see normalize_paths' docstring above.
        compare_files=("out.qual.r",),
        normalize_paths=True,
    ),
    Case(
        name="read_quality_sam_text",
        # DIV-0002/0004: SAM-text input, closed for read_quality.py via
        # open_alignments.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="read_quality.py",
        rust_bin="read_quality",
        py_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
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
        name="RNA_fragment_size_insufficient_fragments",
        # -n above the fixture's real fragment count exercises the
        # count < ncut branch (mean/median/std left as upstream's zero
        # placeholders rather than computed) -- not covered by the
        # equals_cigar case above, which always clears its own -n 1.
        ensure_fixture=ensure_regression_fixtures,
        py_script="RNA_fragment_size.py",
        rust_bin="RNA_fragment_size",
        py_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_rna_equals.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-n",
            "100",
            "-o",
            str(scratch_dir / "out.tsv"),
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            _regression_fixture("regression_rna_equals.bam"),
            "-r",
            _regression_fixture("regression_single_exon.bed12"),
            "-n",
            "100",
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
        name="geneBody_coverage_float_format",
        ensure_fixture=ensure_genebody_coverage_float_fixture,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda scratch_dir: [
            "-i",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_float" / "reads.bam"),
            "-r",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_float" / "model.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_float" / "reads.bam"),
            "-r",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_float" / "model.bed12"),
            "-o",
            str(scratch_dir / "out"),
            "--skip-plot",
        ],
        compare_stream="none",
        compare_files=("out.geneBodyCoverage.txt",),
        # Use exact byte comparison (not numeric) to verify float formatting
        # (e.g., "15.0" in Rust output matches Python's "15.0", not "15")
    ),
    Case(
        name="genebody_coverage_edge_cases",
        ensure_fixture=ensure_genebody_coverage_edge_fixture,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda scratch_dir: [
            "-i",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_edge" / "edge.bam"),
            "-r",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_edge" / "edge.bed12"),
            "-o",
            str(scratch_dir / "gb"),
            "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_edge" / "edge.bam"),
            "-r",
            str(REPO_ROOT / "verification" / "fixtures" / "genebody_coverage_edge" / "edge.bed12"),
            "-o",
            str(scratch_dir / "gb"),
            "--skip-plot",
        ],
        compare_stream="none",
        compare_files=("gb.geneBodyCoverage.txt",),
        # Verifies D (deletion) and N (skip) CIGAR spans mark positions as
        # visited (int 0) not unvisited (float 0.0), even when no reads pass filters
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
        name="tin_subtract_background",
        # -s/--subtract-background branch, not exercised by
        # tin_pair_overlap above (only the default no-background-
        # subtraction path). Verified via live diff before formalizing.
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
            "-s",
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
            "-s",
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
        # not "0.0" -- see crates/commands/src/clipping_profile.rs. Also
        # covers the later DIV-0005 fix: missing --skip-plot/--rscript,
        # missing "Load BAM file ...  Done" progress, and upstream's own
        # "Totoal reads used: N" typo (preserved exactly, not "Total").
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="clipping_profile.py",
        rust_bin="clipping_profile",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.clipping_profile.xls", "out.clipping_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="clipping_profile_pe",
        # Same fixture, PE branch: exercises the "Totoal read-1 used"/
        # "Totoal read-2 used" TWO-SEPARATE-LINES message shape (the
        # port previously combined these into one invented line).
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="clipping_profile.py",
        rust_bin="clipping_profile",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "PE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "PE", "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.clipping_profile.xls", "out.clipping_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="clipping_profile_sam_text",
        # DIV-0002/0004: SAM-text input, closed for clipping_profile.py
        # via open_alignments.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="clipping_profile.py",
        rust_bin="clipping_profile",
        py_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.clipping_profile.xls", "out.clipping_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="insertion_profile_basic",
        # Same DIV-0005 fix pattern as clipping_profile.py (identical
        # code structure upstream): missing --skip-plot/--rscript,
        # missing "Load BAM file ...  Done" progress, upstream's
        # "Totoal reads used: N" typo preserved exactly.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="insertion_profile.py",
        rust_bin="insertion_profile",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.insertion_profile.xls", "out.insertion_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="insertion_profile_sam_text",
        # DIV-0002/0004: SAM-text input, closed for insertion_profile.py
        # via open_alignments.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="insertion_profile.py",
        rust_bin="insertion_profile",
        py_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "SE", "--skip-plot",
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.insertion_profile.xls", "out.insertion_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="insertion_profile_pe",
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="insertion_profile.py",
        rust_bin="insertion_profile",
        py_args=lambda scratch_dir: [
            "-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "PE", "--skip-plot",
        ],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "PE", "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.insertion_profile.xls", "out.insertion_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="bam2fq_single_end",
        # Found by this case (before it was formalized): the CLI's stderr
        # was a raw `{counts:?}` Debug dump instead of upstream's real
        # "Convert BAM/SAM file into FASTQ format ...", "Done", "read
        # count: N" lines (ParseBAM.bam2fq() in qcmodule/SAM.py) -- fixed,
        # now comparing stderr exactly instead of skipping it.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="bam2fq.py",
        rust_bin="bam2fq",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.fastq",),
    ),
    Case(
        name="bam2fq_sam_text",
        # DIV-0002/0004: SAM-text input, closed for bam2fq.py via
        # open_alignments.
        ensure_fixture=ensure_bam_stat_sam_fixture,
        py_script="bam2fq.py",
        rust_bin="bam2fq",
        py_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "-s"],
        rust_args=lambda scratch_dir: ["-i", _bam_stat_sam_fixture_path(), "-o", str(scratch_dir / "out"), "-s"],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.fastq",),
    ),
    Case(
        name="bam2fq_compress",
        # DIV-0007: -c/--compress, closed. Uses flate2 at compression
        # level 9 (matches Python's gzip.open default compresslevel).
        # The .gz container's own bytes are never comparable (Python's
        # gzip.open stamps the current wall time into the header on
        # every run) -- gzip_files decompresses both sides first so only
        # the actual FASTQ content is compared.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="bam2fq.py",
        rust_bin="bam2fq",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "-c"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-o", str(scratch_dir / "out"), "-s", "-c"],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.fastq.gz",),
        gzip_files=("out.fastq.gz",),
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
        name="infer_experiment_sam_text",
        # DIV-0002/0004: SAM-text input, closed for infer_experiment.py via
        # open_alignments.
        ensure_fixture=lambda: (ensure_regression_fixtures(), ensure_bam_stat_sam_fixture()),
        py_script="infer_experiment.py",
        rust_bin="infer_experiment",
        py_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
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
        # oracle-verified case, not a dummy. DIV-0005 fix: this command
        # DOES have --skip-plot/--rscript upstream (an earlier version
        # of this comment claimed otherwise -- wrong; the Rust port
        # simply hadn't wired the flag up yet at the time). Now passed
        # on both sides, so stdout is directly comparable again instead
        # of needing to dodge Rscript's own "null device"/"1" console
        # noise.
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="deletion_profile.py",
        rust_bin="deletion_profile",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stdout",
        stream_format="exact",
        compare_files=("out.deletion_profile.txt", "out.deletion_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="deletion_profile_with_deletion",
        # Exercises the actual deletion-counting logic (a genuine 2-base
        # deletion via CIGAR 5M2D5M) -- the other case only covers the
        # zero-deletions early-exit path.
        ensure_fixture=ensure_deletion_profile_fixture,
        py_script="deletion_profile.py",
        rust_bin="deletion_profile",
        py_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "deletion_profile_basic.bam"),
            "-l", "10", "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "deletion_profile_basic.bam"),
            "-l", "10", "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        compare_stream="stdout",
        stream_format="exact",
        compare_files=("out.deletion_profile.txt", "out.deletion_profile.r"),
        normalize_paths=True,
    ),
    Case(
        name="mismatch_profile_no_mismatches",
        # DIV-0005 fix: added --skip-plot/--rscript + real Rscript
        # invocation for the has-mismatches path (this specific case
        # stays on the "No mismatches found" early exit, which bypasses
        # Rscript on BOTH sides via upstream's own bare `sys.exit()`
        # inside mismatchProfile -- --skip-plot is passed anyway for
        # consistency/robustness, not because this path needs it).
        ensure_fixture=ensure_bam_stat_fixture,
        py_script="mismatch_profile.py",
        rust_bin="mismatch_profile",
        py_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out"), "--skip-plot"],
        rust_args=lambda scratch_dir: ["-i", _nvc_fixture_path(), "-l", "20", "-o", str(scratch_dir / "out"), "--skip-plot"],
        compare_stream="stdout",
        stream_format="exact",
        compare_files=("out.mismatch_profile.xls", "out.mismatch_profile.r"),
        # Upstream's mismatchProfile calls sys.exit() on the no-mismatches path,
        # before the R-script writer runs, so BOTH arms leave a zero-byte .r file.
        # The stream comparison is the substantive check here: both arms printed
        # the same "No mismatches found" text.
        allow_empty_files=("out.mismatch_profile.r",),
    ),
    Case(
        name="mismatch_profile_with_mismatch",
        # Exercises the actual mismatch-counting logic (a genuine A2C
        # mismatch at read position 5, real MD/NM tags) -- the other
        # case only covers the zero-mismatches early-exit path.
        ensure_fixture=ensure_mismatch_profile_fixture,
        py_script="mismatch_profile.py",
        rust_bin="mismatch_profile",
        py_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "mismatch_profile_basic.bam"),
            "-l", "10", "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "mismatch_profile_basic.bam"),
            "-l", "10", "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.mismatch_profile.xls", "out.mismatch_profile.r"),
        normalize_paths=True,
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
        # An alignment with no junctions has no rows to tabulate, and both
        # arms exit before writing the R script, so both artifacts are
        # legitimately zero bytes. The stream comparison is the substantive
        # check: both arms printed the same zero-junction report.
        allow_empty_files=("out.junction.xls", "out.junction_plot.r"),
    ),
    Case(
        name="junction_annotation_sam_text",
        # DIV-0002/0004: SAM-text input, closed for junction_annotation.py
        # via open_alignments.
        ensure_fixture=lambda: (ensure_regression_fixtures(), ensure_bam_stat_sam_fixture()),
        py_script="junction_annotation.py",
        rust_bin="junction_annotation",
        py_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        compare_stream="both",
        stream_format="exact",
        compare_files=("out.junction.xls", "out.junction_plot.r"),
        # An alignment with no junctions has no rows to tabulate, and both
        # arms exit before writing the R script, so both artifacts are
        # legitimately zero bytes. The stream comparison is the substantive
        # check: both arms printed the same zero-junction report.
        allow_empty_files=("out.junction.xls", "out.junction_plot.r"),
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
        name="inner_distance_sam_text",
        # DIV-0002/0004: SAM-text input, closed for inner_distance.py via
        # open_alignments.
        ensure_fixture=lambda: (ensure_regression_fixtures(), ensure_bam_stat_sam_fixture()),
        py_script="inner_distance.py",
        rust_bin="inner_distance",
        py_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
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
        name="read_distribution_sam_text",
        # DIV-0002/0004: SAM-text input, closed for read_distribution.py
        # via open_alignments.
        ensure_fixture=lambda: (ensure_regression_fixtures(), ensure_bam_stat_sam_fixture()),
        py_script="read_distribution.py",
        rust_bin="read_distribution",
        py_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(), "-r", _regression_fixture("regression_single_exon.bed12"),
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
        name="divide_bam_index",
        # DIV-0006: --index, closed. Uses regression_overlap_pair.bam
        # (already known strictly coordinate-sortable/indexable by real
        # htslib -- confirmed by its own pre-existing .bai) rather than
        # bam_stat_basic.bam, which genuinely fails to index under real
        # samtools/htslib: it has a NO_COOR (unmapped) read in the MIDDLE
        # of the file rather than trailing at the end, violating BAM's
        # coordinate-sort convention that htslib's indexer requires
        # (confirmed via a live `pysam.index()` probe -- a pre-existing
        # fixture limitation unrelated to this port's own code). .bai
        # files are expect_files, not compare_files: noodles' and
        # htslib's binning-index encoders produce different bytes for
        # the same logical index (confirmed non-byte-comparable via a
        # live diff), but both are independently verified real,
        # functionally-correct, htslib-readable indexes via a
        # pysam.fetch() probe during development.
        ensure_fixture=ensure_regression_fixtures,
        py_script="divide_bam.py",
        rust_bin="divide_bam",
        py_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_overlap_pair.bam"),
            "-n", "2", "-o", str(scratch_dir / "out"), "--seed", "42", "--index",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_overlap_pair.bam"),
            "-n", "2", "-o", str(scratch_dir / "out"), "--seed", "42", "--index",
        ],
        compare_stream="stderr",
        stream_format="exact",
        expect_files=("out_0.bam.bai", "out_1.bam.bai"),
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
        name="read_hexamer_refs_and_file_output",
        # Previously untested: -r/-g (optional reference-genome/
        # reference-transcript columns, `collect_inputs`) and -o (file
        # output instead of stdout, including the "Created: <path>"
        # stderr line). Deliberately passes the SAME fixture file for
        # both -r and -g to also exercise unique_display_name's
        # collision fallback (bare filename already taken by -r, so -g's
        # column name falls back to the full path) -- previously only
        # unit-tested in isolation, never live-diffed against real
        # upstream's own `unique_display_name`.
        ensure_fixture=ensure_hexamer_fixtures,
        py_script="read_hexamer.py",
        rust_bin="read_hexamer",
        py_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_hexamer_reads.fa"),
            "-r", _regression_fixture("regression_hexamer_ref.fa"),
            "-g", _regression_fixture("regression_hexamer_ref.fa"),
            "-o", str(scratch_dir / "out.tsv"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _regression_fixture("regression_hexamer_reads.fa"),
            "-r", _regression_fixture("regression_hexamer_ref.fa"),
            "-g", _regression_fixture("regression_hexamer_ref.fa"),
            "-o", str(scratch_dir / "out.tsv"),
        ],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.tsv",),
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
    Case(
        name="junction_saturation_sam_text",
        # DIV-0002/0004: SAM-text input, closed for junction_saturation.py
        # via open_alignments. Reuses bam_stat_basic.sam (no spliced
        # reads against this BED) rather than a dedicated .sam conversion
        # of regression_fpkm_fetch_span.bam -- the SAM-text round-trip
        # risk is data-independent (already exercised by 8 other
        # commands' cases), so this only needs to prove the code path
        # runs and matches real upstream, not re-exercise the splice
        # counting logic itself (already covered by junction_saturation_
        # basic on the .bam fixture).
        ensure_fixture=lambda: (ensure_regression_fixtures(), ensure_bam_stat_sam_fixture()),
        py_script="junction_saturation.py",
        rust_bin="junction_saturation",
        py_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(),
            "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _bam_stat_sam_fixture_path(),
            "-r", _regression_fixture("regression_single_exon.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        compare_stream="both",
        stream_format="exact",
        compare_files=("out.junctionSaturation_plot.r",),
        normalize_paths=True,
    ),
    Case(
        name="sc_bamstat_basic",
        # stdout only: upstream's stderr is a fully timestamped
        # `logging` trail (DIV-0019, open -- byte-exact comparison is
        # impossible in principle regardless of message content because
        # of the timestamp prefix). The stdout report -- the actual data
        # output -- IS verified byte-identical here, exercising the
        # RE-tag if/elif/else dead-code bug, sense/antisense/other,
        # chrM detection, and CB/UMI presence tallies all at once. Both
        # sides run with cwd=their own scratch dir since upstream writes
        # `All_reads_uniqID.txt`/`confident_reads_uniqID.txt` into the
        # current directory as an undocumented side effect (DIV-0018,
        # accepted -- not replicated by the Rust port).
        ensure_fixture=ensure_sc_bamstat_fixture,
        py_script="sc_bamStat.py",
        rust_bin="sc_bamStat",
        py_args=lambda scratch_dir: ["-i", _regression_fixture("sc_bamstat_basic.bam")],
        rust_args=lambda scratch_dir: ["-i", _regression_fixture("sc_bamstat_basic.bam")],
        compare_stream="stdout",
        stream_format="exact",
    ),
    Case(
        name="sc_bamstat_chrm_custom",
        # Tests --chrM-id flag with a non-existent contig name, exercising
        # the mitochondrial-read detection code path with a different result
        # (the single chrM read in the fixture is not recognized as
        # mitochondrial, so mitochondrial counts are 0).
        ensure_fixture=ensure_sc_bamstat_fixture,
        py_script="sc_bamStat.py",
        rust_bin="sc_bamStat",
        py_args=lambda scratch_dir: ["-i", _regression_fixture("sc_bamstat_basic.bam"), "--chrM-id", "invalidname"],
        rust_args=lambda scratch_dir: ["-i", _regression_fixture("sc_bamstat_basic.bam"), "--chrM-id", "invalidname"],
        compare_stream="stdout",
        stream_format="exact",
    ),
    Case(
        name="sc_editmatrix_basic",
        # Found by this case (before it was formalized): the edit-count
        # CSV's float/int cell dtype is a real `pandas.DataFrame.
        # from_dict(...).fillna(0)` quirk, but it's decided PER COLUMN
        # here (no transpose in this command's pipeline) -- a column
        # stays plain-integer only if that position had an entry for
        # every substitution key seen anywhere. The previous code always
        # appended ".0" to every cell; this fixture's UMI matrix has
        # exactly one (position, substitution) entry -- a trivially
        # dense single column -- so it now renders "1", not "1.0".
        # --skip-heatmap avoids the Rscript/pheatmap dependency entirely
        # (pure data-file comparison, not the plotting path).
        ensure_fixture=ensure_sc_editmatrix_fixture,
        py_script="sc_editMatrix.py",
        rust_bin="sc_editMatrix",
        py_args=lambda scratch_dir: ["-i", _regression_fixture("sc_editmatrix_basic.bam"), "-o", str(scratch_dir / "out"), "--skip-heatmap"],
        rust_args=lambda scratch_dir: ["-i", _regression_fixture("sc_editmatrix_basic.bam"), "-o", str(scratch_dir / "out"), "--skip-heatmap"],
        compare_stream="none",
        compare_files=("out.CB_edits_count.csv", "out.CB_freq.tsv", "out.UMI_edits_count.csv", "out.UMI_freq.tsv"),
    ),
    Case(
        name="sc_editmatrix_limit",
        # Previously untested: --limit (scbam.barcode_edits' own
        # "total_alignments >= limit" early-stop). sc_editmatrix_basic.bam
        # has 5 records (r1..r5); --limit 3 stops after r3, so r4/r5's
        # edits (including r5's UMI 3-way edit) must NOT appear in any
        # output file.
        ensure_fixture=ensure_sc_editmatrix_fixture,
        py_script="sc_editMatrix.py",
        rust_bin="sc_editMatrix",
        py_args=lambda scratch_dir: ["-i", _regression_fixture("sc_editmatrix_basic.bam"), "-o", str(scratch_dir / "out"), "--skip-heatmap", "--limit", "3"],
        rust_args=lambda scratch_dir: ["-i", _regression_fixture("sc_editmatrix_basic.bam"), "-o", str(scratch_dir / "out"), "--skip-heatmap", "--limit", "3"],
        compare_stream="none",
        compare_files=("out.CB_edits_count.csv", "out.CB_freq.tsv", "out.UMI_edits_count.csv", "out.UMI_freq.tsv"),
    ),
    Case(
        name="sc_seqqual_basic",
        # Found by this case: unlike sc_editMatrix.py, this command's
        # pipeline DOES transpose the matrix before `to_csv`, which
        # forces pandas to upcast EVERY column to a common dtype when
        # any column needed NaN-filling -- so the float/int decision is
        # GLOBAL (one flag for the whole matrix), not per column. This
        # fixture is constructed so every read cycle observes the exact
        # same quality-score set (fully dense) -- the whole count matrix
        # should stay plain-integer, which the previous code got wrong
        # (always appended ".0" unconditionally).
        ensure_fixture=ensure_sc_seqqual_fixture,
        py_script="sc_seqQual.py",
        rust_bin="sc_seqQual",
        py_args=lambda scratch_dir: ["-i", _regression_fixture("regression_sc_seqqual.fq"), "-o", str(scratch_dir / "out"), "--skip-heatmap"],
        rust_args=lambda scratch_dir: ["-i", _regression_fixture("regression_sc_seqqual.fq"), "-o", str(scratch_dir / "out"), "--skip-heatmap"],
        compare_stream="none",
        compare_files=("out.qual_count.csv", "out.qual_percent.csv"),
    ),
    Case(
        name="sc_seqqual_gzip_input",
        # Compressed FASTQ input support: open_text_input dispatches
        # .fq.gz to a gzip decoder, matching upstream's
        # qcmodule.ireader.nopen. Same fixture content as
        # sc_seqqual_basic, just gzip-compressed.
        ensure_fixture=ensure_sc_seqqual_gz_fixture,
        py_script="sc_seqQual.py",
        rust_bin="sc_seqQual",
        py_args=lambda scratch_dir: ["-i", _regression_fixture("regression_sc_seqqual.fq.gz"), "-o", str(scratch_dir / "out"), "--skip-heatmap"],
        rust_args=lambda scratch_dir: ["-i", _regression_fixture("regression_sc_seqqual.fq.gz"), "-o", str(scratch_dir / "out"), "--skip-heatmap"],
        compare_stream="none",
        compare_files=("out.qual_count.csv", "out.qual_percent.csv"),
    ),
    Case(
        name="sc_seqlogo_basic",
        # Same GLOBAL-dtype pandas quirk as sc_seqQual.py (seq2countMat
        # also transposes once before returning) -- this fixture rotates
        # a single edited base through every position across 3 short
        # sequences so every position observes both bases (fully dense),
        # exercising the same previously-broken "always append .0" bug.
        # Asymmetric exit contract (DIV-0016): the pinned upstream
        # environment itself cannot render logos -- its own
        # logomaker+pandas crash (`TypeError: Invalid value
        # '[-0.5 -0.5 -0.5]' for dtype 'int64'`) fires after the CSV is
        # written, before the logo step -- so python is expected to exit
        # 1 (its own bug), while this port is expected to exit 0 AND
        # write the two real PDF logos (rust_expect_files; their content
        # is structurally verified by the render crate's unit tests, and
        # the DNA-sequence data behind them is the byte-identical CSV).
        # The CSV is written by both sides before the crash, so it
        # remains a normal byte-identical comparison file.
        ensure_fixture=ensure_sc_seqlogo_fixture,
        py_script="sc_seqLogo.py",
        rust_bin="sc_seqLogo",
        py_args=lambda scratch_dir: ["-i", _regression_fixture("regression_sc_seqlogo.fa"), "-o", str(scratch_dir / "out"), "--iformat", "fa"],
        rust_args=lambda scratch_dir: ["-i", _regression_fixture("regression_sc_seqlogo.fa"), "-o", str(scratch_dir / "out"), "--iformat", "fa"],
        compare_stream="none",
        compare_files=("out.count_matrix.csv",),
        py_expected_exit=1,
        rust_expected_exit=0,
        rust_expect_files=("out.logo.pdf", "out.logo.mean_centered.pdf"),
    ),
    Case(
        name="fpkm_uq_basic",
        # The real `htseq-count` isn't in this sandbox's toolchain, so
        # both sides run against a committed mock script (--htseq-count
        # points at it directly, no PATH/subprocess-discovery magic
        # needed) that ignores its real arguments and prints a fixed
        # count table -- what's under test is calculate_fpkm's own
        # post-processing, not a real per-alignment scan. Found by this
        # case: several stderr progress lines were completely missing
        # (only the summary-table stats were ported, not the earlier
        # "Read gene information file ..."/"Total genes: N"/"Total
        # protein-coding genes: N" trio from read_gene_information, nor
        # the later un-timestamped "Read gene count file to calculate
        # FPKM and FPKM-UQ: ..." line). Stream comparison isn't used
        # here (upstream's `printlog` lines carry a LOCAL-time timestamp
        # this port deliberately renders in UTC instead, an accepted,
        # disclosed simplification -- see crates/cli/src/bin/fpkm_uq.rs
        # module docs -- so byte-exact stderr comparison isn't
        # meaningful); the two DATA files are what actually matter and
        # both are verified byte-identical.
        ensure_fixture=ensure_fpkm_uq_fixtures,
        py_script="FPKM-UQ.py",
        rust_bin="FPKM_UQ",
        py_args=lambda scratch_dir: [
            "--bam", _regression_fixture("regression_fpkm_uq_dummy.bam"),
            "--gtf", _regression_fixture("regression_fpkm_uq_dummy.gtf"),
            "--info", _regression_fixture("regression_fpkm_uq_genes.info.txt"),
            "-o", str(scratch_dir / "out"),
            "--htseq-count", _regression_fixture("mock_htseq_count.sh"),
        ],
        rust_args=lambda scratch_dir: [
            "--bam", _regression_fixture("regression_fpkm_uq_dummy.bam"),
            "--gtf", _regression_fixture("regression_fpkm_uq_dummy.gtf"),
            "--info", _regression_fixture("regression_fpkm_uq_genes.info.txt"),
            "-o", str(scratch_dir / "out"),
            "--htseq-count", _regression_fixture("mock_htseq_count.sh"),
        ],
        compare_stream="none",
        compare_files=("out.FPKM-UQ.txt", "out.htseq.counts.txt"),
    ),
    Case(
        name="fpkm_uq_with_log2",
        # Tests --log2 flag, which transforms FPKM values to log2(FPKM + 1).
        # Same fixtures and mock htseq-count as fpkm_uq_basic, exercising
        # the log-scale branch of calculate_fpkm's output formatting.
        ensure_fixture=ensure_fpkm_uq_fixtures,
        py_script="FPKM-UQ.py",
        rust_bin="FPKM_UQ",
        py_args=lambda scratch_dir: [
            "--bam", _regression_fixture("regression_fpkm_uq_dummy.bam"),
            "--gtf", _regression_fixture("regression_fpkm_uq_dummy.gtf"),
            "--info", _regression_fixture("regression_fpkm_uq_genes.info.txt"),
            "-o", str(scratch_dir / "out"),
            "--htseq-count", _regression_fixture("mock_htseq_count.sh"),
            "--log2",
        ],
        rust_args=lambda scratch_dir: [
            "--bam", _regression_fixture("regression_fpkm_uq_dummy.bam"),
            "--gtf", _regression_fixture("regression_fpkm_uq_dummy.gtf"),
            "--info", _regression_fixture("regression_fpkm_uq_genes.info.txt"),
            "-o", str(scratch_dir / "out"),
            "--htseq-count", _regression_fixture("mock_htseq_count.sh"),
            "--log2",
        ],
        compare_stream="none",
        compare_files=("out.FPKM-UQ.txt", "out.htseq.counts.txt"),
    ),
    Case(
        name="bam2wig_basic",
        # Found by this case: (1) "Skip multi-hits: {bool}" printed
        # Rust's lowercase "false" instead of Python's capitalized
        # "False" (str(bool) in an f-string); (2) "Total WIG sum: N"
        # rendered as bare "100" via Rust's default f64 Display instead
        # of Python's float str() "100.0" (fixed via python_str_float).
        # Not exercised by this specific default-args case (no --wigsum
        # here, see the module's own doc comment for a --wigsum probe),
        # but both are real, permanent regression guards now that the
        # code paths exist. compare_stream is "none": stderr's final
        # line is a best-effort `wigToBigWig` invocation whose exact
        # failure text is inherently environment-dependent (present vs.
        # absent on PATH) -- the .wig file is the actual data output and
        # is what's verified here.
        ensure_fixture=ensure_track_fixtures,
        py_script="bam2wig.py",
        rust_bin="bam2wig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_reads.bam"),
            "-s", _track_fixture("track_chrom.sizes"),
            "-o", str(scratch_dir / "out"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_reads.bam"),
            "-s", _track_fixture("track_chrom.sizes"),
            "-o", str(scratch_dir / "out"),
        ],
        compare_stream="none",
        compare_files=("out.wig",),
    ),
    Case(
        name="bam2wig_with_wigsum",
        # -t/--wigsum: exercises calWigSum's OWN independent per-
        # chromosome scan (run before the main bamTowig phase), which
        # duplicates the "Processing <chrom> ..." progress line -- a
        # real code path the default-args case above deliberately
        # doesn't exercise (see its own comment). Verified via live
        # diff: all content matches except the final wigToBigWig-not-
        # found line, which is inherently environment-dependent (shell
        # error text vs. this port's own message) -- same reason
        # bam2wig_basic uses compare_stream="none" rather than "stderr".
        ensure_fixture=ensure_track_fixtures,
        py_script="bam2wig.py",
        rust_bin="bam2wig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_reads.bam"),
            "-s", _track_fixture("track_chrom.sizes"),
            "-o", str(scratch_dir / "out"), "-t", "100",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_reads.bam"),
            "-s", _track_fixture("track_chrom.sizes"),
            "-o", str(scratch_dir / "out"), "-t", "100",
        ],
        compare_stream="none",
        compare_files=("out.wig",),
    ),
    Case(
        name="normalize_bigwig_genome",
        # Found by this case: several stderr progress lines were
        # missing entirely ("Get chromosome sizes from BigWig header
        # ...", "Normalizing BigWig file ...", "Writing <chrom> ...")
        # -- same "silent CLI" bug class found repeatedly this session.
        # Fixing "Writing <chrom> ..." required splitting the compute
        # function into two phases (calculate_wigsum/
        # render_normalized_body) so the CLI could interleave its
        # summary prints at the exact point upstream does, between
        # computing the normalization factor and writing chromosome
        # bodies -- previously the whole body was built silently before
        # any of those lines printed, a real (if here unobservable
        # since only one chromosome exists) print-ORDER bug, not just a
        # missing-line one. Default args (no --refgene): exercises
        # calculate_genome_wigsum's own per-chromosome scan.
        ensure_fixture=ensure_track_fixtures,
        py_script="normalize_bigwig.py",
        rust_bin="normalize_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-o", str(scratch_dir / "out.bgr"), "-t", "1000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-o", str(scratch_dir / "out.bgr"), "-t", "1000",
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.bgr",),
        normalize_paths=True,
    ),
    Case(
        name="normalize_bigwig_refgene",
        # Same fix as normalize_bigwig_genome, exercising the OTHER
        # branch instead: calculate_exonic_wigsum's 3 progress lines
        # ("Extract exons from...", "Merge overlapping exons ...",
        # "Calculate WIG sum covered by... only") plus -f wig output.
        ensure_fixture=ensure_track_fixtures,
        py_script="normalize_bigwig.py",
        rust_bin="normalize_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-o", str(scratch_dir / "out.wig"),
            "-f", "wig", "-r", _track_fixture("track_model.bed12"),
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-o", str(scratch_dir / "out.wig"),
            "-f", "wig", "-r", _track_fixture("track_model.bed12"),
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.wig",),
        normalize_paths=True,
    ),
    Case(
        name="overlay_bigwig_geometric_mean",
        # Found by this case: Python's f"{nan:.2f}" renders lowercase
        # "nan", but Rust's default `{:.2}` renders "NaN". A NaN
        # geometricMean result (sqrt of a negative product) IS reachable
        # and printed by upstream -- not filtered by the `value != 0.0`
        # check, since NaN never equals anything, including 0.0 --  and
        # the two track_signal*.bw fixtures happen to have an interval
        # where signal1 is negative and signal2 is positive (or vice
        # versa), exercising exactly this path. Add/Subtract/Product/
        # Max/Min/Average all matched upstream already (verified
        # manually before formalizing this case); Division is a
        # confirmed-broken upstream action (DIV-0015, both sides hard
        # error, not covered by a passing case). compare_stream is
        # "none" here (not "stderr" like the other overlay/normalize
        # cases): computing a NaN also makes numpy itself emit an
        # interpreter-internal `RuntimeWarning: invalid value
        # encountered in sqrt` to stderr, which is not something RSeQC
        # prints and isn't part of this port's contract to replicate
        # (its exact text is a numpy-version implementation detail).
        # The .wig file -- the actual data output, including the
        # lowercase "nan" this case exists to test -- is what matters
        # and is verified byte-identical.
        ensure_fixture=ensure_track_fixtures,
        py_script="overlay_bigwig.py",
        rust_bin="overlay_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "geometricMean", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "geometricMean", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        compare_stream="none",
        compare_files=("out.wig",),
    ),
    Case(
        name="overlay_bigwig_add",
        # Formalizes an action previously only spot-checked manually
        # (see overlay_bigwig_geometric_mean's own comment) into a real
        # harness case, so a future regression in the ordinary
        # non-NaN arithmetic path is actually caught, not just the NaN
        # edge case geometricMean happens to exercise.
        ensure_fixture=ensure_track_fixtures,
        py_script="overlay_bigwig.py",
        rust_bin="overlay_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Add", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Add", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.wig",),
    ),
    Case(
        name="overlay_bigwig_average",
        # Same formalization as overlay_bigwig_add, covering the other
        # arithmetic family (a division-by-2, not just addition).
        ensure_fixture=ensure_track_fixtures,
        py_script="overlay_bigwig.py",
        rust_bin="overlay_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Average", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Average", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.wig",),
    ),
    Case(
        name="overlay_bigwig_subtract",
        # Tests Subtract action: BigWig 1 - BigWig 2.
        ensure_fixture=ensure_track_fixtures,
        py_script="overlay_bigwig.py",
        rust_bin="overlay_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Subtract", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Subtract", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.wig",),
    ),
    Case(
        name="overlay_bigwig_product",
        # Tests Product action: BigWig 1 * BigWig 2.
        ensure_fixture=ensure_track_fixtures,
        py_script="overlay_bigwig.py",
        rust_bin="overlay_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Product", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Product", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.wig",),
    ),
    Case(
        name="overlay_bigwig_max",
        # Tests Max action: maximum of corresponding values.
        ensure_fixture=ensure_track_fixtures,
        py_script="overlay_bigwig.py",
        rust_bin="overlay_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Max", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Max", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.wig",),
    ),
    Case(
        name="overlay_bigwig_min",
        # Tests Min action: minimum of corresponding values.
        ensure_fixture=ensure_track_fixtures,
        py_script="overlay_bigwig.py",
        rust_bin="overlay_bigwig",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Min", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-j", _track_fixture("track_signal2.bw"),
            "-a", "Min", "-o", str(scratch_dir / "out.wig"), "-c", "100000",
        ],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.wig",),
    ),
    Case(
        name="rpkm_saturation_basic",
        # Found by this case: (1) RPKM_saturation.rs's percentile-
        # resampling population was rebuilt INDEPENDENTLY each
        # iteration instead of being CUMULATIVE across iterations --
        # upstream's `ranges`/`ranges_plus`/`ranges_minus` dicts are
        # declared ONCE before the percentile loop and never cleared
        # (same accumulation pattern as junction_saturation.py's
        # `uniqSpliceSites`), so a later percentile's RPKM reflects ALL
        # points sampled so far, not just that iteration's own slice.
        # Getting this wrong produced a completely different, incorrect
        # saturation curve -- a real scientific-correctness bug, not a
        # formatting one. (2) "Load BAM file ... " used `eprintln!`
        # (extra unwanted newline) instead of `eprint!`, and was missing
        # its second space before "Done". (3) The entire per-percentile
        # progress-message pipeline ("sampling N% (...) fragments ...",
        # "assign reads to transcripts in <refbed> ...", a trailing
        # blank line) was completely missing. (4) The CLI printed three
        # "Created ..." lines upstream's main() never prints at all.
        #
        # This fixture (make_rpkm_saturation_fixture.py) deliberately
        # has exactly ONE qualifying alignment (one exon block), so the
        # WHOLE saturation table is independent of random.shuffle's
        # order -- verified deterministic across multiple independent
        # Rust reruns before relying on it here. A larger population
        # would only guarantee the FINAL (100%) column is order-
        # invariant (see the module's own doc comment), not the whole
        # file, since the percentile ranges are cumulative.
        ensure_fixture=ensure_rpkm_saturation_fixture,
        py_script="RPKM_saturation.py",
        rust_bin="RPKM_saturation",
        py_args=lambda scratch_dir: [
            "-i", _regression_fixture("rpkm_saturation_basic.bam"),
            "-r", _regression_fixture("rpkm_saturation_model.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _regression_fixture("rpkm_saturation_basic.bam"),
            "-r", _regression_fixture("rpkm_saturation_model.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.eRPKM.xls", "out.rawCount.xls", "out.saturation.r"),
        normalize_paths=True,
    ),
    Case(
        name="rpkm_saturation_sam_text",
        # DIV-0002/0004: SAM-text input, closed for RPKM_saturation.py via
        # open_alignments. Must reuse the dedicated one-alignment fixture
        # (not bam_stat_basic.sam, which has 10 alignments and is NOT
        # shuffle-order-independent -- confirmed by this case genuinely
        # FAILING intermittently the first time it was tried with that
        # fixture, a fixture-choice mistake, not a code bug).
        ensure_fixture=ensure_rpkm_saturation_sam_fixture,
        py_script="RPKM_saturation.py",
        rust_bin="RPKM_saturation",
        py_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "rpkm_saturation_basic.sam"),
            "-r", _regression_fixture("rpkm_saturation_model.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", str(REPO_ROOT / "verification" / "fixtures" / "rpkm_saturation_basic.sam"),
            "-r", _regression_fixture("rpkm_saturation_model.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.eRPKM.xls", "out.rawCount.xls", "out.saturation.r"),
        normalize_paths=True,
    ),
    Case(
        name="genebody_coverage2_basic",
        # Found by this case: the per-gene progress line
        # ("\t<n> genes finished\r", end=' ') was missing its trailing
        # space from `end=' '` -- the literal's own `\r` was already
        # correctly ported, just not the separate end-of-print space
        # that follows it. Exercises both track_model.bed12 transcripts
        # (a single-exon 100bp-boundary one and a two-exon one whose
        # first exon alone clears the legacy 100bp filter, see DIV-0014
        # in make_track_fixtures.py's own docstring) against real
        # signal data.
        ensure_fixture=ensure_track_fixtures,
        py_script="geneBody_coverage2.py",
        rust_bin="geneBody_coverage2",
        py_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-r", _track_fixture("track_model.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        rust_args=lambda scratch_dir: [
            "-i", _track_fixture("track_signal.bw"), "-r", _track_fixture("track_model.bed12"),
            "-o", str(scratch_dir / "out"), "--skip-plot",
        ],
        compare_stream="stderr",
        stream_format="exact",
        compare_files=("out.geneBodyCoverage.txt", "out.geneBodyCoverage_plot.r"),
        normalize_paths=True,
    ),
    Case(
        name="bam2wig_synthetic_stdout",
        # Found by verification/synthetic_sweep.py: upstream prints
        # "Run wigToBigWig <wig> <sizes> <bw> " to STDOUT and then runs the
        # tool through `subprocess.call(..., shell=True)`, so a missing tool
        # yields the shell's own "not found" line rather than "Failed to
        # call" (qcmodule/SAM.py:2620-2623). The port printed nothing on
        # stdout and its own message on stderr.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="bam2wig.py",
        rust_bin="bam2wig",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-s", _synthetic("chrom.sizes"), "-o", str(d / "out")],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-s", _synthetic("chrom.sizes"), "-o", str(d / "out")],
        compare_stream="both",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.wig",),
    ),
    Case(
        name="fpkm_count_synthetic_stderr",
        # Found by verification/synthetic_sweep.py: the progress report
        # differed from upstream FPKM_count.py:468-613 -- "...from <bed>..."
        # (no space), "Counting total fragment ...  Done" on one line,
        # float-formatted totals ("78.0", not "78") and the per-transcript
        # "\r<n> transcripts finished" counter.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="FPKM_count.py",
        rust_bin="FPKM_count",
        py_args=lambda d: ["-i", _synthetic("pe_placed.bam"), "-r", _synthetic("model.bed12"), "-o", str(d / "out")],
        rust_args=lambda d: ["-i", _synthetic("pe_placed.bam"), "-r", _synthetic("model.bed12"), "-o", str(d / "out")],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        compare_files=("out.FPKM.xls",),
    ),
    Case(
        name="fpkm_uq_synthetic_missing_gene",
        # Found by verification/synthetic_sweep.py: for a count-file gene
        # absent from --info, upstream warns "Warning: <id> is absent from
        # <info path>; skipped" while writing the table, i.e. AFTER the
        # summary lines (FPKM-UQ.py:346-353); the port printed a generic
        # "info file" text before them.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="FPKM-UQ.py",
        rust_bin="FPKM_UQ",
        py_args=lambda d: ["--bam", _synthetic("pe.bam"), "--gtf", _synthetic("model.gtf"),
                           "--info", _synthetic("genes.info.txt"), "-o", str(d / "out"),
                           "--htseq-count", _synthetic("mock_htseq_count.sh")],
        rust_args=lambda d: ["--bam", _synthetic("pe.bam"), "--gtf", _synthetic("model.gtf"),
                             "--info", _synthetic("genes.info.txt"), "-o", str(d / "out"),
                             "--htseq-count", _synthetic("mock_htseq_count.sh")],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        strip_log_prefixes=True,
        compare_files=("out.FPKM-UQ.txt", "out.htseq.counts.txt"),
    ),
    Case(
        name="fpkm_uq_error_prefix",
        # Found by verification/synthetic_sweep.py: 26 upstream scripts
        # report runtime errors via `parser.exit(1, f"{parser.prog}: error:
        # {exc}\n")` (e.g. FPKM-UQ.py:448); the port printed a bare
        # "error: ...". Here the committed mock's gene IDs are all absent
        # from the synthetic info file -> "no protein-coding gene counts".
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="FPKM-UQ.py",
        rust_bin="FPKM_UQ",
        py_args=lambda d: ["--bam", _synthetic("pe.bam"), "--gtf", _synthetic("model.gtf"),
                           "--info", _synthetic("genes.info.txt"), "-o", str(d / "out"),
                           "--htseq-count", _regression_fixture("mock_htseq_count.sh")],
        rust_args=lambda d: ["--bam", _synthetic("pe.bam"), "--gtf", _synthetic("model.gtf"),
                             "--info", _synthetic("genes.info.txt"), "-o", str(d / "out"),
                             "--htseq-count", _regression_fixture("mock_htseq_count.sh")],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        strip_log_prefixes=True,
        expected_exit_code=1,
    ),
    Case(
        name="read_NVC_missing_output_dir",
        # Upstream's validate_args refuses an output prefix whose parent directory
        # does not exist, before reading any input, via parser.error() -> exit 2.
        # 19 of the port's binaries had no such check: they read the whole
        # alignment, computed every metric, and only then failed on the first output
        # open with `No such file or directory (os error 2)` and exit 1 -- an error
        # naming neither the directory nor the flag, indistinguishable from a
        # missing input, arrived at after all the work. That is the audit's "no
        # silent metric loss" failure: the metrics are computed then discarded
        # without ever being reported.
        #
        # Compared on the final error line and the exit status; the usage block
        # above the message is argparse's on one side and clap's on the other
        # (DIV-0025).
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="read_NVC.py",
        rust_bin="read_NVC",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-o", str(d / "absent_dir" / "out"),
                           "--skip-plot"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-o", str(d / "absent_dir" / "out"),
                             "--skip-plot"],
        compare_stream="stderr",
        stream_format="error_line",
        normalize_paths=True,
        py_expected_exit=2,
        rust_expected_exit=2,
    ),
    Case(
        name="read_GC_missing_output_dir",
        # Same upstream validate_args rule, second command, so the contract is not
        # recorded from a single sample.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="read_GC.py",
        rust_bin="read_GC",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-o", str(d / "absent_dir" / "out"),
                           "--skip-plot"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-o", str(d / "absent_dir" / "out"),
                             "--skip-plot"],
        compare_stream="stderr",
        stream_format="error_line",
        normalize_paths=True,
        py_expected_exit=2,
        rust_expected_exit=2,
    ),
    Case(
        name="geneBody_coverage_missing_output_dir",
        # Same rule, and this command additionally reads a gene model, so it is the
        # case where doing the work first was most obviously wasteful.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"),
                           "-o", str(d / "absent_dir" / "out"), "--skip-plot"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"),
                             "-o", str(d / "absent_dir" / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="error_line",
        normalize_paths=True,
        py_expected_exit=2,
        rust_expected_exit=2,
    ),
    Case(
        name="junction_annotation_missing_output_dir",
        # Same rule for a command that writes two files from one prefix, so the
        # refusal is proven to precede both writes.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="junction_annotation.py",
        rust_bin="junction_annotation",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"),
                           "-o", str(d / "absent_dir" / "out")],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"),
                             "-o", str(d / "absent_dir" / "out")],
        compare_stream="stderr",
        stream_format="error_line",
        normalize_paths=True,
        py_expected_exit=2,
        rust_expected_exit=2,
    ),
    Case(
        name="bam2wig_missing_output_dir_is_upstream_inconsistent",
        # The counterpart, and the reason the check is NOT applied everywhere:
        # bam2wig.py and bam2fq.py are the only two upstream scripts that take an
        # output prefix and do not validate its parent. Upstream reaches the work and
        # fails on the first output open with exit 1, so a port that "fixed" this
        # would disagree with upstream on an input upstream accepts (DIV-0026). Only
        # the exit status is compared: upstream's message is Python's OSError
        # repr, `[Errno 2] No such file or directory: '.../out.wig'`, and the port's
        # is Rust's io::Error Display, so the text differs by construction.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="bam2wig.py",
        rust_bin="bam2wig",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-s", _synthetic("chrom.sizes"),
                           "-o", str(d / "absent_dir" / "out")],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-s", _synthetic("chrom.sizes"),
                             "-o", str(d / "absent_dir" / "out")],
        compare_stream="none",
        py_expected_exit=1,
        rust_expected_exit=1,
    ),
    Case(
        name="genebody_coverage_depth_cap_divergence",
        # DIV-0024. The main synthetic fixture has depth around 40, so it cannot
        # reach pysam's default max_depth of 8000 and cannot test it: this suite
        # passed 90/90 while geneBody_coverage disagreed with upstream on 76 of 100
        # bins of a real 8.2M-record alignment. The benchmark harness's structural
        # gate is what caught it.
        #
        # The port reimplements max_depth as a per-position budget over the VISITED
        # set (crates/commands/src/tin.rs:274-296). Upstream's budget is over the
        # pileup buffer. The difference is invisible until the cap binds AND a read
        # is in a deletion, which is why this fixture is 8,100 copies of one
        # 20M20D20M read: the cap binds, the reads become is_del through the middle,
        # and the port reports 0 at positions the clean reads demonstrably cover.
        #
        # Asserted as a divergence WITH A FLOOR rather than tolerated: this case
        # passes only while the two files still differ by at least the recorded
        # amount. When the cap semantics are fixed the files become identical, this
        # case fails, and it says which ledger entry to retire. A divergence that
        # silently changes shape also fails here.
        #
        # Measured on this fixture: 50 of 100 bins differ, largest single-cell
        # difference 20 reads at bin 39, and the port is LOW at every one of them
        # (11..20 reads where clean reads demonstrably cover). The floor of 10 sits
        # between "identical" and the observed value, so a fix trips the case and
        # fixture jitter does not.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda d: ["-i", _synthetic("depth_cap.bam"), "-r", _synthetic("depth_cap.bed12"),
                           "-o", str(d / "out"), "--skip-plot"],
        rust_args=lambda d: ["-i", _synthetic("depth_cap.bam"), "-r", _synthetic("depth_cap.bed12"),
                             "-o", str(d / "out"), "--skip-plot"],
        compare_stream="none",
        divergent_files=(("out.geneBodyCoverage.txt", 10.0),),
    ),
    Case(
        name="genebody_coverage_synthetic_pileup",
        # Found by verification/synthetic_sweep.py: upstream's
        # samfile.pileup() (geneBody_coverage.py:187) runs with pysam's
        # default "samtools" stepper, which drops orphan reads (paired but
        # not properly paired, ignore_orphans=True), and with
        # ignore_overlaps=True, under which htslib rewrites overlapping
        # mates' base qualities (sam.c tweak_overlap_quality: summed for
        # matching bases, 0.8x for mismatches) before min_base_quality=13
        # applies. The port counted orphans and approximated the overlap
        # rule by "higher quality mate wins".
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"), "-o", str(d / "out"),
                           "--skip-plot"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"), "-o", str(d / "out"),
                             "--skip-plot"],
        compare_stream="none",
        compare_files=("out.geneBodyCoverage.txt",),
    ),
    Case(
        name="tin_synthetic_pileup",
        # Same pileup-default root cause as genebody_coverage_synthetic_pileup,
        # via tin.py's genebody_coverage() pileup (per-transcript TIN).
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="tin.py",
        rust_bin="tin",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"), "-o", str(d), "-c", "2"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"), "-o", str(d), "-c", "2"],
        compare_stream="none",
        compare_files=("pe.tin.xls",),
    ),
    Case(
        name="genebody_coverage_synthetic_skewness",
        # Found by verification/synthetic_sweep.py: the per-sample skewness
        # printed to stderr (geneBody_coverage.py:68-77) uses np.std(ddof=1)
        # / np.mean -- numpy's pairwise summation -- and np.float64 ** 3
        # (C pow); the port summed left-to-right and used powi(3), differing
        # in the last digits. Multi-sample so several skewness lines print.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="geneBody_coverage.py",
        rust_bin="geneBody_coverage",
        py_args=lambda d: ["-i", _synthetic("pe.bam") + "," + _synthetic("se.bam") + "," + _synthetic("sc.bam"),
                           "-r", _synthetic("model.bed12"), "-o", str(d / "out"), "--skip-plot"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam") + "," + _synthetic("se.bam") + "," + _synthetic("sc.bam"),
                             "-r", _synthetic("model.bed12"), "-o", str(d / "out"), "--skip-plot"],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        strip_log_prefixes=True,
        compare_files=("out.geneBodyCoverage.txt", "out.geneBodyCoverage.r"),
    ),
    Case(
        name="tin_synthetic_summary_numpy",
        # Same root cause for tin.py's summary (tin.py:541-543: np.mean,
        # np.std over the per-transcript TINs).
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="tin.py",
        rust_bin="tin",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"), "-o", str(d), "-c", "1"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", _synthetic("model.bed12"), "-o", str(d), "-c", "1"],
        compare_stream="none",
        compare_files=("pe.tin.xls", "pe.summary.txt"),
    ),
    Case(
        name="tin_synthetic_multi_input_logging",
        # Found by verification/synthetic_sweep.py: tin.py resolves -i via
        # getBamFiles.get_bam_files (comma list / directory / list file,
        # tin.py:580-599) and logs its INFO progress unconditionally; the
        # port accepted a single path only (a comma list failed with "No such
        # file or directory") and printed progress only with --verbose.
        ensure_fixture=ensure_synthetic_fixtures,
        py_script="tin.py",
        rust_bin="tin",
        py_args=lambda d: ["-i", _synthetic("pe.bam") + "," + _synthetic("se.bam"), "-r", _synthetic("model.bed12"),
                           "-o", str(d), "-c", "1"],
        rust_args=lambda d: ["-i", _synthetic("pe.bam") + "," + _synthetic("se.bam"), "-r", _synthetic("model.bed12"),
                             "-o", str(d), "-c", "1"],
        compare_stream="stderr",
        stream_format="exact",
        normalize_paths=True,
        strip_log_prefixes=True,
        compare_files=("pe.tin.xls", "pe.summary.txt", "se.tin.xls", "se.summary.txt"),
    ),

    # ------------------------------------------------------------------
    # Negative / degenerate-input branch cases.
    #
    # PUBLICATION_PLAN.md 3 item 5 requires error, option, and overwrite
    # branches per command; until 2026-09-30 the 84-case matrix was almost
    # entirely happy-path, which is precisely the gap that let the
    # sliding-window `tin` regression ship. These are the first of that
    # class. Each was checked against real upstream before being written
    # down: where the two sides legitimately differ, the difference is
    # recorded as a per-side exit override (and, if behavioural, in
    # `compatibility/divergences.yaml`) rather than asserted away.
    # ------------------------------------------------------------------
    Case(
        name="bam_stat_empty_bam",
        # Zero records, valid header+index. Upstream prints its full
        # all-zero report and exits 0; so does the port. Pins that the
        # streaming reader handles a legitimately empty file rather than
        # treating EOF as an error.
        ensure_fixture=ensure_malformed_bam_fixtures,
        py_script="bam_stat.py",
        rust_bin="bam_stat",
        py_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "empty.bam")],
        rust_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "empty.bam")],
        compare_stream="stdout",
        stream_format="exact",
        required_labels=("Total records", "Unmapped reads", "Read-1"),
    ),
    Case(
        name="bam_stat_not_a_bam",
        # Truncated at the header: a valid gzip member containing no BAM
        # header. Both sides must fail, and must not be allowed to "agree"
        # on a success: `expected_exit_code=1` plus
        # the harness's identical-exit assertion means a case like this
        # cannot pass by both sides erroring in different, mutually
        # satisfying ways -- only by both failing as upstream does.
        # The messages differ (upstream: "file does not contain alignment
        # data"; port: the underlying IO error), so only the exit status
        # and the stdout stream are compared.
        ensure_fixture=ensure_malformed_bam_fixtures,
        py_script="bam_stat.py",
        rust_bin="bam_stat",
        py_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "junk.bam")],
        rust_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "junk.bam")],
        compare_stream="none",
        expected_exit_code=1,
    ),
    Case(
        name="read_GC_empty_bam",
        # Same empty-but-valid input. Upstream exits 0 here; so does the
        # port, and both write a degenerate report.
        ensure_fixture=ensure_malformed_bam_fixtures,
        py_script="read_GC.py",
        rust_bin="read_GC",
        py_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "empty.bam"),
                           "--skip-plot", "-o", str(d / "gc")],
        rust_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "empty.bam"),
                             "--skip-plot", "-o", str(d / "gc")],
        compare_stream="none",
        compare_files=("gc.GC.xls",),
    ),
    Case(
        name="read_NVC_empty_bam_diverges_from_upstream",
        # DIV-0023. Upstream CRASHES here with an unhandled
        # `UnboundLocalError` (`RNA_read` is only bound inside the
        # per-read loop, so zero records leaves it unbound before the
        # post-loop code dereferences it). The port exits 0 with a
        # well-formed empty report.
        #
        # This is recorded as a DIVERGENCE, not asserted as agreement: the
        # per-side exit overrides make the harness check each side against
        # its own documented value and skip the identical-exit assertion,
        # so a future change that makes either side behave differently in
        # an unrecorded way fails here.
        ensure_fixture=ensure_malformed_bam_fixtures,
        py_script="read_NVC.py",
        rust_bin="read_NVC",
        py_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "empty.bam"),
                           "--skip-plot", "-o", str(d / "nvc")],
        rust_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "empty.bam"),
                             "--skip-plot", "-o", str(d / "nvc")],
        compare_stream="none",
        py_expected_exit=1,
        rust_expected_exit=0,
    ),
    Case(
        name="tin_missing_bai_sidecar",
        # tin.py resolves `-i` through getBamFiles, which requires a `.bai`
        # next to each BAM. Removing it must produce the same "no BAM files
        # found" outcome and the same non-zero exit on both sides -- the
        # silent-success risk here is the port accepting an unindexed file
        # and then quietly scoring everything 0.0, which is the exact shape
        # of the regression fixed in 6cd93e6.
        ensure_fixture=ensure_malformed_bam_fixtures,
        py_script="tin.py",
        rust_bin="tin",
        py_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "noidx.bam"),
                           "-r", _synthetic("model.bed12"), "-o", str(d)],
        rust_args=lambda d: ["-i", str(REPO_ROOT / "verification" / "fixtures" / "noidx.bam"),
                             "-r", _synthetic("model.bed12"), "-o", str(d)],
        compare_stream="stderr",
        stream_format="exact",
        strip_log_prefixes=True,
        expected_exit_code=1,
    ),
    Case(
        name="tin_malformed_bed12_lines",
        # Every line is too short / non-integer. Upstream wraps each line in
        # a broad `except` and skips it with a stderr note, so it succeeds
        # with an empty transcript set; the port does the same, and both
        # must produce the same (header-only) report rather than one side
        # aborting on the first bad line.
        ensure_fixture=ensure_malformed_bam_fixtures,
        py_script="tin.py",
        rust_bin="tin",
        py_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", str(REPO_ROOT / "verification" / "fixtures" / "malformed.bed12"),
                           "-o", str(d)],
        rust_args=lambda d: ["-i", _synthetic("pe.bam"), "-r", str(REPO_ROOT / "verification" / "fixtures" / "malformed.bed12"),
                             "-o", str(d)],
        compare_stream="none",
        compare_files=("pe.tin.xls", "pe.summary.txt"),
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

    require_real_panel()
    if REAL_DATA_DIR:
        print(f"real-data panel: {REAL_DATA_DIR}")
        for k, v in sorted(REAL_EQUIVALENTS.items()):
            print(f"  {_synthetic(k) if _real_path(k) else '(no real equivalent)':>0s} {k} <- {v}")
        print()

    all_ok = True
    skipped = 0
    for case in cases:
        gaps = _real_data_gaps(case)
        if gaps:
            skipped += 1
            print(f"=== {case.name} ===")
            print(f"  SKIPPED: no real-data equivalent for {', '.join(gaps)} "
                  f"(the real panel has no such sample; running this on synthetic "
                  f"input would make a 'real-data' run validate nothing real)")
            continue
        all_ok &= run_case(case)
        print()

    if all_ok:
        print(f"All {len(cases)} case(s) PASSED")
        return 0
    print("One or more cases FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
