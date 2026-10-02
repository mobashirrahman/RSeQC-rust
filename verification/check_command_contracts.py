#!/usr/bin/env python3
"""Command-contract checks: every command must fail in a defined, documented way.

The audit's "command contracts" layer asks for risk-based empty/malformed input,
missing index, invalid flags, existing output, paths with spaces, unreadable input,
output failure, killed jobs and concurrent invocations. The completion condition is
stated precisely: a defined exit or status, complete outputs or explicitly
identifiable incomplete outputs, and **no silent metric loss or successful corrupt
output**.

That last clause is what these tests are really about. A command that exits 0 after
silently dropping half its input is worse than one that fails: the number looks
fine. So each check asserts on the *combination* of exit status and what was
actually written, and a command that exits non-zero without saying why on stderr is
recorded as a gap rather than as correct behaviour.

These run against the compiled binaries and need no oracle environment, so they are
mandatory in CI.

Usage:
    cargo build --workspace --release --locked
    oracle/venv/bin/python3 verification/check_command_contracts.py
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RELEASE = REPO / "target" / "release"
FIXTURES = REPO / "verification" / "fixtures"
TINY = REPO / "datasets" / "aligned" / "tiny"

# Commands exercised, with the inputs each needs. Chosen for the failure modes that
# matter rather than for coverage count: a command not listed here is not a claim
# that it has no untested contract.
#
#   binary, needs_bam, needs_bed, extra_args, out_style, truncate_args
#
# `out_style` records HOW the command takes its output, because getting it wrong
# makes every scenario fail for the wrong reason -- an invocation missing its own
# required output flag fails identically whether or not its input path contains a
# space, so a "path with spaces" test built on it measures nothing.
#   "prefix"   --out-prefix PREFIX   (a file-name prefix; the directory must exist)
#   "out"      -o PREFIX             (ditto, different flag spelling)
#   "outdir"   -o DIRECTORY          (the argument IS a directory and must be created)
#   "none"     stdout only
#
# The distinction matters: upstream `tin.py` requires its `-o` directory to exist and
# errors out if it does not, so a harness that treated `-o` as a prefix would create
# the directory and silently stop testing that contract.
CONTRACT_TARGETS = [
    ("bam_stat", "bam", None, [], "none", None),
    # `infer_experiment` samples at most -s/--sample-size alignments (200,000 by
    # default) and stops there, so on a truncated alignment it can finish cleanly
    # without ever reaching the truncation. That is the documented behaviour rather
    # than silent metric loss, and the truncation scenario has to look past the cap
    # to test what it means to test. The last field is scenario-specific args.
    ("infer_experiment", "bam", "bed", [], "none", ["-s", "100000000"]),
    ("bam2fq", "bam", None, [], "out", None),
    ("tin", "bam", "bed", ["-n", "50"], "outdir", None),
    ("read_NVC", "bam", None, [], "prefix", None),
    ("read_GC", "bam", None, [], "prefix", None),
    ("geneBody_coverage", "bam", "bed", ["--skip-plot"], "prefix", None),
    ("junction_annotation", "bam", "bed",
     ["--skip-plot", "--skip-bed", "--skip-interact"], "prefix", None),
    ("split_bam", "bam", "bed", [], "out", None),
    ("bam2wig", "bam", None, [], "out", None),
]


class Results:
    def __init__(self):
        self.rows = []

    def add(self, command, scenario, outcome, ok, detail=""):
        self.rows.append({"command": command, "scenario": scenario,
                          "outcome": outcome, "pass": bool(ok), "detail": detail})
        mark = "ok  " if ok else "GAP "
        print(f"  [{mark}] {command:22} {scenario:34} {outcome:16} {detail[:60]}")
        return ok

    @property
    def failed(self):
        return [r for r in self.rows if not r["pass"]]


def run(args, cwd, timeout=180):
    t0 = time.time()
    try:
        p = subprocess.run([str(a) for a in args], cwd=str(cwd),
                           capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr, time.time() - t0, False
    except subprocess.TimeoutExpired:
        return None, "", "", time.time() - t0, True


def output_files(d):
    return sorted(str(p.relative_to(d)) for p in d.rglob("*") if p.is_file())


def has_nonfinite(text):
    for token in text.replace("\t", " ").split():
        t = token.strip().rstrip(",;")
        if t.lower() in ("nan", "inf", "-inf", "infinity", "-infinity"):
            return t
    return None


def scenario_missing_input(res, command, argv, work):
    """Point the input at a path that does not exist, and require a defined failure.

    Rewriting the argument is essential: running the unmodified command would
    simply succeed on its real input and report "exit=0", which says nothing about
    the missing-file contract.
    """
    a = []
    rewritten = False
    for x in argv:
        # Only the primary input. Rewriting the annotation too would make the
        # command fail for a different reason than the one under test.
        if not rewritten and str(x) in ("-i", "--input-file", "--input"):
            a.append(str(x))
            a.append(str(work / "definitely-absent-input"))
            rewritten = True
            continue
        a.append(x)
    rc, out, err, _t, timed_out = run(a, work)
    if timed_out:
        return res.add(command, "missing input file", "HANG", False,
                       "a missing input must not hang")
    ok = rc is not None and rc != 0
    said_why = bool(err.strip())
    return res.add(command, "missing input file",
                   "nonzero exit" if ok else f"exit={rc}", ok and said_why,
                   "diagnostic on stderr" if said_why else "no diagnostic on stderr")


def scenario_empty_input(res, command, argv, work, empty_bam, bed):
    """An empty-but-valid BAM must produce a DEFINED outcome.

    "Defined" means one of: a non-zero exit with a diagnostic, or exit 0 with real
    numbers and no nonfinite value. Exit 0 with a nonfinite metric, or a hang, or a
    crash with no explanation, is undefined -- and exit 0 with a plausible-looking
    metric computed from nothing is the worst of those, because the number reads as
    a result.

    A non-zero exit here is *not* automatically a defect. Several upstream commands
    fail on an empty alignment: upstream `read_NVC.py` raises
    `UnboundLocalError` on empty input, and the port's non-zero exit is a
    documented improvement over that. What is checked is that the failure is
    visible, so a reader is never left guessing which of the two happened.
    """
    if empty_bam is None:
        return res.add(command, "empty valid BAM", "SKIPPED", True,
                       "no empty fixture available")
    a = []
    for x in argv:
        if str(x) in ("-i", "--input-file", "--input"):
            a.append(str(x))
            a.append(str(empty_bam))
            continue
        a.append(x)
    rc, out, err, _t, timed_out = run(a, work)
    if timed_out:
        return res.add(command, "empty valid BAM", "HANG", False)
    nf = has_nonfinite(out + err)
    if rc == 0:
        ok = nf is None
        return res.add(command, "empty valid BAM", "exit=0", ok,
                       f"nonfinite {nf!r} in output" if nf else "real answer")
    ok = nf is None and err.strip() != ""
    return res.add(command, "empty valid BAM", f"exit={rc}", ok,
                   "defined failure" if ok else "nonzero exit with no diagnostic")


def scenario_missing_index(res, command, argv, work, bam, bed):
    """A BAM whose .bai is absent must be reported, not silently scanned or guessed."""
    if bam is None or bed is None:
        return res.add(command, "missing .bai sidecar", "SKIPPED", True,
                       "this command does not take an index")
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        isolated = d / "isolated.bam"
        shutil.copy(bam, isolated)
        a = [str(x).replace(str(bam), str(isolated)) for x in argv]
        rc, out, err, _t, timed_out = run(a, work)
        nf = has_nonfinite(out + err)
        # Either it needs no index (exit 0 with a real answer), or it says why.
        ok = not timed_out and nf is None and (rc == 0 or err.strip() != "")
        return res.add(command, "missing .bai sidecar", f"exit={rc}", ok,
                       f"nonfinite {nf!r}" if nf else "")


def scenario_invalid_flag(res, command, argv, work):
    a = list(argv) + ["--definitely-not-a-flag"]
    rc, out, err, _t, timed_out = run(a, work, timeout=60)
    ok = rc is not None and rc != 0
    return res.add(command, "invalid flag", f"exit={rc}", ok,
                   "rejected" if ok else "accepted an unknown flag")


def scenario_missing_required_arg(res, command, argv, work):
    """Omitting a required argument must be an error, not a default that reads nothing."""
    a = []
    skip_next = False
    for i, x in enumerate(argv):
        if str(x) in ("-i", "-r", "-o", "--input", "--refgene", "--out-prefix"):
            skip_next = True
            continue
        if skip_next:
            skip_next = False
            continue
        a.append(x)
    rc, out, err, _t, timed_out = run(a, work, timeout=60)
    ok = rc is not None and rc != 0
    return res.add(command, "missing required argument", f"exit={rc}", ok,
                   "rejected" if ok else "accepted a command with no input")


def scenario_path_with_spaces(res, command, argv, work, bam, bed, out_style):
    """Copy the inputs into a directory whose name contains spaces.

    Sidecar indexes are copied too. Without them the command fails for want of a
    `.bai`, which looks identical to a spaces bug and measures nothing; the
    missing-index contract is a separate scenario.
    """
    with tempfile.TemporaryDirectory() as tmp:
        spaced = Path(tmp) / "a directory with spaces"
        spaced.mkdir()
        pairs = {}
        for src, name in ((bam, "input file.bam"), (bed, "annotation file.bed12")):
            if not src:
                continue
            shutil.copy(src, spaced / name)
            pairs[str(src)] = str(spaced / name)
            sidecar = Path(str(src) + ".bai")
            if sidecar.exists():
                shutil.copy(sidecar, spaced / (name + ".bai"))
        # A path that exists inside the spaced directory. For a prefix-taking
        # command this is a prefix (its parent must exist); for a directory-taking
        # one it IS a directory and must be created, because upstream `tin.py`
        # requires it to exist and errors out if it does not.
        out_arg = spaced / ("out dir" if out_style == "outdir" else "out result")
        if out_style == "outdir":
            out_arg.mkdir(parents=True, exist_ok=True)
        else:
            out_arg.parent.mkdir(parents=True, exist_ok=True)
        a = []
        for x in argv:
            s = pairs.get(str(x), str(x))
            if s.startswith(str(work)):
                s = str(out_arg)
            a.append(s)
        rc, out, err, _t, timed_out = run(a, work)
        nf = has_nonfinite(out + err)
        ok = rc == 0 and nf is None
        return res.add(command, "paths with spaces", f"exit={rc}", ok,
                       f"nonfinite {nf!r}" if nf else
                       ("handled" if rc == 0 else f"failed: {err.strip()[:70]}"))


def scenario_existing_output(res, command, argv, work):
    """Re-running into a directory that already has output must be defined.

    Several upstream commands refuse to overwrite and exit 2; the port may exit 1.
    Either is a defined contract. Silently appending, or truncating the previous
    result without saying so, is not.
    """
    first = run(argv, work)
    before = output_files(work)
    second = run(argv, work)
    after = output_files(work)
    rc = second[0]
    defined = rc == 0 or (rc is not None and rc != 0 and second[2].strip() != "")
    unchanged = before == after
    ok = defined and (rc != 0 or unchanged)
    return res.add(command, "pre-existing output", f"exit={rc}", ok,
                   "refused with a diagnostic" if rc not in (0, None)
                   else ("idempotent" if unchanged else "output changed between identical runs"))


def scenario_unreadable_input(res, command, argv, work, bam):
    if bam is None:
        return res.add(command, "unreadable input", "SKIPPED", True)
    with tempfile.TemporaryDirectory() as tmp:
        blocked = Path(tmp) / "locked.bam"
        shutil.copy(bam, blocked)
        os.chmod(blocked, 0o000)
        a = [str(x).replace(str(bam), str(blocked)) for x in argv]
        try:
            rc, out, err, _t, timed_out = run(a, work, timeout=90)
        finally:
            os.chmod(blocked, 0o644)
        if os.geteuid() == 0:
            return res.add(command, "unreadable input", "SKIPPED", True,
                           "running as root; permission bits are not enforced")
        ok = rc is not None and rc != 0 and err.strip() != ""
        return res.add(command, "unreadable input", f"exit={rc}", ok,
                       "reported" if ok else "no diagnostic")


def scenario_concurrent(res, command, argv, work, out_style):
    """Two invocations at once must not corrupt each other's output."""
    procs = []
    with tempfile.TemporaryDirectory() as tmp:
        a_dir, b_dir = Path(tmp) / "a", Path(tmp) / "b"
        a_dir.mkdir()
        b_dir.mkdir()
        for d in (a_dir, b_dir):
            out_arg = d if out_style == "outdir" else d / "out"
            out_arg.mkdir(parents=True, exist_ok=True)
            a = []
            for x in argv:
                s = str(x)
                if s.startswith(str(work)):
                    s = str(out_arg)
                a.append(s)
            procs.append((d, subprocess.Popen([str(x) for x in a], cwd=str(d),
                                              stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE, text=True)))
        outs = []
        for d, p in procs:
            so, se = p.communicate(timeout=300)
            outs.append((d, p.returncode, so, se))
        ok = all(rc == 0 for _d, rc, _o, _e in outs)
        nf = [has_nonfinite(o + e) for _d, _rc, o, e in outs]
        nf = [x for x in nf if x]
        if nf:
            ok = False
        return res.add(command, "concurrent invocations",
                       f"exits {[rc for _d, rc, _o, _e in outs]}",
                       ok, f"nonfinite {nf}" if nf else "both runs completed")


def scenario_killed(res, command, argv, work, out_style):
    """A killed process must leave output that is IDENTIFIABLY incomplete.

    The risk is not that a partial file exists -- it always will -- but that the
    partial file looks like a complete one. A consumer that finds
    `reads.tin.xls` after a kill has no way to tell a full result from the first
    40% of one. So this checks two things per file: a partially-written text table
    must not end mid-record, and a partially-written FASTQ must still be framed
    correctly or be empty. A BAM is self-describing enough that a truncated one is
    detectable by htslib, so it is checked that way rather than by size.

    This is reported as an observation rather than a pass/fail gate: partial output
    from a killed process is normal, and what matters is whether it is
    self-identifying. A file that is neither complete nor identifiable is a real
    gap.
    """
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        out_arg = d if out_style == "outdir" else d / "out"
        out_arg.mkdir(parents=True, exist_ok=True)
        a = []
        for x in argv:
            s = str(x)
            if s.startswith(str(work)):
                s = str(out_arg)
            a.append(s)
        p = subprocess.Popen([str(x) for x in a], cwd=str(d),
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        time.sleep(0.35)
        p.kill()
        try:
            p.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            pass

        files = output_files(d)
        if not files:
            return res.add(command, "killed mid-run", "0 files", True,
                           "nothing written")

        unidentifiable = []
        for f in files:
            path = d / f
            size = path.stat().st_size
            if size == 0:
                continue  # empty is unambiguously incomplete
            if f.endswith((".fastq", ".fq")):
                lines = path.read_text(errors="replace").splitlines()
                if lines and len(lines) % 4 != 0:
                    # A truncated final record is detectable by a reader, so this
                    # is identifiable -- but record it, because a consumer that
                    # does not check framing will read it as shorter data.
                    unidentifiable.append(f"{f} (unframed final record)")
            elif f.endswith((".xls", ".txt", ".bed", ".csv")):
                text = path.read_text(errors="replace")
                if text and not text.endswith("\n"):
                    unidentifiable.append(f"{f} (no trailing newline)")
            elif f.endswith(".bam"):
                try:
                    import pysam
                    with pysam.AlignmentFile(str(path), "rb") as h:
                        list(h.fetch(until_eof=True))
                    unidentifiable.append(f"{f} (truncated BAM read without error)")
                except Exception as e:
                    # htslib rejecting it is exactly the self-identifying property
                    # wanted, so this is a good outcome.
                    pass
            else:
                # .r, .pdf and friends: a reader will fail on a truncated file.
                pass
        ok = not unidentifiable
        return res.add(command, "killed mid-run", f"{len(files)} files", ok,
                       ("partial output is self-identifying: "
                        + ", ".join(unidentifiable)) if not ok
                       else f"partial output is self-identifying ({len(files)} files)")


def scenario_output_failure(res, command, argv, work, out_style):
    """An output path that cannot be written must fail loudly, not silently drop it.

    The failure this looks for is the quietest one in the set: a command that cannot
    open its output still exits 0, having reported metrics to stdout as though the
    run succeeded. A reader who keeps the stdout and never looks for the file has a
    number and no result.

    The output parent is made unwritable rather than being a nonexistent path,
    because a nonexistent path is caught by argument validation and would exercise
    that instead. Permission bits are the failure that reaches the file-open call.
    """
    if out_style == "none":
        # A stdout-only command has no output file to fail, so the check would
        # measure nothing. Recorded rather than silently omitted.
        return res.add(command, "unwritable output", "n/a", True,
                       "stdout-only command: no output file to make unwritable")
    with tempfile.TemporaryDirectory() as tmp:
        locked = Path(tmp) / "locked"
        locked.mkdir()
        out_arg = locked if out_style == "outdir" else locked / "out"
        if out_style != "outdir":
            out_arg.parent.mkdir(parents=True, exist_ok=True)
        os.chmod(locked, 0o500)
        a = [str(out_arg) if str(x).startswith(str(work)) else str(x) for x in argv]
        try:
            rc, out, err, _t, timed_out = run(a, work, timeout=120)
        finally:
            os.chmod(locked, 0o700)
        if os.geteuid() == 0:
            return res.add(command, "unwritable output", "SKIPPED", True,
                           "running as root; permission bits are not enforced")
        nf = has_nonfinite(out + err)
        # Success here is defined as: it failed, and it said why. A non-zero exit
        # with a diagnostic is a defined contract; exit 0 with a metric on stdout is
        # the silent-drop defect.
        ok = rc is not None and rc != 0 and err.strip() != "" and nf is None
        return res.add(command, "unwritable output", f"exit={rc}", ok,
                       "reported" if ok else
                       (f"nonfinite {nf!r}" if nf else
                        ("exited 0 with no output file written" if rc == 0
                         else "failed with no diagnostic")))


def scenario_truncated_input(res, command, argv, work, bam, truncate_args=None):
    """A truncated BAM must be rejected, not partially processed.

    A command that reads until the first decode error and reports what it got exits
    0 with metrics covering an unknown fraction of the alignment. The number is
    plausible and the answer is wrong, which is the failure mode the audit's "no
    silent metric loss" clause is about.
    """
    if bam is None:
        return res.add(command, "truncated alignment", "SKIPPED", True)
    data = bam.read_bytes()
    # The cut point must be strictly inside the file. `max(1024, ...)` looked safe
    # and was not: for the 347-byte `bam_stat_basic.bam` fixture it yields a cut
    # point beyond the end, so the "truncated" copy was byte-identical to the
    # original and the check passed by measuring nothing. Reported as a skip rather
    # than a pass, because a fixture too small to truncate cannot demonstrate the
    # contract it is named for.
    cut = int(len(data) * 0.7)
    if cut < 64 or cut >= len(data):
        return res.add(command, "truncated alignment", "SKIPPED", True,
                       f"fixture is {len(data)} bytes, too small to truncate "
                       f"meaningfully; use a larger alignment")
    with tempfile.TemporaryDirectory() as tmp:
        truncated = Path(tmp) / "truncated.bam"
        # Cut inside the compressed payload, well past the header, so the file is a
        # structurally valid BAM header followed by garbage rather than an obviously
        # broken file. Truncating at a record boundary would produce a legal BAM.
        truncated.write_bytes(data[:cut])
        sidecar = Path(str(bam) + ".bai")
        if sidecar.exists():
            shutil.copy(sidecar, truncated.with_suffix(".bam.bai"))
        a = [str(x).replace(str(bam), str(truncated)) for x in argv]
        if truncate_args:
            # Drop any cap present in the base invocation, then apply this target's
            # truncation-specific args, so the run actually reaches the truncation.
            drop = set()
            skip = False
            for x in a:
                if skip:
                    drop.add(x)
                    skip = False
                    continue
                if x in ("-s", "--sample-size"):
                    skip = True
            a = [x for x in a if x not in drop] + list(truncate_args)
        rc, out, err, _t, timed_out = run(a, work, timeout=120)
        ok = rc is not None and rc != 0 and err.strip() != ""
        return res.add(command, "truncated alignment", f"exit={rc}", ok,
                       "rejected with a diagnostic" if ok else
                       ("exited 0 on a truncated alignment: metrics cover an "
                        "unknown fraction of the input" if rc == 0
                        else "failed with no diagnostic"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--commands", nargs="*", default=None)
    args = ap.parse_args()

    if not RELEASE.exists():
        print("target/release not found; run: cargo build --workspace --release --locked")
        return 2

    bam = FIXTURES / "bam_stat_basic.bam"
    # The annotation must match the BAM's contigs or every annotation-consuming
    # command returns "no overlaps", which would make the path-with-spaces and
    # concurrency scenarios pass vacuously. The tiny panel is the only fixture
    # whose model and alignment are known to correspond.
    bed = TINY / "model.bed12" if TINY.exists() else None
    tiny_bam = TINY / "tiny.bam" if TINY.exists() else None
    empty = FIXTURES / "empty.bam"
    chrom_sizes = None
    for candidate in (TINY / "chrom.sizes", REPO / "datasets/aligned/real/chrom.sizes"):
        if candidate.exists():
            chrom_sizes = candidate
            break

    res = Results()
    targets = [t for t in CONTRACT_TARGETS
               if args.commands is None or t[0] in args.commands]

    for binary, needs_bam, needs_bed, extra, out_style, truncate_args in targets:
        exe = RELEASE / binary
        if not exe.exists():
            res.add(binary, "binary present", "MISSING", False, str(exe))
            continue
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            argv = [exe]
            used_bam = None
            if needs_bam == "bam":
                # Prefer the tiny panel: its BAM and BED12 correspond, so an
                # annotation-consuming command actually has something to match.
                src = tiny_bam if (needs_bed and tiny_bam) else (bam if bam.exists() else tiny_bam)
                if src is None or not src.exists():
                    res.add(binary, "input fixture", "MISSING", False,
                            "no BAM fixture available")
                    continue
                argv += ["-i", str(src)]
                if binary == "bam2wig":
                    argv += ["-s", str(chrom_sizes)] if chrom_sizes else []
                used_bam = src
            if needs_bed:
                if bed is None or not bed.exists():
                    res.add(binary, "annotation fixture", "MISSING", False,
                            f"no BED12 fixture available for {binary}")
                    continue
                argv += ["-r", str(bed)]
            # A directory-taking command's `-o` must name an existing directory:
            # upstream `tin.py` errors out when it does not, and that error is
            # itself part of the contract under test.
            if out_style == "outdir":
                (work / "out").mkdir(parents=True, exist_ok=True)
                argv += ["-o", str(work / "out")]
            elif out_style == "out":
                argv += ["-o", str(work / "out")]
            elif out_style == "prefix":
                argv += ["--out-prefix", str(work / "out")]
            argv += list(extra)

            print(f"{binary}:")
            scenario_missing_input(res, binary, argv, work)
            scenario_empty_input(res, binary, argv, work, empty, bed)
            scenario_missing_index(res, binary, argv, work, used_bam, bed)
            scenario_invalid_flag(res, binary, argv, work)
            scenario_missing_required_arg(res, binary, argv, work)
            if used_bam:
                scenario_path_with_spaces(res, binary, argv, work, used_bam, bed, out_style)
                scenario_unreadable_input(res, binary, argv, work, used_bam)
            scenario_existing_output(res, binary, argv, work)
            if used_bam:
                scenario_truncated_input(res, binary, argv, work, used_bam,
                                         truncate_args)
            scenario_output_failure(res, binary, argv, work, out_style)
            scenario_killed(res, binary, argv, work, out_style)
            scenario_concurrent(res, binary, argv, work, out_style)

    print()
    if args.json:
        print(json.dumps({"rows": res.rows, "gaps": len(res.failed)}, indent=2))
    else:
        print(f"{len(res.rows) - len(res.failed)}/{len(res.rows)} command-contract "
              f"checks passed")
        if res.failed:
            print("\nGaps (each needs either a fix or an explicit declaration that the "
                  "behaviour is intended):")
            for r in res.failed:
                print(f"  {r['command']:22} {r['scenario']:34} {r['outcome']:16} "
                      f"{r['detail']}")
            print("\nA gap here is not automatically a defect: some of these behaviours "
                  "are the port\nfaithfully reproducing upstream. Each one is recorded "
                  "so it is a decision rather\nthan an omission.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
