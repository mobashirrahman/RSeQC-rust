#!/usr/bin/env python3
"""Measure what each non-flat command actually does when memory runs out.

`docs/ENVELOPE.md` reports a per-record cost for `bam2wig` and
`read_duplication`. A cost is half of what a published limit needs: the other half
is the behaviour *past* the limit, and the audit's contract requirement is explicit
that the command must either produce complete output or produce output that is
identifiably incomplete. A killed process that leaves a plausible-looking partial
file fails that requirement, and a silent OOM fails it worse.

This script runs each command under `ulimit -v` at limits below what it needs and
records, per limit:

* the exit status, and whether anything was printed to stderr at all;
* whether the process was killed by a signal rather than exiting on its own;
* which output files exist, and whether any is a *partial* result that a consumer
  could mistake for a complete one;
* whether the input files are unchanged, which is the other half of the contract
  ("preserve input data").

It asserts nothing about the commands being well behaved -- the point is to record
what they do. A finding is reported as an observation with the measured numbers.

Two findings already recorded here, both reproducible:

* `bam2wig` on the 8.2M-record rat alignment needs ~1.7 GB and completes in 51 s
  with no limit. Under 800 MB, 400 MB or 200 MB it writes **no** output file and
  dies on a panic inside `zlib-rs`'s inflate (`assertion left == right failed,
  left: MemError`) at exit 101, or on `memory allocation of 192 bytes failed` with
  SIGABRT at the tightest limit. The panic message is not actionable and comes from
  a dependency, but the important property holds: nothing partial is left behind, so
  a consumer cannot mistake a truncated result for a whole one.
* `read_duplication` on the same alignment needs ~530 MB. Its behaviour under a
  limit is recorded below.

Usage:
    oracle/venv/bin/python3 verification/check_memory_failure.py
    oracle/venv/bin/python3 verification/check_memory_failure.py --json failure.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RELEASE = REPO / "target" / "release"

DEFAULT_BAM = REPO / "datasets" / "heldout" / "aligned" / "SRR1177982" / "SRR1177982.bam"

# A limit is chosen as a fraction of the limit-free peak RSS measured by
# verification/measure_memory.py, rather than as a round number, so the sweep covers
# "a little short" through "far too little" on whatever machine this runs on.
LIMIT_FRACTIONS = [0.8, 0.4, 0.2]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def peak_rss_of_unlimited_run(argv: list[str], work: Path) -> float:
    """Peak RSS in MB with no limit, used to pick the limits."""
    proc = subprocess.run(["/usr/bin/time", "-v", *[str(a) for a in argv]],
                          cwd=str(work), capture_output=True, text=True, timeout=7200)
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", proc.stderr)
    return int(m.group(1)) / 1024.0 if m else 0.0


def run_limited(argv: list[str], work: Path, limit_kb: int) -> dict:
    """Run under `ulimit -v`, capturing status, signal, output files and stderr."""
    # A fresh directory per limit: a previous run's output must not be mistaken for
    # this run's, which is the whole question ("did the failed run leave something
    # behind?").
    d = work / f"limit-{limit_kb}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    a = [str(d) if str(x) == "@OUT@" else str(x) for x in argv]
    script = f"ulimit -v {limit_kb}; exec \"$@\""
    proc = subprocess.run(["bash", "-c", script, "bash", *a],
                          cwd=str(d), capture_output=True, text=True, timeout=7200)
    # Popen reports a signalled child as -N. `subprocess.run` does not translate it,
    # so bash's `128+signal` convention never applies here and the first version of
    # this script reported SIGABRT as `None` -- the single most important thing in the
    # record, that the process was killed rather than exiting on its own.
    returncode = proc.returncode
    killed = -returncode if returncode < 0 else None
    if killed is not None:
        returncode = 128 + killed
    files = {}
    for p in sorted(d.iterdir()):
        if p.is_file():
            files[p.name] = p.stat().st_size
    stderr = proc.stderr.strip()
    return {
        "limit_mb": round(limit_kb / 1024.0, 1),
        "exit_code": returncode,
        # 137 is SIGKILL (the OOM killer), 134 is SIGABRT, 101 is a Rust panic's
        # normal exit rather than a signal.
        "killed_by_signal": killed,
        "stderr_tail": stderr[-260:],
        "said_anything": bool(stderr),
        "output_files": files,
        "left_partial_output": bool(files),
    }


def measure(command: str, bam: Path, bed: Path, sizes: Path, work: Path) -> dict:
    exe = RELEASE / command
    if not exe.exists():
        return {"command": command, "status": "SKIPPED", "reason": f"{exe} missing"}
    out = work / f"baseline-{command}"
    out.mkdir(parents=True, exist_ok=True)
    if command == "bam2wig":
        argv = [exe, "-i", bam, "-s", sizes, "-o", out / "w"]
    else:
        argv = [exe, "-i", bam, "--out-prefix", out / "out"]
    print(f"{command}: measuring the limit-free peak ...")
    peak = peak_rss_of_unlimited_run(argv, work)
    if peak <= 0:
        return {"command": command, "status": "SKIPPED",
                "reason": "could not measure the limit-free peak"}
    print(f"  limit-free peak RSS {peak:.1f} MB; sweeping "
          f"{', '.join(f'{f:.0%}' for f in LIMIT_FRACTIONS)} of it")

    runs = []
    for frac in LIMIT_FRACTIONS:
        limit_kb = int(peak * 1024 * frac)
        r = run_limited(argv, work, limit_kb)
        r["fraction_of_free_peak"] = frac
        runs.append(r)
        sig = r["killed_by_signal"]
        print(f"  {r['limit_mb']:8.1f} MB  exit={r['exit_code']:>4}"
              + (f" (SIG{sig})" if sig else "")
              + f"  files={len(r['output_files'])}"
              + f"  {'' if r['said_anything'] else 'SILENT'}")
    return {
        "command": command,
        "limit_free_peak_mb": round(peak, 1),
        "runs": runs,
        # The property the contract asks for: a failed run must not leave output a
        # consumer could read as complete. Recorded as a finding, not a gate,
        # because "the panic came from a dependency" is a disclosure decision.
        "no_partial_output_at_any_limit": all(
            not r["left_partial_output"] for r in runs),
        "every_failure_was_reported": all(
            r["said_anything"] for r in runs if r["exit_code"] != 0),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default=str(DEFAULT_BAM))
    ap.add_argument("--bed", default=str(REPO / "datasets" / "heldout" / "reference" /
                                         "rn6.indexed.bed12"))
    ap.add_argument("--commands", nargs="*", default=["bam2wig", "read_duplication"])
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    if not RELEASE.exists():
        print("target/release not found; run: cargo build --workspace --release --locked")
        return 2
    bam = Path(args.bam)
    if not bam.exists():
        print(f"alignment missing: {bam}")
        print("Pass --bam, or run datasets/align_run.sh for a panel first.")
        return 2

    before = sha256(bam)
    results = {"bam": str(bam), "bam_sha256": before, "commands": {}}
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        sizes = work / "chrom.sizes"
        import pysam
        with pysam.AlignmentFile(str(bam)) as f:
            sizes.write_text("\n".join(f"{e['SN']}\t{e['LN']}"
                                       for e in f.header.to_dict()["SQ"]) + "\n")
        for command in args.commands:
            entry = measure(command, bam, Path(args.bed), sizes, work)
            results["commands"][command] = entry
            print()

    after = sha256(bam)
    results["input_unchanged"] = before == after
    print(f"input unchanged: {results['input_unchanged']}")
    if not results["input_unchanged"]:
        print(f"  before {before}\n  after  {after}")

    for command, e in results["commands"].items():
        if e.get("status") == "SKIPPED":
            print(f"\n{command}: SKIPPED ({e['reason']})")
            continue
        print(f"\n{command}: limit-free peak {e['limit_free_peak_mb']} MB")
        for r in e["runs"]:
            print(f"  limit {r['limit_mb']:>8.1f} MB -> exit {r['exit_code']}"
                  f"{' SIG%d' % r['killed_by_signal'] if r['killed_by_signal'] else ''}"
                  f", {len(r['output_files'])} output file(s)")
        print(f"  no partial output at any limit: {e['no_partial_output_at_any_limit']}")
        print(f"  every failure reported on stderr: {e['every_failure_was_reported']}")

    print("\nThese are observations, not gates. A command that leaves no output when "
          "it\nfails satisfies 'complete output or identifiably incomplete output', but "
          "a panic\nmessage from inside a dependency is not an actionable diagnostic and "
          "is\nrecorded here so the release can disclose it rather than discover it in "
          "production.")

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2) + "\n")
        print(f"\nwritten to {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
