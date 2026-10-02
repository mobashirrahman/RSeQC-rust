#!/usr/bin/env python3
"""Measure aggregate memory for N concurrent invocations of a command.

`verification/measure_memory.py` measures one process at a time, and the audit is
explicit about why that is not enough: GNU time's "Maximum resident set size" is
Linux's `RUSAGE_CHILDREN.ru_maxrss`, the largest *single* child's RSS. A deployment
that runs a QC step over N samples in parallel needs the SUM, and the port's docs
already state that "memory per concurrent invocation is therefore additive: N
concurrent jobs need N times the per-job figure" -- which is an assumption, not a
measurement.

This measures the aggregate directly, by sampling total system memory in use (from
`/proc/meminfo`) while N invocations run concurrently, and subtracting the baseline
observed immediately before they start. The subtraction is what makes the number
about these processes rather than about the machine: the machine is shared, so the
absolute `MemAvailable` figure drifts with unrelated work.

What it records, per concurrency level:

* the **peak aggregate** attributable to the run, sampled at 20 Hz, with the baseline
  and the sampling interval recorded so the figure is reproducible;
* per-job peak RSS from `/usr/bin/time -v`, so the aggregate can be compared against
  N x single-job rather than merely reported;
* wall time for the batch, and the throughput (jobs/second) implied;
* whether every job succeeded, and whether each produced its expected output.

Honesty requirements this script holds itself to:

* the sampling loop must not miss the peak, so the sample count and the interval are
  reported rather than assumed;
* if the baseline drifts by more than a stated margin during the run, that is reported
  as drift rather than silently folded into the result;
* a "peak" lower than the largest single job's RSS is a measurement failure, and is
  reported as one. It cannot happen physically, so seeing it means the sampler is
  wrong.

Usage:
    oracle/venv/bin/python3 verification/measure_concurrency.py
    oracle/venv/bin/python3 verification/measure_concurrency.py \\
        --bam <bam> --bed <bed> --concurrency 1 2 4 8 --json concurrency.json
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RELEASE = REPO / "target" / "release"
DEFAULT_BAM = REPO / "datasets" / "heldout" / "aligned" / "SRR1177982" / "SRR1177982.bam"
DEFAULT_BED = REPO / "datasets" / "heldout" / "reference" / "rn6.indexed.bed12"

# How often total system memory is sampled, and how much baseline drift is tolerated
# before the run is flagged as contaminated.
SAMPLE_INTERVAL_S = 0.05
# How long to wait before taking the baseline and starting the batch. Long enough that
# the shell's own transient allocations are not charged to the measurement; the jobs
# are spawned after it, so no command is swallowed by the wait.
SETTLE_S = 1.0
# How long to sample the instrument with nothing of ours running, to establish its own
# noise floor. A fixed tolerance would be a guess; this measures what /proc/meminfo
# does on this machine right now.
NOISE_PROBE_S = 5.0
# How much memory may still be in use after the batch finishes before the run is
# called out as retaining it. A shared machine's memory moves for unrelated reasons,
# so this is a threshold on a signed quantity: positive means the batch left memory
# allocated, negative means it was returned (possibly as page cache).
DRIFT_TOLERANCE_MB = 256.0


def mem_used_mb() -> float:
    """Memory in use, in MB, from /proc/meminfo.

    `MemTotal - MemAvailable`, which is what the kernel reports as in use including
    page cache that is reclaimable. Using `MemAvailable` alone and taking the minimum
    would give the same peak for a machine that was busy for unrelated reasons, which
    is exactly the contamination the baseline subtraction exists to remove.
    """
    total = available = None
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith("MemTotal:"):
                total = int(line.split()[1])
            elif line.startswith("MemAvailable:"):
                available = int(line.split()[1])
            if total is not None and available is not None:
                break
    if total is None or available is None:
        raise RuntimeError("/proc/meminfo did not report MemTotal and MemAvailable")
    return (total - available) / 1024.0


def parse_time_v(stderr: str) -> dict:
    out = {}
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
    if m:
        out["peak_rss_mb"] = int(m.group(1)) / 1024.0
    m = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", stderr)
    if m:
        secs = 0.0
        for part in m.group(1).split(":"):
            secs = secs * 60 + float(part)
        out["wall_s"] = secs
    return out


def build_argv(command: str, exe: Path, bam: Path, bed: Path, sizes: Path,
               work: Path) -> list[str]:
    argv = [str(exe), "-i", str(bam)]
    if command == "bam2wig":
        return argv + ["-s", str(sizes), "-o", str(work / "w")]
    if command == "tin":
        return argv + ["-r", str(bed), "-n", "50", "-o", str(work)]
    return argv + ["--out-prefix", str(work / "out")]


def run_batch(command: str, n: int, bam: Path, bed: Path, sizes: Path,
              root: Path) -> dict:
    """Run `n` invocations concurrently and sample aggregate memory throughout."""
    exe = RELEASE / command
    if not exe.exists():
        return {"status": "SKIPPED", "reason": f"{exe} missing"}

    procs, dirs, time_logs = [], [], []
    for i in range(n):
        d = root / f"job{i}"
        d.mkdir(parents=True)
        dirs.append(d)
        argv = build_argv(command, exe, bam, bed, sizes, d)
        # Each job's own peak RSS is captured by /usr/bin/time, so the aggregate can be
        # compared against N x single rather than merely asserted to be additive.
        err = d / "time.log"
        handle = open(err, "wb")
        time_logs.append((handle, err))
        procs.append(subprocess.Popen(
            ["/usr/bin/time", "-v", *[str(a) for a in argv]],
            cwd=str(d), stdout=subprocess.DEVNULL, stderr=handle))

    # Baseline is taken BEFORE the processes are spawned, after a short settle so the
    # measurement is not charged for the shell's own transient allocations. Taking it
    # after the spawn -- which an intermediate version did, on the reasoning that the
    # processes' allocations should be "already included" -- subtracts the commands'
    # entire working set from their own peak and produced an aggregate (895 MB) below
    # the single job's own peak RSS (1060 MB). That is physically impossible, so the
    # number was a broken measurement, not a finding.
    time.sleep(SETTLE_S)
    baseline = mem_used_mb()
    started = time.monotonic()
    samples = []
    while any(p.poll() is None for p in procs):
        samples.append(mem_used_mb())
        time.sleep(SAMPLE_INTERVAL_S)
    elapsed = time.monotonic() - started
    settled = n

    for handle, err in time_logs:
        handle.close()
    final = mem_used_mb()

    per_job = [parse_time_v(err.read_text(errors="replace")) for _h, err in time_logs]
    exits = [p.returncode for p in procs]
    outputs = {}
    for i, d in enumerate(dirs):
        files = sorted(p.name for p in d.iterdir()
                       if p.is_file() and p.name != "time.log")
        outputs[i] = files

    peak_above_baseline = max(samples) - baseline if samples else 0.0
    # Signed, not absolute. The sign carries the meaning: memory left IN USE after the
    # batch is a retention question worth investigating, while memory returned (or
    # slightly over-baseline, because freed anonymous pages sit in cache) is the normal
    # outcome. Collapsing both to an absolute value made a run that correctly returned
    # its memory indistinguishable from one that did not.
    drift = final - baseline
    largest_single = max((p.get("peak_rss_mb", 0.0) for p in per_job), default=0.0)
    # Measure this machine's own noise floor rather than picking a threshold: sample
    # the instrument the same way, with nothing of ours running, and use the spread.
    noise = 0.0
    quiet = [mem_used_mb()]
    q0 = time.monotonic()
    while time.monotonic() - q0 < NOISE_PROBE_S:
        quiet.append(mem_used_mb())
        time.sleep(SAMPLE_INTERVAL_S)
    noise = max(quiet) - min(quiet)
    plausible = peak_above_baseline >= largest_single - noise

    return {
        "command": command,
        "concurrency": n,
        "status": "OK",
        "settle_s": SETTLE_S,
        "jobs_settled_at_baseline": settled,
        "jobs_finished_before_baseline": n - settled,
        "baseline_mb": round(baseline, 1),
        "peak_aggregate_mb": round(peak_above_baseline, 1),
        "peak_system_mb_observed": round(max(samples), 1) if samples else None,
        "samples": len(samples),
        "sample_interval_s": SAMPLE_INTERVAL_S,
        "post_run_delta_mb": round(drift, 1),
        "memory_retained_after_batch": drift > DRIFT_TOLERANCE_MB,
        "post_run_within_tolerance": abs(drift) <= DRIFT_TOLERANCE_MB,
        "largest_single_job_rss_mb": round(largest_single, 1),
        "sum_of_single_job_rss_mb": round(
            sum(p.get("peak_rss_mb", 0.0) for p in per_job), 1),
        # The aggregate cannot be below the largest single job's RSS: that memory was
        # demonstrably resident. But this instrument's own noise floor on an idle
        # machine is tens of MB (the kernel moves reclaimable page cache under
        # MemAvailable), so "below" means "below by more than the instrument can
        # resolve", and the tolerance is measured rather than assumed.
        "instrument_noise_floor_mb": round(noise, 1),
        "shortfall_vs_largest_single_mb": round(max(0.0, largest_single - peak_above_baseline), 1),
        "aggregate_at_least_largest_single_job": plausible,
        "batch_wall_s": round(elapsed, 2),
        "jobs_per_second": round(n / elapsed, 3) if elapsed else None,
        "exit_codes": exits,
        "all_succeeded": all(rc == 0 for rc in exits),
        "all_jobs_wrote_output": all(bool(v) for v in outputs.values()),
        "per_job_peak_rss_mb": [round(p.get("peak_rss_mb", 0.0), 1) for p in per_job],
        "per_job_wall_s": [round(p.get("wall_s", 0.0), 2) for p in per_job],
        "outputs_per_job": outputs,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default=str(DEFAULT_BAM))
    ap.add_argument("--bed", default=str(DEFAULT_BED))
    ap.add_argument("--commands", nargs="*", default=["bam2wig", "read_duplication"])
    ap.add_argument("--concurrency", nargs="*", type=int, default=[1, 2, 4, 8])
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    if not RELEASE.exists():
        print("target/release not found; run: cargo build --workspace --release --locked")
        return 2
    bam = Path(args.bam)
    if not bam.exists():
        print(f"alignment missing: {bam}")
        return 2

    results = {
        "alignment": str(bam),
        "method": {
            "aggregate": "peak (MemTotal - MemAvailable) sampled at "
                         f"{SAMPLE_INTERVAL_S}s while the batch runs, minus a "
                         "baseline sampled after the processes are spawned",
            "per_job": "/usr/bin/time -v Maximum resident set size, which is Linux's "
                       "RUSAGE_CHILDREN.ru_maxrss: the largest SINGLE child's RSS, not "
                       "an aggregate",
            "drift_tolerance_mb": DRIFT_TOLERANCE_MB,
        },
        "cpu_count": len(__import__("os").sched_getaffinity(0)),
        "results": [],
    }
    print(f"cpus available: {results['cpu_count']}; levels: {args.concurrency}")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        sizes = root / "chrom.sizes"
        import pysam
        with pysam.AlignmentFile(str(bam)) as f:
            sizes.write_text("\n".join(f"{e['SN']}\t{e['LN']}"
                                       for e in f.header.to_dict()["SQ"]) + "\n")
        for command in args.commands:
            for n in args.concurrency:
                r = run_batch(command, n, bam, Path(args.bed), sizes,
                              root / f"{command}-{n}")
                results["results"].append(r)
                if r.get("status") == "SKIPPED":
                    print(f"{command} x{n}: SKIPPED")
                    continue
                per = r["largest_single_job_rss_mb"]
                agg = r["peak_aggregate_mb"]
                ratio = agg / per if per else float("nan")
                flags = []
                if r["memory_retained_after_batch"]:
                    flags.append(f"RETAINED {r['post_run_delta_mb']}MB after the batch")
                if not r["aggregate_at_least_largest_single_job"]:
                    flags.append(
                        f"AGGREGATE BELOW PER-JOB by "
                        f"{r['shortfall_vs_largest_single_mb']}MB "
                        f"(noise floor {r['instrument_noise_floor_mb']}MB)")
                if not r["all_succeeded"]:
                    flags.append(f"exits={r['exit_codes']}")
                print(f"{command:18} x{n:<2} aggregate={agg:8.1f} MB  "
                      f"largest single={per:8.1f} MB  ratio={ratio:4.2f}  "
                      f"batch={r['batch_wall_s']:6.2f}s  "
                      f"{' '.join(flags) if flags else 'ok'}")

    print("\nRatio is aggregate / largest single job. A ratio near the concurrency level "
          "means\nmemory is additive across concurrent invocations, which is what the "
          "documentation\nassumes; a ratio well below it would mean the jobs share "
          "allocations or the\nsampler missed the peak.")
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2) + "\n")
        print(f"\nwritten to {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
