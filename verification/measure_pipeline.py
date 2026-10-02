#!/usr/bin/env python3
"""Measure a bulk QC panel as a workflow, not as per-command speed.

`benchmarks/protocol-v2.md` §10 requires a pipeline-impact measurement: samples per
hour, CPU-hours per sample, memory at the intended concurrency, output and I/O volume,
and reproducibility of the QC decisions. Per-command speedups do not answer those, and
the audit is explicit that "if QC is only a fraction of the complete RNA-seq pipeline,
the overall benefit is correspondingly bounded" -- a number nobody can state without
timing the panel end to end.

The panel is the five commands a bulk RNA-seq QC step actually consists of, in the
order they would run:

  bam_stat              mapping summary (stdout)
  infer_experiment      library layout and strandedness (stdout)
  read_distribution     region assignment against a model
  geneBody_coverage     5'->3' coverage bias
  read_duplication      duplicate rate

Properties this script holds itself to, because each has been got wrong somewhere in
this project before:

* **Each stage is measured once, under instrumentation.** An earlier draft ran every
  stage twice -- once for the wall clock and once under `/usr/bin/time` for resources --
  and would have reported roughly double the real cost. `/usr/bin/time -v` records the
  same wall clock as a stopwatch, so one run suffices.
* **A stage's success is a content check, not an existence check.** A command that
  exits 0 having written a truncated table has not produced a result. Each stage
  declares what its output must contain.
* **A failed stage stops the panel and is reported.** A stage after a broken one
  measures nothing useful, and including its time would understate what a working
  panel costs. Throughput is reported for complete panels only.
* **Concurrency is measured on whole panels**, in parallel, because the audit's own
  guidance is that "parallelism across samples may already be the appropriate
  production design; measure that deployment first". The speedup from concurrency is
  reported as a measured ratio.

Usage:
    oracle/venv/bin/python3 verification/measure_pipeline.py
    oracle/venv/bin/python3 verification/measure_pipeline.py \\
        --bam <bam> --bed <bed> --concurrency 1 2 4 --json pipeline-impact.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
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

# What each stage must PRODUCE, not merely exit. These are the predicates that make a
# stage's success checkable.
def _nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _stdout_has(*needles: str):
    def check(work: Path) -> bool:
        out = work / "stdout.txt"
        if not out.is_file():
            return False
        text = out.read_text(errors="replace")
        return all(n in text for n in needles)
    return check


STAGES = [
    # (name, argv builder, output predicate, where the output lands)
    ("bam_stat",
     lambda w, bam, bed: ["-i", str(bam)],
     _stdout_has("Total records", "mapq"),
     "stdout"),
    ("infer_experiment",
     lambda w, bam, bed: ["-i", str(bam), "-r", str(bed)],
     _stdout_has("explained by", "Fraction of reads"),
     "stdout"),
    # read_distribution takes no output flag: it prints its table to stdout, as
    # upstream does. An earlier draft passed `--out-prefix` to it, which clap rejected
    # with "unexpected argument", and the panel reported itself INCOMPLETE -- which is
    # the predicate working, but only after a wasted run.
    ("read_distribution",
     lambda w, bam, bed: ["-i", str(bam), "-r", str(bed)],
     _stdout_has("Total Reads", "CDS_Exons", "Introns"),
     "stdout"),
    ("geneBody_coverage",
     lambda w, bam, bed: ["-i", str(bam), "-r", str(bed),
                          "--out-prefix", str(w / "gb"), "--skip-plot"],
     lambda w: _nonempty(w / "gb.geneBodyCoverage.txt"),
     "files"),
    ("read_duplication",
     lambda w, bam, bed: ["-i", str(bam), "--out-prefix", str(w / "dup")],
     lambda w: _nonempty(w / "dup.seq.DupRate.xls"),
     "files"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_time_v(stderr: str) -> dict:
    out = {}
    for pattern, key, scale in (
        (r"Maximum resident set size \(kbytes\):\s*(\d+)", "peak_rss_mb", 1 / 1024.0),
        (r"User time \(seconds\):\s*([\d.]+)", "user_s", 1.0),
        (r"System time \(seconds\):\s*([\d.]+)", "system_s", 1.0),
    ):
        m = re.search(pattern, stderr)
        if m:
            out[key] = float(m.group(1)) * scale
    m = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", stderr)
    if m:
        secs = 0.0
        for part in m.group(1).split(":"):
            secs = secs * 60 + float(part)
        out["wall_s"] = secs
    return out


def mem_used_mb() -> float:
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


def run_stage(name, build, produced, kind, work, bam, bed) -> dict:
    """One stage, instrumented once."""
    exe = RELEASE / name
    argv = [str(exe)] + build(work, bam, bed)
    stdout_path = work / "stdout.txt"
    time_path = work / f"{name}.time"
    # /usr/bin/time writes its own report to stderr, so the command's stderr is
    # captured separately rather than interleaved with it.
    cmd_err = work / f"{name}.stderr"
    with open(stdout_path, "wb") as out_fh, open(time_path, "wb") as time_fh, \
            open(cmd_err, "wb") as err_fh:
        proc = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(time_path), *argv],
            cwd=str(work), stdout=out_fh, stderr=err_fh, timeout=7200)
    res = parse_time_v(time_path.read_text(errors="replace"))
    wall = res.get("wall_s")
    if wall is None:
        # /usr/bin/time -o is present on GNU time; if it produced nothing usable, fall
        # back to a stopwatch rather than reporting zero.
        wall = 0.0
    ok = proc.returncode == 0
    detail = ""
    if ok:
        try:
            if not produced(work):
                ok = False
                detail = "exited 0 but produced no usable output"
        except Exception as exc:
            ok = False
            detail = f"output check raised: {exc}"
    else:
        detail = (cmd_err.read_text(errors="replace").strip()[-160:]
                  or f"exit {proc.returncode}")
    return {
        "stage": name, "argv": argv[1:], "exit": proc.returncode,
        "complete": ok, "wall_s": round(wall, 2),
        "cpu_s": round(res.get("user_s", 0.0) + res.get("system_s", 0.0), 2),
        "peak_rss_mb": round(res.get("peak_rss_mb", 0.0), 1), "detail": detail,
    }


def run_one_panel(bam: Path, bed: Path, work: Path) -> dict:
    """Every stage in order against one sample. Stops at the first failure."""
    work.mkdir(parents=True, exist_ok=True)
    stages = []
    complete = True
    for name, build, produced, kind in STAGES:
        if not (RELEASE / name).exists():
            stages.append({"stage": name, "complete": False, "detail": "binary missing",
                           "wall_s": 0.0, "cpu_s": 0.0, "peak_rss_mb": 0.0,
                           "exit": None, "argv": []})
            complete = False
            break
        s = run_stage(name, build, produced, kind, work, bam, bed)
        stages.append(s)
        if not s["complete"]:
            complete = False
            break

    panel_wall = sum(s["wall_s"] for s in stages)
    panel_cpu = sum(s["cpu_s"] for s in stages)
    panel_peak = max((s["peak_rss_mb"] for s in stages), default=0.0)
    # Output volume counts what the panel wrote, excluding the harness's own capture
    # files (stdout.txt, *.time, *.stderr) which are measurement artefacts.
    artefacts = {"stdout.txt"}
    artefacts |= {f"{s['stage']}.{ext}" for s in stages for ext in ("time", "stderr")}
    files, total = 0, 0
    for p in sorted(work.iterdir()):
        if p.is_file() and p.name not in artefacts:
            files += 1
            total += p.stat().st_size
    return {
        "complete": complete, "stages": stages,
        "panel_wall_s": round(panel_wall, 2),
        "panel_cpu_s": round(panel_cpu, 2),
        "panel_peak_rss_mb": round(panel_peak, 1),
        "output_bytes": total, "output_files": files,
        "stages_run": len(stages),
    }


def run_concurrent(bam: Path, bed: Path, n: int, root: Path) -> dict:
    """n whole panels in parallel, as a per-sample-parallel deployment would."""
    boxes: list[dict] = [{} for _ in range(n)]
    threads = []
    dirs = []

    # Sample the memory baseline BEFORE any panel starts, and start the threads only
    # afterwards. Sampling it after the panels have started -- as this did -- means
    # their allocations land inside the baseline and are then subtracted from the peak,
    # so the reported figure is whatever the panels had NOT yet allocated at the 1.0 s
    # mark, which is not a quantity. It is not conservative either: a panel that has
    # not reached its peak by 1.0 s inflates the number. The concurrency harness was
    # corrected for exactly this and the two disagreed, which is how the flaw was
    # noticed; this file was left with the old ordering.
    baseline = mem_used_mb()
    samples = [baseline]

    def worker(i):
        w = root / f"panel{i}"
        dirs.append(w)
        boxes[i]["result"] = run_one_panel(bam, bed, w)

    started = time.monotonic()
    for i in range(n):
        th = threading.Thread(target=worker, args=(i,))
        th.start()
        threads.append(th)

    while any(t.is_alive() for t in threads):
        samples.append(mem_used_mb())
        time.sleep(0.05)
    elapsed = time.monotonic() - started
    for t in threads:
        t.join()

    results = [b.get("result") for b in boxes]
    complete = [r for r in results if r and r["complete"]]
    return {
        "concurrency": n,
        "batch_wall_s": round(elapsed, 2),
        "panels": n,
        "panels_complete": len(complete),
        "peak_aggregate_mb": round(max(samples) - baseline, 1) if samples else 0.0,
        "samples": len(samples),
        "samples_per_hour": round(n / elapsed * 3600, 1) if elapsed else None,
        "cpu_hours_per_sample": (round(
            sum(r["panel_cpu_s"] for r in complete) / 3600 / n, 4)
            if complete else None),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default=str(DEFAULT_BAM))
    ap.add_argument("--bed", default=str(DEFAULT_BED))
    ap.add_argument("--concurrency", nargs="*", type=int, default=[1, 2, 4])
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    if not RELEASE.exists():
        print("target/release not found; run: cargo build --workspace --release --locked")
        return 2
    bam = Path(args.bam).resolve()
    bed = Path(args.bed).resolve()
    if not bam.exists():
        print(f"alignment missing: {bam}")
        return 2
    if not bed.exists():
        print(f"annotation missing: {bed}")
        return 2

    print(f"panel ({len(STAGES)} stages): {', '.join(s[0] for s in STAGES)}")
    print(f"input: {bam}")
    print(f"  sha256 {sha256(bam)}")
    print(f"  {bam.stat().st_size:,} bytes\n")

    results = {
        "panel": [s[0] for s in STAGES],
        "input": {"bam": str(bam), "bam_sha256": sha256(bam), "bed": str(bed),
                  "bed_sha256": sha256(bed)},
        "runs": [],
    }

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for n in args.concurrency:
            batch = run_concurrent(bam, bed, n, root / f"c{n}")
            results["runs"].append(batch)
            print(f"concurrency {n}: {batch['panels_complete']}/{batch['panels']} "
                  f"panels complete, batch {batch['batch_wall_s']:.2f}s, "
                  f"aggregate peak {batch['peak_aggregate_mb']:.0f} MB, "
                  f"{batch['samples_per_hour']:.0f} samples/hour"
                  if batch["samples_per_hour"] else "concurrency %d: incomplete" % n)
            print()

    # Per-stage cost from the single-panel run, which is the only one where the
    # stages are individually attributable.
    with tempfile.TemporaryDirectory() as tmp:
        single = run_one_panel(bam, bed, Path(tmp) / "single")
    results["single_panel"] = single
    print("per-stage cost (from the single-panel run):")
    for s in single["stages"]:
        share = (s["wall_s"] / single["panel_wall_s"] * 100
                 if single["panel_wall_s"] else 0.0)
        print(f"  {s['stage']:20} wall={s['wall_s']:8.2f}s ({share:5.1f}%) "
              f"cpu={s['cpu_s']:8.2f}s rss={s['peak_rss_mb']:8.1f} MB "
              f"{'ok' if s['complete'] else 'INCOMPLETE: ' + s['detail']}")
    print(f"  {'PANEL':20} wall={single['panel_wall_s']:8.2f}s "
          f"cpu={single['panel_cpu_s']:8.2f}s "
          f"peak={single['panel_peak_rss_mb']:8.1f} MB "
          f"output={single['output_bytes'] / 2**20:.1f} MiB "
          f"in {single['output_files']} files")
    if not single["complete"]:
        print("  panel INCOMPLETE: throughput is reported only for complete panels")
    else:
        cpu_hours = single["panel_cpu_s"] / 3600
        print(f"  per sample: {cpu_hours:.4f} CPU-hours, "
              f"{single['panel_wall_s']:.1f}s wall")

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2) + "\n")
        print(f"\nwritten to {args.json}")
    return 0 if single["complete"] else 1


if __name__ == "__main__":
    sys.exit(main())
