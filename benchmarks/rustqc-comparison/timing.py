#!/usr/bin/env python3
"""Card R5: timing harness for the three-arm comparison.

Configurations: U-seq, P-seq, P-par, R-1, R-8 (see RUSTQC_COMPARISON_PLAN.md).
Measurement via /usr/bin/time -v per process; P-par aggregate via
/proc/meminfo sampling (verification/measure_concurrency.py method:
baseline taken BEFORE anything starts).

Usage:
    python3 timing.py --workload W-rat-3c --reps 1 --reps-upstream 1 --seed 0 --outdir raw/ [--timeout 3600]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import re
import signal
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "verification"))
from equivalence import check_presence  # noqa: E402

WORKLOADS = HERE / "workloads.json"
TOOLS_JSON = HERE / "tools" / "tools.json"
UPSTREAM_SCRIPTS = REPO / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHONPATH = str(REPO / "oracle" / "upstream-src" / "src")
ORACLE_PYTHON = REPO / "oracle" / "venv" / "bin" / "python3"
RELEASE = REPO / "target" / "release"
RUSTQC = HERE / "tools" / "rustqc"
TIME_BIN = "/usr/bin/time"

PINNED_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0",
}

CONFIGS = ["U-seq", "P-seq", "P-par", "R-1", "R-8"]


def git_info():
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(REPO),
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        head = "unknown"
    try:
        dirty = bool(subprocess.run(["git", "status", "--short"], cwd=str(REPO),
                                    capture_output=True, text=True).stdout.strip())
    except Exception:
        dirty = None
    return head, dirty


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_time_v(stderr: str) -> dict:
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)", stderr)
    if m:
        parts = [float(x) for x in m.group(1).split(":")]
        secs = 0.0
        for p in parts:
            secs = secs * 60 + p
        out["time_wall_s"] = secs
    m = re.search(r"User time \(seconds\):\s*([\d.]+)", stderr)
    if m:
        out["user_s"] = float(m.group(1))
    m = re.search(r"System time \(seconds\):\s*([\d.]+)", stderr)
    if m:
        out["sys_s"] = float(m.group(1))
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
    if m:
        out["peak_rss_mb"] = int(m.group(1)) / 1024.0
    return out


def _kill_group(pgid, grace_s=2.0):
    for sig, wait in ((signal.SIGTERM, grace_s), (signal.SIGKILL, 1.0)):
        if not pgid:
            return
        try:
            os.killpg(pgid, sig)
        except (ProcessLookupError, PermissionError):
            return
        deadline = time.monotonic() + wait
        while time.monotonic() < deadline:
            try:
                os.killpg(pgid, 0)
            except ProcessLookupError:
                return
            except PermissionError:
                return
            time.sleep(0.05)


def _group_alive(pgid):
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def run_timed(argv, cwd, env, timeout_s):
    """Run argv under /usr/bin/time -v in its own process group.

    Retains the Popen (never killpg(getpgid(0))): on timeout only the
    created group is signalled. Returns (measurements, stdout, stderr,
    exit_code, timed_out).
    """
    cmd = [TIME_BIN, "-v", *argv]
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="timed-cap-") as capdir:
        out_path = Path(capdir) / "stdout"
        err_path = Path(capdir) / "stderr"
        with out_path.open("wb") as out_fh, err_path.open("wb") as err_fh:
            proc = subprocess.Popen(cmd, cwd=str(cwd), env=env,
                                    stdout=out_fh, stderr=err_fh,
                                    start_new_session=True)
            pgid = proc.pid
            timed_out = False
            try:
                proc.wait(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                timed_out = True
                _kill_group(pgid)
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
        wall = time.perf_counter() - t0
        _kill_group(pgid)
        stdout = out_path.read_bytes().decode("utf8", "replace")
        stderr = err_path.read_bytes().decode("utf8", "replace")
    m = parse_time_v(stderr)
    m["wall_s"] = wall
    return m, stdout, stderr, proc.returncode, timed_out


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
    return (total - available) / 1024.0


def tool_argvs(arm: str, bam: str, bed: str, outdir: Path, sample: str):
    """Eight per-tool argv lists for U-seq / P-seq / P-par."""
    if arm == "U":
        py = str(ORACLE_PYTHON)
        return [
            ("bam_stat", [py, str(UPSTREAM_SCRIPTS / "bam_stat.py"), "-i", bam]),
            ("infer_experiment", [py, str(UPSTREAM_SCRIPTS / "infer_experiment.py"), "-i", bam, "-r", bed]),
            ("read_duplication", [py, str(UPSTREAM_SCRIPTS / "read_duplication.py"), "-i", bam, "-o", str(outdir / "dup"), "--skip-plot"]),
            ("read_distribution", [py, str(UPSTREAM_SCRIPTS / "read_distribution.py"), "-i", bam, "-r", bed]),
            ("junction_annotation", [py, str(UPSTREAM_SCRIPTS / "junction_annotation.py"), "-i", bam, "-r", bed, "-o", str(outdir / "ja"), "--skip-plot"]),
            ("junction_saturation", [py, str(UPSTREAM_SCRIPTS / "junction_saturation.py"), "-i", bam, "-r", bed, "-o", str(outdir / "js"), "--skip-plot"]),
            ("inner_distance", [py, str(UPSTREAM_SCRIPTS / "inner_distance.py"), "-i", bam, "-r", bed, "-o", str(outdir / "id"), "--skip-plot"]),
            ("tin", [py, str(UPSTREAM_SCRIPTS / "tin.py"), "-i", bam, "-r", bed, "-n", "100", "-o", str(outdir / "tin")]),
        ]
    else:
        return [
            ("bam_stat", [str(RELEASE / "bam_stat"), "-i", bam]),
            ("infer_experiment", [str(RELEASE / "infer_experiment"), "-i", bam, "-r", bed]),
            ("read_duplication", [str(RELEASE / "read_duplication"), "-i", bam, "-o", str(outdir / "dup"), "--skip-plot"]),
            ("read_distribution", [str(RELEASE / "read_distribution"), "-i", bam, "-r", bed]),
            ("junction_annotation", [str(RELEASE / "junction_annotation"), "-i", bam, "-r", bed, "-o", str(outdir / "ja"), "--skip-plot"]),
            ("junction_saturation", [str(RELEASE / "junction_saturation"), "-i", bam, "-r", bed, "-o", str(outdir / "js"), "--skip-plot"]),
            ("inner_distance", [str(RELEASE / "inner_distance"), "-i", bam, "-r", bed, "-o", str(outdir / "id"), "--skip-plot"]),
            ("tin", [str(RELEASE / "tin"), "-i", bam, "-r", bed, "-n", "100", "-o", str(outdir / "tin")]),
        ]


def run_seq(config: str, bam: str, bed: str, work_root: Path, timeout_s: float):
    """Run eight tools one after another. Returns (per_tool_rows, total)."""
    arm = "U" if config == "U-seq" else "P"
    outdir = work_root / config
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "tin").mkdir(exist_ok=True)
    env = dict(os.environ)
    env.update(PINNED_ENV)
    if arm == "U":
        env["PYTHONPATH"] = ORACLE_PYTHONPATH
    per_tool = []
    t0 = time.perf_counter()
    for name, argv in tool_argvs(arm, bam, bed, outdir, ""):
        load = os.getloadavg()
        m, _out, _err, code, to = run_timed(argv, cwd=outdir, env=env, timeout_s=timeout_s)
        per_tool.append({"tool": name, "argv": argv, "exit": code,
                         "timed_out": to, "loadavg_before": load, **m})
        if to or code != 0:
            # Continue through failures so the total still records all tools;
            # each tool's row carries its own exit/timeout.
            continue
    total_wall = time.perf_counter() - t0
    return per_tool, total_wall


def run_ppar(bam: str, bed: str, work_root: Path, timeout_s: float):
    """Eight port tools started together; batch wall + aggregate memory."""
    outdir = work_root / "P-par"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "tin").mkdir(exist_ok=True)
    env = dict(os.environ)
    env.update(PINNED_ENV)
    # Baseline BEFORE anything starts (measure_concurrency.py method).
    time.sleep(1.0)
    baseline = mem_used_mb()
    started = time.monotonic()
    procs = []
    per_tool_argv = {}
    for name, argv in tool_argvs("P", bam, bed, outdir, ""):
        # Each job under /usr/bin/time -v for its own peak RSS.
        err = outdir / f"{name}.time.log"
        handle = open(err, "wb")
        p = subprocess.Popen([TIME_BIN, "-v", *argv], cwd=str(outdir),
                             stdout=subprocess.DEVNULL, stderr=handle,
                             env=env, start_new_session=True)
        procs.append((name, p, handle, err))
        per_tool_argv[name] = argv
    samples = []
    # Sample aggregate until all finish (or timeout).
    deadline = started + timeout_s
    timed_out = False
    while any(p.poll() is None for _, p, _, _ in procs):
        samples.append(mem_used_mb())
        time.sleep(0.05)
        if time.monotonic() > deadline:
            timed_out = True
            for _, p, _, _ in procs:
                try:
                    _kill_group(p.pid)
                except Exception:
                    pass
            break
    for _, p, handle, _ in procs:
        try:
            p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
        handle.close()
    elapsed = time.monotonic() - started
    for _, p, _, _ in procs:
        _kill_group(p.pid)
    final = mem_used_mb()
    peak_agg = (max(samples) - baseline) if samples else 0.0
    per_tool = []
    for name, p, _, err in procs:
        m = parse_time_v(err.read_text(errors="replace"))
        per_tool.append({"tool": name, "argv": per_tool_argv[name],
                         "exit": p.returncode, "timed_out": timed_out, **m})
    return {"batch_wall_s": elapsed, "baseline_mb": baseline,
            "peak_aggregate_mb": peak_agg, "post_delta_mb": final - baseline,
            "samples": len(samples), "timed_out": timed_out,
            "per_tool": per_tool}


def run_rustqc(bam: str, gtf: str, work_root: Path, threads: int, timeout_s: float):
    outdir = work_root / f"R-{threads}"
    outdir.mkdir(parents=True, exist_ok=True)
    argv = [str(RUSTQC), "rna", bam, "--gtf", gtf, "--paired",
            "--outdir", str(outdir), "--threads", str(threads), "--skip-dup-check"]
    env = dict(os.environ)  # R gets only its own thread flag, not pinned pools
    load = os.getloadavg()
    m, _out, _err, code, to = run_timed(argv, cwd=outdir, env=env, timeout_s=timeout_s)
    return argv, m, code, to, load


def build_schedule(reps: int, reps_upstream: int, seed: int):
    """Blocks; inside each block configs run in random order from seed."""
    rng = random.Random(seed)
    nblocks = max(reps, reps_upstream)
    schedule = []
    for block in range(nblocks):
        cfgs = []
        if block < reps_upstream:
            cfgs.append("U-seq")
        if block < reps:
            cfgs += ["P-seq", "P-par", "R-1", "R-8"]
        rng.shuffle(cfgs)
        for order, cfg in enumerate(cfgs):
            schedule.append({"block": block, "order": order, "config": cfg})
    return schedule


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workload", required=True)
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--reps-upstream", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--timeout", type=float, default=3600.0)
    ap.add_argument("--outdir", required=True, help="raw/ directory for one JSON per run")
    args = ap.parse_args()

    data = json.loads(WORKLOADS.read_text())
    w = data["workloads"][args.workload]
    bam = str(REPO / w["bam"])
    bed = str(REPO / w["bed"])
    gtf = str(REPO / w.get("gtf_subset", w["gtf"]))
    sample = Path(bam).stem
    tools_tag = json.loads(TOOLS_JSON.read_text())["tag"] if TOOLS_JSON.exists() else "unknown"
    head, dirty = git_info()
    try:
        port_sha = sha256_file(RELEASE / "bam_stat")
    except Exception:
        port_sha = "unknown"

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # One unrecorded warm-up read before the first block.
    subprocess.run(["cat", bam], stdout=subprocess.DEVNULL, check=False)

    schedule = build_schedule(args.reps, args.reps_upstream, args.seed)
    print(f"workload {args.workload}: {len(schedule)} runs, seed {args.seed}")
    for item in schedule:
        print(f"  block {item['block']} order {item['order']}: {item['config']}")

    for item in schedule:
        cfg = item["config"]
        block = item["block"]
        with tempfile.TemporaryDirectory(prefix="timing-work-") as work:
            work_root = Path(work)
            row = {"configuration": cfg, "workload": args.workload,
                   "block": block, "order": item["order"], "seed": args.seed,
                   "timeout_s": args.timeout,
                   "git_head": head, "git_dirty": dirty,
                   "port_bam_stat_sha256": port_sha, "rustqc_tag": tools_tag,
                   "loadavg_before_run": os.getloadavg()}
            if cfg in ("U-seq", "P-seq"):
                per_tool, total_wall = run_seq(cfg, bam, bed, work_root, args.timeout)
                row["per_tool"] = per_tool
                row["total_wall_s"] = total_wall
                row["total_cpu_s"] = sum(t.get("user_s", 0) + t.get("sys_s", 0) for t in per_tool)
                row["argv"] = [t["argv"] for t in per_tool]
                row["exit"] = 0 if all(t["exit"] == 0 for t in per_tool) else 1
                row["timed_out"] = any(t["timed_out"] for t in per_tool)
                # presence check on this run's outputs
                u_dir = work_root / cfg
                # check_presence expects (u,p,r) triple; for seq runs fabricate:
                row["presence_ok"] = all((u_dir / f).is_file() and (u_dir / f).stat().st_size > 0
                                         for f in ["dup.seq.DupRate.xls", "ja.junction.xls",
                                                   "js.junctionSaturation_plot.r",
                                                   "id.inner_distance.txt"])
            elif cfg == "P-par":
                res = run_ppar(bam, bed, work_root, args.timeout)
                row.update(res)
                row["argv"] = [t["argv"] for t in res["per_tool"]]
                row["exit"] = 0 if all(t["exit"] == 0 for t in res["per_tool"]) else 1
                row["presence_ok"] = all(((work_root / "P-par") / f).is_file()
                                         for f in ["dup.seq.DupRate.xls", "ja.junction.xls"])
            else:
                threads = 1 if cfg == "R-1" else 8
                argv, m, code, to, load = run_rustqc(bam, gtf, work_root, threads, args.timeout)
                row["argv"] = argv
                row.update(m)
                row["exit"] = code
                row["timed_out"] = to
                row["loadavg_before_run"] = load
                # presence: R arm files
                r_sub = work_root / f"R-{threads}" / "rseqc"
                row["presence_ok"] = (r_sub / "bam_stat").exists()
            fname = f"{args.workload}.{cfg}.b{block}.json"
            (outdir / fname).write_text(json.dumps(row, indent=2) + "\n")
            print(f"wrote {fname}: exit={row.get('exit')} timed_out={row.get('timed_out')} presence={row.get('presence_ok')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
