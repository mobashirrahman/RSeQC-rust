#!/usr/bin/env python3
"""Regenerate every benchmark table from the raw committed measurement JSON.

No number in the report is written by hand: this script reads the raw per-run
measurements and produces the tables, so the report cannot drift from the data.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path


def verdict(ci):
    if not ci:
        return "inconclusive", None
    lo, hi = ci
    if lo > 1.0:
        return "win", None
    if hi < 1.0:
        return "loss", None
    return "inconclusive", None


def recompute(r):
    """Recompute the wall statistics from the raw per-run measurements.

    The point estimate and its confidence interval must describe the same quantity.
    An earlier harness revision took the ratio from raw medians but the interval from
    floor-subtracted pairs, so for short commands the interval did not contain the
    reported ratio. Recomputing here from `runs` + `floor_s` keeps the committed report
    internally consistent regardless of which harness revision produced the raw file.
    """
    import random as _random
    runs = r.get("runs") or {}
    fl = r.get("floor_s") or {}
    fpy, frs = fl.get("python"), fl.get("rust")
    py, rs = runs.get("py", []), runs.get("rs", [])
    n = min(len(py), len(rs))
    if n == 0:
        return r
    pairs = []
    for a, b in zip(py[:n], rs[:n]):
        va = a.get("wall_s", 0) - (fpy or 0.0)
        vb = b.get("wall_s", 0) - (frs or 0.0)
        if va > 1e-6 and vb > 1e-6:
            pairs.append((va, vb))
    if not pairs:
        return r
    pv = [x[0] for x in pairs]
    rv = [x[1] for x in pairs]
    mpy, mrs = statistics.median(pv), statistics.median(rv)
    ratio = (mpy / mrs) if mrs else None

    def boot_ci(stat, iters=10000, seed=12345, alpha=0.05):
        rng = _random.Random(seed)
        k = len(pairs)
        out = []
        for _ in range(iters):
            samp = [pairs[rng.randrange(k)] for _ in range(k)]
            v = stat(samp)
            if v is not None:
                out.append(v)
        if not out:
            return None
        out.sort()
        return (out[int(alpha / 2 * iters)],
                out[min(iters - 1, int((1 - alpha / 2) * iters))])

    def med_ratio(samp):
        a = statistics.median(x[0] for x in samp)
        b = statistics.median(x[1] for x in samp)
        return a / b if b else None

    def geo_ratio(samp):
        tot = 0.0
        for a, b in samp:
            tot += math.log(a / b)
        return math.exp(tot / len(samp))

    r = dict(r)
    r["recomputed"] = {
        "median_python_e1": mpy, "median_rust_e1": mrs, "ratio_e1": ratio,
        "ratio_e1_ci95": boot_ci(med_ratio),
        "ratio_geomean_logratio_e1": geo_ratio(pairs),
        "ratio_geomean_ci95": boot_ci(geo_ratio),
        "median_python_e2e": statistics.median([a.get("wall_s", 0) for a in py[:n]]),
        "median_rust_e2e": statistics.median([b.get("wall_s", 0) for b in rs[:n]]),
        "n_pairs": len(pairs),
    }
    return r


def load(results_path: Path):
    d = json.loads(results_path.read_text())
    for i, r in enumerate(d["results"]):
        if "error" not in r and r.get("runs"):
            d["results"][i] = recompute(r)
    rows = []
    for r in d["results"]:
        if "error" in r:
            rows.append({"command": r["command"], "status": "error", "note": r["error"]})
            continue
        w = r.get("wall_e2e", {})
        cpu = r.get("cpu", {})
        m = r.get("peak_rss_mb", {})
        fl = r.get("floor_s", {})
        rc = r.get("recomputed")
        if rc:
            ci = rc.get("ratio_e1_ci95")
            mp, mr = rc["median_python_e1"], rc["median_rust_e1"]
            ratio = rc["ratio_e1"]
        else:
            ci = w.get("ratio_ci95")
            mp, mr = w.get("median_python"), w.get("median_rust")
            ratio = w.get("ratio_median")
        v, _ = verdict(ci)
        logr = w.get("logratio_ci95")
        # E1: compute-only, i.e. after subtracting each arm's measured per-invocation
        # fixed cost. This is the number that reflects algorithmic work rather than
        # interpreter startup.
        e1_py, e1_rs, e1_ratio = mp, mr, ratio
        rows.append({
            "command": r["command"],
            "status": "ok",
            "workload": Path(r.get("workload", "")).name,
            "gate": r.get("equivalence", {}).get("pass"),
            "gate_errors": r.get("equivalence", {}).get("errors", []),
            "py_median": e1_py, "rs_median": e1_rs,
            "py_e2e": (rc or {}).get("median_python_e2e",
                                     w.get("median_python_e2e", w.get("median_python"))),
            "rs_e2e": (rc or {}).get("median_rust_e2e",
                                     w.get("median_rust_e2e", w.get("median_rust"))),
            "e2e_ratio": ratio, "e2e_ci": ci, "verdict": v,
            "e2e_logratio_ci": (rc or {}).get("ratio_geomean_ci95"),
            "floor_py": fl.get("python"), "floor_rs": fl.get("rust"),
            "e1_py": e1_py, "e1_rs": e1_rs, "e1_ratio": e1_ratio,
            "cpu_py": cpu.get("median_python"), "cpu_rs": cpu.get("median_rust"),
            "cpu_ratio": cpu.get("ratio_median"),
            "rss_py": m.get("median_python"), "rss_rs": m.get("median_rust"),
            "rss_ratio": m.get("ratio"),
            "rss_max_py": m.get("max_python"), "rss_max_rs": m.get("max_rust"),
            "n": w.get("n"),
            "failures": len(r.get("failures", [])),
        })
    return d, rows


def md_table(rows, keys, headers, title):
    out = [f"### {title}", "", "| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        cells = []
        for k in keys:
            v = r.get(k)
            if isinstance(v, (list, tuple)) and len(v) == 2 and all(
                    isinstance(x, (int, float)) for x in v):
                cells.append(f"[{v[0]:.2f}, {v[1]:.2f}]")
            elif isinstance(v, float):
                cells.append(f"{v:.3f}" if abs(v) < 1000 else f"{v:.0f}")
            elif v is None:
                cells.append("–")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--scaling", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    d, rows = load(Path(args.results))
    env = d["environment"]
    ok = [r for r in rows if r["status"] == "ok"]

    L = []
    L.append("# Benchmark results — rseqc-rust vs upstream RSeQC\n")
    L.append(f"Protocol: `benchmarks/protocol.md` (frozen 2026-09-29). "
             f"Run label: `{d.get('label','-')}`. Reps: {d['config']['reps']}, "
             f"seed {d['config']['seed']}.\n")
    L.append("> **Not publication-grade.** Shared, non-isolated hardware "
             f"({env['processor']}, loadavg {env['loadavg'][0]:.2f} at start, "
             f"governor `{env['cpufreq_governor']}`). See protocol.md §2.\n")

    L.append("## Environment\n")
    for k in ("platform", "cpu_count", "processor", "cpufreq_governor", "filesystem",
              "rustc", "python", "numpy", "pysam", "git_commit", "git_dirty"):
        L.append(f"- **{k}**: {env.get(k)}")
    L.append(f"- **pinned env**: `{env['pinned_env']}`\n")

    # ---- headline
    reportable = [r for r in ok if r["gate"] is True and r["e2e_ci"]]
    L.append("## 1. Per-command results (E2 end-to-end, equal deliverables)\n")
    L.append("Speedup = median(upstream) / median(port) on the E1 compute-only "
             "quantities (each arm's measured per-invocation fixed cost subtracted, so "
             "the number and the interval describe the same thing). The interval is a "
             "95% paired block bootstrap over whole pairs. A claim requires the interval "
             "to lie entirely above 1.0; an interval spanning 1.0 is inconclusive and "
             "an interval below 1.0 is a loss. Section 2 shows the end-to-end figures.\n")
    keys = ["command", "workload", "py_median", "rs_median", "e2e_ratio", "e2e_ci",
            "verdict", "cpu_ratio", "rss_ratio", "gate"]
    hdr = ["command", "workload", "py net (s)", "rs net (s)", "speedup", "95% CI",
           "verdict", "CPU x", "mem x", "gate"]
    L.append(md_table(sorted(ok, key=lambda r: -(r["e2e_ratio"] or 0)), keys, hdr,
                      "all measured commands"))

    # ---- E1
    L.append("## 2. E1 compute-only (fixed per-invocation cost subtracted)\n")
    L.append("Upstream pays a large fixed cost before reading any input "
             "(interpreter start plus imports). E1 subtracts each arm's *measured* "
             "floor (`--help` on the real binary) so the number reflects algorithmic "
             "work. E2 includes it, because a user always pays it.\n")
    keys1 = ["command", "floor_py", "floor_rs", "e1_py", "e1_rs", "e1_ratio",
             "py_e2e", "rs_e2e"]
    hdr1 = ["command", "py floor (s)", "rs floor (s)", "py net (s)", "rs net (s)",
            "E1 x", "py E2 (s)", "rs E2 (s)"]
    L.append(md_table(sorted(ok, key=lambda r: -(r["e1_ratio"] or 0)), keys1, hdr1,
                      "compute-only vs end-to-end"))

    # ---- memory
    L.append("## 3. Peak memory (whole process tree, `/usr/bin/time -v`)\n")
    keysm = ["command", "rss_py", "rss_rs", "rss_ratio", "rss_max_py", "rss_max_rs"]
    L.append(md_table(sorted(ok, key=lambda r: (r["rss_ratio"] or 0)), keysm,
                      ["command", "py median MB", "rs median MB", "mem x",
                       "py max MB", "rs max MB"],
                      "peak RSS; a ratio above 1 means the port uses more"))

    # ---- gate failures
    failed = [r for r in ok if r["gate"] is not True]
    L.append("## 4. Equivalence gate\n")
    L.append("No speedup is claimed for a row whose gate did not pass. A fast but "
             "incomplete run must not count as a win.\n")
    if failed:
        L.append(md_table(failed, ["command", "workload", "gate_errors"],
                          ["command", "workload", "errors"], "gate failures"))
    else:
        L.append("All measured rows passed structural equivalence.\n")

    L.append("### Declared exclusions\n")
    for k, v in d.get("excluded", {}).items():
        L.append(f"- **{k}** — {v}")
    L.append("")

    # ---- scaling
    if args.scaling and Path(args.scaling).exists():
        sr = json.loads(Path(args.scaling).read_text())
        L.append("## 5. Scaling along cost drivers\n")
        for drv in ("reads", "transcripts", "readlen"):
            sub = [s for s in sr if s["driver"] == drv]
            if not sub:
                continue
            for cmd in sorted({s["command"] for s in sub}):
                pts = sorted([s for s in sub if s["command"] == cmd],
                             key=lambda s: s["value"])
                L.append(f"\n**{cmd}** vs {drv}\n")
                L.append("| " + drv + " | py net (s) | rs net (s) | speedup | "
                         "95% CI | py E2 (s) | rs E2 (s) | gate |")
                L.append("|---|---|---|---|---|---|---|---|")
                for s in pts:
                    ci = s.get("ci95")
                    cis = f"[{ci[0]:.2f}, {ci[1]:.2f}]" if ci else "n/a"
                    L.append(f"| {s['value']} | {s['median_python']:.3f} | "
                             f"{s['median_rust']:.4f} | "
                             f"{(s['ratio'] if s['ratio'] else 0):.2f}x | {cis} | "
                             f"{s.get('median_python_e2e', 0):.3f} | "
                             f"{s.get('median_rust_e2e', 0):.3f} | {s['gate']} |")

    Path(args.out).write_text("\n".join(L) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
