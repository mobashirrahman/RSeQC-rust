#!/usr/bin/env python3
"""Card R4: three-arm equivalence (U vs P, U vs R) for the eight RSeQC tools.

Reads workloads.json and output-map.json, runs the three arms on one
workload, and for each tool compares P to U and R to U.

Reuses verification/comparators.py; no new numeric comparator.

Verdicts per file, in order:
  byte-identical, numerically-identical, within-tolerance, differs,
  missing, not-comparable.

Usage:
    python3 equivalence.py --workload W-rat-3c --outdir /tmp/eq [--tolerance 0] [--json out.json]
    python3 equivalence.py --workload W-rat-3c --check-presence --u-dir <d> --p-dir <d> --r-dir <d>
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "verification"))
from comparators import compare_numeric_table, compare_text, is_nonfinite  # noqa: E402

HERE = Path(__file__).resolve().parent
WORKLOADS = HERE / "workloads.json"
OUTPUT_MAP = HERE / "output-map.json"
UPSTREAM_SCRIPTS = REPO / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHONPATH = str(REPO / "oracle" / "upstream-src" / "src")
ORACLE_PYTHON = REPO / "oracle" / "venv" / "bin" / "python3"
RELEASE = REPO / "target" / "release"
RUSTQC = HERE / "tools" / "rustqc"

TOOLS = ["bam_stat", "infer_experiment", "read_duplication", "read_distribution",
         "junction_annotation", "junction_saturation", "inner_distance", "tin"]

# junction_saturation shuffles without a seed upstream (see
# compatibility/divergences.yaml); two upstream runs on identical input
# differ, so byte/numeric equality is the wrong gate.
NOT_COMPARABLE = {
    "junction_saturation": "upstream shuffles alignments without an RNG seed; two upstream runs differ",
}

PINNED_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0",
}


class _Proc:
    def __init__(self, returncode, stdout, stderr, timed_out=False):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.timed_out = timed_out


def run(argv, cwd, env_extra=None, timeout=3600):
    """Run argv; a timeout is a recorded result, not a crash."""
    env = dict(os.environ)
    env.update(PINNED_ENV)
    if env_extra:
        env.update(env_extra)
    try:
        proc = subprocess.run(argv, cwd=str(cwd), env=env, capture_output=True,
                              text=True, timeout=timeout)
        return _Proc(proc.returncode, proc.stdout, proc.stderr, False)
    except subprocess.TimeoutExpired as exc:
        out = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        err = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        err += f"\n[equivalence timeout after {timeout:g}s]"
        return _Proc(None, out, err, True)


def verdict_file(a: Path, b: Path, tolerance: float = 0.0):
    """Verdict for two files. Returns (verdict, detail)."""
    a_exists = a.is_file()
    b_exists = b.is_file()
    if not a_exists or not b_exists:
        return ("missing",
                f"present: first={a_exists} second={b_exists} ({a.name})")
    ab = a.read_bytes()
    bb = b.read_bytes()
    if len(ab) == 0 and len(bb) == 0:
        return ("differs", "both files are empty (0 bytes); empty is not agreement")
    if len(ab) == 0 or len(bb) == 0:
        return ("differs", f"one file empty: {len(ab)} vs {len(bb)} bytes")
    # Nonfinite is never agreement, even when both sides are byte-identical.
    for token in ab.decode("utf8", "replace").replace("\t", " ").split():
        if is_nonfinite(token):
            return ("differs", f"nonfinite value {token!r} in first output")
    for token in bb.decode("utf8", "replace").replace("\t", " ").split():
        if is_nonfinite(token):
            return ("differs", f"nonfinite value {token!r} in second output")
    if ab == bb:
        return ("byte-identical", f"{len(ab)} bytes identical")
    ok, _ = compare_text(a, b, rtol=0.0, atol=0.0)
    if ok:
        return ("numerically-identical",
                "same structure, numbers equal after parsing; formatting/whitespace only")
    if tolerance and tolerance > 0:
        ok_tol, _ = compare_text(a, b, rtol=tolerance, atol=tolerance)
        if ok_tol:
            la, lr = largest_diff(a, b)
            return ("within-tolerance",
                    f"largest absolute difference {la:g}, largest relative {lr:g}, tolerance {tolerance:g}")
    # differs: first 20 differing lines + largest difference
    first, la, lr = diff_detail(a, b, n=20)
    detail = f"largest absolute {la:g}, largest relative {lr:g}; first differences:\n" + "\n".join(first[:20])
    return ("differs", detail)


def verdict_text(a_text: str, b_text: str, tolerance: float = 0.0, name="stream"):
    """Verdict for two text blobs (e.g. captured stdout vs R file content)."""
    with tempfile.TemporaryDirectory() as td:
        pa = Path(td) / "a.txt"
        pb = Path(td) / "b.txt"
        pa.write_text(a_text)
        pb.write_text(b_text)
        # Empty-stream guard: two empty streams are not agreement.
        if not a_text.strip() and not b_text.strip():
            return ("missing", f"{name}: both streams empty; nothing was compared")
        if not a_text.strip() or not b_text.strip():
            return ("missing", f"{name}: one stream empty ({len(a_text)} vs {len(b_text)} chars)")
        return verdict_file(pa, pb, tolerance)


def largest_diff(a: Path, b: Path):
    la_lines = a.read_text(errors="replace").splitlines()
    lb_lines = b.read_text(errors="replace").splitlines()
    best_abs = 0.0
    best_rel = 0.0
    for x, y in zip(la_lines, lb_lines):
        for u, v in zip(x.split(), y.split()):
            try:
                fu, fv = float(u), float(v)
            except ValueError:
                continue
            if fu != fu or fv != fv:
                continue
            d = abs(fu - fv)
            best_abs = max(best_abs, d)
            denom = max(abs(fu), abs(fv), 1e-9)
            best_rel = max(best_rel, d / denom)
    return best_abs, best_rel


def diff_detail(a: Path, b: Path, n=20):
    la_lines = a.read_text(errors="replace").splitlines()
    lb_lines = b.read_text(errors="replace").splitlines()
    first = []
    if len(la_lines) != len(lb_lines):
        first.append(f"line count {len(la_lines)} vs {len(lb_lines)}")
    for i, (x, y) in enumerate(zip(la_lines, lb_lines)):
        if x != y:
            first.append(f"line {i+1}: {x[:100]!r} vs {y[:100]!r}")
            if len(first) >= n:
                break
    la, lr = largest_diff(a, b)
    return first, la, lr


def normalize_r_script(text: str, work_dir: str) -> str:
    # R plotting scripts embed absolute output paths; normalise before compare.
    return text.replace(work_dir, "<OUT>")


def run_arms(workload_id: str, scratch: Path, timeout_per_tool: int = 3600):
    data = json.loads(WORKLOADS.read_text())
    w = data["workloads"][workload_id]
    bam = str(REPO / w["bam"])
    bed = str(REPO / w["bed"])
    gtf_subset = str(REPO / w.get("gtf_subset", w["gtf"]))
    sample = Path(bam).stem  # SRR1177982 / ERR10229623

    u_dir = scratch / "U"
    p_dir = scratch / "P"
    r_dir = scratch / "R"
    for d in (u_dir, p_dir, r_dir):
        d.mkdir(parents=True, exist_ok=True)

    results = {"workload": workload_id, "arms": {}, "timeouts": {}}

    def upy(script, args, out_prefix=None):
        argv = [str(ORACLE_PYTHON), str(UPSTREAM_SCRIPTS / script)] + args
        proc = run(argv, cwd=u_dir, env_extra={"PYTHONPATH": ORACLE_PYTHONPATH},
                   timeout=timeout_per_tool)
        if proc.timed_out:
            results["timeouts"][f"U:{script}"] = timeout_per_tool
        return proc

    def rport(binary, args):
        argv = [str(RELEASE / binary)] + args
        proc = run(argv, cwd=p_dir, timeout=timeout_per_tool)
        if proc.timed_out:
            results["timeouts"][f"P:{binary}"] = timeout_per_tool
        return proc

    # --- bam_stat (stdout) ---
    u = upy("bam_stat.py", ["-i", bam])
    (u_dir / "bam_stat.stdout").write_text(u.stdout)
    (u_dir / "bam_stat.stderr").write_text(u.stderr)
    (u_dir / "bam_stat.exit").write_text(str(u.returncode))
    p = rport("bam_stat", ["-i", bam])
    (p_dir / "bam_stat.stdout").write_text(p.stdout)
    (p_dir / "bam_stat.stderr").write_text(p.stderr)
    (p_dir / "bam_stat.exit").write_text(str(p.returncode))

    # --- infer_experiment (stdout) ---
    u = upy("infer_experiment.py", ["-i", bam, "-r", bed])
    (u_dir / "infer_experiment.stdout").write_text(u.stdout)
    (u_dir / "infer_experiment.stderr").write_text(u.stderr)
    (u_dir / "infer_experiment.exit").write_text(str(u.returncode))
    p = rport("infer_experiment", ["-i", bam, "-r", bed])
    (p_dir / "infer_experiment.stdout").write_text(p.stdout)
    (p_dir / "infer_experiment.stderr").write_text(p.stderr)
    (p_dir / "infer_experiment.exit").write_text(str(p.returncode))

    # --- read_duplication (files) ---
    u = upy("read_duplication.py", ["-i", bam, "-o", str(u_dir / "dup"), "--skip-plot"])
    (u_dir / "dup.stdout").write_text(u.stdout)
    (u_dir / "dup.stderr").write_text(u.stderr)
    p = rport("read_duplication", ["-i", bam, "-o", str(p_dir / "dup"), "--skip-plot"])
    (p_dir / "dup.stdout").write_text(p.stdout)
    (p_dir / "dup.stderr").write_text(p.stderr)

    # --- read_distribution (stdout) ---
    u = upy("read_distribution.py", ["-i", bam, "-r", bed])
    (u_dir / "read_distribution.stdout").write_text(u.stdout)
    (u_dir / "read_distribution.stderr").write_text(u.stderr)
    (u_dir / "read_distribution.exit").write_text(str(u.returncode))
    p = rport("read_distribution", ["-i", bam, "-r", bed])
    (p_dir / "read_distribution.stdout").write_text(p.stdout)
    (p_dir / "read_distribution.stderr").write_text(p.stderr)
    (p_dir / "read_distribution.exit").write_text(str(p.returncode))

    # --- junction_annotation (files + stdout total line) ---
    u = upy("junction_annotation.py", ["-i", bam, "-r", bed, "-o", str(u_dir / "ja"), "--skip-plot"])
    (u_dir / "ja.stdout").write_text(u.stdout)
    (u_dir / "ja.stderr").write_text(u.stderr)
    p = rport("junction_annotation", ["-i", bam, "-r", bed, "-o", str(p_dir / "ja"), "--skip-plot"])
    (p_dir / "ja.stdout").write_text(p.stdout)
    (p_dir / "ja.stderr").write_text(p.stderr)

    # --- junction_saturation (files + stderr curve) ---
    u = upy("junction_saturation.py", ["-i", bam, "-r", bed, "-o", str(u_dir / "js"), "--skip-plot"])
    (u_dir / "js.stdout").write_text(u.stdout)
    (u_dir / "js.stderr").write_text(u.stderr)
    p = rport("junction_saturation", ["-i", bam, "-r", bed, "-o", str(p_dir / "js"), "--skip-plot"])
    (p_dir / "js.stdout").write_text(p.stdout)
    (p_dir / "js.stderr").write_text(p.stderr)

    # --- inner_distance (files) ---
    u = upy("inner_distance.py", ["-i", bam, "-r", bed, "-o", str(u_dir / "id"), "--skip-plot"])
    (u_dir / "id.stdout").write_text(u.stdout)
    (u_dir / "id.stderr").write_text(u.stderr)
    p = rport("inner_distance", ["-i", bam, "-r", bed, "-o", str(p_dir / "id"), "--skip-plot"])
    (p_dir / "id.stdout").write_text(p.stdout)
    (p_dir / "id.stderr").write_text(p.stderr)

    # --- tin (files) ---
    (u_dir / "tin").mkdir(exist_ok=True)
    (p_dir / "tin").mkdir(exist_ok=True)
    u = upy("tin.py", ["-i", bam, "-r", bed, "-n", "100", "-o", str(u_dir / "tin")])
    (u_dir / "tin.stdout").write_text(u.stdout)
    (u_dir / "tin.stderr").write_text(u.stderr)
    p = rport("tin", ["-i", bam, "-r", bed, "-n", "100", "-o", str(p_dir / "tin")])
    (p_dir / "tin.stdout").write_text(p.stdout)
    (p_dir / "tin.stderr").write_text(p.stderr)

    # --- RustQC single pass ---
    r_argv = [str(RUSTQC), "rna", bam, "--gtf", gtf_subset, "--paired",
              "--outdir", str(r_dir), "--threads", "1", "--skip-dup-check"]
    (r_dir / "argv.txt").write_text(" ".join(r_argv))
    r = run(r_argv, cwd=r_dir, timeout=timeout_per_tool)
    (r_dir / "stdout.txt").write_text(r.stdout)
    (r_dir / "stderr.txt").write_text(r.stderr)
    (r_dir / "exit.txt").write_text(str(r.returncode))
    if r.timed_out:
        results["timeouts"]["R:rustqc"] = timeout_per_tool
    results["rustqc_argv"] = r_argv
    results["rustqc_exit"] = r.returncode
    return {"u_dir": u_dir, "p_dir": p_dir, "r_dir": r_dir,
            "sample": sample, "results": results}


def compare_arms(u_dir: Path, p_dir: Path, r_dir: Path, sample: str,
                 tolerance: float = 0.0, timeouts: dict | None = None):
    """Compare already-run arms. Returns {tool: {'P': (verdict, detail, file), 'R': ...}}.

    A tool with a recorded timeout keeps its partial outputs on disk, but a
    partial file is not a result: any verdict involving it is `missing` with
    the timeout as the reason, not `differs`.
    """
    table = {}
    timeouts = timeouts or {}

    def timed_out(*names):
        return any(k in timeouts for k in names)

    def cmp_pair(a: Path, b: Path):
        return verdict_file(a, b, tolerance)

    # bam_stat: U stdout vs P stdout; U stdout vs R file
    v_p, d_p = verdict_text((u_dir / "bam_stat.stdout").read_text() if (u_dir / "bam_stat.stdout").exists() else "",
                            (p_dir / "bam_stat.stdout").read_text() if (p_dir / "bam_stat.stdout").exists() else "",
                            tolerance, "bam_stat stdout")
    r_files = list((r_dir / "rseqc" / "bam_stat").glob("*.bam_stat.txt")) if (r_dir / "rseqc" / "bam_stat").exists() else []
    if r_files:
        v_r, d_r = verdict_text((u_dir / "bam_stat.stdout").read_text(),
                                r_files[0].read_text(), tolerance, "bam_stat U-stdout vs R-file")
        r_detail = f"{r_files[0].name}: {d_r}"
    else:
        v_r, r_detail = "missing", "R bam_stat file absent"
    table["bam_stat"] = {"P": (v_p, d_p, "stdout"), "R": (v_r, r_detail, r_files[0].name if r_files else "-")}

    # infer_experiment
    v_p, d_p = verdict_text((u_dir / "infer_experiment.stdout").read_text() if (u_dir / "infer_experiment.stdout").exists() else "",
                            (p_dir / "infer_experiment.stdout").read_text() if (p_dir / "infer_experiment.stdout").exists() else "",
                            tolerance, "infer_experiment stdout")
    r_files = list((r_dir / "rseqc" / "infer_experiment").glob("*.infer_experiment.txt")) if (r_dir / "rseqc" / "infer_experiment").exists() else []
    if r_files:
        v_r, d_r = verdict_text((u_dir / "infer_experiment.stdout").read_text(), r_files[0].read_text(),
                                tolerance, "infer_experiment U-stdout vs R-file")
        r_detail = f"{r_files[0].name}: {d_r}"
    else:
        v_r, r_detail = "missing", "R infer_experiment file absent"
    table["infer_experiment"] = {"P": (v_p, d_p, "stdout"), "R": (v_r, r_detail, r_files[0].name if r_files else "-")}

    # read_duplication: two xls files
    for which in ["dup.seq.DupRate.xls", "dup.pos.DupRate.xls"]:
        a = u_dir / which.replace("dup.", "dup.")
        # U files are u_dir/dup.seq... ; P files p_dir/dup.seq...
        a = u_dir / which
        b = p_dir / which
        # byte/numeric compare with path normalisation not needed (no paths inside)
        v, d = cmp_pair(a, b)
        table.setdefault("read_duplication", {}).setdefault("P_files", {})[which] = (v, d)
    # collapse P to worst verdict across its files
    p_verts = [v for v, _ in table["read_duplication"]["P_files"].values()]
    table["read_duplication"]["P"] = (worst(p_verts), "; ".join(f"{k}={v}" for k, (v, _) in table["read_duplication"]["P_files"].items()), "2 xls")
    # R vs U
    r_seq = list((r_dir / "rseqc" / "read_duplication").glob("*.seq.DupRate.xls")) if (r_dir / "rseqc" / "read_duplication").exists() else []
    r_pos = list((r_dir / "rseqc" / "read_duplication").glob("*.pos.DupRate.xls")) if (r_dir / "rseqc" / "read_duplication").exists() else []
    r_verts = []
    r_details = []
    for which, rf in [("seq", r_seq[0] if r_seq else None), ("pos", r_pos[0] if r_pos else None)]:
        uf = u_dir / f"dup.{which}.DupRate.xls"
        if rf is None:
            r_verts.append("missing")
            r_details.append(f"{which}: R file absent")
        else:
            v, d = cmp_pair(uf, rf)
            r_verts.append(v)
            r_details.append(f"{rf.name}: {v}: {d[:200]}")
    table["read_duplication"]["R"] = (worst(r_verts), " | ".join(r_details), "2 xls")
    # extra R layout files
    table["read_duplication"]["R_extra"] = sorted(p.name for p in (r_dir / "rseqc" / "read_duplication").glob("*")) if (r_dir / "rseqc" / "read_duplication").exists() else []

    # read_distribution
    v_p, d_p = verdict_text((u_dir / "read_distribution.stdout").read_text() if (u_dir / "read_distribution.stdout").exists() else "",
                            (p_dir / "read_distribution.stdout").read_text() if (p_dir / "read_distribution.stdout").exists() else "",
                            tolerance, "read_distribution stdout")
    r_files = list((r_dir / "rseqc" / "read_distribution").glob("*.read_distribution.txt")) if (r_dir / "rseqc" / "read_distribution").exists() else []
    if r_files:
        v_r, d_r = verdict_text((u_dir / "read_distribution.stdout").read_text(), r_files[0].read_text(),
                                tolerance, "read_distribution U-stdout vs R-file")
        r_detail = f"{r_files[0].name}: {d_r}"
    else:
        v_r, r_detail = "missing", "R read_distribution file absent"
    table["read_distribution"] = {"P": (v_p, d_p, "stdout"), "R": (v_r, r_detail, r_files[0].name if r_files else "-")}

    # junction_annotation: .xls primary + .bed + Interact.bed
    for suffix in ["ja.junction.xls", "ja.junction.bed", "ja.junction.Interact.bed"]:
        a = u_dir / suffix
        b = p_dir / suffix
        v, d = cmp_pair(a, b)
        table.setdefault("junction_annotation", {}).setdefault("P_files", {})[suffix] = (v, d)
    p_verts = [v for v, _ in table["junction_annotation"]["P_files"].values()]
    table["junction_annotation"]["P"] = (worst(p_verts), "; ".join(f"{k}={v}" for k, (v, _) in table["junction_annotation"]["P_files"].items()), "xls+bed+interact")
    # R vs U on .xls/.bed/Interact
    r_map = {"ja.junction.xls": "*.junction.xls", "ja.junction.bed": "*.junction.bed",
             "ja.junction.Interact.bed": "*.junction.Interact.bed"}
    r_verts, r_details = [], []
    for usuffix, pat in r_map.items():
        uf = u_dir / usuffix
        cands = list((r_dir / "rseqc" / "junction_annotation").glob(pat)) if (r_dir / "rseqc" / "junction_annotation").exists() else []
        # exclude .r (matched by *.junction.xls? no)
        cands = [c for c in cands if c.suffix in (".xls", ".bed")]
        if not cands:
            r_verts.append("missing")
            r_details.append(f"{usuffix}: R file absent")
        else:
            v, d = cmp_pair(uf, cands[0])
            r_verts.append(v)
            r_details.append(f"{cands[0].name}: {v}: {d[:200]}")
    table["junction_annotation"]["R"] = (worst(r_verts), " | ".join(r_details), "xls+bed+interact")
    table["junction_annotation"]["R_extra"] = sorted(p.name for p in (r_dir / "rseqc" / "junction_annotation").glob("*")) if (r_dir / "rseqc" / "junction_annotation").exists() else []

    # junction_saturation: not-comparable by design
    reason = NOT_COMPARABLE["junction_saturation"]
    table["junction_saturation"] = {"P": ("not-comparable", reason, "js.junctionSaturation_plot.r + stderr curve"),
                                    "R": ("not-comparable", reason, "R summary + plot")}

    # inner_distance: .txt + _freq.txt
    for suffix in ["id.inner_distance.txt", "id.inner_distance_freq.txt"]:
        a = u_dir / suffix
        b = p_dir / suffix
        v, d = cmp_pair(a, b)
        table.setdefault("inner_distance", {}).setdefault("P_files", {})[suffix] = (v, d)
    p_verts = [v for v, _ in table["inner_distance"]["P_files"].values()]
    table["inner_distance"]["P"] = (worst(p_verts), "; ".join(f"{k}={v}" for k, (v, _) in table["inner_distance"]["P_files"].items()), "txt+freq")
    r_verts, r_details = [], []
    for usuffix, pat in [("id.inner_distance.txt", "*.inner_distance.txt"),
                         ("id.inner_distance_freq.txt", "*.inner_distance_freq.txt")]:
        uf = u_dir / usuffix
        cands = list((r_dir / "rseqc" / "inner_distance").glob(pat)) if (r_dir / "rseqc" / "inner_distance").exists() else []
        cands = [c for c in cands if "_freq" in c.name or "freq" not in pat or True]
        # disambiguate: for .txt pattern both files match; pick exact kind
        if usuffix.endswith("_freq.txt"):
            cands = [c for c in cands if "freq" in c.name]
        else:
            cands = [c for c in cands if "freq" not in c.name and "mean" not in c.name and "summary" not in c.name]
        if not cands:
            r_verts.append("missing")
            r_details.append(f"{usuffix}: R file absent")
        else:
            v, d = cmp_pair(uf, cands[0])
            # inner_distance.txt contains read names in column 1; sample prefix
            # differs by arm only if outdirs differ -- normalise sample names?
            r_verts.append(v)
            r_details.append(f"{cands[0].name}: {v}: {d[:200]}")
    table["inner_distance"]["R"] = (worst(r_verts), " | ".join(r_details), "txt+freq")
    table["inner_distance"]["R_extra"] = sorted(p.name for p in (r_dir / "rseqc" / "inner_distance").glob("*")) if (r_dir / "rseqc" / "inner_distance").exists() else []

    # tin: <sample>.tin.xls + <sample>.summary.txt
    base = sample + ".tin.xls"
    sbase = sample + ".summary.txt"
    if timed_out("U:tin.py"):
        msg = "U tin timed out; U outputs partial, no verdict"
        table["tin"] = {"P_files": {base: ("missing", msg), sbase: ("missing", msg)}}
        table["tin"]["P"] = ("missing", msg, "tin.xls+summary")
        table["tin"]["R"] = ("missing", msg + " (R files listed, not judged)", "tin.xls+summary")
        table["tin"]["R_extra"] = sorted(p.name for p in (r_dir / "rseqc" / "tin").glob("*")) if (r_dir / "rseqc" / "tin").exists() else []
        return table
    a = u_dir / "tin" / base
    b = p_dir / "tin" / base
    v1, d1 = cmp_pair(a, b)
    a2 = u_dir / "tin" / sbase
    b2 = p_dir / "tin" / sbase
    v2, d2 = cmp_pair(a2, b2)
    table["tin"] = {"P_files": {base: (v1, d1), sbase: (v2, d2)}}
    table["tin"]["P"] = (worst([v1, v2]), f"{base}={v1}; {sbase}={v2}: {(d1 if v1 not in ('byte-identical','numerically-identical') else d2)[:300]}", "tin.xls+summary")
    r_tin = list((r_dir / "rseqc" / "tin").glob("*.tin.xls")) if (r_dir / "rseqc" / "tin").exists() else []
    r_sum = list((r_dir / "rseqc" / "tin").glob("*.summary.txt")) if (r_dir / "rseqc" / "tin").exists() else []
    r_verts, r_details = [], []
    for uf, cands, label in [(a, r_tin, "tin.xls"), (a2, r_sum, "summary")]:
        if not cands:
            r_verts.append("missing")
            r_details.append(f"{label}: R file absent")
        else:
            v, d = cmp_pair(uf, cands[0])
            r_verts.append(v)
            r_details.append(f"{cands[0].name}: {v}: {d[:200]}")
    table["tin"]["R"] = (worst(r_verts), " | ".join(r_details), "tin.xls+summary")

    return table


ORDER = ["byte-identical", "numerically-identical", "within-tolerance",
         "not-comparable", "differs", "missing"]


def worst(verdicts):
    for v in reversed(ORDER):
        if v in verdicts:
            return v
    return verdicts[0] if verdicts else "missing"


def check_presence(u_dir: Path, p_dir: Path, r_dir: Path, sample: str) -> bool:
    """File-presence check used by timing.py: every expected artifact exists non-empty."""
    checks = [
        (u_dir / "bam_stat.stdout"), (p_dir / "bam_stat.stdout"),
        (u_dir / "infer_experiment.stdout"), (p_dir / "infer_experiment.stdout"),
        (u_dir / "dup.seq.DupRate.xls"), (p_dir / "dup.seq.DupRate.xls"),
        (u_dir / "dup.pos.DupRate.xls"), (p_dir / "dup.pos.DupRate.xls"),
        (u_dir / "read_distribution.stdout"), (p_dir / "read_distribution.stdout"),
        (u_dir / "ja.junction.xls"), (p_dir / "ja.junction.xls"),
        (u_dir / "js.junctionSaturation_plot.r"), (p_dir / "js.junctionSaturation_plot.r"),
        (u_dir / "id.inner_distance.txt"), (p_dir / "id.inner_distance.txt"),
        (u_dir / "id.inner_distance_freq.txt"), (p_dir / "id.inner_distance_freq.txt"),
        (u_dir / "tin" / f"{sample}.tin.xls"), (p_dir / "tin" / f"{sample}.tin.xls"),
    ]
    for p in checks:
        if not p.is_file() or p.stat().st_size == 0:
            return False
    # R arm: one file per tool
    for sub, pat in [("bam_stat", "*.bam_stat.txt"), ("infer_experiment", "*.infer_experiment.txt"),
                     ("read_duplication", "*.seq.DupRate.xls"), ("read_distribution", "*.read_distribution.txt"),
                     ("junction_annotation", "*.junction.xls"), ("junction_saturation", "*.junctionSaturation_plot.r"),
                     ("inner_distance", "*.inner_distance.txt"), ("tin", "*.tin.xls")]:
        cands = list((r_dir / "rseqc" / sub).glob(pat)) if (r_dir / "rseqc" / sub).exists() else []
        if not cands or not any(c.stat().st_size > 0 for c in cands):
            return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workload", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--tolerance", type=float, default=0.0)
    ap.add_argument("--json", default=None)
    ap.add_argument("--check-presence", action="store_true")
    ap.add_argument("--u-dir", default=None)
    ap.add_argument("--p-dir", default=None)
    ap.add_argument("--r-dir", default=None)
    args = ap.parse_args()

    if args.check_presence:
        data = json.loads(WORKLOADS.read_text())
        sample = Path(data["workloads"][args.workload]["bam"]).stem
        ok = check_presence(Path(args.u_dir), Path(args.p_dir), Path(args.r_dir), sample)
        print("presence PASS" if ok else "presence FAIL")
        return 0 if ok else 1

    scratch = Path(args.outdir) if args.outdir else Path(tempfile.mkdtemp(prefix="eq-"))
    scratch.mkdir(parents=True, exist_ok=True)
    arms = run_arms(args.workload, scratch)
    table = compare_arms(arms["u_dir"], arms["p_dir"], arms["r_dir"], arms["sample"],
                         args.tolerance, arms["results"].get("timeouts", {}))

    print(f"workload {args.workload} (tolerance={args.tolerance})")
    print(f"{'tool':<20} {'P vs U':<22} {'R vs U':<22} file(s)")
    for t in TOOLS:
        row = table.get(t, {})
        pv = row.get("P", ("-", "", ""))[0]
        rv = row.get("R", ("-", "", ""))[0]
        fn = row.get("P", ("", "", ""))[2] if len(row.get("P", ())) > 2 else ""
        print(f"{t:<20} {pv:<22} {rv:<22} {fn}")
        for arm in ("P", "R"):
            if arm in row and isinstance(row[arm], tuple):
                print(f"  {arm}: {row[arm][1][:500]}")
    # R extras
    print("R-only layout files (listed, not judged):")
    for t in TOOLS:
        extra = table.get(t, {}).get("R_extra")
        if extra:
            print(f"  {t}: {', '.join(extra[:10])}" + (" ..." if len(extra) > 10 else ""))

    serial = {t: {"P": {"verdict": v, "detail": d}, "R": {"verdict": v2, "detail": d2}}
              for t, ((v, d, _), (v2, d2, _2)) in
              ((t, (table[t]["P"], table[t]["R"])) for t in TOOLS if "P" in table.get(t, {}) and "R" in table.get(t, {}))}
    out = {"workload": args.workload, "tolerance": args.tolerance,
           "u_dir": str(arms["u_dir"]), "p_dir": str(arms["p_dir"]), "r_dir": str(arms["r_dir"]),
           "timeouts": arms["results"].get("timeouts", {}),
           "rustqc_exit": arms["results"].get("rustqc_exit"),
           "verdicts": serial}
    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
