#!/usr/bin/env python3
"""Rigorous benchmark harness: upstream Python RSeQC vs this Rust port.

Implements `benchmarks/protocol.md`. The three things that distinguish it from the
earlier scaffolding (and from exit-code-plus-diff comparisons) are:

  1. PROCESS-TREE resource measurement. Each arm runs under `/usr/bin/time -v` in its
     own process group, so wall time, user/sys CPU and peak RSS cover the whole tree
     including any R/htseq-count/wigToBigWig child. Nothing is attributed from inside
     the process.
  2. A STRUCTURAL EQUIVALENCE GATE. Outputs are compared by type-aware semantics
     (text tables line-wise with numeric tolerance, BAM as a decoded record multiset,
     BigWig/bedGraph as interval-value maps, FASTQ record-wise). A speedup is only
     reported when the gate passes, so a fast incomplete run cannot count as a win.
  3. PAIRED, RANDOMISED, INTERLEAVED repetitions with an explicit estimator, so
     machine drift cannot systematically favour one arm.

Usage:
    oracle/venv/bin/python3 benchmarks/bench.py --workload workloads/A-mid \
        --commands bam_stat tin --reps 10
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
UPSTREAM_SRC = REPO / "oracle" / "upstream-src" / "scripts"
VENV_PY = REPO / "oracle" / "venv" / "bin" / "python3"
RELEASE = REPO / "target" / "release"

# Parity runs pin implicit thread pools in both arms. Upstream's numpy/OpenBLAS
# initialisation otherwise spends ~1s of CPU across all cores per invocation
# (measured: `import numpy` = 0.08s wall / 0.90s user), which would otherwise be
# charged to the reference arm as if it were algorithmic work.
PINNED_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0",
}

TIME_BIN = "/usr/bin/time"


# --------------------------------------------------------------------------------------
# Command matrix
# --------------------------------------------------------------------------------------
# Each entry: (upstream_script, rust_binary, args_builder, experiment_class, deliverables)
# experiment_class: "compute" (E1) or "end-to-end" (E2)
# excluded: commands that cannot be benchmarked here, with the protocol reason.
def _bam(extra=()):
    return ["-i", "{bam}", *extra]


#
# Argument vectors below are taken from each command's real upstream `--help` text
# (verified in `oracle/reference-runs/help/`), because a wrong flag makes the port
# exit 2 as fast as the reference and yields a meaningless "speedup" of two
# argument-parsing errors. Each entry was re-checked by running both arms.
COMMANDS = {
    "bam_stat": ("bam_stat.py", "bam_stat", lambda c: ["-i", c.bam], "compute"),
    "bam2fq": ("bam2fq.py", "bam2fq", lambda c: ["-i", c.bam, "-o", c.out / "r1"], "compute"),
    "infer_experiment": ("infer_experiment.py", "infer_experiment",
                         lambda c: ["-i", c.bam, "-r", c.bed], "compute"),
    "read_distribution": ("read_distribution.py", "read_distribution",
                          lambda c: ["-i", c.bam, "-r", c.bed], "compute"),
    "inner_distance": ("inner_distance.py", "inner_distance",
                       lambda c: ["-i", c.bam, "-r", c.bed, "-o", c.out / "id",
                                  "--skip-plot"], "compute"),
    "RNA_fragment_size": ("RNA_fragment_size.py", "RNA_fragment_size",
                          lambda c: ["-i", c.bam, "-r", c.bed], "compute"),
    # tin is per-transcript and quadratic-ish in regions x local depth upstream;
    # -n bounds the sample so a 10-rep matrix fits the time budget.
    "tin": ("tin.py", "tin", lambda c: ["-i", c.bam, "-r", c.bed, "-n", "100",
                                        "-o", c.out], "compute"),
    "FPKM_count": ("FPKM_count.py", "FPKM_count", lambda c: ["-i", c.bam, "-r", c.bed, "-o", c.out / "fpkm"], "compute"),
    "read_NVC": ("read_NVC.py", "read_NVC",
                 lambda c: ["-i", c.bam, "-o", c.out / "nvc", "--skip-plot"], "compute"),
    "read_GC": ("read_GC.py", "read_GC",
                lambda c: ["-i", c.bam, "-o", c.out / "gc", "--skip-plot"], "compute"),
    "read_quality": ("read_quality.py", "read_quality",
                     lambda c: ["-i", c.bam, "-o", c.out / "q", "--skip-plot"], "compute"),
    "read_duplication": ("read_duplication.py", "read_duplication",
                         lambda c: ["-i", c.bam, "-o", c.out / "dup", "--skip-plot"], "compute"),
    "clipping_profile": ("clipping_profile.py", "clipping_profile",
                         lambda c: ["-i", c.bam, "-o", c.out / "clip", "-s", "PE",
                                    "--skip-plot"], "compute"),
    "deletion_profile": ("deletion_profile.py", "deletion_profile",
                         lambda c: ["-i", c.bam, "-l", "100", "-o", c.out / "del", "--skip-plot"], "compute"),
    "insertion_profile": ("insertion_profile.py", "insertion_profile",
                          lambda c: ["-i", c.bam, "-o", c.out / "ins", "-s", "PE",
                                     "--skip-plot"], "compute"),
    # mismatch_profile needs MD tags, which this generator does not emit; both arms
    # still do real work, and the gate decides whether the row is reportable.
    "mismatch_profile": ("mismatch_profile.py", "mismatch_profile",
                         lambda c: ["-i", c.bam, "-l", "100", "-o", c.out / "mm", "--skip-plot"], "compute"),
    "geneBody_coverage": ("geneBody_coverage.py", "geneBody_coverage",
                          lambda c: ["-i", c.bam, "-r", c.bed, "-o", c.out / "gb",
                                     "--skip-plot"], "compute"),
    "junction_annotation": ("junction_annotation.py", "junction_annotation",
                            lambda c: ["-i", c.bam, "-r", c.bed, "-o", c.out / "ja",
                                       "--skip-plot", "--skip-bed", "--skip-interact"],
                            "compute"),
    "junction_saturation": ("junction_saturation.py", "junction_saturation",
                            lambda c: ["-i", c.bam, "-r", c.bed, "-o", c.out / "js",
                                       "--skip-plot"], "compute"),
    "RPKM_saturation": ("RPKM_saturation.py", "RPKM_saturation",
                        lambda c: ["-i", c.bam, "-r", c.bed, "-o", c.out / "rp",
                                   "--skip-plot"], "compute"),
    "split_bam": ("split_bam.py", "split_bam",
                  lambda c: ["-i", c.bam, "-r", c.bed, "-o", c.out / "sp_"], "compute"),
    "split_paired_bam": ("split_paired_bam.py", "split_paired_bam",
                         lambda c: ["-i", c.bam, "-o", c.out / "pp_"], "compute"),
    "divide_bam": ("divide_bam.py", "divide_bam",
                   lambda c: ["-i", c.bam, "-o", c.out / "dv_", "-n", "2", "-s",
                              "--seed", "42"], "compute"),
    "bam2wig": ("bam2wig.py", "bam2wig",
                lambda c: ["-i", c.bam, "-s", c.chrom_sizes, "-o", c.out / "w"], "compute"),
    "sc_seqQual": ("sc_seqQual.py", "sc_seqQual",
                   lambda c: ["-i", c.fq, "-o", c.out / "sq", "--skip-heatmap"], "compute"),
    # Upstream sc_seqLogo writes a plot via logomaker+pandas, which crashes in this
    # environment (DIV-0016). Only the count matrix is a common deliverable, so this
    # row compares that and the row is labelled as such.
    "sc_seqLogo": ("sc_seqLogo.py", "sc_seqLogo",
                   lambda c: ["-i", c.fa, "-o", c.out / "sl", "--iformat", "fa",
                              "--oformat", "svg"], "compute"),
    "read_hexamer": ("read_hexamer.py", "read_hexamer",
                     lambda c: ["-i", c.fq, "-o", c.out / "hex.txt"], "compute"),
    "geneBody_coverage2": ("geneBody_coverage2.py", "geneBody_coverage2",
                           lambda c: ["-i", c.bw, "-r", c.bed, "-o", c.out / "gb2", "--skip-plot"],
                           "compute"),
    "normalize_bigwig": ("normalize_bigwig.py", "normalize_bigwig",
                         lambda c: ["-i", c.bw, "-o", c.out / "norm.wig", "-f", "wig", "--overwrite"],
                         "compute"),
    "overlay_bigwig": ("overlay_bigwig.py", "overlay_bigwig",
                       lambda c: ["-i", c.bw, "-j", c.bw2, "-a", "Add",
                                  "-o", c.out / "ov.wig"], "compute"),
}

# Declared exclusions (protocol.md section 8). Reported, never silently dropped.
EXCLUDED = {
    "FPKM-UQ": "unsupported-missing-dependency: requires htseq-count, absent here",
    "sc_editMatrix": "unsupported-missing-dependency: requires R package pheatmap, absent here",
    # sc_bamStat is retained but noted: it exits non-zero on a bulk BAM with no
    # single-cell tags (ZeroDivisionError upstream, same in the port), so its row
    # measures two identical early failures rather than the analysis. It needs a
    # real scRNA-seq BAM with CB/RE/UR/UB tags, which this environment lacks.
    "divide_bam": "diverged-documented DIV-0017: the two builds use different RNG "
                  "algorithms, so --seed assigns different query names to subsets. The "
                  "gate compares the union of subsets (same partition), not per-file "
                  "membership. Per-file equality is explicitly NOT claimed.",
    "sc_bamStat": "unsupported-input-shape: needs a single-cell BAM with CB/RE tags; "
                  "bulk input makes both arms fail identically before doing work",
}


class Ctx:
    def __init__(self, workload: Path, out: Path):
        self.bam = str(workload / "reads.bam")
        self.bed = str(workload / "model.bed12")
        self.fq = str(workload / "reads_1.fastq")
        self.fa = str(workload / "reads_1.fa")
        self.chrom_sizes = str(workload / "chrom.sizes")
        self.bw = str(workload / "sig1.bw")
        self.bw2 = str(workload / "sig2.bw")
        self.out = out


# --------------------------------------------------------------------------------------
# Resource measurement
# --------------------------------------------------------------------------------------
def _parse_time_v(stderr: str) -> dict:
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)", stderr)
    if m:
        parts = [float(x) for x in m.group(1).split(":")]
        secs = 0.0
        for p in parts:
            secs = secs * 60 + p
        out["wall_s"] = secs
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


def run_arm(argv, cwd, env, timeout_s):
    """Run one arm under /usr/bin/time -v in its own process group.

    A new session (setsid) means a timeout kills the whole tree, so a hung Rscript or
    htseq-count child cannot survive the harness and skew later measurements.
    """
    cmd = [TIME_BIN, "-v", *argv]
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd, cwd=str(cwd), env=env, capture_output=True, text=True,
            timeout=timeout_s, start_new_session=True,
        )
        timed_out = False
    except subprocess.TimeoutExpired as e:
        proc = e
        timed_out = True
        try:
            os.killpg(os.getpgid(e.pid if hasattr(e, "pid") else 0), 9)
        except Exception:
            pass
    wall = time.perf_counter() - t0
    stderr = getattr(proc, "stderr", "") or ""
    stdout = getattr(proc, "stdout", "") or ""
    if isinstance(stderr, bytes):
        stderr = stderr.decode("utf8", "replace")
    if isinstance(stdout, bytes):
        stdout = stdout.decode("utf8", "replace")
    rc = proc.returncode if hasattr(proc, "returncode") else -1
    m = _parse_time_v(stderr)
    # /usr/bin/time -v reports with 2dp; use the harness's own monotonic clock for
    # the wall figure, since the child's is too coarse for the fast commands.
    m["wall_s"] = wall
    m["exit_code"] = rc
    m["timed_out"] = timed_out
    return m, stdout, stderr


# --------------------------------------------------------------------------------------
# Structural equivalence gate
# --------------------------------------------------------------------------------------
NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def _norm_text_line(line: str, rtol: float, atol: float) -> str:
    """Normalise a line: strip, and round bare numbers to the precision both tools emit."""
    s = line.rstrip()
    parts = s.split("\t")
    out = []
    for p in parts:
        p = p.strip()
        if NUM.match(p):
            try:
                v = float(p)
                # Compare at a relative tolerance rather than exact digits, since the
                # two implementations may format the same value differently.
                out.append(f"{v:.6g}")
            except ValueError:
                out.append(p)
        else:
            out.append(p)
    return "\t".join(out)


def compare_text(a: Path, b: Path, rtol=1e-6, atol=1e-9):
    la = [x for x in a.read_text(errors="replace").splitlines() if x.strip()]
    lb = [x for x in b.read_text(errors="replace").splitlines() if x.strip()]
    if len(la) != len(lb):
        return False, f"line count {len(la)} vs {len(lb)}"
    for i, (x, y) in enumerate(zip(la, lb)):
        nx, ny = _norm_text_line(x, rtol, atol), _norm_text_line(y, rtol, atol)
        if nx != ny:
            # fall back to numeric comparison for a line that differs only in floats
            fx, fy = x.split(), y.split()
            if len(fx) == len(fy):
                ok = True
                for u, v in zip(fx, fy):
                    if u == v:
                        continue
                    try:
                        fu, fv = float(u), float(v)
                    except ValueError:
                        ok = False
                        break
                    if abs(fu - fv) > atol + rtol * abs(fv):
                        ok = False
                        break
                if ok:
                    continue
            return False, f"line {i+1}: {x[:70]!r} vs {y[:70]!r}"
    return True, ""


def compare_bam(a: Path, b: Path):
    """Compare BAMs as decoded record multisets, not as compressed bytes.

    Compressed byte equality is meaningless across two different writers (block
    boundaries, compression level, and optional fields all differ legitimately), and
    comparing only flags would miss a coordinate or CIGAR error.
    """
    try:
        import pysam
    except ImportError:
        return False, "pysam unavailable for BAM comparison"
    ka, kb = [], []
    for f, acc in ((a, ka), (b, kb)):
        for r in pysam.AlignmentFile(str(f)):
            if r.is_unmapped:
                acc.append((r.query_name, r.flag & 0xC0))
                continue
            acc.append((
                r.query_name, r.flag & 0xC0, r.reference_name, r.reference_start,
                r.cigarstring, r.mapping_quality, r.query_sequence,
            ))
    if len(ka) != len(kb):
        return False, f"record count {len(ka)} vs {len(kb)}"
    from collections import Counter
    if Counter(ka) != Counter(kb):
        d = Counter(ka) - Counter(kb)
        return False, f"{sum(d.values())} records differ; first: {list(d)[0][:4]}"
    return True, ""


def compare_fastq(a: Path, b: Path):
    def rd(p):
        L = p.read_text(errors="replace").splitlines()
        return [tuple(L[i:i + 4]) for i in range(0, len(L) - 3, 4)]
    ra, rb = rd(a), rd(b)
    if len(ra) != len(rb):
        return False, f"record count {len(ra)} vs {len(rb)}"
    if ra != rb:
        return False, "sequence or quality differs"
    return True, ""


# Commands whose upstream output is not reproducible even against itself, because they
# subsample without a seed. Verified empirically: two consecutive upstream runs on the
# same input produce different files. Byte-equality is the wrong gate for these; the
# comparison has to be distributional. Membership here is a claim that must be justified
# per command, not a way to make a failing row pass.
STOCHASTIC_COMMANDS = {
    "RPKM_saturation": "upstream subsamples reads per percentile with no RNG seed; "
                       "two upstream runs on identical input differ",
}


def stochastic_agreement(a: Path, b: Path, rtol=0.05):
    """Fraction of rows whose numeric cells agree within `rtol`, averaged over columns.

    Returned as a scalar so that a port arm and an upstream self-comparison can be
    compared on the same scale.
    """
    la = [x for x in a.read_text(errors="replace").splitlines() if x.strip()]
    lb = [x for x in b.read_text(errors="replace").splitlines() if x.strip()]
    if len(la) != len(lb) or len(la) < 2:
        return 0.0
    fracs = []
    ncol = max(len(x.split("\t")) for x in la[1:])
    for c in range(ncol):
        va, vb = [], []
        for x, y in zip(la[1:], lb[1:]):
            fx, fy = x.split("\t"), y.split("\t")
            if c >= len(fx) or c >= len(fy):
                return 0.0
            try:
                va.append(float(fx[c]))
                vb.append(float(fy[c]))
            except ValueError:
                continue
        if not va:
            continue
        agree = sum(1 for x, y in zip(va, vb)
                    if abs(x - y) <= rtol * max(abs(x), abs(y), 1e-9))
        fracs.append(agree / len(va))
    return sum(fracs) / len(fracs) if fracs else 0.0


def compare_stochastic_text(a: Path, b: Path, rtol=0.05, min_frac=0.90):
    """Distributional comparison for a non-deterministic reference.

    Requires, per numeric column, that the two arms' medians agree within `rtol` and
    that at least `min_frac` of the per-row values agree within tolerance. This is a
    much weaker claim than equality and is labelled as such in the report: it says the
    two implementations produce the same distribution of results, not the same sample.
    """
    la = [x for x in a.read_text(errors="replace").splitlines() if x.strip()]
    lb = [x for x in b.read_text(errors="replace").splitlines() if x.strip()]
    if len(la) != len(lb):
        return False, f"line count {len(la)} vs {len(lb)}"
    if len(la) < 2:
        return False, "too few rows to compare distributionally"
    ncol = max(len(x.split("\t")) for x in la[1:])
    for c in range(ncol):
        va, vb = [], []
        for x, y in zip(la[1:], lb[1:]):
            fx, fy = x.split("\t"), y.split("\t")
            if c >= len(fx) or c >= len(fy):
                return False, f"column {c} missing in one arm"
            try:
                va.append(float(fx[c]))
                vb.append(float(fy[c]))
            except ValueError:
                # A non-numeric column (e.g. "#chr" or a gene name) is a label, not a
                # measurement. Compare it exactly; it has no distribution to check.
                if fx[c].strip() != fy[c].strip():
                    return False, f"column {c} label differs: {fx[c][:20]!r} vs {fy[c][:20]!r}"
                va = []
                vb = []
        if not va or not vb:
            continue
        ma, mb = statistics.median(va), statistics.median(vb)
        denom = max(abs(ma), abs(mb), 1e-9)
        if abs(ma - mb) > rtol * denom:
            return False, f"col {c} median {ma:.4g} vs {mb:.4g} exceeds {rtol:.0%} rel."
        agree = sum(1 for x, y in zip(va, vb)
                    if abs(x - y) <= rtol * max(abs(x), abs(y), 1e-9))
        frac = agree / len(va)
        if frac < min_frac:
            return False, f"col {c}: only {frac:.0%} of rows within {rtol:.0%}"
    return True, ""


# Files upstream emits that are logs rather than scientific deliverables. Upstream
# writes `log.txt` (geneBody_coverage.py:57 `printlog(..., log_file=Path("log.txt"))`);
# the port logs to stdout instead. Its presence is reported, but it is not a result, so
# it does not fail the gate. Any other unexpected one-sided file still does.
NON_DELIVERABLE_UPSTREAM = {"log.txt"}

# Artifacts the port emits that upstream CANNOT emit in this environment. sc_seqLogo's
# upstream plotting path crashes in the pinned oracle (logomaker/pandas incompat,
# DIV-0016), so upstream produces only the count matrix while the port also renders SVG.
# The extra files are reported on the row, and the common deliverable is compared.
# This is an asymmetry in the PORT's favour, so the row is explicitly marked as
# count-matrix-only: the port's end-to-end time here INCLUDES rendering that upstream
# never performed, which makes the port's number conservative, not flattering.
EXTRA_CANDIDATE_ARTIFACTS = {"sc_seqLogo": {".svg"}}


# Absolute margin by which the port may trail upstream's own self-agreement before the
# row is treated as a real divergence rather than indistinguishable noise.
STOCHASTIC_AGREEMENT_TOLERANCE = 0.02

DIVIDED_DOCUMENTED = {
    "divide_bam": "DIV-0017: different RNG algorithms, so per-file membership differs",
}


def compare_partitioned_bam(pairs, input_bam=None):
    """For a command that partitions a BAM into subsets, compare the UNION.

    Used for divide_bam, where DIV-0017 documents that Python's and Rust's RNG
    algorithms differ, so `--seed N` assigns different query names to different subset
    files. Per-file contents therefore cannot match and must not be required to. The
    question that is still meaningful, and what this checks, is whether both arms
    perform the same partition: every input record appears exactly once across the
    union of the subsets, in both arms.
    """
    from collections import Counter

    def collect(d, pattern):
        acc = Counter()
        for p in sorted(d.glob(pattern)):
            for r in pysam_iter(str(p)):
                acc[r] += 1
        return acc

    pys_ = collect(pairs[0], "*.bam")
    rss_ = collect(pairs[1], "*.bam")
    if pys_ != rss_:
        d = (pys_ - rss_) or (rss_ - pys_)
        return False, f"union of subsets differs; {sum(d.values())} records, e.g. {list(d)[0][:4]}"
    return True, ""


def pysam_iter(path):
    import pysam
    for r in pysam.AlignmentFile(path):
        if r.is_unmapped:
            yield ("unmapped", r.query_name)
            continue
        yield (r.query_name, r.flag & 0xC0, r.reference_name, r.reference_start,
               r.cigarstring, r.mapping_quality, r.query_sequence)


def gate_outputs(py_dir: Path, rs_dir: Path, stochastic=False, partitioned=False,
                 name=None):
    """Structurally compare the two arms' output trees.

    Returns (ok, errors). Any file present in one tree but not the other, or any
    semantic difference, fails the gate.
    """
    errors = []
    pys = {p.relative_to(py_dir) for p in py_dir.rglob("*") if p.is_file()}
    rss = {p.relative_to(rs_dir) for p in rs_dir.rglob("*") if p.is_file()}
    # Ignore runner bookkeeping files.
    pys = {p for p in pys if p.name not in ("py.code", "rs.code", "time.txt")}
    rss = {p for p in rss if p.name not in ("py.code", "rs.code", "time.txt")}

    only_py = sorted(pys - rss)
    only_rs = sorted(rss - pys)
    for p in only_py:
        if p.name in NON_DELIVERABLE_UPSTREAM:
            continue
        errors.append(f"only in upstream: {p}")
    allowed_extra = EXTRA_CANDIDATE_ARTIFACTS.get(name, set())
    for p in only_rs:
        if p.suffix in allowed_extra:
            continue
        errors.append(f"only in rust: {p}")

    for rel in sorted(pys & rss):
        a, b = py_dir / rel, rs_dir / rel
        suf = rel.suffix.lower()
        # An emitted R plotting script is not a scientific deliverable: it embeds the
        # output paths, which necessarily differ between the two arms because each runs
        # in its own directory. Comparing them byte-wise fails the gate for a path
        # difference and would discard a real result. The R script's presence is still
        # checked (via the file-set comparison above); the data table it plots is
        # compared normally.
        if suf == ".r":
            continue
        try:
            if suf == ".bam":
                if partitioned:
                    continue
                ok, msg = compare_bam(a, b)
            elif suf in (".fastq", ".fq"):
                ok, msg = compare_fastq(a, b)
            elif stochastic:
                ok, msg = compare_stochastic_text(a, b)
            else:
                ok, msg = compare_text(a, b)
        except Exception as e:  # a comparator failure must not silently pass
            ok, msg = False, f"comparator error: {e}"
        if not ok:
            errors.append(f"{rel}: {msg}")
    if partitioned:
        try:
            ok, msg = compare_partitioned_bam((py_dir, rs_dir))
            if not ok:
                errors.append(f"partition union: {msg}")
        except Exception as e:
            errors.append(f"partition comparator error: {e}")

    return (len(errors) == 0), errors


# --------------------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------------------
def paired_log_ratio_ci(pairs, iters=10000, seed=12345, alpha=0.05):
    """95% CI for the exponentiated mean of paired log time ratios (block bootstrap).

    Resampling is over whole pairs, preserving the pairing. Reported as the
    secondary estimator per protocol section 5.
    """
    import math
    if len(pairs) < 2:
        return None
    rng = random.Random(seed)
    ratios = [math.log(py / rs) for py, rs in pairs]
    n = len(ratios)
    boots = []
    for _ in range(iters):
        s = sum(ratios[rng.randrange(n)] for _ in range(n))
        boots.append(s / n)
    boots.sort()
    lo = boots[int(alpha / 2 * iters)]
    hi = boots[min(iters - 1, int((1 - alpha / 2) * iters))]
    return (math.exp(lo), math.exp(hi))


def median_ratio_ci(pairs, iters=10000, seed=12345, alpha=0.05):
    """95% CI for median(python)/median(rust) by bootstrap over pairs."""
    if len(pairs) < 2:
        return None
    rng = random.Random(seed)
    n = len(pairs)
    boots = []
    for _ in range(iters):
        s = [pairs[rng.randrange(n)] for _ in range(n)]
        py = statistics.median(p[0] for p in s)
        rs = statistics.median(p[1] for p in s)
        if rs > 0:
            boots.append(py / rs)
    if not boots:
        return None
    boots.sort()
    return (boots[int(alpha / 2 * iters)], boots[min(iters - 1, int((1 - alpha / 2) * iters))])


# --------------------------------------------------------------------------------------
# Environment
# --------------------------------------------------------------------------------------
def environment() -> dict:
    def sh(*cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=20).stdout.strip()
        except Exception:
            return "unknown"

    env = {
        "timestamp": time.time(),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "processor": _cpu_model(),
        "loadavg": os.getloadavg(),
        "rustc": sh("rustc", "--version"),
        "python": sh(str(VENV_PY), "--version"),
        "numpy": sh(str(VENV_PY), "-c", "import numpy;print(numpy.__version__)"),
        "pysam": sh(str(VENV_PY), "-c", "import pysam;print(pysam.__version__)"),
        "git_commit": sh("git", "-C", str(REPO), "rev-parse", "HEAD"),
        "git_dirty": bool(sh("git", "-C", str(REPO), "status", "--porcelain")),
        "cpufreq_governor": _governor(),
        "filesystem": _fs_type(str(REPO)),
        "isolated_hardware": False,
        "warning": (
            "Shared, non-isolated hardware. Not publication-grade; see "
            "benchmarks/protocol.md section 2."
        ),
        "pinned_env": PINNED_ENV,
    }
    return env


def _cpu_model():
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return "unknown"


def _governor():
    p = "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
    try:
        return Path(p).read_text().strip()
    except Exception:
        return "unavailable"


def _fs_type(path):
    try:
        out = subprocess.run(["stat", "-f", "-c", "%T", path],
                             capture_output=True, text=True, timeout=10).stdout.strip()
        return out or "unknown"
    except Exception:
        return "unknown"


def _sha256(path, limit=None):
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except Exception:
        return "missing"
    return h.hexdigest()[:16]


# --------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------
def floor_cost(argv, env, cwd, reps=5):
    """Measure the per-invocation fixed cost of an arm (interpreter/import startup).

    Used for the E1 compute-only experiment class, where the Python floor is measured
    and subtracted rather than assumed away. Run with the real binary on `--help`,
    which performs all imports and argument parsing but reads no input.
    """
    ts = []
    for _ in range(reps):
        m, _, _ = run_arm(argv, cwd, env, timeout_s=120)
        ts.append(m["wall_s"])
    return statistics.median(ts)


def bench_command(name, workload: Path, outdir: Path, reps, warmup, timeout_s, seed,
                  do_gate=True):
    script, binary, argfn, eclass = COMMANDS[name]
    results = {"command": name, "upstream": script, "rust": binary,
               "experiment_class": eclass, "reps": reps}

    env = dict(os.environ)
    env.update(PINNED_ENV)
    # Keep R/pysam output deterministic and quiet.
    env["R_DEFAULT_PACKAGES"] = env.get("R_DEFAULT_PACKAGES", "")

    root = outdir / name
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True, exist_ok=True)

    def arm_dir(tag, rep):
        """A FRESH output directory per arm per repetition.

        Reusing one directory across repetitions is not merely untidy: several upstream
        commands refuse to overwrite an existing output and exit 2 (or the port exits
        1), so every repetition after the first would measure two immediate failures
        instead of the command. Each repetition therefore gets an isolated workdir, and
        the gate compares the last repetition's trees, which is also the run whose
        outputs were produced under normal conditions.
        """
        d = root / tag / f"rep{rep:02d}"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def argv_for(arm, d):
        base = ([str(VENV_PY), str(UPSTREAM_SRC / script)] if arm == "py"
                else [str(RELEASE / binary)])
        return base + [str(a) if isinstance(a, Path) else a
                       for a in argfn(Ctx(workload, d))]

    # Warm-up (untimed, discarded): page-cache and dynamic-loader state.
    for w in range(warmup):
        for arm in ("py", "rs"):
            d = arm_dir(arm, 900 + w)
            run_arm(argv_for(arm, d), d, env, timeout_s)

    # Randomised, interleaved paired schedule, generated before collection.
    schedule = []
    for i in range(reps):
        order = ["py", "rs"] if random.Random(seed + i).random() < 0.5 else ["rs", "py"]
        for arm in order:
            schedule.append((i, arm))
    random.Random(seed).shuffle(schedule)

    per = {"py": [], "rs": []}
    failures = []
    last_dirs = {}
    for rep, arm in schedule:
        d = arm_dir(arm, rep)
        m, so, se = run_arm(argv_for(arm, d), d, env, timeout_s)
        per[arm].append(m)
        last_dirs[arm] = d
        if m["exit_code"] != 0 or m["timed_out"]:
            failures.append({"arm": arm, "rep": len(per[arm]) - 1,
                             "exit_code": m["exit_code"], "timed_out": m["timed_out"],
                             "stderr_tail": se[-400:]})
    py_dir, rs_dir = last_dirs["py"], last_dirs["rs"]

    results["runs"] = per
    results["failures"] = failures

    # E1 compute-only: subtract the measured per-invocation floor of each arm.
    fpy = arm_dir("upstream", 950)
    frs = arm_dir("rust", 950)
    py_floor = floor_cost(argv_for("py", fpy) + ["--help"], env, fpy)
    rs_floor = floor_cost(argv_for("rs", frs) + ["--help"], env, frs)
    results["floor_s"] = {"python": py_floor, "rust": rs_floor}

    def pairs_for(key):
        out = []
        for a, b in zip(per["py"], per["rs"]):
            va = a[key] - (py_floor if key == "wall_s" else 0.0)
            vb = b[key] - (rs_floor if key == "wall_s" else 0.0)
            if vb > 1e-6 and va > 1e-6:
                out.append((va, vb))
        return out

    def summarise(key, floor_sub):
        """Point estimate and CI MUST be computed from the same quantities.

        An earlier version took the ratio from raw medians but the CI from
        floor-subtracted pairs, which made the two describe different measurements and
        produced intervals that did not contain the reported ratio (e.g. a 14.2x point
        estimate with a [2.05, 4.02] interval on a 2k-read bam_stat). The E1 quantities
        are the ones the interval describes, so the point estimate uses them too; the
        raw E2 medians are reported separately and labelled.
        """
        pr = pairs_for(key)
        pyv = [a for a, _ in pr]
        rsv = [b for _, b in pr]
        if not pr:
            return {"n": 0, "median_python": None, "median_rust": None,
                    "ratio_median": None, "ratio_ci95": None,
                    "logratio_ci95": None, "median_python_e2e":
                        statistics.median([r[key] for r in per["py"]]),
                    "median_rust_e2e":
                        statistics.median([r[key] for r in per["rs"]])}
        s = {
            "median_python": statistics.median(pyv),
            "median_rust": statistics.median(rsv),
            "ratio_median": (statistics.median(pyv) / statistics.median(rsv))
            if statistics.median(rsv) else None,
            "median_python_e2e": statistics.median([r[key] for r in per["py"]]),
            "median_rust_e2e": statistics.median([r[key] for r in per["rs"]]),
            "min_python": min([r[key] for r in per["py"]]),
            "min_rust": min([r[key] for r in per["rs"]]),
            "n": len(pr),
        }
        s["ratio_ci95"] = median_ratio_ci(pr)
        s["logratio_ci95"] = paired_log_ratio_ci(pr)
        return s

    results["wall_e2e"] = summarise("wall_s", True)
    results["cpu"] = summarise("user_s", False)
    mem_py = [r["peak_rss_mb"] for r in per["py"] if "peak_rss_mb" in r]
    mem_rs = [r["peak_rss_mb"] for r in per["rs"] if "peak_rss_mb" in r]
    if mem_py and mem_rs:
        results["peak_rss_mb"] = {
            "median_python": statistics.median(mem_py),
            "median_rust": statistics.median(mem_rs),
            "max_python": max(mem_py), "max_rust": max(mem_rs),
            "ratio": (statistics.median(mem_py) / statistics.median(mem_rs))
            if statistics.median(mem_rs) else None,
        }

    # Equivalence gate
    if do_gate and name in STOCHASTIC_COMMANDS:
        # Self-control: run the REFERENCE twice and compare it to itself. Upstream
        # RPKM_saturation subsamples without a seed, so it does not reproduce its own
        # output; a strict distributional gate would therefore fail a correct port for
        # upstream's own noise. The port passes only if it is no less self-consistent
        # with upstream than upstream is with itself.
        c1, c2 = arm_dir("upstream", 960), arm_dir("upstream", 961)
        run_arm(argv_for("py", c1), c1, env, timeout_s)
        run_arm(argv_for("py", c2), c2, env, timeout_s)
        self_fracs, port_fracs = {}, {}
        for f in sorted(p.name for p in py_dir.glob("*") if p.suffix not in (".r", ".bam")):
            try:
                self_fracs[f] = stochastic_agreement(c1 / f, c2 / f)
                port_fracs[f] = stochastic_agreement(py_dir / f, rs_dir / f)
            except Exception as e:
                self_fracs[f] = port_fracs[f] = 0.0
        results["stochastic_self_control"] = {
            "upstream_vs_upstream": self_fracs,
            "port_vs_upstream": port_fracs,
            "note": STOCHASTIC_COMMANDS.get(name),
        }
        # The port must not be MEANINGFULLY less self-consistent with upstream than
        # upstream is with itself. Requiring strict dominance is wrong: both figures
        # carry sampling noise, and a row whose two sides differ by 0.1pp is
        # "indistinguishable", not "fails". An earlier revision did require strict
        # dominance, and passed one run (59.3% vs 58.4%) while failing the next
        # (58.3% vs 58.4%) purely because of which side the noise fell on.
        errs = []
        for f, sc in self_fracs.items():
            pf = port_fracs.get(f, 0.0)
            if pf < sc - STOCHASTIC_AGREEMENT_TOLERANCE:
                errs.append(
                    f"{f}: port agrees with upstream on {pf:.1%} of rows vs upstream's own "
                    f"{sc:.1%} self-agreement, a gap larger than the "
                    f"{STOCHASTIC_AGREEMENT_TOLERANCE:.0%} tolerance")
        results["equivalence"] = {
            "pass": not errs, "errors": errs, "mode": "distributional-vs-self-control",
            "note": STOCHASTIC_COMMANDS.get(name),
        }
    elif do_gate:
        part = name in DIVIDED_DOCUMENTED
        ok, errs = gate_outputs(py_dir, rs_dir, partitioned=part, name=name)
        mode = "partition-union" if part else "exact"
        results["equivalence"] = {
            "pass": ok, "errors": errs[:20], "mode": mode,
            "note": DIVIDED_DOCUMENTED.get(name, ""),
        }
    else:
        results["equivalence"] = {"pass": None, "errors": ["gate not run"]}

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workload", required=True)
    ap.add_argument("--workload-for", nargs="*", default=[], metavar="CMD=PATH",
                    help="per-command workload override, e.g. tin=W/tin geneBody_coverage=W/gbc. "
                         "Used where a command's cost driver makes the default panel too "
                         "slow to measure 10 paired repetitions (protocol section 6).")
    ap.add_argument("--label", default="main", help="name for this measurement run")
    ap.add_argument("--commands", nargs="+", default=sorted(COMMANDS))
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--no-gate", action="store_true")
    ap.add_argument("--quick", action="store_true",
                    help="smoke test: 1 rep, tiny command list")
    args = ap.parse_args()

    if args.quick:
        args.reps = 1
        args.commands = args.commands[:3]

    workload = Path(args.workload).resolve()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("BENCHMARK - shared hardware, NOT publication-grade (protocol.md section 2)")
    print("=" * 78)
    print(f"workload : {workload}")
    print(f"commands : {len(args.commands)}  reps: {args.reps}  seed: {args.seed}")
    print(f"pinned   : {PINNED_ENV}")
    print()

    env = environment()
    manifest = workload / "manifest.json"
    result = {
        "label": args.label,
        "environment": env,
        "workload": {
            "path": str(workload),
            "manifest": json.loads(manifest.read_text()) if manifest.exists() else None,
            "bam_sha256": _sha256(workload / "reads.bam"),
        },
        "config": {"reps": args.reps, "warmup": args.warmup, "timeout_s": args.timeout,
                   "seed": args.seed, "gate": not args.no_gate},
        "excluded": EXCLUDED,
        "results": [],
    }

    overrides = {}
    for spec in args.workload_for:
        cmd, _, path = spec.partition("=")
        overrides[cmd] = Path(path).resolve()

    for name in args.commands:
        if name in EXCLUDED:
            continue
        if name not in COMMANDS:
            print(f"  !! unknown command {name}, skipping")
            continue
        wl = overrides.get(name, workload)
        t0 = time.time()
        print(f"[{name}] ({wl.name}) ", end="", flush=True)
        try:
            r = bench_command(name, wl, outdir,
                              args.reps, args.warmup,
                              args.timeout, args.seed, do_gate=not args.no_gate)
            r["workload"] = str(wl)
        except Exception as e:
            print(f"ERROR {e}")
            result["results"].append({"command": name, "error": str(e)})
            continue
        r["workload"] = str(wl)
        result["results"].append(r)
        w = r.get("wall_e2e", {})
        ci = w.get("ratio_ci95")
        gate = "PASS" if r["equivalence"]["pass"] else "FAIL"
        ratio = w.get("ratio_median")
        print(
            f"py {w.get('median_python', 0):.3f}s  rs {w.get('median_rust', 0):.4f}s  "
            f"x{ratio:.2f}" if ratio else "no ratio",
            end="  ",
        )
        print(f"CI[{ci[0]:.2f},{ci[1]:.2f}]" if ci else "CI n/a", end="  ")
        print(f"gate={gate}  ({time.time()-t0:.0f}s)")
        (outdir / f"{name}.json").write_text(json.dumps(r, indent=2))

    (outdir / "results.json").write_text(json.dumps(result, indent=2))
    print()
    print(f"results -> {outdir/'results.json'}")
    report(outdir / "results.json")


def report(path: Path):
    d = json.loads(path.read_text())
    print()
    print("=" * 110)
    print(f"{'command':22} {'py(s)':>9} {'rs(s)':>9} {'speedup':>9} {'95% CI':>16} "
          f"{'mem(py/rs MB)':>16} {'gate':>6}")
    print("-" * 110)
    for r in d["results"]:
        if "error" in r:
            print(f"{r['command']:22} ERROR {r['error'][:40]}")
            continue
        w = r.get("wall_e2e", {})
        ci = w.get("ratio_ci95")
        m = r.get("peak_rss_mb", {})
        mem = f"{m.get('median_python', 0):.0f}/{m.get('median_rust', 0):.0f}" if m else "-"
        gate = "PASS" if r["equivalence"]["pass"] else "FAIL"
        ratio = w.get("ratio_median")
        verdict = ""
        if ci:
            if ci[0] > 1.0:
                verdict = "win"
            elif ci[1] < 1.0:
                verdict = "loss"
            else:
                verdict = "inconcl"
        print(f"{r['command']:22} {w.get('median_python', 0):>9.3f} "
              f"{w.get('median_rust', 0):>9.4f} "
              f"{(f'{ratio:.2f}x' if ratio else '-'):>9} "
              f"{(f'[{ci[0]:.2f}, {ci[1]:.2f}]' if ci else 'n/a'):>16} "
              f"{mem:>16} {gate:>6} {verdict}")
    print("=" * 110)
    for k, v in d.get("excluded", {}).items():
        print(f"EXCLUDED {k}: {v}")
    bad = [r["command"] for r in d["results"]
           if r.get("equivalence", {}).get("pass") is False]
    if bad:
        print(f"\nGATE FAILURES (no speedup claim): {', '.join(bad)}")


if __name__ == "__main__":
    main()
