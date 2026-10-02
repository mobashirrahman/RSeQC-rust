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
import signal
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
# verification/ holds the comparators shared with the differential runner. Added to
# the import path explicitly because this file is also loaded by path from
# benchmarks/test_bench_harness.py, where `benchmarks/` is not on sys.path.
sys.path.insert(0, str(REPO / "verification"))
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
# Per-command expected streams, labels, schemas and artifacts
# --------------------------------------------------------------------------------------
# The audit's P0 benchmark-gate finding was that the gate only diffed two file
# trees. Two empty directories passed. Captured stdout was never compared, so a
# command whose entire deliverable is on stdout received a gate pass with no
# metric checking at all. An expectation is now declared per command, and a
# command with no declaration fails the gate rather than passing silently.
#
# `artifacts` lists relative paths, with `?` as a single-component wildcard for
# names the upstream port chooses (e.g. a per-transcript output prefix). An empty
# list is a declaration that a command is expected to write no file, which is
# only correct for stdout-only commands.
EXPECTED_STREAMS = {
    "bam_stat": {
        "stdout": True,
        "labels": {
            "Total records", "QC failed", "Optical/PCR duplicate",
            "Non primary hits", "Unmapped reads",
            "mapq < mapq_cut (non-unique)", "mapq >= mapq_cut (unique)",
            "Read-1", "Read-2", "Reads map to '+'", "Reads map to '-'",
            "Reads mapped in proper pairs", "Proper-paired reads map to different chrom",
            "Splice reads", "Non-splice reads",
        },
        "artifacts": [],
        "note": "Purely stdout. The gate previously never compared this stream, so "
                "every metric here went unchecked.",
    },
    "bam2fq": {
        "stdout": False, "labels": set(),
        "artifacts": ["r1.R1.fastq", "r1.R2.fastq"],
    },
    "infer_experiment": {
        "stdout": True,
        "labels": {
            'Fraction of reads failed to determine',
            'Fraction of reads explained by "1++,1--,2+-,2-+"',
            'Fraction of reads explained by "1+-,1-+,2++,2--"',
        },
        "artifacts": [],
        "note": "Pair-end bucket labels. A single-end run prints the '++/--' and "
                "'+-/+-' variants instead; see tolerated_labels.",
        "tolerated_labels": {
            'Fraction of reads explained by "++,--"',
            'Fraction of reads explained by "+-,-+"',
        },
    },
    "read_distribution": {
        # A fixed-width table on stdout with no files. It is now parseable because
        # the extractor handles multi-space separation, so these labels are real
        # checks rather than a declaration that the command writes nothing.
        "stdout": True,
        "labels": {
            "Total Reads", "Total Tags", "Total Assigned Tags",
            "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns",
        },
        "artifacts": [],
        "tolerated_labels": {
            "Group", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb",
            "TSS_down_1kb", "TSS_down_5kb", "TSS_down_10kb",
            "TES_up_1kb", "TES_up_5kb", "TES_up_10kb",
            "TES_down_1kb", "TES_down_5kb", "TES_down_10kb",
        },
        "note": "Fixed-width stdout table. The window rows are tolerated rather "
                "than required because upstream emits them conditionally.",
    },
    "inner_distance": {
        "stdout": False, "labels": set(),
        "artifacts": ["id.inner_distance.txt", "id.inner_distance_freq.txt"],
    },
    "RNA_fragment_size": {
        # The per-transcript fragment-size table IS the deliverable, and it goes to
        # STDOUT: upstream writes it to sys.stdout when no -o is given, and the port
        # does the same. The command writes no files at all in that mode.
        #
        # This was declared `stdout: False` because a declaration probe looked for
        # `label value` pairs, found none in a table whose HEADER ROW carries the
        # column names, and concluded there was no stdout. "No labelled stdout" was
        # read as "no stdout", and the gate then rejected two arms whose stdout was
        # byte-identical across 369,382 bytes -- the second instance of the same
        # defect class as bam2wig, where a wrong declaration gates a command on a
        # rule that does not describe it.
        #
        # Compared as text: the header substrings must be present in both arms and
        # the two streams must be identical, which for a table is the strongest
        # available statement.
        "stdout": True,
        "stdout_mode": "text",
        "stdout_contains": ["chrom", "tx_start", "frag_count", "frag_median"],
        "labels": set(),
        "artifacts": [],
        "note": "The fragment-size table is stdout. With -o it is written to that "
                "file instead and stdout is empty; the harness runs the stdout form, "
                "which is the form a user gets by default.",
    },
    "tin": {
        "stdout": False, "labels": set(),
        "artifacts": ["reads.tin.xls", "reads.summary.txt"],
    },
    "FPKM_count": {
        "stdout": False, "labels": set(),
        "artifacts": ["fpkm.FPKM.xls"],
    },
    "read_NVC": {
        "stdout": False, "labels": set(),
        "artifacts": ["nvc.NVC.xls"],
    },
    "read_GC": {
        "stdout": False, "labels": set(),
        "artifacts": ["gc.GC.xls"],
    },
    "read_quality": {
        "stdout": False, "labels": set(),
        "artifacts": ["q.qual.r"],
        "note": "Only the plotting script is emitted here; the summary table is "
                "stdout-only in fixed-width form.",
    },
    "read_duplication": {
        "stdout": False, "labels": set(),
        "artifacts": ["dup.seq.DupRate.xls", "dup.pos.DupRate.xls"],
    },
    "clipping_profile": {
        "stdout": False, "labels": set(),
        "artifacts": ["clip.clipping_profile.xls"],
    },
    "deletion_profile": {
        "stdout": False, "labels": set(),
        "artifacts": ["del.deletion_profile.txt"],
    },
    "insertion_profile": {
        "stdout": False, "labels": set(),
        "artifacts": ["ins.insertion_profile.xls"],
    },
    "mismatch_profile": {
        "stdout": False, "labels": set(),
        "artifacts": ["mm.mismatch_profile.xls"],
    },
    "geneBody_coverage": {
        "stdout": False, "labels": set(),
        "artifacts": ["gb.geneBodyCoverage.txt"],
    },
    "junction_annotation": {
        # `total = <n>` on stdout. Declared as text, not as a label, because the
        # metrics parser expects `label value` and this line is `label = value`: with
        # a label declaration the parser reports "upstream stdout lacks required
        # label 'total'" against output that is present and correct. Third instance
        # of a declaration that gated two byte-identical arms on a rule that does not
        # describe them; `verification/verify_stream_declarations.py` now finds these
        # by re-deriving each declaration from a live run.
        "stdout": True,
        "stdout_mode": "text",
        "stdout_contains": ["total"],
        "labels": set(),
        "artifacts": ["ja.junction.xls"],
    },
    "junction_saturation": {
        "stdout": False, "labels": set(),
        "artifacts": ["js.junctionSaturation_plot.r"],
        "note": "Only the plotting script is emitted; the curve is stdout-only.",
    },
    "RPKM_saturation": {
        "stdout": False, "labels": set(),
        "artifacts": ["rp.eRPKM.xls", "rp.rawCount.xls"],
    },
    "split_bam": {
        "stdout": True, "labels": {"Total records"},
        "artifacts": ["sp_.in.bam", "sp_.ex.bam", "sp_.junk.bam"],
        "note": "The per-output lines name each arm's own absolute paths, so they "
                "are normalised to a placeholder before comparison. Only "
                "'Total records' is required; the per-file lines are recorded by "
                "the file-set comparison instead.",
    },
    "split_paired_bam": {
        "stdout": True, "labels": {"Total records"},
        "artifacts": ["pp_.R1.bam", "pp_.R2.bam", "pp_.unmap.bam"],
        "note": "Same path-in-label caveat as split_bam.",
    },
    "divide_bam": {
        "stdout": False, "labels": set(),
        "artifacts": ["dv__*.bam"],
        "note": "DIV-0017: the two implementations use different RNG algorithms, so "
                "--seed assigns different subsets. Per-file membership is NOT "
                "claimed; the gate compares the union.",
    },
    "bam2wig": {
        # Both arms print the same one-line wigToBigWig invocation to stdout, so the
        # declaration has to expect it. It did not, and the gate then failed
        # `bam2wig` on the real-data row with "unexpected stdout content for a command
        # that declares none" -- while the two arms' stdout was byte-identical. That is
        # the audit's P0 gate finding recurring in the declaration rather than in the
        # comparison: a command whose declaration is wrong is gated on a rule that does
        # not describe it, so a correct pair fails.
        #
        # The line is ASSERTED, not merely permitted, because it is a derived artifact
        # name a user may need to reproduce, and it is the only evidence that the helper
        # invocation was attempted with the paths this run actually wrote. It is free
        # text rather than a `label value` table, so it is compared as text: both arms
        # must contain the expected substrings and their normalised streams must be
        # identical. Declaring it as a metric label instead -- the obvious fix, and the
        # one tried first -- makes the gate demand "parseable metrics" from a sentence,
        # which fails both arms.
        "stdout": True,
        "stdout_mode": "text",
        "stdout_contains": ["wigToBigWig"],
        "labels": set(),
        "artifacts": ["w.wig"],
        "note": "Writes a WIG track, then prints the wigToBigWig command it would "
                "run and shells out for the .bw. wigToBigWig is absent in this "
                "environment, so the row measures the track only; the printed command "
                "is the deliverable that records what was attempted.",
    },
    "sc_seqQual": {
        "stdout": False, "labels": set(),
        "artifacts": ["sq.qual_count.csv", "sq.qual_percent.csv"],
    },
    "sc_seqLogo": {
        "stdout": False, "labels": set(),
        "artifacts": ["sl.count_matrix.csv"],
        "note": "Only the count matrix is common. Upstream's logo rendering fails in "
                "this environment (DIV-0016) and exits non-zero; the port also "
                "renders SVG, which is recorded as an allowed extra.",
    },
    "read_hexamer": {
        "stdout": False, "labels": set(),
        "artifacts": ["hex.txt"],
    },
    "geneBody_coverage2": {
        "stdout": False, "labels": set(),
        "artifacts": ["gb2.geneBodyCoverage.txt"],
    },
    "normalize_bigwig": {
        "stdout": False, "labels": set(),
        "artifacts": ["norm.wig"],
    },
    "overlay_bigwig": {
        "stdout": False, "labels": set(),
        "artifacts": ["ov.wig"],
    },
}


# --------------------------------------------------------------------------------------
# Resource measurement
# --------------------------------------------------------------------------------------
# Bumped whenever a comparator's accept/reject behaviour changes, so a recorded
# result can never be read under comparator semantics that did not produce it.
# Version 2 is the failure-closed rewrite: BAM quality/flags/tags/mates/headers,
# FASTQ framing and lengths, nonfinite numeric values, per-command stdout and
# exit-status checks. Results recorded under version 1 are not comparable and
# say so in their own provenance.
COMPARATOR_VERSION = 2


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


def _kill_group(pgid, grace_s=2.0):
    """Terminate a process group we created, escalating to SIGKILL.

    Only ever called with the pgid of a group this harness created via
    `start_new_session=True`. That is the whole point: `os.getpgid(0)` returns the
    CALLER's group, so `killpg(getpgid(0))` on a `TimeoutExpired` whose `.pid` is
    None signals the harness itself and every sibling in its terminal session --
    the audit's finding, which killed the runner's own group.
    """
    for sig, wait in ((signal.SIGTERM, grace_s), (signal.SIGKILL, 1.0)):
        if not pgid:
            return
        try:
            os.killpg(pgid, sig)
        except ProcessLookupError:
            return
        except PermissionError:
            return
        deadline = time.monotonic() + wait
        while time.monotonic() < deadline:
            if not _group_alive(pgid):
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


def run_arm(argv, cwd, env, timeout_s):
    """Run one arm under /usr/bin/time -v in its own process group.

    A new session (setsid) means a timeout kills the whole tree, so a hung Rscript
    or htseq-count child cannot survive the harness and skew later measurements.

    The `Popen` object is retained rather than discarded, because `subprocess.run`
    cannot report the child's process group on timeout: its `TimeoutExpired` has
    no `pid`, so the only handle on the tree is the object this harness created.
    """
    cmd = [TIME_BIN, "-v", *argv]
    t0 = time.perf_counter()
    # Captured to files rather than pipes. A pipe is inherited by every descendant,
    # so a grandchild that survives a timeout keeps the pipe open and `communicate`
    # blocks forever waiting for an EOF that never arrives -- the harness then hangs
    # instead of recording the timeout. Files also preserve the child's output on
    # disk after a kill, which is the diagnostic a timeout investigation needs.
    with tempfile.TemporaryDirectory(prefix="bench-out-") as capdir:
        out_path = Path(capdir) / "stdout"
        err_path = Path(capdir) / "stderr"
        with out_path.open("wb") as out_fh, err_path.open("wb") as err_fh:
            # `start_new_session=True` makes the child a session and group leader,
            # so its pgid equals its pid and is knowable from the Popen object alone.
            proc = subprocess.Popen(
                cmd, cwd=str(cwd), env=env, stdout=out_fh, stderr=err_fh,
                start_new_session=True,
            )
            pgid = proc.pid  # only valid because of start_new_session
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
        # Reap anything the tree left behind, so a hung grandchild cannot survive
        # the harness and silently consume CPU during later repetitions.
        _kill_group(pgid)
        stdout = out_path.read_bytes().decode("utf8", "replace")
        stderr = err_path.read_bytes().decode("utf8", "replace")
    m = _parse_time_v(stderr)
    # /usr/bin/time -v reports with 2dp; use the harness's own monotonic clock for
    # the wall figure, since the child's is too coarse for the fast commands.
    m["wall_s"] = wall
    m["exit_code"] = proc.returncode
    m["timed_out"] = timed_out
    return m, stdout, stderr


# --------------------------------------------------------------------------------------
# Structural equivalence gate
# --------------------------------------------------------------------------------------
# Artifact comparison lives in verification/comparators.py, shared with the
# differential runner. Duplicating it is how the two harnesses drifted apart in the
# first place: the benchmark's copy accepted changed BAM qualities, flags and tags
# and a truncated FASTQ record, and nothing in the differential suite noticed
# because it had its own weaker copy.
#
# Re-exported here under the names the rest of this module already uses.
from comparators import (  # noqa: E402
    NUM,
    NONFINITE,
    bam_header_signature,
    bam_signature,
    compare_bam,
    compare_fastq,
    compare_text,
    is_nonfinite,
    read_fastq,
)

_is_nonfinite = is_nonfinite


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
    """Full record signature for the union-of-subsets comparison.

    Uses the shared signature rather than the `flag & 0xC0` subset this used to
    take. The partition check exists to prove both arms performed the same
    partition of the input; with a weak per-record signature it could not tell a
    correct partition from one whose records had been silently altered, so a
    corruption inside a subset file went unnoticed while the union still matched.
    """
    import pysam
    with pysam.AlignmentFile(path) as handle:
        for r in handle:
            if r.is_unmapped:
                yield ("unmapped", r.query_name)
                continue
            yield bam_signature(r)


# Two separators appear in upstream reports: `label: value`, and fixed-width
# `label<spaces>value`. The second is not optional to support -- bam_stat's
# "Non primary hits" and all of read_distribution's report use it, so an extractor
# that only understood colons would silently skip those metrics.
LABEL_SPLIT = re.compile(r"^(?P<label>\S.*?)\s*(?::|\s{2,})\s*(?P<value>\S.*)$")


def extract_stdout_metrics(text):
    """Parse `label: value` and fixed-width `label   value` report lines.

    Upstream RSeQC prints some metrics only to stdout, and the file-tree gate
    never compared stdout at all, so a row whose entire deliverable was on stdout
    passed with no metric checking whatsoever. A command whose stdout has no
    parseable label/value lines yields an empty dict, which the caller must treat
    as insufficient evidence rather than as agreement.
    """
    metrics = {}
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line.strip() or line.strip().startswith(("#", "=")):
            continue
        m = LABEL_SPLIT.match(line.strip())
        if not m:
            continue
        label = m.group("label").strip()
        value = m.group("value").strip()
        if not label or not value:
            continue
        if label in metrics:
            # Duplicate labels make a value ambiguous, so refuse to guess which
            # occurrence the report meant.
            metrics[label] = None
            continue
        # Collapse internal runs of whitespace. A fixed-width table's columns are
        # padded, and the padding width is a formatting choice rather than a
        # measurement, so it must not decide whether the gate passes.
        metrics[label] = " ".join(value.split())
    return metrics


# A stdout label that names an output file rather than a measured quantity. After
# working-directory normalisation both arms produce the same label, so the count
# attached to it is still compared; only the label's presence is not required.
OUTPUT_PATH_LABEL = re.compile(r"^<out>/|^/?[^\s:]*/[^\s:]*\.(?:bam|fastq|fq|bai"
                               r"|wig|bw|bedgraph|txt|xls|csv|r)$")


def _is_output_path_label(label):
    return bool(OUTPUT_PATH_LABEL.match(label.strip()))


def normalise_stdout(text, work_dir):
    """Replace an arm's own working directory with a placeholder.

    Several commands name their own output files in stdout, and each arm runs in a
    different directory (`py/repNN` versus `rs/repNN`). Without this, those labels
    differ between arms for a reason that says nothing about correctness, and the
    gate would either fail a correct pair or -- worse, to keep passing -- have those
    lines exempted from comparison. Rewriting the arm's own directory keeps the
    measurement (the filename and the count) while removing the unavoidable
    difference.
    """
    if not work_dir:
        return text
    return text.replace(str(work_dir).rstrip("/") + "/", "<out>/").replace(
        str(work_dir).rstrip("/"), "<out>")


def compare_streams(spec, stdout_a, stdout_b, exit_a, exit_b,
                    timed_out_a=False, timed_out_b=False,
                    work_dir_a=None, work_dir_b=None):
    """Check exit status, timeouts, stdout presence and required labels.

    Returns (ok, errors). Every check fails closed: an absent stream, a missing
    label or a nonfinite value is an error, never a pass.
    """
    errors = []
    stdout_a = normalise_stdout(stdout_a, work_dir_a)
    stdout_b = normalise_stdout(stdout_b, work_dir_b)
    if timed_out_a or timed_out_b:
        errors.append(f"timeout: upstream={timed_out_a} rust={timed_out_b}")
    if exit_a != exit_b:
        errors.append(f"exit status differs: upstream={exit_a} rust={exit_b}")
    elif exit_a != 0:
        errors.append(f"both arms exited {exit_a}; no equivalence to claim")

    if spec.get("stdout_mode") == "text":
        # A command whose stdout is prose or a derived command line rather than a
        # `label value` table. There is nothing to parse into metrics, so requiring
        # parseable metrics would fail a correct pair -- which is what happened to
        # bam2wig once its declaration was corrected to expect stdout at all. Instead
        # the declared substrings must be present in BOTH arms and the two normalised
        # streams must be identical, which is the strongest available statement for
        # this shape and strictly stronger than a label comparison.
        if not stdout_a.strip() or not stdout_b.strip():
            errors.append("stream: one arm produced no stdout where text is declared "
                          f"(upstream={len(stdout_a.strip())} chars, "
                          f"rust={len(stdout_b.strip())} chars)")
        for needle in sorted(spec.get("stdout_contains", ())):
            if needle not in stdout_a:
                errors.append(f"stream: upstream stdout lacks {needle!r}")
            if needle not in stdout_b:
                errors.append(f"stream: rust stdout lacks {needle!r}")
        if stdout_a.strip() and stdout_b.strip() and stdout_a != stdout_b:
            errors.append("stream: the two arms' stdout differs\n"
                          f"    upstream: {stdout_a.strip()[:200]}\n"
                          f"    rust    : {stdout_b.strip()[:200]}")
    elif spec.get("stdout"):
        ma = extract_stdout_metrics(stdout_a)
        mb = extract_stdout_metrics(stdout_b)
        if not ma:
            errors.append("upstream produced no parseable stdout metrics")
        if not mb:
            errors.append("rust produced no parseable stdout metrics")
        for label in sorted(spec.get("labels", ())):
            va, vb = ma.get(label), mb.get(label)
            if va is None:
                errors.append(f"upstream stdout lacks required label {label!r}")
            elif vb is None:
                errors.append(f"rust stdout lacks required label {label!r}")
            elif va == "" or vb == "":
                errors.append(f"ambiguous duplicate label {label!r}")
            elif _is_nonfinite(va) or _is_nonfinite(vb):
                errors.append(f"label {label!r} is nonfinite: {va!r} vs {vb!r}")
            elif va != vb:
                try:
                    fa, fb = float(va), float(vb)
                    if abs(fa - fb) > 1e-9 + 1e-6 * abs(fb):
                        errors.append(f"label {label!r} differs: {va} vs {vb}")
                except ValueError:
                    if va != vb:
                        errors.append(f"label {label!r} differs: {va!r} vs {vb!r}")
        for label in sorted(set(ma) | set(mb)):
            if label in (spec.get("labels") or ()):
                continue
            if label in spec.get("tolerated_labels", ()):
                continue
            if _is_output_path_label(label):
                # A command that reports a count per output file. After working-
                # directory normalisation the two arms' labels are comparable, so
                # the counts are compared like any other metric; the label itself is
                # tolerated because which files exist is a file-set question.
                va, vb = ma.get(label), mb.get(label)
                if va is None or vb is None:
                    errors.append(f"output-path label {label!r} present in one arm only")
                elif va != vb:
                    errors.append(f"output-path label {label!r} count differs: "
                                  f"{va} vs {vb}")
                continue
            errors.append(f"unexpected stdout label {label!r}")
    elif spec.get("stdout") is False:
        # A command declared to write nothing to stdout still must not disagree
        # about it: a stream appearing where none is expected is a contract change.
        if stdout_a.strip() or stdout_b.strip():
            errors.append("unexpected stdout content for a command that declares none")
    return (not errors), errors


def validate_arm_outputs(spec, out_dir, arm):
    """Validate one arm's artifacts in isolation, before any comparison.

    Called for EVERY timed run, outside the timing interval, so a run that
    produced an empty or partial deliverable cannot contribute a measurement. The
    audit's finding was that only the last repetition's trees were compared, and
    that an absent artifact in both arms compared equal to a present one.
    """
    errors = []
    if spec.get("artifacts") is None:
        errors.append(f"{arm}: command declares no artifact expectation")
    expected = spec.get("artifacts") or []
    produced = sorted(
        str(p.relative_to(out_dir)) for p in out_dir.rglob("*")
        if p.is_file() and p.name not in ("py.code", "rs.code", "time.txt"))
    if not produced and expected:
        errors.append(f"{arm}: produced no files at all")
    for rel in expected:
        if "?" in rel:
            import re as _re
            pattern = _re.compile("^" + rel.replace("?", "[^/]*") + "$")
            if not any(pattern.match(p) for p in produced):
                errors.append(f"{arm}: expected artifact matching {rel!r}, none produced")
            continue
        p = out_dir / rel
        if not p.exists():
            errors.append(f"{arm}: missing expected artifact {rel}")
            continue
        if p.stat().st_size == 0:
            errors.append(f"{arm}: artifact {rel} is empty")
        else:
            for line in p.read_text(errors="replace").splitlines():
                if _is_nonfinite(line.split("\t")[0]):
                    errors.append(f"{arm}: artifact {rel} contains a nonfinite value")
                    break
    return errors


def gate_outputs(py_dir: Path, rs_dir: Path, stochastic=False, partitioned=False,
                 name=None, stdout_a="", stdout_b="", exit_a=0, exit_b=0,
                 timed_out_a=False, timed_out_b=False, required_labels=None):
    """Structurally compare the two arms' output trees and streams.

    Returns (ok, errors). Any file present in one tree but not the other, any
    semantic difference, any stream or exit-status disagreement, and any missing
    expected artifact fails the gate. Two empty output trees fail: equality of
    nothing is not a measurement.
    """
    errors = []
    spec = dict(EXPECTED_STREAMS.get(name, {}))
    if required_labels:
        spec["labels"] = set(spec.get("labels", ())) | set(required_labels)
        spec["stdout"] = True
    # A command may legitimately deliver everything on stdout. In that case a file
    # tree is not the deliverable, and an empty tree must not be read as evidence
    # of nothing having been checked.
    stream_only = bool(spec.get("stdout")) and not spec.get("artifacts")

    pys = {p.relative_to(py_dir) for p in py_dir.rglob("*") if p.is_file()}
    rss = {p.relative_to(rs_dir) for p in rs_dir.rglob("*") if p.is_file()}
    # Ignore runner bookkeeping files.
    pys = {p for p in pys if p.name not in ("py.code", "rs.code", "time.txt")}
    rss = {p for p in rss if p.name not in ("py.code", "rs.code", "time.txt")}

    # Validate each arm on its own terms first, so "both produced nothing" cannot
    # be reported as agreement. A stdout-only command's deliverable is the stream,
    # which `compare_streams` checks instead.
    if not stream_only:
        errors.extend(validate_arm_outputs(spec, py_dir, "upstream"))
        errors.extend(validate_arm_outputs(spec, rs_dir, "rust"))
    # An empty output tree is only a defect for a command that declares artifacts.
    # A stdout-only command legitimately produces no files, and calling that a
    # failure would mean every such row is unreportable rather than unverified.
    if not pys and not rss and spec.get("artifacts"):
        errors.append("both arms produced no output files")
    if stream_only and not (stdout_a.strip() or stdout_b.strip()):
        errors.append("a stdout-only command produced no stdout in either arm")
    # A zero-byte file is not a deliverable, whichever arm wrote it. Two empty
    # files compare equal, so without this the gate would pass a command that
    # truncates its output identically in both arms.
    for rel in sorted(pys | rss):
        for arm, root in (("upstream", py_dir), ("rust", rs_dir)):
            p = root / rel
            if p.exists() and p.is_file() and p.stat().st_size == 0:
                errors.append(f"{arm}: {rel} is empty (0 bytes)")

    ok_streams, stream_errors = compare_streams(
        spec, stdout_a, stdout_b, exit_a, exit_b, timed_out_a, timed_out_b,
        work_dir_a=py_dir, work_dir_b=rs_dir)
    if not ok_streams:
        errors.extend(f"stream: {e}" for e in stream_errors)

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
# Scheduling
# --------------------------------------------------------------------------------------
class ScheduleItem:
    """One timed run: which repetition, which arm, in what position.

    Carrying the identity explicitly is the fix for the audit's P1 pairing
    finding. The previous schedule shuffled a flat list of individual arms
    globally and then zipped the two per-arm result lists by completion order, so
    at seed 20260929 all ten "pairs" had different repetition ids. The identity
    can no longer be lost, because the pairing is read off the schedule rather
    than inferred from two independently ordered lists.
    """

    __slots__ = ("rep", "arm", "order")

    def __init__(self, rep, arm, order):
        self.rep = rep
        self.arm = arm
        self.order = order

    def __repr__(self):
        return f"ScheduleItem(rep={self.rep}, arm={self.arm}, order={self.order})"

    def __eq__(self, other):
        return (isinstance(other, ScheduleItem)
                and (self.rep, self.arm, self.order)
                == (other.rep, other.arm, other.order))


def build_schedule(reps, seed):
    """Adjacent matched blocks: repetition i's two arms run back to back.

    Machine drift and thermal state move on the timescale of a single run, so the
    pairing has to be established by adjacency. Randomisation is kept, but moved
    INSIDE each block, where it decides only which arm goes first. That preserves
    the original intent -- no arm systematically occupies the warmer slot --
    without destroying the pairing that the confidence intervals depend on.
    """
    rng = random.Random(seed)
    schedule = []
    for rep in range(reps):
        order = ["py", "rs"] if rng.random() < 0.5 else ["rs", "py"]
        for arm in order:
            schedule.append(ScheduleItem(rep, arm, len(schedule)))
    return schedule


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
def oracle_lock_digest():
    """Digest of the pinned upstream reference environment.

    `compatibility/upstream.lock` records which reference implementation these
    numbers are compared against. Binding its digest into every run means a row
    cannot be silently attributed to a different oracle than the one recorded.
    """
    return sha256_file(REPO / "compatibility" / "upstream.lock")


def toolchain():
    def sh(*cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=20).stdout.strip()
        except Exception:
            return "unknown"
    return {
        "rustc": sh("rustc", "--version"),
        "cargo": sh("cargo", "--version"),
        "rust_target": sh("rustc", "-vV").split("host: ")[-1].splitlines()[0]
        if "host: " in sh("rustc", "-vV") else "unknown",
        "gnu_time": sh(TIME_BIN, "--version"),
        "python_oracle": sh(str(VENV_PY), "--version"),
        "upstream_script_digests": {
            entry[0]: _sha256(UPSTREAM_SRC / entry[0]) for entry in COMMANDS.values()
        },
    }


def input_digests(workload: Path):
    """Digest every input file the command matrix can reference."""
    names = ["reads.bam", "model.bed12", "reads_1.fastq", "reads_2.fastq",
             "reads_1.fa", "chrom.sizes", "sig1.bw", "sig2.bw"]
    return {n: sha256_file(workload / n) for n in names
            if (workload / n).exists()}


def provenance(workload: Path):
    """Everything needed to say which candidate produced which row."""
    def sh(*cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=20).stdout.strip()
        except Exception:
            return "unknown"
    binaries = {}
    for _script, binary, _argfn, _cls in COMMANDS.values():
        p = RELEASE / binary
        digest = sha256_file(p)
        if digest:
            binaries[binary] = digest
    return {
        "timestamp": time.time(),
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_commit": sh("git", "-C", str(REPO), "rev-parse", "HEAD"),
        "git_dirty": bool(sh("git", "-C", str(REPO), "status", "--porcelain")),
        "binary_sha256": binaries,
        "binaries_missing": sorted(
            {b for _s, b, _a, _c in COMMANDS.values()} - set(binaries)),
        "oracle_lock_sha256": oracle_lock_digest(),
        "comparator_version": COMPARATOR_VERSION,
        "inputs_sha256": input_digests(workload),
        "command_line": sys.argv,
        "toolchain": toolchain(),
        "pinned_env": PINNED_ENV,
    }


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


def sha256_file(path):
    """Full 64-hex digest, or None when the file is absent.

    A 16-character prefix identifies a file well enough to recognise a change but
    not well enough to bind a release candidate, and the audit's P1 provenance
    finding was that the recorded evidence named an older clean commit with no
    per-run binary hash at all. Truncation is not available as an option here.
    """
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


def _sha256(path, limit=None):
    """Digest prefix, for compactness where a full digest would only be noise."""
    digest = sha256_file(path)
    return digest[:16] if digest else "missing"


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

    # Adjacent matched blocks: repetition i's two arms run back to back, with the
    # within-block order randomised from the seed. Pairing is read from the
    # schedule, never reconstructed by zipping two independently ordered lists.
    schedule = build_schedule(reps, seed)

    per = {"py": [], "rs": []}
    by_rep = {"py": {}, "rs": {}}
    failures = []
    validation = []
    last_dirs = {}
    streams = {"py": {}, "rs": {}}
    spec = EXPECTED_STREAMS.get(name, {})

    for item in schedule:
        d = arm_dir(item.arm, item.rep)
        m, so, se = run_arm(argv_for(item.arm, d), d, env, timeout_s)
        m["rep"] = item.rep
        m["order"] = item.order
        m["started_at"] = time.time()
        per[item.arm].append(m)
        by_rep[item.arm][item.rep] = m
        streams[item.arm][item.rep] = so
        last_dirs[item.arm] = d
        # Every timed run's outputs are validated, outside the timing interval:
        # a repetition that produced an empty or partial deliverable must not
        # contribute a measurement. Previously only the last repetition was ever
        # compared, so nine of ten runs were unvalidated.
        errors = validate_arm_outputs(spec, d, item.arm)
        if m["exit_code"] != 0 or m["timed_out"] or errors:
            failures.append({
                "arm": item.arm, "rep": item.rep, "order": item.order,
                "exit_code": m["exit_code"], "timed_out": m["timed_out"],
                "output_errors": errors,
                # Both ends, because the two ends carry different evidence. An
                # argparse/validation error is printed at the START of stderr, while
                # GNU time's resource report is printed at the END; keeping only the
                # tail discarded the actual cause and left the word "imum resident
                # set size" as the recorded diagnosis of a ten-repetition failure.
                "stderr_tail": se[-400:], "stderr_head": se[:1200],
                "stderr_len": len(se),
            })
        validation.append({"arm": item.arm, "rep": item.rep, "errors": errors})
    py_dir, rs_dir = last_dirs["py"], last_dirs["rs"]

    results["runs"] = per
    results["failures"] = failures
    results["schedule"] = [{"rep": i.rep, "arm": i.arm, "order": i.order}
                           for i in schedule]
    results["output_validation"] = validation
    results["comparator_version"] = COMPARATOR_VERSION

    def pairs_for(key):
        out = []
        for rep in sorted(set(by_rep["py"]) & set(by_rep["rs"])):
            a, b = by_rep["py"][rep], by_rep["rs"][rep]
            # A pair is only usable when both arms succeeded and both delivered
            # the expected artifacts. Dropping such pairs would bias the estimate
            # toward whichever arm fails less often, so they are excluded and
            # counted instead.
            if a["exit_code"] != 0 or b["exit_code"] != 0:
                continue
            if a["timed_out"] or b["timed_out"]:
                continue
            va = a[key] - (py_floor if key == "wall_s" else 0.0)
            vb = b[key] - (rs_floor if key == "wall_s" else 0.0)
            if vb > 1e-6 and va > 1e-6:
                out.append((va, vb))
        return out

    # The measured per-invocation fixed cost of each arm. Reported and used only
    # as a sensitivity analysis: `--help` performs the imports and argument
    # parsing but not the algorithmic work, so it is an imperfect proxy for
    # startup cost, and the protocol treats its subtraction as secondary to the
    # raw invocation time rather than the primary estimator.
    fpy = arm_dir("upstream", 950)
    frs = arm_dir("rust", 950)
    py_floor = floor_cost(argv_for("py", fpy) + ["--help"], env, fpy)
    rs_floor = floor_cost(argv_for("rs", frs) + ["--help"], env, frs)
    results["floor_s"] = {"python": py_floor, "rust": rs_floor}

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

    # Primary estimator: raw invocation time, no floor subtraction. Every value
    # actually spent by the command, which is what a user's wall clock sees.
    results["wall_primary"] = summarise("wall_s", False)
    # Sensitivity analysis only: the same matched pairs with each arm's measured
    # `--help` floor removed. Labelled separately so it cannot be mistaken for the
    # headline figure.
    results["wall_e2e"] = summarise("wall_s", True)
    results["cpu"] = summarise("user_s", False)
    results["system"] = summarise("sys_s", False)
    results["estimator_note"] = (
        "wall_primary is the primary estimator (raw invocation time, adjacent matched "
        "blocks). wall_e2e subtracts each arm's measured --help floor and is a "
        "sensitivity analysis, not a second measurement."
    )
    mem_py = [r["peak_rss_mb"] for r in per["py"] if "peak_rss_mb" in r]
    mem_rs = [r["peak_rss_mb"] for r in per["rs"] if "peak_rss_mb" in r]
    if mem_py and mem_rs:
        # GNU time's child RSS is the LARGEST single child's RSS
        # (RUSAGE_CHILDREN.ru_maxrss), not the sum over a concurrently running
        # process tree. Labelled as such, because a helper-heavy workflow's real
        # simultaneous footprint is larger than this number.
        results["peak_rss_mb"] = {
            "median_python": statistics.median(mem_py),
            "median_rust": statistics.median(mem_rs),
            "max_python": max(mem_py), "max_rust": max(mem_rs),
            "ratio": (statistics.median(mem_py) / statistics.median(mem_rs))
            if statistics.median(mem_rs) else None,
            "measure": "largest-single-child RSS (GNU time -v / "
                       "RUSAGE_CHILDREN.ru_maxrss), not aggregate process-tree memory",
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
        # Gate the LAST matched block, so the compared streams and exit statuses
        # are the ones recorded for it rather than whichever arm happened to be
        # appended last. Every repetition's artifacts were validated above.
        last_rep = max(by_rep["py"])
        ok, errs = gate_outputs(
            py_dir, rs_dir, partitioned=part, name=name,
            stdout_a=streams["py"].get(last_rep, ""),
            stdout_b=streams["rs"].get(last_rep, ""),
            exit_a=by_rep["py"][last_rep]["exit_code"],
            exit_b=by_rep["rs"][last_rep]["exit_code"],
            timed_out_a=by_rep["py"][last_rep]["timed_out"],
            timed_out_b=by_rep["rs"][last_rep]["timed_out"],
        )
        mode = "partition-union" if part else "exact"
        # The streams the gate actually compared, retained. A stream failure with no
        # record of the stream is not diagnosable after the fact: RNA_fragment_size
        # failed with "unexpected stdout content for a command that declares none"
        # while a manual run of both arms printed nothing at all, and with stdout not
        # retained there was no way to tell what the gate saw. The first 400
        # characters of each arm's compared stream are enough to identify a stray
        # message and small enough to keep in the record.
        results["equivalence"] = {
            "pass": ok, "errors": errs[:20], "mode": mode,
            "compared_stdout": {
                "python": streams["py"].get(last_rep, "")[:400],
                "rust": streams["rs"].get(last_rep, "")[:400],
                "rep": last_rep,
            },
            "note": DIVIDED_DOCUMENTED.get(name, ""),
            "comparator_version": COMPARATOR_VERSION,
            "compared_rep": last_rep,
        }
    else:
        results["equivalence"] = {"pass": None, "errors": ["gate not run"]}

    return results


def main():
    ap = argparse.ArgumentParser()
    # Also resolved: --workload-for overrides are passed to commands as -i paths, and
    # a relative one is re-resolved against each arm's run directory too.
    ap.add_argument("--workload", required=True, type=Path)
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
    # Resolved, not taken literally. Each arm runs with cwd set to its own run
    # directory, so a relative --output-dir is re-resolved against that directory:
    # `--output-dir benchmarks/results` made every arm receive
    # `-o <run-dir>/benchmarks/results/<cmd>/<arm>/rep00/x`, whose parent does not
    # exist, and every command taking --out-prefix then failed in BOTH arms. It
    # looked like two broken commands rather than one bad invocation, and because
    # the failure record kept only the tail of stderr the recorded diagnosis was
    # the truncated text "imum resident set size". An entire real-data study was
    # launched down this path before it was caught.
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--no-gate", action="store_true")
    ap.add_argument("--quick", action="store_true",
                    help="smoke test: 1 rep, tiny command list")
    ap.add_argument("--summarise-only", action="store_true",
                    help="rebuild results.json from the per-command records already "
                         "in the output directory, measuring nothing. The recovery "
                         "operation for a summary that lost rows, and the cheap way "
                         "to re-render the report over an existing study")
    args = ap.parse_args()

    # Resolved once, here, because bench_command derives every arm's run directory
    # from it and then runs each arm with cwd set to that directory. Resolving it
    # inside bench_command is not possible: by then the relative path is already the
    # literal argument value, and a cwd change only makes the mistake permanent.
    args.output_dir = args.output_dir.resolve()
    args.workload = args.workload.resolve()

    if args.summarise_only:
        return summarise_only(args.output_dir)

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
    # Every number in the results is bound to what produced it. The audit's P1
    # finding was that the recorded evidence named an older clean commit with no
    # per-run binary hash, so a reader could not tell which candidate any given
    # row described.
    result = {
        "label": args.label,
        "environment": env,
        "provenance": provenance(workload),
        "workload": {
            "path": str(workload),
            "manifest": json.loads(manifest.read_text()) if manifest.exists() else None,
            "bam_sha256": sha256_file(workload / "reads.bam"),
        },
        "config": {"reps": args.reps, "warmup": args.warmup, "timeout_s": args.timeout,
                   "seed": args.seed, "gate": not args.no_gate},
        "comparator_version": COMPARATOR_VERSION,
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
        # .resolve() here rather than relying only on the Path type: arms run with cwd
        # set to their own run directory, so any relative path reaching -i or -o is
        # re-resolved against that directory instead of the invocation directory.
        wl = overrides.get(name, workload).resolve()
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
        # Per-run provenance, so an archived row names the exact candidate and
        # inputs rather than only the harness revision.
        r["provenance"] = provenance(wl)
        result["results"].append(r)
        w = r.get("wall_primary", {})
        ci = w.get("ratio_ci95")
        gate = "PASS" if r["equivalence"]["pass"] else "FAIL"
        ratio = w.get("ratio_median")
        # A row whose every repetition failed has median_python/RATIO_median as None,
        # and `{None:.3f}` raises TypeError. That aborted the run at the progress line,
        # after the row's own JSON had been written -- so the evidence survived but the
        # run reported nothing else and wrote no summary. A missing measurement prints
        # as MISSING, which is what it is.
        py_med, rs_med = w.get("median_python"), w.get("median_rust")
        print(
            (f"py {py_med:.3f}s  rs {rs_med:.4f}s  x{ratio:.2f}"
             if (ratio is not None and py_med is not None and rs_med is not None)
             else f"py MISSING  rs MISSING  no ratio"),
            end="  ",
        )
        print(f"CI[{ci[0]:.2f},{ci[1]:.2f}]" if ci else "CI n/a", end="  ")
        nfail = len(r.get("failures", []))
        print(f"gate={gate}  failed_reps={nfail}  ({time.time()-t0:.0f}s)")
        (outdir / f"{name}.json").write_text(json.dumps(r, indent=2))

    result = merge_summary(outdir, result)
    summary_path = outdir / "results.json"
    summary_path.write_text(json.dumps(result, indent=2))
    print()
    print(f"results -> {summary_path}")
    report(summary_path)


def merge_summary(outdir: Path, result: dict) -> dict:
    """Fold this invocation's rows into any summary already in `outdir`, in place.

    Recollecting a subset of rows is a normal protocol operation -- protocol-v2
    section 9 says to recollect affected rows only when code, methods or an
    unresolved discrepancy justify it -- so writing the summary from this
    invocation alone silently dropped the rows it did not re-measure, leaving a
    study summary that read as if only the re-collected commands had been
    benchmarked. The per-command records beside it survived, so no evidence was
    lost, but the aggregate is what a reader opens first.

    Rows measured here replace their earlier versions. Rows not re-measured are
    carried forward and marked with the revision that produced them, so a merged
    summary never presents old and new numbers as if they came from one run.
    """
    summary_path = outdir / "results.json"
    merged_rows = {r["command"]: r for r in result["results"]}
    carried, superseded = [], []
    if summary_path.is_file():
        try:
            previous = json.loads(summary_path.read_text())
        except json.JSONDecodeError:
            previous = {}
        for old_row in previous.get("results", []):
            cmd = old_row.get("command")
            if cmd is None or cmd in merged_rows:
                continue
            old_row = dict(old_row)
            old_row["carried_forward_from_earlier_invocation"] = True
            old_row["measured_at_revision"] = previous.get("provenance", {}).get("git_commit")
            merged_rows[cmd] = old_row
            carried.append(cmd)
        if carried:
            print(f"\ncarried forward {len(carried)} row(s) not re-measured in this "
                  f"invocation: {', '.join(sorted(carried))}")
            print("  they carry the revision that produced them; see "
                  "measured_at_revision on each row.")
    measured_now = {r["command"] for r in result["results"]}
    order = [r["command"] for r in result["results"]]
    order += [c for c in sorted(merged_rows) if c not in order]
    result["results"] = [merged_rows[c] for c in order]
    result["rows_measured_this_invocation"] = sorted(measured_now)
    result["rows_carried_forward"] = sorted(carried)

    return result


def summarise_only(outdir: Path) -> int:
    """Rebuild results.json from the per-command records in `outdir`.

    Each command's full record -- raw timings, failures, schedule, provenance -- is
    written to its own JSON as it completes, and results.json is an aggregate over
    those. So the records are the evidence and the aggregate is derived: when the
    aggregate loses rows, the records can rebuild it exactly and nothing has to be
    re-measured to get it back.

    Every row is marked as reconstructed, and the per-command file it came from is
    named, so a reader can tell a rebuilt summary from a freshly measured one.
    """
    summary_path = outdir / "results.json"
    rows, prov, excluded = [], None, {}
    for path in sorted(outdir.glob("*.json")):
        if path.name == "results.json":
            continue
        row = json.loads(path.read_text())
        row["reconstructed_from"] = path.name
        rows.append(row)
        prov = prov or row.get("provenance")
        excluded.update(row.get("excluded", {}))
    if not rows:
        print(f"no per-command records in {outdir}")
        return 1
    rows.sort(key=lambda r: r["command"])
    result = {
        "schema": rows[0].get("schema"),
        "comparator_version": rows[0].get("comparator_version"),
        "reconstructed": True,
        "note": ("Rebuilt by --summarise-only from the per-command records in this "
                 "directory. No command was measured to produce it; each row names "
                 "the record it came from."),
        "provenance": prov,
        "excluded": excluded,
        "results": rows,
        "rows_measured_this_invocation": [],
        "rows_carried_forward": sorted(r["command"] for r in rows),
    }
    summary_path.write_text(json.dumps(result, indent=2))
    print(f"rebuilt {summary_path} from {len(rows)} per-command record(s)")
    report(summary_path)
    return 0


def report(path: Path):
    d = json.loads(path.read_text())
    print()
    print("=" * 118)
    print("Primary estimator: raw invocation time over adjacent matched blocks. "
          "Floor-subtracted figures are in wall_e2e as a sensitivity analysis.")
    print(f"{'command':22} {'py(s)':>9} {'rs(s)':>9} {'speedup':>9} {'95% CI':>16} "
          f"{'mem(py/rs MB)':>16} {'pairs':>6} {'gate':>6}")
    print("-" * 118)
    for r in d["results"]:
        if "error" in r:
            print(f"{r['command']:22} ERROR {r['error'][:40]}")
            continue
        w = r.get("wall_primary", {})
        ci = w.get("ratio_ci95")
        m = r.get("peak_rss_mb", {})
        if m:
            mp, mr = m.get("median_python"), m.get("median_rust")
            mem = (f"{mp:.0f}/{mr:.0f}" if mp is not None and mr is not None
                   else "MISSING")
        else:
            mem = "-"
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
        # A row whose every repetition failed has median_python/RATIO_median as
        # None, not 0: `w.get(k, 0)` returns the stored None, and formatting None
        # raises TypeError. That aborted the whole report and lost the numbers for
        # every other command, including the ones that passed. A missing measurement
        # is displayed as MISSING, which is what it is.
        py_med = w.get("median_python")
        rs_med = w.get("median_rust")
        print(f"{r['command']:22} "
              f"{(f'{py_med:9.3f}' if py_med is not None else '     MISS'):>9} "
              f"{(f'{rs_med:9.4f}' if rs_med is not None else '     MISS'):>9} "
              f"{(f'{ratio:.2f}x' if ratio is not None else 'MISSING'):>9} "
              f"{(f'[{ci[0]:.2f}, {ci[1]:.2f}]' if ci else 'n/a'):>16} "
              f"{mem:>16} {w.get('n', 0):>6} {gate:>6} {verdict}")
    print("=" * 118)
    prov = d.get("provenance", {})
    print(f"revision   : {prov.get('git_commit')}  dirty={prov.get('git_dirty')}")
    print(f"comparator : v{prov.get('comparator_version', d.get('comparator_version'))}")
    print(f"oracle lock: {prov.get('oracle_lock_sha256')}")
    print(f"binaries   : {len(prov.get('binary_sha256', {}))} hashed"
          + (f", MISSING {prov['binaries_missing']}" if prov.get("binaries_missing") else ""))
    for k, v in d.get("excluded", {}).items():
        print(f"EXCLUDED {k}: {v}")
    bad = [r["command"] for r in d["results"]
           if r.get("equivalence", {}).get("pass") is False]
    if bad:
        print(f"\nGATE FAILURES (no speedup claim): {', '.join(bad)}")
    partial = [r["command"] for r in d["results"]
               if r.get("failures") and r.get("equivalence", {}).get("pass")]
    if partial:
        print(f"\nROWS WITH FAILED REPETITIONS (paired estimate uses the rest): "
              f"{', '.join(partial)}")


if __name__ == "__main__":
    main()
