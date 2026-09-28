#!/usr/bin/env python3
"""Synthetic sweep: byte-for-byte upstream-vs-Rust comparison of all 33 commands.

For every (seed, size) workload produced by verification/synthetic_data.py,
each case below runs the upstream script (oracle/upstream-src/scripts) and
the Rust binary (target/release) in SEPARATE temporary run directories used
as each side's CWD, then compares:

* the exit status;
* stdout and stderr, after replacing each side's own run-directory path with
  ``<RUN>`` and any ``YYYY-MM-DD HH:MM:SS`` timestamp with ``<TS>``;
* every file either side wrote, byte-for-byte after the same run-directory
  normalisation -- never numerically.  BAM outputs are compared as decoded
  header text + SAM record text via pysam (container/BGZF bytes are an
  encoder detail); ``.gz`` outputs are compared decompressed (gzip stamps a
  wall-clock mtime); ``.bai`` indexes only need to exist on both sides.

Genuinely non-comparable artefacts are skipped with the reason recorded in
the case (see SKIP_FILES and the per-case ``notes``).

Usage:
    oracle/venv/bin/python3 verification/synthetic_sweep.py [--seeds 1,2,3] [--sizes 150,1500]
        [--only REGEX] [--jobs N] [--keep DIR]
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import dataclasses
import fnmatch
import gzip
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
ORACLE_PYTHON = REPO_ROOT / "oracle" / "venv" / "bin" / "python3"
ORACLE_SCRIPTS = REPO_ROOT / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHONPATH = str(REPO_ROOT / "oracle" / "upstream-src" / "src")
RUST_BIN_DIR = REPO_ROOT / "target" / "release"
DATA_GEN = REPO_ROOT / "verification" / "synthetic_data.py"
COMMITTED_MOCK_HTSEQ = REPO_ROOT / "verification" / "fixtures" / "mock_htseq_count.sh"

TIMESTAMP_RE = re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(?:[.,]\d+)?")
# Timestamp-bearing log prefixes, stripped from stderr lines before comparing:
# Python logging's "%(asctime)s [%(levelname)s] " (DIV-0019) and the
# printlog() "@ <timestamp>: " of geneBody_coverage.py/FPKM-UQ.py.
LOG_PREFIX_RE = re.compile(rb"^(?:<TS> \[[A-Z]+\] |@ <TS>: )", re.M)
# Interpreter-level Python warnings (e.g. numpy's "RuntimeWarning: invalid
# value encountered in sqrt") name the installed module path and source line
# of the upstream package; they are not program output the port can mirror.
PY_WARNING_RE = re.compile(rb"^\S+\.py:\d+: \w*Warning: .*\n(?:  .*\n)?", re.M)

# Files never compared, with the reason.
SKIP_FILES = {
    "log.txt": "DIV-0022: upstream printlog() appends timestamped lines to ./log.txt; not replicated",
    "All_reads_uniqID.txt": "DIV-0018: sc_bamStat CWD pollution, not replicated",
    "confident_reads_uniqID.txt": "DIV-0018: sc_bamStat CWD pollution, not replicated",
}

STRAND_PE = "1++,1--,2+-,2-+"
STRAND_SE = "++,--"


@dataclasses.dataclass
class Case:
    name: str
    script: str  # upstream script file name
    binary: str  # Rust binary name
    args: list[str]  # "{D}" = data dir, "{O}" = run dir
    # "records" (default) compares decoded BAM records per file; "multiset"
    # compares the union of records across all output BAMs (divide_bam,
    # DIV-0017: RNG-driven assignment, content must still be identical).
    bam_mode: str = "records"
    py_exit: int | None = None  # per-side exit expectations (None: must merely match)
    rust_exit: int | None = None
    compare_stderr: bool = True
    stdout_mask: str | None = None  # regex replaced on both stdouts (RNG-dependent numbers only)
    skip_globs: tuple[str, ...] = ()  # extra per-case file skips
    notes: str = ""


def C(name, script, binary, *args, **kw) -> Case:
    return Case(name, script, binary, list(args), **kw)


def build_cases() -> list[Case]:
    cs: list[Case] = []
    add = cs.append
    # ---------------------------------------------------------------- bam_stat
    for bam in ("pe", "se", "sc"):
        add(C(f"bam_stat_{bam}", "bam_stat.py", "bam_stat", "-i", f"{{D}}/{bam}.bam"))
    add(C("bam_stat_pe_q0", "bam_stat.py", "bam_stat", "-i", "{D}/pe.bam", "-q", "0"))
    add(C("bam_stat_pe_q255", "bam_stat.py", "bam_stat", "-i", "{D}/pe.bam", "-q", "255"))
    # ------------------------------------------------------------------ bam2fq
    add(C("bam2fq_pe", "bam2fq.py", "bam2fq", "-i", "{D}/pe.bam", "-o", "{O}/out"))
    add(C("bam2fq_pe_single", "bam2fq.py", "bam2fq", "-i", "{D}/pe.bam", "-o", "{O}/out", "-s"))
    add(C("bam2fq_se_single", "bam2fq.py", "bam2fq", "-i", "{D}/se.bam", "-o", "{O}/out", "-s"))
    add(C("bam2fq_pe_compress", "bam2fq.py", "bam2fq", "-i", "{D}/pe.bam", "-o", "{O}/out", "-c"))
    # ----------------------------------------------------------------- bam2wig
    # wigToBigWig is absent here, so the final "Failed to call" line is part
    # of the compared stderr on both sides.
    base = ["-s", "{D}/chrom.sizes", "-o", "{O}/out"]
    add(C("bam2wig_pe", "bam2wig.py", "bam2wig", "-i", "{D}/pe.bam", *base))
    add(C("bam2wig_pe_wigsum", "bam2wig.py", "bam2wig", "-i", "{D}/pe.bam", *base, "-t", "1000000"))
    add(C("bam2wig_pe_skipmulti", "bam2wig.py", "bam2wig", "-i", "{D}/pe.bam", *base, "-u", "-q", "20"))
    add(C("bam2wig_pe_strand", "bam2wig.py", "bam2wig", "-i", "{D}/pe.bam", *base, "-d", STRAND_PE, "-t", "5000000"))
    add(C("bam2wig_se_strand", "bam2wig.py", "bam2wig", "-i", "{D}/se.bam", *base, "-d", STRAND_SE, "-u"))
    # -------------------------------------------------------- clipping_profile
    for lay, bam in (("PE", "pe"), ("SE", "se")):
        add(C(f"clipping_profile_{bam}", "clipping_profile.py", "clipping_profile",
              "-i", f"{{D}}/{bam}.bam", "-o", "{O}/out", "-s", lay, "--skip-plot"))
    add(C("clipping_profile_pe_q0", "clipping_profile.py", "clipping_profile",
          "-i", "{D}/pe.bam", "-o", "{O}/out", "-s", "PE", "-q", "0", "--skip-plot"))
    add(C("clipping_profile_se_as_pe", "clipping_profile.py", "clipping_profile",
          "-i", "{D}/se.bam", "-o", "{O}/out", "-s", "PE", "--skip-plot"))
    # ------------------------------------------------ deletion/mismatch/insert
    for tool in ("deletion_profile", "mismatch_profile"):
        add(C(f"{tool}_pe", f"{tool}.py", tool, "-i", "{D}/pe.bam", "-l", "50", "-o", "{O}/out", "--skip-plot"))
        add(C(f"{tool}_se_n20_q0", f"{tool}.py", tool, "-i", "{D}/se.bam", "-l", "50", "-o", "{O}/out",
              "-n", "20", "-q", "0", "--skip-plot"))
        add(C(f"{tool}_pe_l40", f"{tool}.py", tool, "-i", "{D}/pe.bam", "-l", "40", "-o", "{O}/out", "--skip-plot"))
    for lay, bam in (("PE", "pe"), ("SE", "se")):
        add(C(f"insertion_profile_{bam}", "insertion_profile.py", "insertion_profile",
              "-i", f"{{D}}/{bam}.bam", "-o", "{O}/out", "-s", lay, "--skip-plot"))
    add(C("insertion_profile_pe_q0", "insertion_profile.py", "insertion_profile",
          "-i", "{D}/pe.bam", "-o", "{O}/out", "-s", "PE", "-q", "0", "--skip-plot"))
    # --------------------------------------------------------------- divide_bam
    dv = "DIV-0017: record->file assignment is RNG-driven; compared as a multiset across outputs"
    add(C("divide_bam_pe", "divide_bam.py", "divide_bam", "-i", "{D}/pe.bam", "-n", "3", "-o", "{O}/out",
          bam_mode="multiset", stdout_mask=r"(?<=\t)\d+$", notes=dv))
    add(C("divide_bam_pe_skipunmap_index", "divide_bam.py", "divide_bam", "-i", "{D}/pe.bam", "-n", "2",
          "-o", "{O}/out", "-s", "--seed", "7", "--index", bam_mode="multiset", stdout_mask=r"(?<=\t)\d+$",
          notes=dv))
    # --------------------------------------------------------------- FPKM_count
    fb = ["-o", "{O}/out", "-r", "{D}/model.bed12"]
    add(C("fpkm_count_pe", "FPKM_count.py", "FPKM_count", "-i", "{D}/pe_placed.bam", *fb))
    add(C("fpkm_count_pe_strand_u_e", "FPKM_count.py", "FPKM_count", "-i", "{D}/pe.bam", *fb,
          "-d", STRAND_PE, "-u", "-e"))
    add(C("fpkm_count_pe_q0_s05", "FPKM_count.py", "FPKM_count", "-i", "{D}/pe_placed.bam", *fb, "-q", "0",
          "-s", "0.5"))
    add(C("fpkm_count_pe_s0_u", "FPKM_count.py", "FPKM_count", "-i", "{D}/pe.bam", *fb, "-s", "0", "-u"))
    add(C("fpkm_count_se_strand", "FPKM_count.py", "FPKM_count", "-i", "{D}/se_placed.bam", *fb, "-d", STRAND_SE))
    add(C("fpkm_count_unplaced_unmapped", "FPKM_count.py", "FPKM_count", "-i", "{D}/pe.bam", *fb,
          py_exit=1, rust_exit=0, compare_stderr=False, skip_globs=("out.FPKM.xls",),
          notes="DIV-0023: upstream crashes (getrname(-1) -> None.upper()) on coordinate-less unmapped reads "
                "without -u; the port skips them"))
    # ------------------------------------------------------------------ FPKM-UQ
    ub = ["--bam", "{D}/pe.bam", "--gtf", "{D}/model.gtf", "--info", "{D}/genes.info.txt", "-o", "{O}/out"]
    add(C("fpkm_uq_workload", "FPKM-UQ.py", "FPKM_UQ", *ub, "--htseq-count", "{D}/mock_htseq_count.sh"))
    add(C("fpkm_uq_workload_log2", "FPKM-UQ.py", "FPKM_UQ", *ub, "--htseq-count", "{D}/mock_htseq_count.sh", "--log2"))
    add(C("fpkm_uq_committed_mock", "FPKM-UQ.py", "FPKM_UQ", *ub, "--htseq-count", str(COMMITTED_MOCK_HTSEQ)))
    add(C("fpkm_uq_print_command", "FPKM-UQ.py", "FPKM_UQ", *ub, "--htseq-count", "{D}/mock_htseq_count.sh",
          "--print-htseq-command"))
    # ------------------------------------------------------- geneBody_coverage2
    for bw in ("sig1", "sig2"):
        add(C(f"genebody_coverage2_{bw}", "geneBody_coverage2.py", "geneBody_coverage2",
              "-i", f"{{D}}/{bw}.bw", "-r", "{D}/model.bed12", "-o", "{O}/out", "--skip-plot"))
    add(C("genebody_coverage2_png", "geneBody_coverage2.py", "geneBody_coverage2",
          "-i", "{D}/sig1.bw", "-r", "{D}/model.bed12", "-o", "{O}/out", "-t", "png", "--skip-plot"))
    # -------------------------------------------------------- geneBody_coverage
    gb = ["-r", "{D}/model.bed12", "-o", "{O}/out", "--skip-plot"]
    add(C("genebody_coverage_pe", "geneBody_coverage.py", "geneBody_coverage", "-i", "{D}/pe.bam", *gb))
    add(C("genebody_coverage_multi", "geneBody_coverage.py", "geneBody_coverage",
          "-i", "{D}/pe.bam,{D}/se.bam,{D}/sc.bam", *gb))
    add(C("genebody_coverage_listfile_l300_png", "geneBody_coverage.py", "geneBody_coverage",
          "-i", "{D}/bams.txt", *gb, "-l", "300", "-f", "png"))
    # --------------------------------------------------------- infer_experiment
    for bam in ("pe", "se"):
        add(C(f"infer_experiment_{bam}", "infer_experiment.py", "infer_experiment",
              "-i", f"{{D}}/{bam}.bam", "-r", "{D}/model.bed12"))
    add(C("infer_experiment_pe_s50_q0", "infer_experiment.py", "infer_experiment",
          "-i", "{D}/pe.bam", "-r", "{D}/model.bed12", "-s", "50", "-q", "0"))
    # ----------------------------------------------------------- inner_distance
    ib = ["-o", "{O}/out", "-r", "{D}/model.bed12", "--skip-plot"]
    add(C("inner_distance_pe", "inner_distance.py", "inner_distance", "-i", "{D}/pe.bam", *ib))
    add(C("inner_distance_pe_bins_q0", "inner_distance.py", "inner_distance", "-i", "{D}/pe.bam", *ib,
          "-k", "150", "-l", "-100", "-u", "300", "-s", "7", "-q", "0"))
    add(C("inner_distance_se", "inner_distance.py", "inner_distance", "-i", "{D}/se.bam", *ib))
    # ------------------------------------------------------ junction_annotation
    jb = ["-r", "{D}/model.bed12", "-o", "{O}/out", "--skip-plot"]
    add(C("junction_annotation_pe", "junction_annotation.py", "junction_annotation", "-i", "{D}/pe.bam", *jb))
    add(C("junction_annotation_se_m30_q0", "junction_annotation.py", "junction_annotation",
          "-i", "{D}/se.bam", *jb, "-m", "30", "-q", "0"))
    add(C("junction_annotation_pe_skipbed_skipinteract", "junction_annotation.py", "junction_annotation",
          "-i", "{D}/pe.bam", *jb, "--skip-bed", "--skip-interact"))
    # ------------------------------------------------------ junction_saturation
    # sat_junction.bam: every read passing the filters carries the same
    # intron, so the random.shuffle() of the junction list is a no-op.
    sb = ["-o", "{O}/out", "-r", "{D}/model.bed12", "--skip-plot"]
    add(C("junction_saturation_default", "junction_saturation.py", "junction_saturation",
          "-i", "{D}/sat_junction.bam", *sb))
    add(C("junction_saturation_steps", "junction_saturation.py", "junction_saturation",
          "-i", "{D}/sat_junction.bam", *sb, "-l", "10", "-u", "90", "-s", "20", "-v", "3", "-m", "100"))
    # -------------------------------------------------------- normalize_bigwig
    add(C("normalize_bigwig_bgr", "normalize_bigwig.py", "normalize_bigwig",
          "-i", "{D}/sig1.bw", "-o", "{O}/out.bgr"))
    add(C("normalize_bigwig_wig_t_c", "normalize_bigwig.py", "normalize_bigwig",
          "-i", "{D}/sig2.bw", "-o", "{O}/out.wig", "-f", "wig", "-t", "12345678", "-c", "777"))
    add(C("normalize_bigwig_refgene", "normalize_bigwig.py", "normalize_bigwig",
          "-i", "{D}/sig1.bw", "-o", "{O}/out.wig", "-f", "wig", "-r", "{D}/model.bed12"))
    add(C("normalize_bigwig_refgene_bgr_sig2", "normalize_bigwig.py", "normalize_bigwig",
          "-i", "{D}/sig2.bw", "-o", "{O}/out.bgr", "-r", "{D}/model.bed12", "-c", "5000"))
    # ----------------------------------------------------------- overlay_bigwig
    for act in ("Add", "Average", "Max", "Min", "Product", "Subtract", "geometricMean"):
        add(C(f"overlay_bigwig_{act}", "overlay_bigwig.py", "overlay_bigwig",
              "-i", "{D}/sig1.bw", "-j", "{D}/sig2.bw", "-a", act, "-o", "{O}/out.wig"))
    add(C("overlay_bigwig_Subtract_chunk", "overlay_bigwig.py", "overlay_bigwig",
          "-i", "{D}/sig2.bw", "-j", "{D}/sig1.bw", "-a", "Subtract", "-o", "{O}/out.wig", "-c", "1000"))
    add(C("overlay_bigwig_Division", "overlay_bigwig.py", "overlay_bigwig",
          "-i", "{D}/sig1.bw", "-j", "{D}/sig2.bw", "-a", "Division", "-o", "{O}/out.wig",
          compare_stderr=False, skip_globs=("out.wig",),
          notes="DIV-0015: upstream crashes with a Python traceback (ndarray.__div__) after opening out.wig; "
                "only the failure exit is compared"))
    # -------------------------------------------------------- read_distribution
    for bam in ("pe", "se", "sc"):
        add(C(f"read_distribution_{bam}", "read_distribution.py", "read_distribution",
              "-i", f"{{D}}/{bam}.bam", "-r", "{D}/model.bed12"))
    # ---------------------------------------------------------- read_duplication
    add(C("read_duplication_pe", "read_duplication.py", "read_duplication", "-i", "{D}/pe.bam", "-o", "{O}/out",
          "--skip-plot"))
    add(C("read_duplication_se_u5_q0", "read_duplication.py", "read_duplication", "-i", "{D}/se.bam",
          "-o", "{O}/out", "-u", "5", "-q", "0", "--skip-plot"))
    # ------------------------------------------------------------------ read_GC
    add(C("read_GC_pe", "read_GC.py", "read_GC", "-i", "{D}/pe.bam", "-o", "{O}/out", "--skip-plot"))
    add(C("read_GC_se_q0", "read_GC.py", "read_GC", "-i", "{D}/se.bam", "-o", "{O}/out", "-q", "0", "--skip-plot"))
    # ------------------------------------------------------------- read_hexamer
    add(C("read_hexamer_fa", "read_hexamer.py", "read_hexamer", "-i", "{D}/reads.fa"))
    add(C("read_hexamer_fq", "read_hexamer.py", "read_hexamer", "-i", "{D}/reads.fq"))
    add(C("read_hexamer_all_refs", "read_hexamer.py", "read_hexamer", "-i", "{D}/reads.fa,{D}/barcodes.fq",
          "-r", "{D}/genome.fa", "-g", "{D}/mrna.fa", "-o", "{O}/hex.txt"))
    add(C("read_hexamer_skip_missing", "read_hexamer.py", "read_hexamer", "-i", "{D}/reads.fq,{D}/nope.fa",
          "--skip-missing", "-g", "{D}/mrna.fa"))
    # ----------------------------------------------------------------- read_NVC
    add(C("read_NVC_pe", "read_NVC.py", "read_NVC", "-i", "{D}/pe.bam", "-o", "{O}/out", "--skip-plot"))
    add(C("read_NVC_se_x_q0", "read_NVC.py", "read_NVC", "-i", "{D}/se.bam", "-o", "{O}/out", "-x", "-q", "0",
          "--skip-plot"))
    # ------------------------------------------------------------- read_quality
    add(C("read_quality_pe", "read_quality.py", "read_quality", "-i", "{D}/pe.bam", "-o", "{O}/out", "--skip-plot"))
    add(C("read_quality_se_r1_q0", "read_quality.py", "read_quality", "-i", "{D}/se.bam", "-o", "{O}/out",
          "-r", "1", "-q", "0", "--skip-plot"))
    # -------------------------------------------------------- RNA_fragment_size
    add(C("rna_fragment_size_pe", "RNA_fragment_size.py", "RNA_fragment_size", "-i", "{D}/pe.bam",
          "-r", "{D}/model.bed12"))
    add(C("rna_fragment_size_pe_q0_n1_o", "RNA_fragment_size.py", "RNA_fragment_size", "-i", "{D}/pe.bam",
          "-r", "{D}/model.bed12", "-q", "0", "-n", "1", "-o", "{O}/frag.txt"))
    # ---------------------------------------------------------- RPKM_saturation
    # sat_rpkm.bam: every qualifying read yields the same exon-block key, so
    # the random.shuffle() of the block lists is a no-op.
    rb = ["-o", "{O}/out", "-r", "{D}/model.bed12", "--skip-plot"]
    add(C("rpkm_saturation_default", "RPKM_saturation.py", "RPKM_saturation", "-i", "{D}/sat_rpkm.bam", *rb))
    add(C("rpkm_saturation_strand_steps", "RPKM_saturation.py", "RPKM_saturation", "-i", "{D}/sat_rpkm.bam", *rb,
          "-d", STRAND_PE, "-l", "10", "-u", "80", "-s", "15", "-c", "5"))
    # -------------------------------------------------------------- sc_bamStat
    add(C("sc_bamstat_default", "sc_bamStat.py", "sc_bamStat", "-i", "{D}/sc.bam"))
    add(C("sc_bamstat_tags", "sc_bamStat.py", "sc_bamStat", "-i", "{D}/sc.bam", "--cb-tag", "CR",
          "--umi-tag", "UR", "--chrM-id", "chr10", "--verbose"))
    # ----------------------------------------------------------- sc_editMatrix
    add(C("sc_editmatrix_default", "sc_editMatrix.py", "sc_editMatrix", "-i", "{D}/sc.bam", "-o", "{O}/out",
          "--skip-heatmap"))
    add(C("sc_editmatrix_limit", "sc_editMatrix.py", "sc_editMatrix", "-i", "{D}/sc.bam", "-o", "{O}/out",
          "--skip-heatmap", "--limit", "37", "--verbose"))
    # ------------------------------------------------------------- sc_seqLogo
    logo = dict(skip_globs=("out.logo.*",),
                notes="DIV-0016: logo images are rendered natively (not logomaker/matplotlib) and are not "
                      "byte-comparable; count matrix, stderr and exit compared")
    add(C("sc_seqlogo_fa", "sc_seqLogo.py", "sc_seqLogo", "-i", "{D}/barcodes.fa", "-o", "{O}/out",
          "--iformat", "fa", **logo))
    add(C("sc_seqlogo_fq_gz_excludeN_n", "sc_seqLogo.py", "sc_seqLogo", "-i", "{D}/barcodes.fq.gz", "-o", "{O}/out",
          "--exclude-N", "-n", "25", "--oformat", "svg", "--step-size", "7", **logo))
    # ------------------------------------------------------------- sc_seqQual
    for inp in ("barcodes.fq", "reads.fq", "reads.fq.gz", "reads.fq.bz2"):
        add(C(f"sc_seqqual_{inp.replace('.', '_')}", "sc_seqQual.py", "sc_seqQual", "-i", f"{{D}}/{inp}",
              "-o", "{O}/out", "--skip-heatmap"))
    add(C("sc_seqqual_n", "sc_seqQual.py", "sc_seqQual", "-i", "{D}/reads.fq", "-o", "{O}/out", "-n", "17",
          "--skip-heatmap", "--verbose"))
    # --------------------------------------------------------------- split_bam
    add(C("split_bam_pe", "split_bam.py", "split_bam", "-i", "{D}/pe.bam", "-r", "{D}/model.bed12", "-o", "{O}/out"))
    add(C("split_bam_se_verbose", "split_bam.py", "split_bam", "-i", "{D}/se.bam", "-r", "{D}/model.bed12",
          "-o", "{O}/out", "--verbose"))
    # -------------------------------------------------------- split_paired_bam
    add(C("split_paired_bam_pe", "split_paired_bam.py", "split_paired_bam", "-i", "{D}/pe.bam", "-o", "{O}/out"))
    add(C("split_paired_bam_se_verbose", "split_paired_bam.py", "split_paired_bam", "-i", "{D}/se.bam",
          "-o", "{O}/out", "--verbose"))
    # --------------------------------------------------------------------- tin
    add(C("tin_pe", "tin.py", "tin", "-i", "{D}/pe.bam", "-r", "{D}/model.bed12", "-o", "{O}"))
    add(C("tin_multi_c2_n20_s", "tin.py", "tin", "-i", "{D}/pe.bam,{D}/se.bam", "-r", "{D}/model.bed12",
          "-c", "2", "-n", "20", "-s", "-o", "{O}"))
    add(C("tin_listfile_c1", "tin.py", "tin", "-i", "{D}/bams.txt", "-r", "{D}/model.bed12", "-c", "1",
          "-o", "{O}", "--verbose"))
    return cs


# --------------------------------------------------------------------- runner
def run(argv, cwd: Path, env_extra=None, timeout=600):
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    p = subprocess.Popen(argv, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         start_new_session=True)
    try:
        out, err = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        out, err = p.communicate()
        return None, out, err + b"\n[timeout]"
    return p.returncode, out, err


def norm(data: bytes, run_dir: Path) -> bytes:
    data = data.replace(str(run_dir).encode(), b"<RUN>")
    return TIMESTAMP_RE.sub("<TS>", data.decode("utf-8", "surrogateescape")).encode("utf-8", "surrogateescape")


def bam_records(path: Path, run_dir: Path) -> list[str]:
    import pysam

    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as fh:
        hdr = str(fh.header).replace(str(run_dir), "<RUN>")
        return ["@" + line for line in hdr.split("\n@") if line] + [r.to_string() for r in fh]


def list_files(root: Path) -> set[str]:
    return {str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()}


def first_diff(a: bytes, b: bytes, ctx=160) -> str:
    al, bl = a.split(b"\n"), b.split(b"\n")
    for i, (x, y) in enumerate(zip(al, bl)):
        if x != y:
            return f"line {i + 1}:\n      py  : {x[:ctx]!r}\n      rust: {y[:ctx]!r}"
    return f"line count py={len(al)} rust={len(bl)}; extra: {(al[len(bl):] or bl[len(al):])[:2]!r}"


def compare(case: Case, py, rs, py_dir: Path, rs_dir: Path, data_dir: Path) -> list[str]:
    probs = []
    (pc, po, pe), (rc, ro, re_) = py, rs
    if case.py_exit is not None or case.rust_exit is not None:
        if case.py_exit is not None and pc != case.py_exit:
            probs.append(f"python exit {pc} != expected {case.py_exit}: {pe[-400:]!r}")
        if case.rust_exit is not None and rc != case.rust_exit:
            probs.append(f"rust exit {rc} != expected {case.rust_exit}: {re_[-400:]!r}")
    elif pc != rc:
        probs.append(f"exit py={pc} rust={rc}\n      py stderr tail: {pe[-300:]!r}\n      rust stderr tail: {re_[-300:]!r}")
    npo, nro = norm(po, py_dir), norm(ro, rs_dir)
    # normalize both run_dir and data_dir paths
    npo = npo.replace(str(data_dir).encode(), b"<DATA>")
    nro = nro.replace(str(data_dir).encode(), b"<DATA>")
    if case.stdout_mask:
        npo, nro = (re.sub(case.stdout_mask.encode(), b"<MASKED>", x, flags=re.M) for x in (npo, nro))
    if npo != nro:
        probs.append("stdout " + first_diff(npo, nro))
    if case.compare_stderr:
        npe, nre = norm(pe, py_dir), norm(re_, rs_dir)
        # normalize both run_dir and data_dir paths in stderr
        npe = npe.replace(str(data_dir).encode(), b"<DATA>")
        nre = nre.replace(str(data_dir).encode(), b"<DATA>")
        npe = PY_WARNING_RE.sub(b"", LOG_PREFIX_RE.sub(b"", npe))
        nre = LOG_PREFIX_RE.sub(b"", nre)
        if npe != nre:
            probs.append("stderr " + first_diff(npe, nre))
    pf, rf = list_files(py_dir), list_files(rs_dir)
    skip = set(SKIP_FILES)

    def skipped(f):
        return Path(f).name in skip or any(fnmatch.fnmatch(f, g) for g in case.skip_globs)

    pf = {f for f in pf if not skipped(f)}
    rf = {f for f in rf if not skipped(f)}
    if pf != rf:
        probs.append(f"file sets differ: only-py={sorted(pf - rf)} only-rust={sorted(rf - pf)}")
    if case.bam_mode == "multiset":
        def body(d, files):
            return sorted(x for f in sorted(files) if f.endswith(".bam")
                          for x in bam_records(d / f, d) if not x.startswith("@"))

        pb, rb = body(py_dir, pf), body(rs_dir, rf)
        if pb != rb:
            probs.append(f"BAM record multiset differs: py={len(pb)} rust={len(rb)}")
        for f in sorted(pf & rf):
            if f.endswith(".bam"):
                ph = [x for x in bam_records(py_dir / f, py_dir) if x.startswith("@")]
                rh = [x for x in bam_records(rs_dir / f, rs_dir) if x.startswith("@")]
                if ph != rh:
                    probs.append(f"{f}: BAM header differs py={ph!r} rust={rh!r}")
    for f in sorted(pf & rf):
        a, b = py_dir / f, rs_dir / f
        if f.endswith(".bai"):
            if a.stat().st_size == 0 or b.stat().st_size == 0:
                probs.append(f"{f}: empty index")
            continue
        if f.endswith(".bam"):
            if case.bam_mode == "multiset":
                continue
            pr, rr = bam_records(a, py_dir), bam_records(b, rs_dir)
            if pr != rr:
                for i, (x, y) in enumerate(zip(pr, rr)):
                    if x != y:
                        probs.append(f"{f}: BAM record {i} differs\n      py  : {x[:200]}\n      rust: {y[:200]}")
                        break
                else:
                    probs.append(f"{f}: BAM record count py={len(pr)} rust={len(rr)}")
            continue
        da, db = a.read_bytes(), b.read_bytes()
        if f.endswith(".gz"):
            da, db = gzip.decompress(da), gzip.decompress(db)
        da, db = norm(da, py_dir), norm(db, rs_dir)
        if da != db:
            probs.append(f"{f}: " + first_diff(da, db))
    return probs


def run_case(case: Case, data_dir: Path, keep: Path | None, tag: str):
    py_dir = Path(tempfile.mkdtemp(prefix=f"sw_{case.name}_py_"))
    rs_dir = Path(tempfile.mkdtemp(prefix=f"sw_{case.name}_rs_"))
    try:
        def subst(d):
            return [a.replace("{D}", str(data_dir)).replace("{O}", str(d)) for a in case.args]

        py = run([str(ORACLE_PYTHON), str(ORACLE_SCRIPTS / case.script)] + subst(py_dir), py_dir,
                 {"PYTHONPATH": ORACLE_PYTHONPATH})
        rs = run([str(RUST_BIN_DIR / case.binary)] + subst(rs_dir), rs_dir)
        probs = compare(case, py, rs, py_dir, rs_dir, data_dir)
        if probs and keep is not None:
            dest = keep / tag / case.name
            if dest.exists():
                shutil.rmtree(dest)
            dest.mkdir(parents=True)
            shutil.copytree(py_dir, dest / "py")
            shutil.copytree(rs_dir, dest / "rust")
            for side, (_, o, e) in (("py", py), ("rust", rs)):
                (dest / f"{side}.stdout").write_bytes(o)
                (dest / f"{side}.stderr").write_bytes(e)
            (dest / "cmd.txt").write_text(" ".join(subst(Path("<RUN>"))) + "\n")
        return probs
    finally:
        shutil.rmtree(py_dir, ignore_errors=True)
        shutil.rmtree(rs_dir, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--sizes", default="150,1500")
    ap.add_argument("--only", default=None, help="regex on case names")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--keep", type=Path, default=None, help="copy failing runs here")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    cases = build_cases()
    if args.only:
        cases = [c for c in cases if re.search(args.only, c.name)]
    if args.list:
        for c in cases:
            print(c.name, c.notes)
        print(f"{len(cases)} cases over {len({c.script for c in cases})} commands")
        return 0
    if not RUST_BIN_DIR.is_dir():
        print("build first: cargo build --workspace --release", file=sys.stderr)
        return 2
    seeds = [int(s) for s in args.seeds.split(",")]
    sizes = [int(s) for s in args.sizes.split(",")]
    failures: dict[str, list[str]] = {}
    total = 0
    with tempfile.TemporaryDirectory(prefix="rseqc_sweep_data_") as tmp:
        for size in sizes:
            for seed in seeds:
                tag = f"size{size}_seed{seed}"
                ddir = Path(tmp) / tag
                code, _, err = run([str(ORACLE_PYTHON), str(DATA_GEN), "--seed", str(seed), "--size", str(size),
                                    "--output-dir", str(ddir)], REPO_ROOT)
                if code != 0:
                    print(f"data generation failed for {tag}:\n{err.decode()}", file=sys.stderr)
                    return 2
                with cf.ThreadPoolExecutor(args.jobs) as ex:
                    futs = {ex.submit(run_case, c, ddir, args.keep, tag): c for c in cases}
                    for fut in cf.as_completed(futs):
                        c = futs[fut]
                        total += 1
                        probs = fut.result()
                        if probs:
                            failures.setdefault(c.name, []).append(f"[{tag}] " + "\n    ".join(probs))
                            print(f"FAIL {tag} {c.name}", flush=True)
                n_fail = sum(1 for v in failures.values() for x in v if x.startswith(f"[{tag}]"))
                print(f"{tag}: {len(cases) - n_fail}/{len(cases)} passed", flush=True)
    print()
    for name in sorted(failures):
        print(f"=== {name}")
        for f in failures[name][:3]:
            print("  " + f)
    n_fail = sum(len(v) for v in failures.values())
    print(f"\n{total - n_fail}/{total} runs passed; {len(failures)} failing case(s) "
          f"over {len(cases)} cases x {len(seeds) * len(sizes)} workloads "
          f"({len({c.script for c in cases})} commands)")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
