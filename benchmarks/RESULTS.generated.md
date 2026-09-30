# Benchmark results — rseqc-rust vs upstream RSeQC

Protocol: `benchmarks/protocol.md` (frozen 2026-09-29). Run label: `windowed-all`. Reps: 10, seed 20260929.

> **Not publication-grade.** Shared, non-isolated hardware (AMD Ryzen 7 3700X 8-Core Processor, loadavg 0.86 at start, governor `schedutil`). See protocol.md §2.

## Environment

- **platform**: Linux-6.14.0-37-generic-x86_64-with-glibc2.39
- **cpu_count**: 16
- **processor**: AMD Ryzen 7 3700X 8-Core Processor
- **cpufreq_governor**: schedutil
- **filesystem**: ext2/ext3
- **rustc**: rustc 1.98.1 (48a229cea 2026-09-01)
- **python**: Python 3.13.15
- **numpy**: 2.5.3
- **pysam**: 0.24.1
- **git_commit**: ee7b6e266ecdd47e41f48ff003b8fd5c69f65d8f
- **git_dirty**: True
- **pinned env**: `{'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONHASHSEED': '0'}`

## 1. Per-command results (E2 end-to-end, equal deliverables)

Speedup = median(upstream) / median(port) on the E1 compute-only quantities (each arm's measured per-invocation fixed cost subtracted, so the number and the interval describe the same thing). The interval is a 95% paired block bootstrap over whole pairs. A claim requires the interval to lie entirely above 1.0; an interval spanning 1.0 is inconclusive and an interval below 1.0 is a loss. Section 2 shows the end-to-end figures.

### all measured commands

| command | workload | py net (s) | rs net (s) | speedup | 95% CI | verdict | CPU x | mem x | gate |
|---|---|---|---|---|---|---|---|---|---|
| normalize_bigwig | W-mid | 47.802 | 0.642 | 74.510 | [66.37, 78.93] | win | 69.536 | 7.875 | True |
| geneBody_coverage | W-gbc | 61.928 | 1.329 | 46.611 | [46.38, 46.85] | win | 47.977 | 2.936 | True |
| tin | W-tn | 148.619 | 9.106 | 16.322 | [16.26, 16.42] | win | 16.781 | 2.066 | True |
| read_NVC | W-mid | 14.720 | 1.325 | 11.106 | [10.99, 11.26] | win | 11.267 | 12.356 | True |
| inner_distance | W-mid | 7.975 | 0.730 | 10.926 | [10.84, 11.02] | win | 12.138 | 3.998 | True |
| RNA_fragment_size | W-mid | 54.892 | 5.364 | 10.233 | [10.17, 10.33] | win | 10.322 | 2.181 | True |
| sc_seqLogo | W-mid | 0.325 | 0.040 | 8.021 | [7.70, 8.32] | win | 25.333 | 15.385 | True |
| sc_seqQual | W-mid | 0.427 | 0.054 | 7.909 | [7.26, 8.26] | win | 17.300 | 15.670 | True |
| RPKM_saturation | W-mid | 29.788 | 4.128 | 7.216 | [7.09, 7.35] | win | 7.538 | 1.949 | True |
| geneBody_coverage2 | W-mid | 25.120 | 4.218 | 5.955 | [5.93, 6.03] | win | 16.513 | 12.098 | True |
| overlay_bigwig | W-mid | 0.470 | 0.108 | 4.338 | [4.25, 4.42] | win | 6.375 | 5.328 | True |
| bam_stat | W-mid | 1.801 | 0.426 | 4.227 | [4.20, 4.26] | win | 4.524 | 13.073 | True |
| FPKM_count | W-mid | 33.900 | 8.124 | 4.173 | [4.16, 4.45] | win | 4.215 | 1.577 | True |
| read_hexamer | W-mid | 0.562 | 0.158 | 3.559 | [3.50, 3.62] | win | 4.067 | 3.996 | True |
| clipping_profile | W-mid | 1.585 | 0.457 | 3.472 | [3.44, 3.49] | win | 3.716 | 12.934 | True |
| read_quality | W-mid | 10.037 | 2.912 | 3.447 | [3.42, 3.47] | win | 3.490 | 12.076 | True |
| insertion_profile | W-mid | 1.556 | 0.462 | 3.370 | [3.35, 3.40] | win | 3.578 | 12.920 | True |
| infer_experiment | W-mid | 0.612 | 0.182 | 3.364 | [3.26, 3.40] | win | 3.941 | 7.909 | True |
| read_distribution | W-mid | 4.703 | 1.455 | 3.234 | [3.22, 3.25] | win | 3.237 | 3.056 | True |
| deletion_profile | W-mid | 1.125 | 0.354 | 3.175 | [3.14, 3.20] | win | 3.471 | 14.714 | True |
| split_bam | W-mid | 8.501 | 2.731 | 3.113 | [3.10, 3.12] | win | 3.185 | 2.375 | True |
| junction_annotation | W-mid | 1.984 | 0.713 | 2.782 | [2.76, 2.80] | win | 3.045 | 1.485 | True |
| mismatch_profile | W-mid | 1.337 | 0.496 | 2.693 | [2.66, 2.72] | win | 2.885 | 10.319 | True |
| split_paired_bam | W-mid | 7.977 | 3.286 | 2.427 | [2.42, 2.44] | win | 2.436 | 8.152 | True |
| junction_saturation | W-mid | 2.246 | 0.926 | 2.425 | [2.41, 2.44] | win | 2.562 | 2.147 | True |
| bam2wig | W-mid | 21.563 | 8.967 | 2.405 | [2.38, 2.42] | win | 2.427 | 1.401 | True |
| bam2fq | W-mid | 1.869 | 0.821 | 2.278 | [2.27, 2.30] | win | 2.561 | 8.052 | True |
| read_GC | W-mid | 1.583 | 0.795 | 1.990 | [1.98, 2.00] | win | 2.096 | 12.450 | True |
| read_duplication | W-mid | 2.641 | 1.745 | 1.514 | [1.50, 1.54] | win | 1.583 | 1.465 | True |

## 2. E1 compute-only (fixed per-invocation cost subtracted)

Upstream pays a large fixed cost before reading any input (interpreter start plus imports). E1 subtracts each arm's *measured* floor (`--help` on the real binary) so the number reflects algorithmic work. E2 includes it, because a user always pays it.

### compute-only vs end-to-end

| command | py floor (s) | rs floor (s) | py net (s) | rs net (s) | E1 x | py E2 (s) | rs E2 (s) |
|---|---|---|---|---|---|---|---|
| normalize_bigwig | 0.089 | 0.003 | 47.802 | 0.642 | 74.510 | 47.891 | 0.644 |
| geneBody_coverage | 0.089 | 0.003 | 61.928 | 1.329 | 46.611 | 62.017 | 1.332 |
| tin | 0.092 | 0.003 | 148.619 | 9.106 | 16.322 | 148.710 | 9.108 |
| read_NVC | 0.091 | 0.003 | 14.720 | 1.325 | 11.106 | 14.811 | 1.328 |
| inner_distance | 0.091 | 0.003 | 7.975 | 0.730 | 10.926 | 8.067 | 0.733 |
| RNA_fragment_size | 0.085 | 0.003 | 54.892 | 5.364 | 10.233 | 54.977 | 5.367 |
| sc_seqLogo | 0.509 | 0.003 | 0.325 | 0.040 | 8.021 | 0.834 | 0.043 |
| sc_seqQual | 0.514 | 0.003 | 0.427 | 0.054 | 7.909 | 0.941 | 0.057 |
| RPKM_saturation | 0.094 | 0.003 | 29.788 | 4.128 | 7.216 | 29.882 | 4.131 |
| geneBody_coverage2 | 0.080 | 0.003 | 25.120 | 4.218 | 5.955 | 25.199 | 4.221 |
| overlay_bigwig | 0.087 | 0.003 | 0.470 | 0.108 | 4.338 | 0.556 | 0.111 |
| bam_stat | 0.089 | 0.003 | 1.801 | 0.426 | 4.227 | 1.890 | 0.429 |
| FPKM_count | 0.093 | 0.003 | 33.900 | 8.124 | 4.173 | 33.993 | 8.127 |
| read_hexamer | 0.077 | 0.003 | 0.562 | 0.158 | 3.559 | 0.639 | 0.160 |
| clipping_profile | 0.091 | 0.003 | 1.585 | 0.457 | 3.472 | 1.676 | 0.460 |
| read_quality | 0.094 | 0.003 | 10.037 | 2.912 | 3.447 | 10.131 | 2.915 |
| insertion_profile | 0.093 | 0.003 | 1.556 | 0.462 | 3.370 | 1.649 | 0.464 |
| infer_experiment | 0.090 | 0.003 | 0.612 | 0.182 | 3.364 | 0.702 | 0.185 |
| read_distribution | 0.094 | 0.003 | 4.703 | 1.455 | 3.234 | 4.797 | 1.458 |
| deletion_profile | 0.092 | 0.003 | 1.125 | 0.354 | 3.175 | 1.217 | 0.357 |
| split_bam | 0.091 | 0.003 | 8.501 | 2.731 | 3.113 | 8.593 | 2.734 |
| junction_annotation | 0.092 | 0.003 | 1.984 | 0.713 | 2.782 | 2.076 | 0.716 |
| mismatch_profile | 0.090 | 0.003 | 1.337 | 0.496 | 2.693 | 1.427 | 0.500 |
| split_paired_bam | 0.056 | 0.003 | 7.977 | 3.286 | 2.427 | 8.033 | 3.289 |
| junction_saturation | 0.093 | 0.003 | 2.246 | 0.926 | 2.425 | 2.340 | 0.929 |
| bam2wig | 0.089 | 0.003 | 21.563 | 8.967 | 2.405 | 21.652 | 8.970 |
| bam2fq | 0.089 | 0.002 | 1.869 | 0.821 | 2.278 | 1.958 | 0.823 |
| read_GC | 0.093 | 0.003 | 1.583 | 0.795 | 1.990 | 1.676 | 0.798 |
| read_duplication | 0.091 | 0.003 | 2.641 | 1.745 | 1.514 | 2.732 | 1.748 |

## 3. Peak memory (whole process tree, `/usr/bin/time -v`)

### peak RSS; a ratio above 1 means the port uses more

| command | py median MB | rs median MB | mem x | py max MB | rs max MB |
|---|---|---|---|---|---|
| bam2wig | 1111 | 793.232 | 1.401 | 1112 | 793.578 |
| read_duplication | 265.584 | 181.236 | 1.465 | 265.898 | 181.555 |
| junction_annotation | 58.861 | 39.637 | 1.485 | 59.152 | 39.902 |
| FPKM_count | 92.025 | 58.352 | 1.577 | 92.258 | 59.703 |
| RPKM_saturation | 348.303 | 178.666 | 1.949 | 348.500 | 180.840 |
| tin | 42.639 | 20.639 | 2.066 | 42.836 | 20.926 |
| junction_saturation | 65.938 | 30.713 | 2.147 | 66.164 | 31.055 |
| RNA_fragment_size | 38.303 | 17.562 | 2.181 | 38.516 | 17.699 |
| split_bam | 132.676 | 55.854 | 2.375 | 133.004 | 56.023 |
| geneBody_coverage | 41.084 | 13.994 | 2.936 | 41.336 | 14.207 |
| read_distribution | 147.641 | 48.312 | 3.056 | 147.781 | 49.781 |
| read_hexamer | 32.311 | 8.086 | 3.996 | 32.473 | 8.211 |
| inner_distance | 212.844 | 53.232 | 3.998 | 213.031 | 53.473 |
| overlay_bigwig | 45.572 | 8.553 | 5.328 | 45.730 | 8.660 |
| normalize_bigwig | 66.232 | 8.410 | 7.875 | 66.438 | 8.516 |
| infer_experiment | 41.088 | 5.195 | 7.909 | 41.344 | 5.289 |
| bam2fq | 38.326 | 4.760 | 8.052 | 38.738 | 4.898 |
| split_paired_bam | 25.379 | 3.113 | 8.152 | 25.645 | 3.168 |
| mismatch_profile | 38.574 | 3.738 | 10.319 | 38.898 | 3.965 |
| read_quality | 38.729 | 3.207 | 12.076 | 38.918 | 3.301 |
| geneBody_coverage2 | 34.143 | 2.822 | 12.098 | 34.410 | 3.000 |
| read_NVC | 38.662 | 3.129 | 12.356 | 38.875 | 3.254 |
| read_GC | 38.469 | 3.090 | 12.450 | 38.762 | 3.211 |
| insertion_profile | 38.582 | 2.986 | 12.920 | 38.898 | 3.055 |
| clipping_profile | 38.600 | 2.984 | 12.934 | 38.836 | 3.102 |
| bam_stat | 37.941 | 2.902 | 13.073 | 38.277 | 3.027 |
| deletion_profile | 38.424 | 2.611 | 14.714 | 38.840 | 2.703 |
| sc_seqLogo | 98.832 | 6.424 | 15.385 | 99.250 | 6.570 |
| sc_seqQual | 100.082 | 6.387 | 15.670 | 100.309 | 6.562 |

## 4. Equivalence gate

No speedup is claimed for a row whose gate did not pass. A fast but incomplete run must not count as a win.

All measured rows passed structural equivalence.

### Declared exclusions

- **FPKM-UQ** — unsupported-missing-dependency: requires htseq-count, absent here
- **sc_editMatrix** — unsupported-missing-dependency: requires R package pheatmap, absent here
- **divide_bam** — diverged-documented DIV-0017: the two builds use different RNG algorithms, so --seed assigns different query names to subsets. The gate compares the union of subsets (same partition), not per-file membership. Per-file equality is explicitly NOT claimed.
- **sc_bamStat** — unsupported-input-shape: needs a single-cell BAM with CB/RE tags; bulk input makes both arms fail identically before doing work

## 5. Scaling along cost drivers


**bam_stat** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.011 | 0.0034 | 3.24x | [2.08, 5.00] | 0.098 | 0.006 | True |
| 10000 | 0.048 | 0.0140 | 3.44x | [2.60, 3.87] | 0.136 | 0.017 | True |
| 50000 | 0.231 | 0.0573 | 4.02x | [3.79, 4.17] | 0.320 | 0.060 | True |
| 200000 | 0.899 | 0.2090 | 4.30x | [4.12, 4.37] | 0.988 | 0.212 | True |
| 800000 | 3.476 | 0.7569 | 4.59x | [4.57, 4.66] | 3.564 | 0.760 | True |

**read_NVC** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.074 | 0.0081 | 9.09x | [6.86, 10.02] | 0.167 | 0.011 | True |
| 10000 | 0.376 | 0.0376 | 10.01x | [9.79, 10.53] | 0.465 | 0.041 | True |
| 50000 | 1.838 | 0.1751 | 10.49x | [10.41, 10.91] | 1.928 | 0.178 | True |
| 200000 | 7.304 | 0.6446 | 11.33x | [11.00, 11.47] | 7.393 | 0.647 | True |
| 800000 | 29.257 | 2.5412 | 11.51x | [11.35, 11.59] | 29.348 | 2.544 | True |

**read_quality** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.049 | 0.0162 | 3.03x | [2.67, 3.35] | 0.140 | 0.019 | True |
| 10000 | 0.253 | 0.0790 | 3.20x | [3.03, 3.38] | 0.344 | 0.082 | True |
| 50000 | 1.266 | 0.3772 | 3.36x | [3.31, 3.60] | 1.357 | 0.380 | True |
| 200000 | 5.011 | 1.4634 | 3.42x | [3.36, 3.44] | 5.101 | 1.466 | True |
| 800000 | 20.085 | 5.7585 | 3.49x | [3.48, 3.55] | 20.176 | 5.761 | True |

**FPKM_count** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 7.263 | 0.5076 | 14.31x | [14.12, 14.54] | 7.359 | 0.511 | True |
| 1000 | 13.198 | 0.8857 | 14.90x | [14.32, 15.03] | 13.289 | 0.889 | True |
| 5000 | 14.146 | 1.6900 | 8.37x | [8.25, 9.20] | 14.238 | 1.693 | True |
| 20000 | 15.428 | 2.2370 | 6.90x | [6.62, 7.15] | 15.520 | 2.240 | True |

**read_distribution** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 1.136 | 0.2338 | 4.86x | [4.77, 5.01] | 1.228 | 0.237 | True |
| 1000 | 1.359 | 0.2601 | 5.22x | [5.10, 5.30] | 1.452 | 0.263 | True |
| 5000 | 1.636 | 0.3808 | 4.30x | [4.24, 4.38] | 1.730 | 0.384 | True |
| 20000 | 1.913 | 0.5075 | 3.77x | [3.75, 3.83] | 2.005 | 0.510 | True |

**read_hexamer** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 0.204 | 0.0600 | 3.41x | [3.38, 3.58] | 0.280 | 0.063 | True |
| 75 | 0.286 | 0.0942 | 3.04x | [2.96, 3.11] | 0.361 | 0.097 | True |
| 100 | 0.377 | 0.1079 | 3.49x | [3.36, 3.63] | 0.453 | 0.111 | True |
| 150 | 0.553 | 0.1802 | 3.07x | [2.90, 3.10] | 0.627 | 0.183 | True |

**read_quality** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 1.345 | 0.3758 | 3.58x | [3.42, 3.72] | 1.438 | 0.379 | True |
| 75 | 1.938 | 0.5536 | 3.50x | [3.42, 3.53] | 2.029 | 0.556 | True |
| 100 | 2.518 | 0.7389 | 3.41x | [3.39, 3.47] | 2.610 | 0.742 | True |
| 150 | 3.705 | 1.0954 | 3.38x | [3.35, 3.43] | 3.795 | 1.098 | True |

## 6. Findings

These are the conclusions the data support. Every number in sections 1-5 is regenerated
from `results-main.json` / `results-scaling.json` by `analyze.py`; none is hand-written.
The tables in this section restate them alongside measurements taken outside the harness
(differential output checks, profiles), and are marked as such.

### 6.1 All 29 measured commands are faster, and all 29 passed the equivalence gate

The first measured run found three commands that were **reproducibly slower in the port**,
each with a confidence interval entirely below 1.0. All three were root-caused and fixed,
and this run re-measures them as wins: `inner_distance` **10.93x**, `bam2fq` **2.28x**,
`infer_experiment` **3.36x** (the last was 0.18x, then 1.09x after its first fix).

| | first run | this run |
|---|---|---|
| gates failed | 1 (`RPKM_saturation`) | **0** |
| commands slower than upstream | 3 | **0** |
| fastest command | `infer_experiment` 0.18x (as a loss) | `normalize_bigwig` 74.5x |
| slowest command | `infer_experiment` 0.18x | `read_duplication` 1.51x |

**The lesson from those three is still the most useful thing in this study.** They shared
a symptom (all were `open_alignments` callers) and I inferred a single shared cause -- the
eager whole-file decode. The profile refuted it: `bam_stat` performs the identical decode
in 0.50 s while `infer_experiment` took 4.07 s. The user/sys split then showed three
*different* signatures -- CPU-bound, syscall-bound, and algorithmic -- and three unrelated
fixes. **A correlation across two metrics is evidence for a hypothesis, not proof of one.**

### 6.2 A single speedup number is meaningless; the curve is the result

`bam_stat` across a 400x range of input sizes (section 5):

| reads | compute-only (E1) | end-to-end (E2) |
|---|---|---|
| 2,000 | 3.24x | 16.44x |
| 10,000 | 3.44x | 8.18x |
| 50,000 | 4.02x | 5.31x |
| 200,000 | 4.30x | 4.66x |
| 800,000 | 4.59x | 4.69x |

Same code, same machine, **opposite trends**: the end-to-end figure *falls* by 3.5x
while the compute-only figure *rises* toward ~4.6x. The reason is the reference's
~0.09 s fixed interpreter/import cost, which is a constant that dominates small inputs
and vanishes on large ones. At 2,000 reads the port's own arithmetic is only 3.2x faster;
the other 13x is CPython startup that any real user pays but that has nothing to do with
RSeQC.

A benchmark run at 1,000 reads would have reported ~69x for this command -- almost entirely
CPython startup. **Any claim about this port must state its workload size, and must say
whether it includes interpreter startup.**

### 6.3 Upstream's CPU-time column is mostly BLAS thread initialisation

Upstream RSeQC has no `multiprocessing` and no explicit threading (verified across all 33
scripts and `qcmodule`). Nevertheless `import numpy` alone costs **0.08 s wall but 0.90 s
of CPU** on this 16-thread machine, because OpenBLAS starts a thread pool per process;
with `OPENBLAS_NUM_THREADS=1` that falls to 0.05 s.

So an unpinned CPU measurement of upstream is ~90% thread-pool start-up, and the resulting
wall-vs-CPU gap looks like parallelism that does not exist. All runs pin
`OPENBLAS_NUM_THREADS=1` and friends in **both** arms (protocol §4.1). This trap would
silently corrupt any naive benchmark of this comparison.

### 6.4 Memory: the port is lighter than upstream on all 29 commands

The first run found the port using **4-8x more** memory than upstream on 14 of 29 commands,
with peak RSS growing linearly at ~226 bytes/record while upstream stayed flat at ~39 MB.
At 50M read pairs -- an ordinary human RNA-seq BAM -- that extrapolated to **~23 GB in the
port against ~39 MB upstream**: a tool that could not open a normal dataset. Cause:
`open_alignments` returning `Vec<io::Result<Record>>`, decoding the whole file up front.

Three separate whole-file structures had to go, and all three are fixed:

| what | mechanism now | peak RSS before -> after | upstream |
|---|---|---|---|
| `open_alignments` (14 commands) | streaming `AlignmentRecords`; BAM decodes into one reusable buffer, SAM in 4096-record batches | `bam_stat` at 800k: 364 MB -> **2.9 MB** | 37.9 MB |
| `tin`'s own read index | sliding window over transcripts in coordinate order | 366 MB -> **20.3 MB** | 42.3 MB |
| `geneBody_coverage`'s own index | same sliding window | 75.1 MB -> **14.4 MB** | 40.6 MB |

The last two are exact rather than approximate, and for slightly different reasons.
`tin` scores one transcript at a time, so only reads that can still reach an unscored
transcript need to be resident: a read with `end <= tx_start` cannot overlap any later
transcript. `geneBody_coverage` aggregates across *all* transcripts at once, so that
argument does not apply -- but its accumulation is a **sum plus an OR**, both commutative,
so visiting transcripts in coordinate order instead of BED order cannot change the result.

**Result: every one of the 29 measured commands is now lighter than upstream**, from
`bam2wig` at 1.40x lighter (793 MB against 1111 MB, the largest absolute figure on either
side) to `sc_seqQual` at 15.6x. Streaming the reader also made things *faster*, not just
smaller, so none of this is a time-for-memory trade.

The memory sweep confirms the reader fix holds across a 400x range rather than at one
point: `bam_stat` peak RSS is flat at 2.9-3.0 MB from 2,000 to 800,000 reads.

CRAM remains whole-file buffered, and the reason is a real API constraint rather than a
preference: `noodles-cram` 0.99 exposes record iteration only as `records(&header)`,
which is **single-use**, and re-entering it on a drained reader returns a spurious
`TryFromIntError` instead of EOF (verified directly). Bounding it would need a
self-referential reader, not justified for the least common input format here.

### 6.5 `RPKM_saturation` is not deterministic upstream, and the port is as close to it as upstream is to itself

Two consecutive upstream runs on identical input differ -- upstream subsamples reads per
percentile with no RNG seed. Measured here: the port agrees with upstream on **58.91%** of
`rp.eRPKM.xls` rows and **59.01%** of `rp.rawCount.xls`; **upstream agrees with itself on
58.35%** and 58.47%. This time the port is *ahead* of upstream's own reproducibility by
0.56 and 0.54 pp. The gate does not require that -- it requires the port to be no more
than 2 pp behind (`CHANGES.md` §1.6), and it would have passed either way, because two
figures differing by 0.5 pp on a 10% disagreement baseline are not distinguishable.

The honest reading is that the port is statistically indistinguishable from the reference
on a comparison nothing deterministic can pass. It is not evidence that the port is
correct -- it is evidence that the port is no *worse* than the thing it replaces, on a
metric where the thing it replaces does not agree with itself.

### 6.6 Preregistered expectations: which held and which did not

Protocol §10 predicted six things. Outcome:

| Prediction | Result |
|---|---|
| E1 speedups far below E2 speedups | **Held.** `bam_stat` at 2,000 reads: 3.24x compute-only vs 16.44x end-to-end. The gap closes as input grows. |
| Speedup falls with input size | **Held for E2** (16.44x -> 4.69x on `bam_stat`). **Refined**: E1 *rises* to 4.59x; the two converge at large input, which is the defensible number. |
| Port loses memory on eager-decode commands | **Held, and it was the most consequential finding.** Now resolved: the port is lighter than upstream on all 29 commands (§6.4). |
| Whole-output buffering shows output-dependent memory growth | **Partly held.** `bam2wig` is the worst absolute case (793 MB port / 1111 MB upstream) but the port wins it. |
| `read_hexamer` scales worse than linear in read length | **Not held.** 3.42x / 3.09x / 3.51x / 3.07x at 50/75/100/150 bp, and `read_quality` is flat at 3.46x -> 3.34x over the same range. The per-base `String` allocation is not a bottleneck at these sizes. Prediction withdrawn. |
| At least one command is slower in the port | **Held, three times** (`infer_experiment`, `bam2fq`, `inner_distance`) -- all three now fixed and re-measured as wins (§6.1). |

The memory sweep also confirms the streaming fix holds across a 400x range rather than at
one point: `bam_stat` peak RSS is flat at 2.9-3.0 MB from 2,000 to 800,000 reads, while
upstream sits flat at ~37.9 MB. Before the fix the port tracked input size linearly at
~226 bytes/record.

---

## 7. What this study does not establish

- **Not publication-grade.** Shared, non-isolated hardware, `schedutil` governor, load
  average 0.86/1.80/2.31 at the start of this run (the machine was quieter than for
  previous runs, which is why absolute times are lower here across the board). Ratios are
  paired and interleaved so machine drift cannot favour one arm, but absolute times are
  not trustworthy to better than roughly +-10%, and the intervals measure run-to-run
  scatter on *this* machine, not the uncertainty of a publication measurement.
- **Run-to-run, the absolute numbers move more than the ratios do.** Comparing this run
  with the previous one on the same workloads and build lineage, `geneBody_coverage` moved
  36.2x -> 46.6x and `tin` 15.2x -> 16.3x. Nothing was changed in either; the machine was
  simply less loaded. Treat any single figure here as the ratio measured on a quiet-ish
  sandbox, not a portable constant.
- **A previous full run was discarded for exactly this reason and re-measured.** An
  orphaned `cargo test` process sat at 100% of one core with 20 GB resident for the entire
  duration of the preceding run, contaminating every row. It was killed and the whole
  suite re-run from scratch (`CHANGES.md` §2.5). The harness records load average at start
  but does not compare load before against load after, so it did not catch this itself --
  a real gap in the instrumentation.
- **No cold-cache series.** Warm-cache only, as declared in protocol §2.
- **No thread scaling.** The port has no first-party parallelism and no `--threads` flag,
  so a 1/2/4/8 sweep is impossible. This is the most valuable untested axis, and the
  streaming rewrite makes it more valuable: decode is no longer hidden behind a large
  allocation, so it is now the dominant remaining cost for the I/O-bound commands
  (`bam2fq` 2.28x, `read_GC` 1.99x, `read_duplication` 1.51x) rather than being masked.
- **Reads are simulated from real transcripts**, not a real sequencer run, and the gene
  model is real but expression is uniform-over-transcripts. Commands whose cost driver is
  duplicate structure (`read_duplication`), error artefacts, or expression skew are
  therefore measured on a favourable input. `read_duplication` has the smallest speedup of
  any winning command (1.51x), consistent with that caveat.
- **CRAM is still whole-file buffered**, so a CRAM input would undo the memory result for
  every command. Nothing here measures a CRAM workload; §6.4 explains why the port cannot
  currently stream it.
- **No single-cell panel.** `sc_bamStat` and `sc_editMatrix` are not benchmarkable here
  (section 4). `sc_seqLogo`/`sc_seqQual` run on bulk-derived FASTQ, not real barcode data.
- **Four commands could not be measured**: `FPKM-UQ` (needs `htseq-count`), `sc_editMatrix`
  (needs R `pheatmap`), `sc_bamStat` (needs CB/RE-tagged input), `divide_bam` (documented
  RNG divergence, DIV-0017).
- **Reproducibility gate != scientific validity.** Every row here means "this workload's
  outputs match". It does not mean the computation is right on real data, which is the T4
  held-out question and remains open.
- **The benchmark's workloads cannot catch everything.** The 14 MB win for
  `geneBody_coverage` and the 20 MB win for `tin` were measured on generated reads that
  are cleanly block-separated by chromosome, and those workloads did *not* catch the
  stream-consumption regression that `verification/run_diff.py` did -- because a bug in
  chromosome-transition handling needs an interleaved record shape that the generator
  does not emit (`CHANGES.md` §2.7). A benchmark whose inputs this project generates is
  testing the project against its own assumptions. The independent fixture set remains the
  stronger correctness signal, and the two are not substitutes.
