# Benchmark results — rseqc-rust vs upstream RSeQC

Protocol: `benchmarks/protocol.md` (frozen 2026-09-29). Run label: `streaming-windowed`. Reps: 10, seed 20260929.

> **Not publication-grade.** Shared, non-isolated hardware (AMD Ryzen 7 3700X 8-Core Processor, loadavg 4.18 at start, governor `schedutil`). See protocol.md §2.

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
- **git_commit**: 1b3f40e6b875db2dca1ee8e604c9518a0e672b79
- **git_dirty**: True
- **pinned env**: `{'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONHASHSEED': '0'}`

## 1. Per-command results (E2 end-to-end, equal deliverables)

Speedup = median(upstream) / median(port) on the E1 compute-only quantities (each arm's measured per-invocation fixed cost subtracted, so the number and the interval describe the same thing). The interval is a 95% paired block bootstrap over whole pairs. A claim requires the interval to lie entirely above 1.0; an interval spanning 1.0 is inconclusive and an interval below 1.0 is a loss. Section 2 shows the end-to-end figures.

### all measured commands

| command | workload | py net (s) | rs net (s) | speedup | 95% CI | verdict | CPU x | mem x | gate |
|---|---|---|---|---|---|---|---|---|---|
| normalize_bigwig | W-mid | 53.123 | 0.711 | 74.716 | [70.72, 77.31] | win | 69.522 | 7.969 | True |
| geneBody_coverage | W-gbc | 68.828 | 1.901 | 36.202 | [36.01, 36.53] | win | 37.538 | 0.541 | True |
| tin | W-tn | 168.429 | 11.061 | 15.227 | [15.11, 15.51] | win | 15.558 | 2.081 | True |
| inner_distance | W-mid | 9.155 | 0.819 | 11.177 | [10.95, 11.46] | win | 12.535 | 4.107 | True |
| read_NVC | W-mid | 16.330 | 1.517 | 10.763 | [10.67, 10.89] | win | 10.977 | 12.588 | True |
| RNA_fragment_size | W-mid | 58.968 | 6.376 | 9.248 | [7.51, 9.35] | win | 9.327 | 2.189 | True |
| sc_seqLogo | W-mid | 0.375 | 0.047 | 7.983 | [7.35, 8.57] | win | 23.375 | 15.362 | True |
| sc_seqQual | W-mid | 0.465 | 0.061 | 7.654 | [7.34, 7.84] | win | 20.900 | 15.602 | True |
| RPKM_saturation | W-mid | 32.621 | 4.699 | 6.943 | [6.86, 7.29] | win | 7.234 | 1.948 | True |
| geneBody_coverage2 | W-mid | 28.422 | 4.645 | 6.119 | [6.07, 6.20] | win | 17.009 | 11.354 | True |
| overlay_bigwig | W-mid | 0.529 | 0.118 | 4.471 | [4.32, 4.51] | win | 6.444 | 5.425 | True |
| bam_stat | W-mid | 2.017 | 0.470 | 4.289 | [4.20, 4.36] | win | 4.622 | 12.929 | True |
| FPKM_count | W-mid | 36.355 | 9.428 | 3.856 | [3.06, 3.96] | win | 3.881 | 1.554 | True |
| read_hexamer | W-mid | 0.650 | 0.180 | 3.603 | [3.50, 3.65] | win | 4.147 | 4.081 | True |
| clipping_profile | W-mid | 1.748 | 0.500 | 3.492 | [3.45, 3.54] | win | 3.771 | 12.780 | True |
| infer_experiment | W-mid | 0.692 | 0.202 | 3.430 | [3.38, 3.54] | win | 4.000 | 7.934 | True |
| insertion_profile | W-mid | 1.754 | 0.518 | 3.387 | [3.34, 3.44] | win | 3.630 | 12.916 | True |
| read_quality | W-mid | 11.206 | 3.311 | 3.384 | [3.37, 3.41] | win | 3.422 | 12.154 | True |
| read_distribution | W-mid | 5.451 | 1.652 | 3.300 | [3.25, 3.33] | win | 3.310 | 3.057 | True |
| split_bam | W-mid | 10.062 | 3.102 | 3.244 | [3.22, 3.29] | win | 3.329 | 2.377 | True |
| deletion_profile | W-mid | 1.248 | 0.396 | 3.154 | [3.09, 3.23] | win | 3.493 | 14.473 | True |
| junction_annotation | W-mid | 2.245 | 0.838 | 2.678 | [2.58, 2.74] | win | 3.033 | 1.487 | True |
| mismatch_profile | W-mid | 1.485 | 0.555 | 2.678 | [2.67, 2.71] | win | 2.925 | 10.053 | True |
| split_paired_bam | W-mid | 9.053 | 3.704 | 2.444 | [2.42, 2.46] | win | 2.459 | 8.039 | True |
| bam2wig | W-mid | 24.198 | 9.993 | 2.421 | [2.40, 2.44] | win | 2.453 | 1.401 | True |
| bam2fq | W-mid | 2.092 | 0.903 | 2.316 | [2.28, 2.36] | win | 2.649 | 8.142 | True |
| junction_saturation | W-mid | 2.809 | 1.319 | 2.130 | [2.11, 2.17] | win | 2.240 | 2.143 | True |
| read_GC | W-mid | 1.749 | 0.889 | 1.968 | [1.96, 2.00] | win | 2.110 | 12.719 | True |
| read_duplication | W-mid | 3.139 | 2.127 | 1.476 | [1.45, 1.48] | win | 1.545 | 1.467 | True |

## 2. E1 compute-only (fixed per-invocation cost subtracted)

Upstream pays a large fixed cost before reading any input (interpreter start plus imports). E1 subtracts each arm's *measured* floor (`--help` on the real binary) so the number reflects algorithmic work. E2 includes it, because a user always pays it.

### compute-only vs end-to-end

| command | py floor (s) | rs floor (s) | py net (s) | rs net (s) | E1 x | py E2 (s) | rs E2 (s) |
|---|---|---|---|---|---|---|---|
| normalize_bigwig | 0.116 | 0.003 | 53.123 | 0.711 | 74.716 | 53.239 | 0.714 |
| geneBody_coverage | 0.110 | 0.003 | 68.828 | 1.901 | 36.202 | 68.939 | 1.904 |
| tin | 0.119 | 0.002 | 168.429 | 11.061 | 15.227 | 168.548 | 11.064 |
| inner_distance | 0.108 | 0.003 | 9.155 | 0.819 | 11.177 | 9.263 | 0.822 |
| read_NVC | 0.117 | 0.002 | 16.330 | 1.517 | 10.763 | 16.447 | 1.520 |
| RNA_fragment_size | 0.102 | 0.002 | 58.968 | 6.376 | 9.248 | 59.070 | 6.378 |
| sc_seqLogo | 0.651 | 0.002 | 0.375 | 0.047 | 7.983 | 1.026 | 0.049 |
| sc_seqQual | 0.672 | 0.003 | 0.465 | 0.061 | 7.654 | 1.137 | 0.063 |
| RPKM_saturation | 0.105 | 0.003 | 32.621 | 4.699 | 6.943 | 32.726 | 4.701 |
| geneBody_coverage2 | 0.096 | 0.002 | 28.422 | 4.645 | 6.119 | 28.517 | 4.647 |
| overlay_bigwig | 0.109 | 0.002 | 0.529 | 0.118 | 4.471 | 0.638 | 0.120 |
| bam_stat | 0.111 | 0.003 | 2.017 | 0.470 | 4.289 | 2.128 | 0.473 |
| FPKM_count | 0.102 | 0.002 | 36.355 | 9.428 | 3.856 | 36.457 | 9.430 |
| read_hexamer | 0.092 | 0.002 | 0.650 | 0.180 | 3.603 | 0.742 | 0.183 |
| clipping_profile | 0.110 | 0.003 | 1.748 | 0.500 | 3.492 | 1.857 | 0.503 |
| infer_experiment | 0.107 | 0.002 | 0.692 | 0.202 | 3.430 | 0.799 | 0.204 |
| insertion_profile | 0.110 | 0.002 | 1.754 | 0.518 | 3.387 | 1.864 | 0.520 |
| read_quality | 0.118 | 0.002 | 11.206 | 3.311 | 3.384 | 11.324 | 3.314 |
| read_distribution | 0.115 | 0.002 | 5.451 | 1.652 | 3.300 | 5.567 | 1.654 |
| split_bam | 0.115 | 0.003 | 10.062 | 3.102 | 3.244 | 10.176 | 3.105 |
| deletion_profile | 0.115 | 0.002 | 1.248 | 0.396 | 3.154 | 1.363 | 0.398 |
| junction_annotation | 0.144 | 0.004 | 2.245 | 0.838 | 2.678 | 2.389 | 0.842 |
| mismatch_profile | 0.115 | 0.003 | 1.485 | 0.555 | 2.678 | 1.600 | 0.558 |
| split_paired_bam | 0.069 | 0.003 | 9.053 | 3.704 | 2.444 | 9.122 | 3.707 |
| bam2wig | 0.111 | 0.002 | 24.198 | 9.993 | 2.421 | 24.308 | 9.996 |
| bam2fq | 0.108 | 0.002 | 2.092 | 0.903 | 2.316 | 2.200 | 0.906 |
| junction_saturation | 0.115 | 0.003 | 2.809 | 1.319 | 2.130 | 2.924 | 1.321 |
| read_GC | 0.117 | 0.003 | 1.749 | 0.889 | 1.968 | 1.866 | 0.891 |
| read_duplication | 0.119 | 0.002 | 3.139 | 2.127 | 1.476 | 3.258 | 2.130 |

## 3. Peak memory (whole process tree, `/usr/bin/time -v`)

### peak RSS; a ratio above 1 means the port uses more

| command | py median MB | rs median MB | mem x | py max MB | rs max MB |
|---|---|---|---|---|---|
| geneBody_coverage | 40.586 | 75.055 | 0.541 | 41.020 | 75.578 |
| bam2wig | 1111 | 793.055 | 1.401 | 1111 | 793.223 |
| read_duplication | 264.982 | 180.580 | 1.467 | 265.297 | 180.840 |
| junction_annotation | 58.152 | 39.102 | 1.487 | 58.418 | 39.781 |
| FPKM_count | 91.676 | 59.004 | 1.554 | 91.797 | 59.418 |
| RPKM_saturation | 347.867 | 178.549 | 1.948 | 348.102 | 181.047 |
| tin | 42.264 | 20.311 | 2.081 | 42.742 | 20.594 |
| junction_saturation | 65.334 | 30.482 | 2.143 | 65.871 | 30.633 |
| RNA_fragment_size | 37.865 | 17.299 | 2.189 | 38.258 | 17.543 |
| split_bam | 132.137 | 55.592 | 2.377 | 132.441 | 56.012 |
| read_distribution | 146.980 | 48.076 | 3.057 | 147.605 | 49.992 |
| read_hexamer | 31.775 | 7.787 | 4.081 | 32.035 | 8.086 |
| inner_distance | 212.189 | 51.662 | 4.107 | 212.641 | 53.164 |
| overlay_bigwig | 44.969 | 8.289 | 5.425 | 45.398 | 8.707 |
| infer_experiment | 40.814 | 5.145 | 7.934 | 40.914 | 5.254 |
| normalize_bigwig | 65.725 | 8.248 | 7.969 | 66.020 | 8.410 |
| split_paired_bam | 24.824 | 3.088 | 8.039 | 25.320 | 3.188 |
| bam2fq | 37.945 | 4.660 | 8.142 | 38.398 | 4.828 |
| mismatch_profile | 38.129 | 3.793 | 10.053 | 38.793 | 3.980 |
| geneBody_coverage2 | 33.707 | 2.969 | 11.354 | 34.000 | 3.004 |
| read_quality | 38.646 | 3.180 | 12.154 | 38.863 | 3.352 |
| read_NVC | 37.984 | 3.018 | 12.588 | 38.441 | 3.082 |
| read_GC | 37.934 | 2.982 | 12.719 | 38.809 | 3.273 |
| clipping_profile | 38.090 | 2.980 | 12.780 | 38.457 | 3.074 |
| insertion_profile | 38.041 | 2.945 | 12.916 | 38.625 | 3.086 |
| bam_stat | 37.727 | 2.918 | 12.929 | 37.996 | 3.051 |
| deletion_profile | 37.963 | 2.623 | 14.473 | 38.609 | 2.684 |
| sc_seqLogo | 97.963 | 6.377 | 15.362 | 98.812 | 6.523 |
| sc_seqQual | 99.312 | 6.365 | 15.602 | 99.789 | 6.480 |

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
| 2000 | 0.012 | 0.0036 | 3.37x | [1.78, 5.22] | 0.114 | 0.006 | True |
| 10000 | 0.047 | 0.0153 | 3.09x | [2.84, 3.91] | 0.151 | 0.018 | True |
| 50000 | 0.253 | 0.0665 | 3.80x | [3.74, 4.25] | 0.364 | 0.069 | True |
| 200000 | 0.991 | 0.2221 | 4.46x | [4.11, 4.61] | 1.098 | 0.224 | True |
| 800000 | 3.889 | 0.8466 | 4.59x | [4.49, 4.70] | 3.998 | 0.849 | True |

**read_NVC** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.081 | 0.0097 | 8.33x | [7.18, 9.44] | 0.188 | 0.012 | True |
| 10000 | 0.404 | 0.0434 | 9.31x | [8.78, 10.12] | 0.511 | 0.046 | True |
| 50000 | 2.025 | 0.1936 | 10.46x | [9.97, 10.95] | 2.133 | 0.196 | True |
| 200000 | 8.125 | 0.7371 | 11.02x | [10.66, 11.09] | 8.239 | 0.740 | True |
| 800000 | 32.518 | 2.9072 | 11.19x | [10.94, 11.30] | 32.639 | 2.910 | True |

**read_quality** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.050 | 0.0198 | 2.53x | [2.10, 3.08] | 0.161 | 0.022 | True |
| 10000 | 0.282 | 0.0897 | 3.14x | [2.99, 3.37] | 0.389 | 0.092 | True |
| 50000 | 1.404 | 0.4141 | 3.39x | [3.15, 3.43] | 1.511 | 0.417 | True |
| 200000 | 5.541 | 1.6271 | 3.41x | [3.32, 3.48] | 5.651 | 1.630 | True |
| 800000 | 22.401 | 6.5292 | 3.43x | [3.40, 3.43] | 22.516 | 6.532 | True |

**FPKM_count** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 8.218 | 0.6136 | 13.39x | [13.20, 13.53] | 8.338 | 0.617 | True |
| 1000 | 14.974 | 1.2234 | 12.24x | [11.96, 12.58] | 15.095 | 1.226 | True |
| 5000 | 16.256 | 2.6114 | 6.23x | [6.17, 6.40] | 16.375 | 2.614 | True |
| 20000 | 17.689 | 3.8404 | 4.61x | [4.23, 4.68] | 17.807 | 3.843 | True |

**read_distribution** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 1.304 | 0.2583 | 5.05x | [4.96, 5.12] | 1.424 | 0.261 | True |
| 1000 | 1.608 | 0.2970 | 5.41x | [5.39, 5.61] | 1.734 | 0.299 | True |
| 5000 | 1.953 | 0.4404 | 4.44x | [4.33, 4.50] | 2.072 | 0.443 | True |
| 20000 | 2.301 | 0.5877 | 3.92x | [3.73, 3.99] | 2.421 | 0.590 | True |

**read_hexamer** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 0.230 | 0.0673 | 3.42x | [3.31, 3.53] | 0.323 | 0.070 | True |
| 75 | 0.333 | 0.1076 | 3.09x | [3.03, 3.16] | 0.428 | 0.111 | True |
| 100 | 0.438 | 0.1247 | 3.51x | [3.37, 3.77] | 0.535 | 0.127 | True |
| 150 | 0.623 | 0.2032 | 3.07x | [2.97, 3.13] | 0.718 | 0.206 | True |

**read_quality** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 1.500 | 0.4337 | 3.46x | [3.43, 3.52] | 1.621 | 0.436 | True |
| 75 | 2.161 | 0.6356 | 3.40x | [3.35, 3.44] | 2.280 | 0.639 | True |
| 100 | 2.814 | 0.8392 | 3.35x | [3.34, 3.39] | 2.933 | 0.842 | True |
| 150 | 4.214 | 1.2615 | 3.34x | [3.32, 3.36] | 4.329 | 1.265 | True |

## 6. Findings

These are the conclusions the data support. Every number in sections 1-5 is regenerated
from `results-main.json` / `results-scaling.json` by `analyze.py`; none is hand-written.
The tables in this section restate them alongside measurements taken outside the harness
(differential output checks, profiles), and are marked as such.

### 6.1 All 29 measured commands are faster, and all 29 passed the equivalence gate

The first measured run found three commands that were **reproducibly slower in the port**,
each with a confidence interval entirely below 1.0. All three were root-caused and fixed,
and this run re-measures them as wins: `inner_distance` **11.18x**, `bam2fq` **2.32x**,
`infer_experiment` **3.43x** (the last was 0.18x, then 1.09x after its first fix).

| | first run | this run |
|---|---|---|
| gates failed | 1 (`RPKM_saturation`) | **0** |
| commands slower than upstream | 3 | **0** |
| slowest command | `infer_experiment` 0.18x | `read_duplication` 1.48x |

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
| 2,000 | 3.37x | 17.61x |
| 10,000 | 3.09x | 8.54x |
| 50,000 | 3.80x | 5.29x |
| 200,000 | 4.46x | 4.89x |
| 800,000 | 4.59x | 4.71x |

Same code, same machine, **opposite trends**: the end-to-end figure *falls* by 3.7x
while the compute-only figure *rises* toward ~4.6x. The reason is the reference's
~0.09 s fixed interpreter/import cost, which is a constant that dominates small inputs
and vanishes on large ones. At 2,000 reads the port's own arithmetic is only 3.4x faster;
the other 14x is CPython startup that any real user pays but that has nothing to do with
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

### 6.4 Memory: the port is now lighter than upstream on 28 of 29 commands

The first run found the port using **4-8x more** memory than upstream on 14 of 29 commands,
with peak RSS growing linearly at ~226 bytes/record while upstream stayed flat at ~39 MB.
At 50M read pairs -- an ordinary human RNA-seq BAM -- that extrapolated to **~23 GB in the
port against ~39 MB upstream**: a tool that could not open a normal dataset. Cause:
`open_alignments` returning `Vec<io::Result<Record>>`, decoding the whole file up front.

That was fixed by returning a streaming `AlignmentRecords` iterator, and `tin`'s own
whole-file index by a sliding window. Both are covered in `CHANGELOG.md`; the second is
new to this run.

`tin` in particular went from **the worst row in the table** (366 MB, 8.5x heavier than
upstream) to **20.3 MB, 2.08x lighter than upstream**, and got faster at the same time
(19.8 s -> 11.1 s), because it no longer builds and sorts a whole-file index:

| | whole-file index | sliding window |
|---|---|---|
| `tin` peak RSS (600k reads) | 366 MB | **20.3 MB** |
| resident reads | 597,048 | mean 3,424 (max 16,438) |
| `tin` wall (E2) | 19.82 s | **11.06 s** |
| upstream | 42.3 MB / 168.4 s | 42.3 MB / 168.4 s |

The single remaining exception is **`geneBody_coverage`** (75.1 MB vs upstream's 40.6 MB).
It routes through the same `build_read_index` helper but computes coverage over *all*
positions of *all* transcripts in one pass, so it has no per-transcript window to stream;
it needs the whole-file index. It is the last command heavier than upstream, and unlike
`tin` it is not a small change -- it would need the same restructuring `tin` just got.

### 6.5 `RPKM_saturation` is not deterministic upstream, and the port is as close to it as upstream is to itself

Two consecutive upstream runs on identical input differ -- upstream subsamples reads per
percentile with no RNG seed. Measured here: the port agrees with upstream on **58.18%** of
`rp.eRPKM.xls` rows and **58.30%** of `rp.rawCount.xls`; **upstream agrees with itself on
58.43%** and 58.54%. The two are within 0.25 pp, inside the 2 pp tolerance declared in
`CHANGES.md` §1.6. This run is the first in which the row's gate passed (the previous run
failed at 58.3% vs 58.4% under the superseded strict-dominance rule). The honest reading
is that the port is indistinguishable from the reference on a comparison that nothing
deterministic can pass -- not that the port has been shown correct.

### 6.6 Preregistered expectations: which held and which did not

Protocol §10 predicted six things. Outcome:

| Prediction | Result |
|---|---|
| E1 speedups far below E2 speedups | **Held.** `bam_stat` at 2,000 reads: 3.37x compute-only vs 17.61x end-to-end. The gap closes as input grows. |
| Speedup falls with input size | **Held for E2** (17.61x -> 4.71x on `bam_stat`). **Refined**: E1 *rises* to 4.59x; the two converge at large input, which is the defensible number. |
| Port loses memory on eager-decode commands | **Held, and it was the most consequential finding.** Now resolved for 28 of 29 commands; only `geneBody_coverage` remains heavier than upstream (§6.4). |
| Whole-output buffering shows output-dependent memory growth | **Partly held.** `bam2wig` is the worst absolute case (793 MB port / 1111 MB upstream) but the port wins it. |
| `read_hexamer` scales worse than linear in read length | **Not held.** 3.42x / 3.09x / 3.51x / 3.07x at 50/75/100/150 bp, and `read_quality` is flat at 3.46x -> 3.34x over the same range. The per-base `String` allocation is not a bottleneck at these sizes. Prediction withdrawn. |
| At least one command is slower in the port | **Held, three times** (`infer_experiment`, `bam2fq`, `inner_distance`) -- all three now fixed and re-measured as wins (§6.1). |

The memory sweep also confirms the streaming fix holds across a 400x range rather than at
one point: `bam_stat` peak RSS is flat at 2.8-3.0 MB from 2,000 to 800,000 reads, while
upstream sits flat at ~37.7 MB. Before the fix the port tracked input size linearly at
~226 bytes/record.

---

## 7. What this study does not establish

- **Not publication-grade.** Shared, non-isolated hardware, `schedutil` governor, ~10 users
  on the machine, load average 4.2/5.0/5.0 at the start of this run. Ratios are paired
  and interleaved so
  machine drift cannot favour one arm, but absolute times are not trustworthy to better
  than roughly +-10%, and the intervals measure run-to-run scatter on *this* machine, not
  the uncertainty of a publication measurement.
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
  (`bam2fq` 2.32x, `read_GC` 1.97x, `read_duplication` 1.48x) rather than being masked.
- **Reads are simulated from real transcripts**, not a real sequencer run, and the gene
  model is real but expression is uniform-over-transcripts. Commands whose cost driver is
  duplicate structure (`read_duplication`), error artefacts, or expression skew are
  therefore measured on a favourable input. `read_duplication` has the smallest speedup of
  any winning command (1.48x), consistent with that caveat.
- **`geneBody_coverage` still uses more memory than upstream** (75.1 vs 40.6 MB, §6.4), so
  the port cannot yet claim to be uniformly lighter.
- **No single-cell panel.** `sc_bamStat` and `sc_editMatrix` are not benchmarkable here
  (section 4). `sc_seqLogo`/`sc_seqQual` run on bulk-derived FASTQ, not real barcode data.
- **Four commands could not be measured**: `FPKM-UQ` (needs `htseq-count`), `sc_editMatrix`
  (needs R `pheatmap`), `sc_bamStat` (needs CB/RE-tagged input), `divide_bam` (documented
  RNG divergence, DIV-0017).
- **Reproducibility gate != scientific validity.** Every row here means "this workload's
  outputs match". It does not mean the computation is right on real data, which is the T4
  held-out question and remains open.
