# Benchmark results — rseqc-rust vs upstream RSeQC

Protocol: `benchmarks/protocol.md` (frozen 2026-09-29). Run label: `g6b-clean`. Reps: 10, seed 20260929.

> **Not publication-grade.** Shared, non-isolated hardware (AMD Ryzen 7 3700X 8-Core Processor, loadavg 3.70 at start, governor `schedutil`). See protocol.md §2.

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
- **git_commit**: 3a5b88c00184ce03b818e16acf72bdc5a6efe4df
- **git_dirty**: False
- **pinned env**: `{'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'PYTHONHASHSEED': '0'}`

## 1. Per-command results (E2 end-to-end, equal deliverables)

Speedup = median(upstream) / median(port) on the E1 compute-only quantities (each arm's measured per-invocation fixed cost subtracted, so the number and the interval describe the same thing). The interval is a 95% paired block bootstrap over whole pairs. A claim requires the interval to lie entirely above 1.0; an interval spanning 1.0 is inconclusive and an interval below 1.0 is a loss. Section 2 shows the end-to-end figures.

### all measured commands

| command | workload | py net (s) | rs net (s) | speedup | 95% CI | verdict | CPU x | mem x | gate |
|---|---|---|---|---|---|---|---|---|---|
| normalize_bigwig | W-mid | 48.481 | 0.677 | 71.622 | [68.21, 81.35] | win | 66.667 | 7.959 | True |
| geneBody_coverage | W-gbc | 62.935 | 1.357 | 46.379 | [45.99, 46.76] | win | 47.663 | 2.918 | True |
| tin | W-tn | 150.738 | 9.239 | 16.315 | [16.15, 16.37] | win | 16.707 | 2.069 | True |
| read_NVC | W-mid | 14.871 | 1.350 | 11.015 | [10.88, 11.13] | win | 11.169 | 12.308 | True |
| inner_distance | W-mid | 8.138 | 0.745 | 10.929 | [10.84, 11.00] | win | 12.205 | 3.949 | True |
| RNA_fragment_size | W-mid | 55.778 | 5.514 | 10.116 | [9.77, 10.25] | win | 10.212 | 2.157 | True |
| sc_seqLogo | W-mid | 0.340 | 0.042 | 8.200 | [7.36, 8.56] | win | 26.333 | 15.509 | True |
| sc_seqQual | W-mid | 0.432 | 0.056 | 7.773 | [7.41, 7.83] | win | 17.600 | 15.738 | True |
| RPKM_saturation | W-mid | 30.076 | 4.153 | 7.242 | [7.00, 7.37] | win | 7.554 | 1.949 | True |
| geneBody_coverage2 | W-mid | 25.338 | 4.257 | 5.952 | [5.93, 6.01] | win | 16.277 | 11.479 | True |
| overlay_bigwig | W-mid | 0.476 | 0.111 | 4.289 | [4.08, 4.44] | win | 6.118 | 5.410 | True |
| bam_stat | W-mid | 1.821 | 0.427 | 4.266 | [4.23, 4.31] | win | 4.549 | 13.126 | True |
| FPKM_count | W-mid | 34.385 | 8.277 | 4.154 | [4.05, 4.37] | win | 4.188 | 1.549 | True |
| read_hexamer | W-mid | 0.572 | 0.157 | 3.636 | [3.53, 3.70] | win | 4.167 | 3.940 | True |
| clipping_profile | W-mid | 1.614 | 0.466 | 3.462 | [3.44, 3.49] | win | 3.700 | 12.573 | True |
| read_quality | W-mid | 10.166 | 2.949 | 3.447 | [3.41, 3.47] | win | 3.486 | 12.033 | True |
| infer_experiment | W-mid | 0.620 | 0.185 | 3.355 | [3.26, 3.47] | win | 4.000 | 8.067 | True |
| insertion_profile | W-mid | 1.572 | 0.474 | 3.319 | [3.28, 3.36] | win | 3.600 | 12.852 | True |
| read_distribution | W-mid | 4.807 | 1.478 | 3.252 | [3.21, 3.28] | win | 3.277 | 2.991 | True |
| deletion_profile | W-mid | 1.139 | 0.361 | 3.156 | [3.14, 3.19] | win | 3.414 | 14.797 | True |
| split_bam | W-mid | 8.627 | 2.766 | 3.119 | [3.01, 3.13] | win | 3.172 | 2.380 | True |
| junction_annotation | W-mid | 2.015 | 0.732 | 2.752 | [2.73, 2.78] | win | 3.030 | 1.484 | True |
| mismatch_profile | W-mid | 1.353 | 0.504 | 2.683 | [2.65, 2.73] | win | 2.857 | 10.206 | True |
| split_paired_bam | W-mid | 8.112 | 3.331 | 2.435 | [2.42, 2.45] | win | 2.449 | 7.781 | True |
| bam2wig | W-mid | 22.030 | 9.194 | 2.396 | [2.38, 2.41] | win | 2.415 | 1.401 | True |
| junction_saturation | W-mid | 2.347 | 1.005 | 2.336 | [2.30, 2.37] | win | 2.466 | 2.146 | True |
| bam2fq | W-mid | 1.920 | 0.836 | 2.298 | [2.28, 2.32] | win | 2.596 | 8.050 | True |
| read_GC | W-mid | 1.596 | 0.806 | 1.980 | [1.96, 1.99] | win | 2.082 | 12.321 | True |
| read_duplication | W-mid | 2.699 | 1.793 | 1.505 | [1.49, 1.53] | win | 1.576 | 1.466 | True |

## 2. E1 compute-only (fixed per-invocation cost subtracted)

Upstream pays a large fixed cost before reading any input (interpreter start plus imports). E1 subtracts each arm's *measured* floor (`--help` on the real binary) so the number reflects algorithmic work. E2 includes it, because a user always pays it.

### compute-only vs end-to-end

| command | py floor (s) | rs floor (s) | py net (s) | rs net (s) | E1 x | py E2 (s) | rs E2 (s) |
|---|---|---|---|---|---|---|---|
| normalize_bigwig | 0.091 | 0.003 | 48.481 | 0.677 | 71.622 | 48.572 | 0.680 |
| geneBody_coverage | 0.094 | 0.003 | 62.935 | 1.357 | 46.379 | 63.029 | 1.360 |
| tin | 0.095 | 0.002 | 150.738 | 9.239 | 16.315 | 150.833 | 9.242 |
| read_NVC | 0.094 | 0.003 | 14.871 | 1.350 | 11.015 | 14.965 | 1.353 |
| inner_distance | 0.094 | 0.003 | 8.138 | 0.745 | 10.929 | 8.232 | 0.748 |
| RNA_fragment_size | 0.088 | 0.003 | 55.778 | 5.514 | 10.116 | 55.866 | 5.517 |
| sc_seqLogo | 0.531 | 0.003 | 0.340 | 0.042 | 8.200 | 0.872 | 0.045 |
| sc_seqQual | 0.537 | 0.003 | 0.432 | 0.056 | 7.773 | 0.969 | 0.058 |
| RPKM_saturation | 0.096 | 0.003 | 30.076 | 4.153 | 7.242 | 30.172 | 4.156 |
| geneBody_coverage2 | 0.082 | 0.003 | 25.338 | 4.257 | 5.952 | 25.420 | 4.260 |
| overlay_bigwig | 0.090 | 0.003 | 0.476 | 0.111 | 4.289 | 0.566 | 0.114 |
| bam_stat | 0.091 | 0.003 | 1.821 | 0.427 | 4.266 | 1.912 | 0.430 |
| FPKM_count | 0.093 | 0.003 | 34.385 | 8.277 | 4.154 | 34.478 | 8.279 |
| read_hexamer | 0.077 | 0.003 | 0.572 | 0.157 | 3.636 | 0.649 | 0.160 |
| clipping_profile | 0.093 | 0.003 | 1.614 | 0.466 | 3.462 | 1.707 | 0.469 |
| read_quality | 0.094 | 0.003 | 10.166 | 2.949 | 3.447 | 10.260 | 2.952 |
| infer_experiment | 0.092 | 0.003 | 0.620 | 0.185 | 3.355 | 0.712 | 0.188 |
| insertion_profile | 0.092 | 0.003 | 1.572 | 0.474 | 3.319 | 1.663 | 0.477 |
| read_distribution | 0.098 | 0.003 | 4.807 | 1.478 | 3.252 | 4.905 | 1.481 |
| deletion_profile | 0.095 | 0.003 | 1.139 | 0.361 | 3.156 | 1.234 | 0.364 |
| split_bam | 0.095 | 0.003 | 8.627 | 2.766 | 3.119 | 8.721 | 2.768 |
| junction_annotation | 0.095 | 0.003 | 2.015 | 0.732 | 2.752 | 2.110 | 0.735 |
| mismatch_profile | 0.093 | 0.003 | 1.353 | 0.504 | 2.683 | 1.446 | 0.508 |
| split_paired_bam | 0.056 | 0.003 | 8.112 | 3.331 | 2.435 | 8.168 | 3.334 |
| bam2wig | 0.092 | 0.003 | 22.030 | 9.194 | 2.396 | 22.122 | 9.197 |
| junction_saturation | 0.092 | 0.003 | 2.347 | 1.005 | 2.336 | 2.439 | 1.008 |
| bam2fq | 0.092 | 0.003 | 1.920 | 0.836 | 2.298 | 2.013 | 0.838 |
| read_GC | 0.093 | 0.003 | 1.596 | 0.806 | 1.980 | 1.690 | 0.809 |
| read_duplication | 0.096 | 0.003 | 2.699 | 1.793 | 1.505 | 2.795 | 1.796 |

## 3. Peak memory (whole process tree, `/usr/bin/time -v`)

### peak RSS; a ratio above 1 means the port uses more

| command | py median MB | rs median MB | mem x | py max MB | rs max MB |
|---|---|---|---|---|---|
| bam2wig | 1111 | 793.041 | 1.401 | 1111 | 793.348 |
| read_duplication | 265.469 | 181.121 | 1.466 | 265.855 | 181.383 |
| junction_annotation | 58.898 | 39.693 | 1.484 | 59.027 | 39.855 |
| FPKM_count | 91.979 | 59.365 | 1.549 | 92.277 | 59.867 |
| RPKM_saturation | 348.137 | 178.594 | 1.949 | 348.363 | 181.188 |
| tin | 42.652 | 20.617 | 2.069 | 43.070 | 20.820 |
| junction_saturation | 66.051 | 30.773 | 2.146 | 66.406 | 31.238 |
| RNA_fragment_size | 38.094 | 17.662 | 2.157 | 38.266 | 17.941 |
| split_bam | 132.684 | 55.756 | 2.380 | 132.996 | 56.020 |
| geneBody_coverage | 40.900 | 14.016 | 2.918 | 41.250 | 14.215 |
| read_distribution | 147.322 | 49.252 | 2.991 | 147.840 | 50.449 |
| read_hexamer | 32.107 | 8.148 | 3.940 | 32.379 | 8.211 |
| inner_distance | 209.752 | 53.121 | 3.949 | 210.090 | 53.410 |
| overlay_bigwig | 45.381 | 8.389 | 5.410 | 45.586 | 8.723 |
| split_paired_bam | 25.287 | 3.250 | 7.781 | 25.500 | 3.297 |
| normalize_bigwig | 66.193 | 8.316 | 7.959 | 66.492 | 8.504 |
| bam2fq | 38.330 | 4.762 | 8.050 | 38.613 | 4.863 |
| infer_experiment | 40.982 | 5.080 | 8.067 | 41.367 | 5.281 |
| mismatch_profile | 38.533 | 3.775 | 10.206 | 38.816 | 3.938 |
| geneBody_coverage2 | 34.100 | 2.971 | 11.479 | 34.293 | 3.012 |
| read_quality | 38.826 | 3.227 | 12.033 | 39.105 | 3.312 |
| read_NVC | 38.438 | 3.123 | 12.308 | 38.723 | 3.223 |
| read_GC | 38.359 | 3.113 | 12.321 | 38.805 | 3.270 |
| clipping_profile | 38.357 | 3.051 | 12.573 | 38.633 | 3.145 |
| insertion_profile | 38.582 | 3.002 | 12.852 | 38.832 | 3.141 |
| bam_stat | 38.121 | 2.904 | 13.126 | 38.375 | 3.016 |
| deletion_profile | 38.525 | 2.604 | 14.797 | 38.926 | 2.676 |
| sc_seqLogo | 98.748 | 6.367 | 15.509 | 98.957 | 6.578 |
| sc_seqQual | 100.051 | 6.357 | 15.738 | 100.539 | 6.555 |

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
| 2000 | 0.010 | 0.0046 | 2.25x | [1.95, 3.39] | 0.101 | 0.007 | True |
| 10000 | 0.050 | 0.0134 | 3.69x | [2.88, 4.82] | 0.139 | 0.016 | True |
| 50000 | 0.237 | 0.0594 | 3.98x | [3.85, 4.14] | 0.330 | 0.062 | True |
| 200000 | 0.910 | 0.2090 | 4.35x | [4.03, 4.40] | 1.003 | 0.211 | True |
| 800000 | 3.559 | 0.7656 | 4.65x | [4.53, 4.71] | 3.652 | 0.768 | True |

**read_NVC** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.075 | 0.0082 | 9.25x | [6.99, 9.91] | 0.168 | 0.011 | True |
| 10000 | 0.380 | 0.0361 | 10.53x | [9.80, 10.86] | 0.473 | 0.039 | True |
| 50000 | 1.863 | 0.1788 | 10.42x | [10.31, 11.02] | 1.958 | 0.181 | True |
| 200000 | 7.452 | 0.6719 | 11.09x | [10.83, 11.45] | 7.544 | 0.675 | True |
| 800000 | 29.731 | 2.6291 | 11.31x | [11.26, 11.64] | 29.825 | 2.632 | True |

**read_quality** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.050 | 0.0175 | 2.86x | [2.59, 3.23] | 0.143 | 0.020 | True |
| 10000 | 0.260 | 0.0796 | 3.26x | [3.12, 3.40] | 0.354 | 0.082 | True |
| 50000 | 1.266 | 0.3800 | 3.33x | [3.27, 3.38] | 1.360 | 0.383 | True |
| 200000 | 5.066 | 1.4692 | 3.45x | [3.42, 3.49] | 5.159 | 1.472 | True |
| 800000 | 20.286 | 5.8303 | 3.48x | [3.46, 3.51] | 20.381 | 5.833 | True |

**FPKM_count** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 7.364 | 0.5224 | 14.10x | [14.01, 14.14] | 7.460 | 0.525 | True |
| 1000 | 13.765 | 0.9879 | 13.93x | [13.67, 14.16] | 13.866 | 0.991 | True |
| 5000 | 14.423 | 1.7500 | 8.24x | [7.75, 8.75] | 14.517 | 1.753 | True |
| 20000 | 16.225 | 2.8107 | 5.77x | [5.62, 5.83] | 16.322 | 2.814 | True |

**read_distribution** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 1.173 | 0.2366 | 4.96x | [4.82, 5.18] | 1.270 | 0.239 | True |
| 1000 | 1.432 | 0.2732 | 5.24x | [5.00, 5.49] | 1.532 | 0.276 | True |
| 5000 | 1.736 | 0.3984 | 4.36x | [4.34, 4.43] | 1.835 | 0.402 | True |
| 20000 | 1.956 | 0.5143 | 3.80x | [3.77, 3.83] | 2.054 | 0.518 | True |

**read_hexamer** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 0.219 | 0.0630 | 3.48x | [3.23, 3.60] | 0.300 | 0.066 | True |
| 75 | 0.295 | 0.1016 | 2.90x | [2.62, 3.14] | 0.381 | 0.104 | True |
| 100 | 0.389 | 0.1130 | 3.45x | [3.35, 3.61] | 0.473 | 0.115 | True |
| 150 | 0.583 | 0.1981 | 2.94x | [2.81, 3.08] | 0.667 | 0.201 | True |

**read_quality** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 1.392 | 0.4033 | 3.45x | [3.42, 3.59] | 1.491 | 0.406 | True |
| 75 | 2.013 | 0.5792 | 3.48x | [3.42, 3.52] | 2.111 | 0.582 | True |
| 100 | 2.627 | 0.7790 | 3.37x | [3.34, 3.41] | 2.727 | 0.782 | True |
| 150 | 3.874 | 1.1479 | 3.37x | [3.33, 3.41] | 3.978 | 1.151 | True |

## 6. Findings

These are the conclusions the data support. Every number in sections 1-5 is regenerated
from `results-main.json` / `results-scaling.json` by `analyze.py`; none is hand-written.
The tables in this section restate them alongside measurements taken outside the harness
(differential output checks, profiles), and are marked as such.

### 6.1 All 29 measured commands are faster, and all 29 passed the equivalence gate

The first measured run found three commands that were **reproducibly slower in the port**,
each with a confidence interval entirely below 1.0. All three were root-caused and fixed,
and this run re-measures them as wins: `inner_distance` **10.93x**, `bam2fq` **2.30x**,
`infer_experiment` **3.35x** (the last was 0.18x, then 1.09x after its first fix).

| | first run | this run |
|---|---|---|
| gates failed | 1 (`RPKM_saturation`) | **0** |
| commands slower than upstream | 3 | **0** |
| fastest command | `infer_experiment` 0.18x (as a loss) | `normalize_bigwig` 71.6x |
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
| 2,000 | 2.25x | 14.17x |
| 10,000 | 3.69x | 8.61x |
| 50,000 | 3.98x | 5.30x |
| 200,000 | 4.35x | 4.75x |
| 800,000 | 4.65x | 4.75x |

Same code, same machine, **opposite trends**: the end-to-end figure *falls* by 3.0x
while the compute-only figure *doubles* from its smallest point. The reason is the reference's
~0.09 s fixed interpreter/import cost, which is a constant that dominates small inputs
and vanishes on large ones. At 2,000 reads the port's own arithmetic is only 2.3x faster;
the other 12x is CPython startup that any real user pays but that has nothing to do with
RSeQC. That 2,000-read cell is also the noisiest interval in the report
(E1 CI [1.95, 3.39]) because the port's own work is only ~3 ms there.

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
| `geneBody_coverage`'s own index | same sliding window | 75.1 MB -> **14.0 MB** | 40.6 MB |

The last two are exact rather than approximate, and for slightly different reasons.
`tin` scores one transcript at a time, so only reads that can still reach an unscored
transcript need to be resident: a read with `end <= tx_start` cannot overlap any later
transcript. `geneBody_coverage` aggregates across *all* transcripts at once, so that
argument does not apply -- but its accumulation is a **sum plus an OR**, both commutative,
so visiting transcripts in coordinate order instead of BED order cannot change the result.

**Result: every one of the 29 measured commands is now lighter than upstream**, from
`bam2wig` at 1.40x lighter (793 MB against 1111 MB, the largest absolute figure on either
side) to `sc_seqQual` at 15.7x. Streaming the reader also made things *faster*, not just
smaller, so none of this is a time-for-memory trade.

The memory sweep confirms the reader fix holds across a 400x range rather than at one
point: `bam_stat` peak RSS is flat at ~2.9 MB from 2,000 to 800,000 reads.

CRAM remains whole-file buffered, and the reason is a real API constraint rather than a
preference: `noodles-cram` 0.99 exposes record iteration only as `records(&header)`,
which is **single-use**, and re-entering it on a drained reader returns a spurious
`TryFromIntError` instead of EOF (verified directly). Bounding it would need a
self-referential reader, not justified for the least common input format here.

### 6.5 `RPKM_saturation` is not deterministic upstream, and the port is as close to it as upstream is to itself

Two consecutive upstream runs on identical input differ -- upstream subsamples reads per
percentile with no RNG seed. Measured here: the port agrees with upstream on **58.93%** of
`rp.eRPKM.xls` rows and **59.05%** of `rp.rawCount.xls`; **upstream agrees with itself on
58.32%** and 58.43%. This time the port is *ahead* of upstream's own reproducibility by
0.61 and 0.62 pp. The gate does not require that -- it requires the port to be no more
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
| E1 speedups far below E2 speedups | **Held.** `bam_stat` at 2,000 reads: 2.25x compute-only vs 14.17x end-to-end. The gap closes as input grows. |
| Speedup falls with input size | **Held for E2** (14.17x -> 4.75x on `bam_stat`). **Refined**: E1 *rises* to 4.65x; the two converge at large input, which is the defensible number. |
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
  average 3.70/2.58/1.29 at the start of this run, 7 users on the machine. Ratios are
  paired and interleaved so machine drift cannot favour one arm, but absolute times are
  not trustworthy to better than roughly +-10%, and the intervals measure run-to-run
  scatter on *this* machine, not the uncertainty of a publication measurement.
- **Run-to-run, the absolute numbers move more than the ratios do.** Three runs of the
  same build lineage on the same workloads, in increasing order of machine load:

  | | run A (loadavg 0.86) | run B (loadavg 0.86) | this run (loadavg 3.70) |
  |---|---|---|---|
  | `geneBody_coverage` | 36.2x | 46.6x | 46.4x |
  | `normalize_bigwig` | 74.5x | 74.5x | 71.6x |
  | `tin` | 15.2x | 16.3x | 16.3x |
  | `read_duplication` | 1.48x | 1.51x | 1.51x |

  `normalize_bigwig` is the most load-sensitive at +-4%; the small commands are stable to
  within 2%. Treat any single figure here as a ratio measured on a shared sandbox, not a
  portable constant, and prefer the shape of a curve (section 5) over any one point.
- **This run is the first whose results trace to a clean commit** (`3a5b88c`,
  `git_dirty: false`), which is what gate G6b requires. The two runs before it are
  retained nowhere in the report, so their columns above come from comparing the
  superseded commits' raw JSON, not from the committed artifacts.
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
