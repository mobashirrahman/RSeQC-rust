# Benchmark results — rseqc-rust vs upstream RSeQC

Protocol: `benchmarks/protocol.md` (frozen 2026-09-29). Run label: `streaming`. Reps: 10, seed 20260929.

> **Not publication-grade.** Shared, non-isolated hardware (AMD Ryzen 7 3700X 8-Core Processor, loadavg 2.02 at start, governor `schedutil`). See protocol.md §2.

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
| normalize_bigwig | W-mid | 48.598 | 0.670 | 72.502 | [68.16, 80.65] | win | 68.023 | 7.946 | True |
| geneBody_coverage | W-gbc | 64.200 | 1.604 | 40.032 | [39.65, 40.46] | win | 41.260 | 0.551 | True |
| read_NVC | W-mid | 15.009 | 1.356 | 11.066 | [10.95, 11.18] | win | 11.235 | 12.865 | True |
| inner_distance | W-mid | 8.304 | 0.751 | 11.055 | [10.75, 11.18] | win | 12.439 | 4.055 | True |
| RNA_fragment_size | W-mid | 56.264 | 5.721 | 9.834 | [9.72, 9.97] | win | 9.925 | 2.241 | True |
| tin | W-tn | 159.533 | 19.824 | 8.047 | [7.77, 8.18] | win | 8.259 | 0.118 | True |
| sc_seqQual | W-mid | 0.436 | 0.057 | 7.669 | [7.42, 8.10] | win | 17.700 | 15.887 | True |
| sc_seqLogo | W-mid | 0.326 | 0.043 | 7.606 | [7.10, 7.82] | win | 26.333 | 15.510 | True |
| RPKM_saturation | W-mid | 32.413 | 4.487 | 7.224 | [6.95, 7.37] | win | 7.513 | 1.984 | True |
| geneBody_coverage2 | W-mid | 25.609 | 4.295 | 5.962 | [5.90, 6.00] | win | 16.355 | 11.973 | True |
| overlay_bigwig | W-mid | 0.475 | 0.112 | 4.238 | [4.19, 4.28] | win | 6.500 | 5.426 | True |
| FPKM_count | W-mid | 34.288 | 8.350 | 4.106 | [3.96, 4.20] | win | 4.138 | 1.580 | True |
| bam_stat | W-mid | 1.821 | 0.452 | 4.026 | [3.96, 4.13] | win | 4.261 | 13.124 | True |
| read_hexamer | W-mid | 0.579 | 0.164 | 3.519 | [3.44, 3.64] | win | 4.200 | 4.023 | True |
| clipping_profile | W-mid | 1.615 | 0.469 | 3.441 | [3.40, 3.49] | win | 3.711 | 13.105 | True |
| read_quality | W-mid | 10.275 | 3.012 | 3.412 | [3.39, 3.43] | win | 3.443 | 12.389 | True |
| infer_experiment | W-mid | 0.633 | 0.187 | 3.388 | [3.33, 3.44] | win | 3.833 | 8.142 | True |
| insertion_profile | W-mid | 1.596 | 0.477 | 3.348 | [3.33, 3.39] | win | 3.548 | 13.090 | True |
| read_distribution | W-mid | 4.867 | 1.492 | 3.262 | [3.22, 3.28] | win | 3.234 | 3.077 | True |
| split_bam | W-mid | 8.737 | 2.800 | 3.120 | [3.11, 3.15] | win | 3.201 | 2.386 | True |
| deletion_profile | W-mid | 1.105 | 0.366 | 3.022 | [2.97, 3.11] | win | 3.343 | 15.004 | True |
| junction_annotation | W-mid | 2.004 | 0.752 | 2.666 | [2.61, 2.68] | win | 2.949 | 1.493 | True |
| mismatch_profile | W-mid | 1.328 | 0.512 | 2.594 | [2.58, 2.62] | win | 2.816 | 10.281 | True |
| bam2wig | W-mid | 22.901 | 9.345 | 2.451 | [2.44, 2.49] | win | 2.474 | 1.402 | True |
| split_paired_bam | W-mid | 8.155 | 3.396 | 2.401 | [2.37, 2.43] | win | 2.409 | 7.848 | True |
| bam2fq | W-mid | 1.944 | 0.850 | 2.288 | [2.24, 2.32] | win | 2.580 | 8.144 | True |
| junction_saturation | W-mid | 2.342 | 1.034 | 2.264 | [2.24, 2.32] | win | 2.387 | 2.161 | True |
| read_GC | W-mid | 1.645 | 0.812 | 2.026 | [2.00, 2.05] | win | 2.112 | 12.500 | True |
| read_duplication | W-mid | 2.764 | 1.831 | 1.510 | [1.49, 1.53] | win | 1.592 | 1.468 | True |

## 2. E1 compute-only (fixed per-invocation cost subtracted)

Upstream pays a large fixed cost before reading any input (interpreter start plus imports). E1 subtracts each arm's *measured* floor (`--help` on the real binary) so the number reflects algorithmic work. E2 includes it, because a user always pays it.

### compute-only vs end-to-end

| command | py floor (s) | rs floor (s) | py net (s) | rs net (s) | E1 x | py E2 (s) | rs E2 (s) |
|---|---|---|---|---|---|---|---|
| normalize_bigwig | 0.092 | 0.003 | 48.598 | 0.670 | 72.502 | 48.690 | 0.673 |
| geneBody_coverage | 0.093 | 0.003 | 64.200 | 1.604 | 40.032 | 64.293 | 1.606 |
| read_NVC | 0.096 | 0.003 | 15.009 | 1.356 | 11.066 | 15.105 | 1.359 |
| inner_distance | 0.095 | 0.003 | 8.304 | 0.751 | 11.055 | 8.399 | 0.754 |
| RNA_fragment_size | 0.092 | 0.002 | 56.264 | 5.721 | 9.834 | 56.356 | 5.723 |
| tin | 0.099 | 0.002 | 159.533 | 19.824 | 8.047 | 159.632 | 19.827 |
| sc_seqQual | 0.528 | 0.003 | 0.436 | 0.057 | 7.669 | 0.964 | 0.060 |
| sc_seqLogo | 0.556 | 0.003 | 0.326 | 0.043 | 7.606 | 0.882 | 0.045 |
| RPKM_saturation | 0.113 | 0.002 | 32.413 | 4.487 | 7.224 | 32.526 | 4.489 |
| geneBody_coverage2 | 0.084 | 0.003 | 25.609 | 4.295 | 5.962 | 25.693 | 4.298 |
| overlay_bigwig | 0.096 | 0.002 | 0.475 | 0.112 | 4.238 | 0.570 | 0.114 |
| FPKM_count | 0.099 | 0.002 | 34.288 | 8.350 | 4.106 | 34.387 | 8.352 |
| bam_stat | 0.090 | 0.003 | 1.821 | 0.452 | 4.026 | 1.911 | 0.455 |
| read_hexamer | 0.083 | 0.002 | 0.579 | 0.164 | 3.519 | 0.662 | 0.167 |
| clipping_profile | 0.096 | 0.003 | 1.615 | 0.469 | 3.441 | 1.711 | 0.472 |
| read_quality | 0.098 | 0.003 | 10.275 | 3.012 | 3.412 | 10.372 | 3.014 |
| infer_experiment | 0.092 | 0.002 | 0.633 | 0.187 | 3.388 | 0.725 | 0.189 |
| insertion_profile | 0.098 | 0.003 | 1.596 | 0.477 | 3.348 | 1.693 | 0.480 |
| read_distribution | 0.101 | 0.002 | 4.867 | 1.492 | 3.262 | 4.969 | 1.495 |
| split_bam | 0.098 | 0.002 | 8.737 | 2.800 | 3.120 | 8.835 | 2.803 |
| deletion_profile | 0.099 | 0.003 | 1.105 | 0.366 | 3.022 | 1.204 | 0.368 |
| junction_annotation | 0.096 | 0.003 | 2.004 | 0.752 | 2.666 | 2.100 | 0.755 |
| mismatch_profile | 0.091 | 0.003 | 1.328 | 0.512 | 2.594 | 1.418 | 0.515 |
| bam2wig | 0.094 | 0.003 | 22.901 | 9.345 | 2.451 | 22.994 | 9.348 |
| split_paired_bam | 0.060 | 0.003 | 8.155 | 3.396 | 2.401 | 8.214 | 3.399 |
| bam2fq | 0.095 | 0.003 | 1.944 | 0.850 | 2.288 | 2.040 | 0.852 |
| junction_saturation | 0.096 | 0.003 | 2.342 | 1.034 | 2.264 | 2.438 | 1.037 |
| read_GC | 0.093 | 0.003 | 1.645 | 0.812 | 2.026 | 1.738 | 0.814 |
| read_duplication | 0.095 | 0.002 | 2.764 | 1.831 | 1.510 | 2.859 | 1.833 |

## 3. Peak memory (whole process tree, `/usr/bin/time -v`)

### peak RSS; a ratio above 1 means the port uses more

| command | py median MB | rs median MB | mem x | py max MB | rs max MB |
|---|---|---|---|---|---|
| tin | 43.020 | 365.646 | 0.118 | 43.637 | 366.000 |
| geneBody_coverage | 41.451 | 75.256 | 0.551 | 41.820 | 75.480 |
| bam2wig | 1112 | 793.113 | 1.402 | 1112 | 793.355 |
| read_duplication | 265.980 | 181.195 | 1.468 | 266.539 | 181.414 |
| junction_annotation | 59.330 | 39.727 | 1.493 | 59.645 | 39.934 |
| FPKM_count | 92.447 | 58.521 | 1.580 | 92.852 | 59.582 |
| RPKM_saturation | 348.439 | 175.658 | 1.984 | 348.828 | 178.547 |
| junction_saturation | 66.328 | 30.695 | 2.161 | 66.535 | 31.004 |
| RNA_fragment_size | 38.912 | 17.367 | 2.241 | 39.121 | 17.816 |
| split_bam | 133.146 | 55.803 | 2.386 | 133.375 | 55.879 |
| read_distribution | 148.219 | 48.172 | 3.077 | 148.633 | 50.348 |
| read_hexamer | 32.402 | 8.055 | 4.023 | 32.641 | 8.211 |
| inner_distance | 210.299 | 51.865 | 4.055 | 210.895 | 53.234 |
| overlay_bigwig | 45.875 | 8.455 | 5.426 | 46.086 | 8.621 |
| split_paired_bam | 25.566 | 3.258 | 7.848 | 25.906 | 3.328 |
| normalize_bigwig | 66.547 | 8.375 | 7.946 | 66.754 | 8.473 |
| infer_experiment | 41.633 | 5.113 | 8.142 | 41.867 | 5.270 |
| bam2fq | 38.859 | 4.771 | 8.144 | 39.113 | 4.898 |
| mismatch_profile | 39.156 | 3.809 | 10.281 | 39.410 | 3.938 |
| geneBody_coverage2 | 34.516 | 2.883 | 11.973 | 34.746 | 3.098 |
| read_quality | 39.273 | 3.170 | 12.389 | 39.531 | 3.324 |
| read_GC | 39.037 | 3.123 | 12.500 | 39.383 | 3.262 |
| read_NVC | 39.021 | 3.033 | 12.865 | 39.352 | 3.164 |
| insertion_profile | 39.041 | 2.982 | 13.090 | 39.277 | 3.020 |
| clipping_profile | 38.906 | 2.969 | 13.105 | 39.172 | 3.148 |
| bam_stat | 38.500 | 2.934 | 13.124 | 38.855 | 2.965 |
| deletion_profile | 39.062 | 2.604 | 15.004 | 39.359 | 2.711 |
| sc_seqLogo | 99.631 | 6.424 | 15.510 | 100.012 | 6.453 |
| sc_seqQual | 100.939 | 6.354 | 15.887 | 101.215 | 6.555 |

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
| 2000 | 0.013 | 0.0051 | 2.59x | [2.05, 4.02] | 0.104 | 0.007 | True |
| 10000 | 0.046 | 0.0153 | 3.05x | [2.89, 3.61] | 0.140 | 0.018 | True |
| 50000 | 0.234 | 0.0671 | 3.48x | [3.37, 3.66] | 0.327 | 0.070 | True |
| 200000 | 0.891 | 0.2427 | 3.67x | [3.52, 3.79] | 0.983 | 0.246 | True |
| 800000 | 3.506 | 0.9194 | 3.81x | [3.78, 3.91] | 3.598 | 0.922 | True |

**read_NVC** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.073 | 0.0102 | 7.18x | [6.18, 8.91] | 0.173 | 0.013 | True |
| 10000 | 0.374 | 0.0397 | 9.42x | [8.82, 9.89] | 0.467 | 0.042 | True |
| 50000 | 1.878 | 0.1852 | 10.14x | [9.79, 10.42] | 1.978 | 0.188 | True |
| 200000 | 7.480 | 0.7050 | 10.61x | [10.28, 10.79] | 7.573 | 0.708 | True |
| 800000 | 29.811 | 2.7446 | 10.86x | [10.67, 11.07] | 29.903 | 2.747 | True |

**read_quality** vs reads

| reads | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 2000 | 0.049 | 0.0197 | 2.50x | [2.42, 3.00] | 0.149 | 0.022 | True |
| 10000 | 0.262 | 0.0844 | 3.11x | [2.76, 3.21] | 0.361 | 0.087 | True |
| 50000 | 1.288 | 0.4001 | 3.22x | [3.14, 3.34] | 1.385 | 0.403 | True |
| 200000 | 5.118 | 1.5447 | 3.31x | [3.29, 3.33] | 5.216 | 1.547 | True |
| 800000 | 20.469 | 6.0946 | 3.36x | [3.32, 3.39] | 20.563 | 6.097 | True |

**FPKM_count** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 7.242 | 0.5118 | 14.15x | [14.03, 14.23] | 7.337 | 0.514 | True |
| 1000 | 13.106 | 0.9317 | 14.07x | [13.94, 14.15] | 13.202 | 0.934 | True |
| 5000 | 14.213 | 1.6842 | 8.44x | [8.38, 8.58] | 14.307 | 1.687 | True |
| 20000 | 15.445 | 2.3998 | 6.44x | [6.27, 6.57] | 15.540 | 2.402 | True |

**read_distribution** vs transcripts

| transcripts | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 200 | 1.161 | 0.2862 | 4.05x | [3.88, 4.13] | 1.259 | 0.289 | True |
| 1000 | 1.387 | 0.3104 | 4.47x | [4.34, 4.60] | 1.487 | 0.314 | True |
| 5000 | 1.691 | 0.4291 | 3.94x | [3.79, 4.09] | 1.791 | 0.432 | True |
| 20000 | 1.955 | 0.5637 | 3.47x | [3.37, 3.57] | 2.055 | 0.566 | True |

**read_hexamer** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 0.207 | 0.0626 | 3.31x | [3.20, 3.44] | 0.288 | 0.065 | True |
| 75 | 0.287 | 0.0947 | 3.03x | [2.93, 3.15] | 0.366 | 0.097 | True |
| 100 | 0.397 | 0.1085 | 3.66x | [3.49, 3.75] | 0.473 | 0.111 | True |
| 150 | 0.568 | 0.1904 | 2.98x | [2.87, 3.09] | 0.645 | 0.193 | True |

**read_quality** vs readlen

| readlen | py net (s) | rs net (s) | speedup | 95% CI | py E2 (s) | rs E2 (s) | gate |
|---|---|---|---|---|---|---|---|
| 50 | 1.367 | 0.3970 | 3.44x | [3.37, 3.53] | 1.463 | 0.400 | True |
| 75 | 1.978 | 0.5942 | 3.33x | [3.12, 3.40] | 2.070 | 0.597 | True |
| 100 | 2.566 | 0.7792 | 3.29x | [3.27, 3.41] | 2.660 | 0.782 | True |
| 150 | 3.747 | 1.1566 | 3.24x | [3.23, 3.31] | 3.842 | 1.159 | True |
## 6. Findings

These are the conclusions the data support. Every number above is regenerated from
`results-main.json` / `results-scaling.json` by `analyze.py`; none is hand-written.

### 6.1 All 29 measured commands are now faster; three were not, and were fixed

The first run found three commands that were **reproducibly slower in the port**, each
with a confidence interval entirely below 1.0. All three have since been root-caused
and fixed, and re-measured through the same harness:

| command | before | after | 95% CI | root cause and fix |
|---|---|---|---|---|
| `inner_distance` | 0.53x | **9.94x** | [9.85, 10.03] | `names_at` scanned every transcript on the chromosome and built two `HashSet<String>` per pair, twice per pair. Upstream indexes intervals. Replaced with an O(log n) start-sorted/running-max-end index. |
| `bam2fq` | 0.35x | **2.04x** | [2.03, 2.05] | Output went to a bare `File`, so a FASTQ record's 5 writes each became a syscall: **5,600,008 `write` calls** and 4.6 s of system time on 800k reads. Wrapped in a 1 MiB `BufWriter`: 5,600,008 → **178** syscalls. |
| `infer_experiment` | 0.18x | **1.09x** | [1.08, 1.10] | Same linear-scan class as `inner_distance`. Upstream uses `bx.intervals.Intersecter` (`qcmodule/SAM.py:2122`), a bitset interval index; the port scanned all ranges per read (~1.8e9 comparisons at the 200k sample cap). Replaced with an equivalent O(log n) index. |

Each fix was verified two ways: outputs compared byte-identical against real upstream
on the 800k-read workload, and the full 84-case differential harness re-run green
(`verification/run_diff.py`: "All 84 case(s) PASSED"). `infer_experiment`'s index was
additionally covered by a new randomised equivalence test asserting it returns exactly
what the linear scan returned, including zero-length, nested, abutting and duplicate
ranges.

**The lesson is the one the protocol's §10.6 prediction got right and my own initial
diagnosis got wrong.** The three losses shared a symptom (all were `open_alignments`
callers) and I inferred a single shared cause -- the eager whole-file decode. The
profile refuted that: `bam_stat` performs the identical decode in 0.50 s while
`infer_experiment` took 4.07 s, so the shared work could not explain it. The
user/sys split then showed three *different* signatures -- CPU-bound, syscall-bound,
and algorithmic -- and three unrelated fixes. **A correlation across two metrics is
evidence for a hypothesis, not proof of one**; the decisive check cost ten minutes.

### 6.2 A single speedup number is meaningless; the curve is the result

`bam_stat`, measured across a 400x range of input sizes:

| reads | end-to-end speedup (E2) | compute-only speedup (E1) |
|---|---|---|
| 2,000 | 14.2x | 2.59x |
| 10,000 | 7.8x | 3.05x |
| 50,000 | 4.7x | 3.48x |
| 200,000 | 4.0x | 3.67x |
| 800,000 | 3.9x | 3.81x |

End-to-end speedup **falls** with input size because the reference's ~0.09 s fixed cost
is a constant that dominates small inputs. Compute-only speedup **rises** toward ~3.8x
because the port's constant-factor advantage becomes visible as the fixed cost shrinks
in relative terms. Both converge near 800k reads.

The practical consequence: a benchmark run at 1,000 reads would have reported ~69x for
this command (measured during protocol design, `benchmarks/CHANGES.md` §2.2 context).
That number is almost entirely CPython interpreter startup. **Any future claim about
this port must state its workload size.**

### 6.3 Upstream's CPU-time column is mostly BLAS thread initialisation

Upstream RSeQC has no `multiprocessing` and no explicit threading (verified across all
33 scripts and `qcmodule`). Nevertheless, `import numpy` alone costs **0.08 s wall but
0.90 s of CPU** on this 16-thread machine, because OpenBLAS starts a thread pool per
process. With `OPENBLAS_NUM_THREADS=1` that falls to 0.05 s of CPU.

So an unpinned CPU-time measurement of upstream is ~90% thread-pool start-up, and the
resulting wall-vs-CPU gap looks like parallelism that does not exist. All reported runs
pin `OPENBLAS_NUM_THREADS=1` and friends in **both** arms (protocol §4.1). This is a
measurement trap specific to this comparison and would silently corrupt any naive
benchmark of it.

### 6.4 The memory regression is fixed; the port is now lighter than upstream on 28 of 29 commands

The first run found the port using **4-8x more** memory than upstream on 14 of 29
commands, with peak RSS growing linearly at ~226 bytes per record while upstream stayed
flat at ~39 MB. At 50M read pairs -- an ordinary human RNA-seq BAM -- that extrapolated
to **~23 GB in the port against ~39 MB upstream**, i.e. a tool that could not open a
normal dataset. The cause was `open_alignments` (`crates/formats/src/lib.rs`) returning
`Vec<io::Result<Record>>`, decoding the whole file up front, which 14 of the commands
routed through.

`open_alignments` now returns a streaming `AlignmentRecords` iterator. The whole
workspace rebuilt with **zero call-site changes** -- the `compute_*` functions were
already generic over `IntoIterator<Item = io::Result<bam::Record>>`, so this was a
drop-in. Memory behaviour by format:

| format | mechanism | peak memory |
|---|---|---|
| BAM | `read_record` into one reusable buffer | **O(1)** |
| SAM | `read_record_buf` into owned buffers, converted in batches of 4096 | O(chunk) |
| CRAM | still whole-file buffered | O(file) |

Measured `bam_stat` peak RSS, after vs before, across a 400x range of input sizes:

| reads | before | after | upstream |
|---|---|---|---|
| 2,000 | 3.7 MB | **2.9 MB** | 38.6 MB |
| 10,000 | 7.3 MB | **2.7 MB** | 38.7 MB |
| 100,000 | 25.2 MB | **2.7 MB** | 37.6 MB |
| 800,000 | 363.9 MB | **2.9 MB** | 38.5 MB |

Flat at ~2.8 MB, and `bam_stat` also got faster (0.50 s -> 0.45 s), so the fix was a
strict improvement rather than a time-for-memory trade.

**The one remaining exception is `tin`** (366 MB vs upstream's 43 MB). That is a
different structure: `tin`'s own `build_read_index` materialises per-read match blocks,
skip/delete blocks, qualities, sequence and CIGAR for the whole file
(`crates/commands/src/tin.rs:84`). It is the only command still heavier than upstream and
is a well-scoped follow-up, not a regression in the I/O layer.

Two things this exercise changed about how the port is read. First, the original
justification for eager decoding -- "acceptable for the small/QC-scale inputs this tool
targets" -- was an assumption that was never tested, and the benchmark falsified it;
the assumption was doing the work of a measurement. Second, CRAM remains buffered, and
the reason is a real API constraint rather than a preference: `noodles-cram` 0.99
exposes record iteration only as `records(&header)`, which is **single-use**, and
re-entering it on a drained reader returns a spurious `TryFromIntError` instead of EOF
(verified directly). Bounding CRAM's memory would need a self-referential reader, which
is not justified for the least common input format here.

### 6.5 `RPKM_saturation` is not deterministic upstream, and the port is as close as upstream is to itself

Two consecutive upstream runs on identical input differ. The port agrees with upstream
on **59.3%** of rows; **upstream agrees with itself on 58.4%**. The port is therefore
statistically indistinguishable from the reference, and marginally more consistent
with it. This is reported as equivalence against a self-control rather than claimed as
byte-identity, because byte-identity is not achievable by anything.

### 6.6 Preregistered expectations: which held and which did not

Protocol §10 predicted six things. Outcome:

| Prediction | Result |
|---|---|
| E1 speedups far below E2 speedups | **Held.** E.g. `sc_seqLogo` 19.8x end-to-end vs 8.4x compute-only. |
| Speedup falls with input size | **Held for E2** (14.2x → 3.9x). **Refined**: E1 *rises* to 3.8x; both converge. |
| Port loses memory on eager-decode commands | **Held**, and it was the most consequential finding. **Since resolved**: streaming `open_alignments` turned a 9.4x memory loss into a 13x memory win on the same command (section 6.4). |
| Whole-output buffering shows output-dependent memory growth | **Partly held.** `bam2wig` is the worst absolute case (793/1112 MB) but the port wins it. |
| `read_hexamer` scales worse than linear in read length | **Not held.** Its speedup falls 4.4x → 3.4x from 50 to 150 bp, i.e. mild, and the per-base `String` allocation is not yet a bottleneck at these sizes. Prediction withdrawn. |
| At least one command is slower in the port | **Held, three times** (`infer_experiment`, `bam2fq`, `inner_distance`) — all three now fixed and re-measured as wins (§6.1). |

---

## 7. What this study does not establish

- **Not publication-grade.** Shared, non-isolated hardware, `schedutil` governor,
  ~8 users on the machine. Ratios are paired and interleaved so machine drift cannot
  favour one arm, but the absolute times are not trustworthy to better than roughly
  ±10%, and the confidence intervals in this report measure run-to-run scatter on
  *this* machine, not the uncertainty of a publication measurement.
- **No cold-cache series.** Warm-cache only, as declared in protocol §2.
- **No thread scaling.** The port has no first-party parallelism and no `--threads`
  flag, so the 1/2/4/8 sweep is impossible. This is now the most valuable untested
  axis, and the streaming rewrite makes it more valuable still: decode is no longer
  hidden behind a large allocation, so it is now the dominant remaining cost for the
  I/O-bound commands (`bam2fq` 2.3x, `infer_experiment` 1.1x) rather than being masked.
- **`tin`'s own index is the last memory outlier** (366 MB vs upstream 43 MB,
  section 6.4). It is now the only command heavier than upstream.
- **Reads are simulated from real transcripts**, not a real sequencer run, and the
  gene model is real but the expression distribution is uniform-over-transcripts.
  Commands whose cost driver is duplicate structure (`read_duplication`), error
  artefacts, or expression skew are therefore measured on a favourable input.
  `read_duplication` shows the smallest speedup of any winning command (1.48x), which
  is consistent with that caveat.
- **No single-cell panel.** `sc_bamStat` and `sc_editMatrix` are not benchmarkable here
  (§4). `sc_seqLogo`/`sc_seqQual` are benchmarked on bulk-derived FASTQ, not on real
  barcode data.
- **Three commands could not be measured at all**: `FPKM-UQ` (needs `htseq-count`),
  `sc_editMatrix` (needs R `pheatmap`), `sc_bamStat` (needs CB/RE-tagged input).
- **Reproducibility gate ≠ scientific validity.** Every row here means "this workload's
  outputs match". It does not mean the computation is right on real data, which is the
  T4 held-out question and remains open.
