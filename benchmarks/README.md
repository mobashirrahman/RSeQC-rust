# Benchmark Runner

This directory contains scaffolding to measure and compare performance between the upstream Python RSeQC implementation and the Rust port.

**Important:** Results from this sandbox are **not** suitable for publication. All measurements here run on shared, non-isolated hardware. See [testing.md](../testing.md) section 12 for the protocol required for publishable performance claims.

## Overview

The benchmark runner implements section 12 of [testing.md](../testing.md), which specifies:

1. **Deterministic workload generation** (`generate_workload.py`)
   - Synthetic paired-end RNA-seq BAM with coordinate sorting and indexing
   - Matching BED12 gene model
   - FASTQ files
   - Seeded reproducibility

2. **Command execution and timing** (`run_benchmarks.py`)
   - Runs each command with both Python and Rust implementations
   - Captures wall time, CPU time, and peak memory
   - Records environment (CPU, kernel, versions)
   - Verifies output compatibility

3. **Analysis** (section 12.2 and 12.3)
   - Raw paired measurements
   - Median wall times and spread
   - Bootstrap confidence intervals (requires post-processing)
   - Separate scaling studies (independent sweeps of read count, depth, transcript count)

## Quick Start

### Generate a small synthetic workload

```bash
cd /path/to/RSeQC-rust

# Create a workload with 1000 paired-end reads
oracle/venv/bin/python3 benchmarks/generate_workload.py \
    --size 1000 \
    --output-dir workloads/test_1000
```

This creates:
- `workloads/test_1000/reads.bam` — coordinate-sorted BAM
- `workloads/test_1000/reads.bam.bai` — BAM index
- `workloads/test_1000/model.bed12` — gene annotation
- `workloads/test_1000/reads_1.fastq`, `reads_2.fastq` — paired FASTQ
- `workloads/test_1000/manifest.json` — seed, version, parameters

### Run benchmarks

```bash
# Run a 3-repetition benchmark on the 1000-read workload
oracle/venv/bin/python3 benchmarks/run_benchmarks.py \
    --workload workloads/test_1000 \
    --commands bam_stat read_distribution geneBody_coverage \
    --reps 3 \
    --output-dir benchmark_results/test_1000
```

This produces:
- `benchmark_results/test_1000/results.json` — raw per-run measurements and environment info
- Console summary with median times and speedups

## File Format Reference

### Manifest (workload)

```json
{
  "generator_version": "1.0.0",
  "seed": 42,
  "num_reads": 1000,
  "files": {
    "bam": "reads.bam",
    "bam_index": "reads.bam.bai",
    "bed12": "model.bed12",
    "fastq_1": "reads_1.fastq",
    "fastq_2": "reads_2.fastq"
  },
  "metadata": {
    "seed": 42,
    "num_reads": 1000,
    "fragment_size": 100,
    "read_length": 100,
    "chrom": "chr1",
    "chrom_length": 10000,
    "genes": [...]
  }
}
```

### Results

Raw results in `results.json`:

```json
{
  "environment": {
    "timestamp": 1234567890.123,
    "platform": "Linux-6.14.0-37-generic-x86_64-with-glibc2.2.5",
    "cpu_count": 16,
    "processor": "Intel(R) Xeon(R)",
    "git_commit": "abc123def456...",
    "rustc_version": "rustc 1.81.0",
    "python_version": "Python 3.10.0",
    "loadavg": [0.5, 0.4, 0.3],
    "warning": "Results from non-isolated hardware..."
  },
  "benchmark_results": [
    {
      "command": "bam_stat",
      "workload": "workloads/test_1000",
      "repetitions": 3,
      "python": {
        "runs": [
          {
            "repetition": 1,
            "exit_code": 0,
            "wall_time_s": 0.123,
            "cpu_time_s": 0.120,
            "peak_rss_mb": 45.2,
            "timed_out": false
          },
          ...
        ]
      },
      "rust": {
        "runs": [...]
      },
      "compatibility": {
        "match": true,
        "errors": []
      }
    }
  ]
}
```

## Scaling Studies

To satisfy section 12.3, run independent sweeps along each cost driver:

```bash
# Scaling by read count
for size in 100 500 1000 5000 10000; do
    oracle/venv/bin/python3 benchmarks/generate_workload.py \
        --size $size \
        --output-dir workloads/scale_reads_$size
    
    oracle/venv/bin/python3 benchmarks/run_benchmarks.py \
        --workload workloads/scale_reads_$size \
        --commands bam_stat read_distribution geneBody_coverage \
        --reps 10 \
        --output-dir benchmark_results/scale_reads_$size
done
```

Aggregate these results and produce plots showing each command's cost curve.

## Publication Requirements (testing.md §12)

A performance claim is publishable only after:

### §12.1 Protocol frozen before benchmarking
- [ ] Each command/mode and input size defined
- [ ] All required outputs specified
- [ ] Setup/indexing costs identified separately
- [ ] Threading/compression settings locked
- [ ] Environment fully recorded (see `environment` above)
- [ ] Wall time, CPU, peak memory measured for whole process tree
- [ ] Upstream setup cost included; Rust indexing cost included

### §12.2 Repetitions and analysis
- [ ] At least 10 paired repetitions for primary deterministic cases
- [ ] Randomized interleaving (not back-to-back)
- [ ] Raw measurements recorded and published
- [ ] Median times and 95% CI computed
- [ ] Speedup estimator clearly defined (e.g., exponentiated mean of log ratios)
- [ ] Bootstrap or resampling methodology named

### §12.3 Scaling studies
- [ ] Independent sweeps of read count, covered bases, depth, transcript count
- [ ] Maximum operating point and failure mode documented
- [ ] Semantic consumers retested after each scale increase

### Hardware isolation
- [ ] Measurements on isolated hardware (not this sandbox)
- [ ] Documented CPU model, core count, NUMA, cache state
- [ ] Power/frequency management locked
- [ ] No concurrent load during measurement
- [ ] Separate data-only, compute-only, and end-to-end measurements

### Output compatibility
- [ ] All commands produce identical outputs to upstream (or documented divergence)
- [ ] Outputs verified across all tested scales

## Gitignore

All benchmark outputs (workloads and results) are gitignored:

```
benchmarks/workloads/
benchmarks/benchmark_results/
```

Commit the runner code, not the data or numbers.

## Next Steps

1. **Extend generator**: Add support for spliced reads, different strand compositions, depth variation, and other features from section 6.2
2. **Add comparators**: Implement BAM/BigWig structural comparison (not just exit codes)
3. **Implement analysis**: Bootstrap confidence intervals and speedup ratio computation
4. **Scale sweeps**: Implement cost-driver sweeps as described in §12.3
5. **Real data validation**: Add support for real dataset benchmarks
6. **Isolated hardware**: Run on isolated hardware per §12.1 for publishable results

## References

- [testing.md §12 — Fair performance and scalability measurements](../testing.md#12-fair-performance-and-scalability-measurements)
- [testing.md §5 — Verification harness requirements](../testing.md#5-make-the-verification-harness-capable-of-failing-correctly)
