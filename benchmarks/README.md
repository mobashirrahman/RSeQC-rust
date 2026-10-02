# Benchmark suite

> **Status.** The version-1 results in [`RESULTS.generated.md`](RESULTS.generated.md)
> are **historical development evidence**, retained for continuity. The 2026-10-01
> readiness audit found comparator false-pass paths, lost repetition pairing, and an
> unsafe timeout process-group fallback in the harness that produced them, so those
> gate passes and paired intervals do not qualify any publication claim. All of those
> defects are now closed and tested; the replacement study is specified in
> [`protocol-v2.md`](protocol-v2.md) and has not yet been run at scale.
>
> See the [audit and executed probes](../docs/READINESS_AUDIT_2026-10-01.md) and
> [`CHANGES.md`](CHANGES.md).

**Version-1 results, as measured.** Results from a full 10-repetition,
gate-checked run are in [`RESULTS.generated.md`](RESULTS.generated.md), regenerated
from [`results-main.json`](results-main.json) and
[`results-scaling.json`](results-scaling.json) by `analyze.py`. The protocol those
runs followed is frozen in [`protocol.md`](protocol.md).

> **Not publication-grade.** Measured on shared, non-isolated hardware, with a
> harness whose comparators and pairing have since been found defective. See
> [`protocol.md` §2](protocol.md), RESULTS §7, and `protocol-v2.md`.

## What is here

| File | Role |
|---|---|
| `protocol-v2.md` | **The next study's frozen protocol.** Repaired comparators, declared per-command expectations, adjacent matched blocks, per-run provenance. |
| `protocol.md` | Version-1 preregistration, superseded. Retained unchanged: its recorded results are only interpretable against the protocol that produced them. |
| `CHANGES.md` | Every post-freeze change, split into protocol clarifications and harness bugs whose results were discarded and re-measured. |
| `generate_workload_real.py` | Tier A/B workload generator: real hg38 chromosomes + real RefSeq BED12, reads simulated from real transcript sequences. Validates every generated BAM against the reference before writing a manifest. |
| `generate_workload.py` | Tier C legacy generator (1 contig, 3 genes). Retained for pinning old numbers; **not used for any claim**. |
| `bench.py` | The harness. Process-tree resources, structural equivalence gate over declared streams and artifacts, adjacent matched-block repetitions, block-bootstrap intervals, per-run provenance. |
| `derive_expected_streams.py` | Runs both arms once per command and records what they actually emit, as a draft for `EXPECTED_STREAMS`. Observes behaviour; does not certify correctness. |
| `test_bench_harness.py` | Harness-credibility tests. Every test asserts that a deliberate defect **fails**, so a comparator that stops rejecting corruption turns the suite red. |
| `scaling.py` | Cost-driver sweeps (read count, transcript count, read length). |
| `analyze.py` | Regenerates every table in RESULTS.generated.md from the raw JSON. |
| `run_benchmarks.py` | **Superseded** by `bench.py`. Kept for reference; it covered 5 commands and compared exit codes plus a raw file diff. |
| `results-main.json` | Version-1 raw per-run measurements. Comparator version 1; not comparable with version-2 rows. |
| `results-scaling.json` | Version-1 raw per-run measurements for each scaling point. |

Comparators are shared with the differential runner
([`verification/comparators.py`](../verification/comparators.py)) so the two
harnesses cannot drift apart again — which is how the benchmark's comparators came
to accept corrupted BAM qualities, flags and tags while the differential suite had
its own weaker copies.

## Headline

> Every figure below was measured under comparator version 1 and a schedule whose
> pairs were not the pairs it claimed. The *relative* orderings are still
> informative; the gate passes and intervals are not evidence of equivalence.

- **All 29 measured commands were faster, and all 29 passed the version-1
  equivalence gate.** The first run found three that were not (`infer_experiment`
  0.18x, `bam2fq` 0.35x, `inner_distance` 0.53x); all three were root-caused, fixed
  and re-measured at 3.35x, 2.30x and 10.93x. See RESULTS §6.1 — notably the three
  shared a *symptom* but had three *different* causes.
- Speedup is **not** a single number. For `bam_stat` the end-to-end figure falls from
  14.2x at 2k reads to 4.8x at 800k while the compute-only figure *rises* from 2.3x to
  4.7x -- opposite trends on the same code, because ~0.09 s of the reference's time is
  fixed interpreter/import cost. Any claim must state its workload size and whether it
  includes interpreter startup.
- The port uses **less** memory than upstream on **all 29** commands. The first run found
  the opposite on 14 commands (whole-file eager decode in `open_alignments`); that reader
  now streams, and the two whole-file per-read indexes left over (`tin`'s and
  `geneBody_coverage`'s) have been replaced with a sliding window (366 MB -> 20 MB and
  75 MB -> 14 MB).
- Upstream burns ~0.9 s of CPU per invocation in OpenBLAS thread-pool start-up with no
  `multiprocessing` anywhere, so thread env vars are pinned in both arms or the CPU
  column is meaningless.

See RESULTS.generated.md §6 for the full findings and §7 for what is not established.

## Overview

The benchmark runner implements section 12 of [testing.md](../testing.md):

1. **Workload generation** (`generate_workload_real.py`, `generate_workload.py`)
   - Real hg38 chromosomes and RefSeq BED12, reads simulated from real transcript
     sequences; validated against the reference before a manifest is written
   - Legacy synthetic generator retained only for pinning old numbers

2. **Command execution and measurement** (`bench.py`)
   - Each command runs in both arms as a real process invocation
   - Captures wall, user and system CPU, and peak memory
   - Records the environment, the toolchain, and per-run provenance including a
     SHA256 of every binary invoked
   - Gates equivalence on declared streams, labels and artifacts, per command

3. **Scheduling and statistics**
   - Adjacent matched blocks: repetition *i*'s two arms run back to back, with the
     within-block order randomised from a recorded seed
   - Block bootstrap over whole pairs, seed recorded
   - Predeclared precision rule for adding repetitions
   - Separate cost-driver sweeps, labelling measured rather than requested counts

The version-1 harness globally shuffled individual arms and then zipped two
independently ordered result lists, so at seed 20260929 all ten "pairs" had
different repetition IDs. The intervals it reported described unpaired
measurements. That is fixed; see [`protocol-v2.md`](protocol-v2.md) §7.1.

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
# Check the oracle environment matches its lock before trusting any result
oracle/venv/bin/python3 verification/check_oracle_env.py

# Run a 3-repetition benchmark on a workload, with the gate enabled
oracle/venv/bin/python3 benchmarks/bench.py \
    --workload workloads/test_1000 \
    --commands bam_stat read_distribution geneBody_coverage \
    --reps 3 \
    --label smoke \
    --output-dir benchmark_results/test_1000
```

This produces:
- `benchmark_results/test_1000/results.json` — raw per-run measurements, per-run
  provenance, gate outcomes and failures
- `benchmark_results/test_1000/<command>.json` — the same for one command
- A console table whose primary column is raw invocation time over adjacent matched
  blocks, with the number of usable pairs and any failed repetitions shown

To check the harness itself:

```bash
oracle/venv/bin/python3 benchmarks/test_bench_harness.py
```

Every test in that suite asserts that a deliberate defect **fails**. If it does not,
a comparator has stopped detecting corruption and every gate pass it has issued is
void.

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
for size in 100 500 1000 5000 10000; do
    oracle/venv/bin/python3 benchmarks/generate_workload.py \
        --size $size \
        --output-dir workloads/scale_reads_$size

    oracle/venv/bin/python3 benchmarks/bench.py \
        --workload workloads/scale_reads_$size \
        --commands bam_stat read_distribution geneBody_coverage \
        --reps 10 \
        --label "scale_reads_$size" \
        --output-dir benchmark_results/scale_reads_$size
done
```

Report **measured** counts, never requested ones. Version 1's "20,000 transcript"
point actually contained 9,179, because the panel was truncated to what was
available; a reader comparing that point against a different generator would have
been comparing different inputs.

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
