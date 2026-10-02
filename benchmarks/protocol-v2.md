# Benchmark protocol, version 2

**Status: FROZEN for the next study. Version 1's results are retained as
historical development evidence and are NOT publication claims.**

- Protocol version: 2.0.0
- Frozen: 2026-10-01
- Supersedes: `benchmarks/protocol.md` version 1.0.0 (2026-09-29), retained
  unchanged in git history; its recorded results stay in
  `benchmarks/results-main.json` and `benchmarks/RESULTS.generated.md`, labelled
  historical.
- Reference implementation: upstream RSeQC 5.0.5, commit
  `59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24`, in the executable environment pinned
  by `compatibility/upstream.lock` and verified by
  `verification/check_oracle_env.py`.
- Candidate implementation: this repository at a clean commit,
  `cargo build --workspace --release --locked`.
- Governing requirements: `testing.md` §12, `docs/PUBLICATION_PLAN.md`.

## Why version 2 exists

Version 1's harness had defects that make its own numbers uninterpretable, and the
2026-10-01 readiness audit found them by executing the harness rather than reading
it. None of these are "the port is slow" findings; all of them are "the measurement
could not detect a difference" findings.

| Defect | Consequence for version 1's numbers | Closure in this version |
|---|---|---|
| The output gate compared two file trees only. Two empty directories passed. | A command that produced nothing at all could receive a gate pass. | §6.2: per-command declared streams, labels, schemas and artifacts; emptiness fails. |
| Captured stdout was never compared. | Every stdout-only metric (`bam_stat`, `infer_experiment`, `read_distribution`, the split commands) was unchecked. | §6.2: stdout is part of the deliverable and is compared. |
| The BAM comparator read `flag & 0xC0` and omitted qualities, mate fields and tags. | Corrupted quality strings, duplicate flags, NM tags and mate positions compared equal. | §5.1: full record signature plus header reference dictionary and sort order. |
| The FASTQ reader stepped by four lines over `len - 3`. | A truncated final record was discarded and compared equal to a complete file. | §5.2: strict framing and sequence/quality length checks. |
| Text comparison normalised `nan`/`inf` like any other token. | A finite metric replaced by NaN compared equal in both directions. | §5.3: nonfinite values are rejected, including symmetric ones. |
| Only the last repetition's outputs were gated. | Nine of ten repetitions per row were unvalidated. | §6.3: every timed run's artifacts are validated, outside the timing interval. |
| The timeout handler called `killpg(getpgid(0))`, and `subprocess.run`'s `TimeoutExpired` carries no PID. | `getpgid(0)` returns the *caller's* group, so a timeout could kill the harness and its whole session. | §7.2: retained `Popen`, target only the group this harness created. |
| The schedule globally shuffled individual arms and results were zipped by completion order. | At seed 20260929 all ten "pairs" had different repetition IDs. The intervals described unpaired measurements. | §7.1: adjacent matched blocks, pairing read from the schedule. |
| Rows recorded an older clean commit, with no per-run binary hash. | A row could not be attributed to the candidate that produced it. | §8: per-run provenance including binary SHA256. |

Because the comparators changed, `COMPARATOR_VERSION` in `benchmarks/bench.py` was
incremented to 2, and every recorded row carries the version that produced it. Rows
recorded under version 1 are not comparable with rows recorded under version 2 and
must not be pooled.

## 1. Claim scope

The unit of analysis is a **(command, experiment class, workload)** triple.
Per-command results are primary. No cross-command aggregate is a headline number;
an aggregate must name its command set, its weighting and its missing entries.

**Supported:** per-command runtime and memory ratios with uncertainty intervals, on
named workloads and hardware; scaling behaviour along each cost driver; a
per-command workflow impact figure; an explicit list of what could not be measured.

**Not supported:** any claim of scientific equivalence (that is T3/T4); any claim
about macOS, Windows or ARM; any thread-scaling claim (the port exposes no thread
flag); any claim that this study's numbers generalise beyond the datasets,
scheduling sessions and hardware named in the manifest.

## 2. Estimator

**Primary estimator: actual command-invocation elapsed time, CPU time and peak
memory, with equal delivered work.** Each arm runs as a real process invocation;
nothing is attributed from inside the process. Reported per command as
`wall_primary`, `cpu`, `system` and `peak_rss_mb`.

Experiment classes stay distinct and are never merged:

| Class | Definition | Primary? |
|---|---|---|
| `compute` | The command's own algorithmic work, whole process including its own startup. | yes |
| `data-only` | Data transformation with no plotting, as with `--skip-plot`. | yes |
| `plotting` | Plot rendering, including any Rscript the arm invokes. | separate |
| `end-to-end` | The complete pipeline a user runs, conversion and plotting included. | separate |

A `compute`-class row is not a full-pipeline measurement, and a report must not
present it as one. The 29 rows recorded under version 1 were all `compute`-class;
they measure neither plotting nor the conversion pipeline.

**Sensitivity analysis only.** Subtracting a `--help` floor is reported as
`wall_e2e` and is never the headline. `--help` performs imports and argument
parsing but no algorithmic work, so it is an imperfect proxy for startup cost, and
subtracting it adds uncertainty rather than removing it. Version 1 used the
floor-subtracted figure as its primary estimator; that inversion is corrected here,
and version 1's `wall_e2e` field becomes the secondary analysis.

## 3. Data

- Real whole-genome RNA-seq, plus controlled stress workloads for cost drivers.
- Report alignment records and read pairs **separately**. They differ by a factor
  that depends on the library, so a single ambiguous "size" makes two studies
  incomparable.
- Full annotation for primary realistic workloads. Reduced transcript panels are
  diagnostic sweeps and are labelled as such.
- Prefer public, independently aligned BAMs with verified provenance where building
  a whole-genome STAR index locally is impractical. A BAM generated by this
  repository's own generator tests the port against this project's assumptions.

## 4. Scale

Three sizes (small, intermediate, production) plus selected high-risk drivers: unique
sequences, covered bases, local depth, transcript and exon counts, introns, output
volume, CRAM buffers.

**Label measured counts, never requested counts.** Version 1's "20,000 transcript"
point actually contained 9,179, because the panel was truncated to what was
available. Requested and measured counts are separate fields, and a point where they
differ names the discrepancy.

## 5. Validity: the comparators

A speedup cannot be earned by producing less, or by producing something different
in a way the comparator cannot see. Every comparator therefore fails closed.

### 5.1 BAM

Compared as decoded records, never as compressed bytes: block boundaries,
compression level and optional-field encoding differ legitimately between writers.

Compared fields: query name, **all** flag bits, reference name, start and end,
CIGAR, mapping quality, sequence, **base qualities**, mate reference, mate
position, template length, reverse and mate-reverse bits, and all tags with their
types. Also compared: the header's reference dictionary (names and lengths) and
declared sort order. `@PG` lines are not compared, since each arm legitimately names
a different program.

`require_sorted=True` compares record order instead of a multiset, and is the
correct mode for a command whose contract is coordinate-sorted output. It is a
per-command parameter, not a global setting: several commands legitimately emit
unsorted subsets, and for those, per-file membership is not a meaningful claim.

### 5.2 FASTQ

Strict framing. A line count that is not a multiple of four is a failure naming a
trailing partial record. Every record must have a header starting `@`, a separator
starting `+`, and equal sequence and quality lengths. Header, sequence, separator and
quality are compared independently, so the failure names which field differs.

### 5.3 Text tables

Line-wise with numeric tolerance. Two rules beyond equality:

- **Nonfinite values are rejected, not compared.** `nan` and `inf` mean a division by
  zero or an empty input reached the output. Two arms both producing NaN is a
  shared defect, not agreement, so the comparison fails even when the tokens match.
- Zero-byte files fail. Two empty files are byte-identical, which is not a result.

### 5.4 Declared expectations

Each command declares, in `EXPECTED_STREAMS`: whether it writes to stdout, which
stdout labels must be present, which additional labels are tolerated, and which
artifacts it must produce. A command with no declaration **fails** the gate rather
than passing by omission. The declarations were derived by observing both arms
(`benchmarks/derive_expected_streams.py`), then reviewed; observing behaviour does
not certify correctness, and a declaration is a claim to be checked, not a fact to
be trusted.

## 6. The gate

A row yields a reportable speedup only if the gate passes. Failing marks the row
failed and claims nothing. A fast incomplete run is not a win.

### 6.1 Ordering

1. Each arm's declared artifacts validated in isolation, so "both produced nothing"
   cannot be reported as agreement.
2. Exit statuses compared; a shared non-zero exit is a failure, not agreement.
3. Timeouts recorded on either side fail the gate.
4. Declared stdout labels compared; missing labels, nonfinite values and unexpected
   labels fail. Labels naming an output file have each arm's working directory
   normalised first, because the two arms necessarily run in different directories.
5. File-set comparison: any file in one tree and not the other fails, except
   upstream's non-deliverable log files and documented port-only artifacts.
6. Content comparison by §5.

### 6.2 Failures are retained

Timeouts, out-of-memory conditions, crashes and incomplete outputs stay in the
table with their exit status and stderr tail. They are never replaced by an
invented runtime. A successful-only subset is never averaged without saying so, and
the paired estimate reports how many pairs survived.

### 6.3 Every timed run is validated

Version 1 gated only the last repetition. Each timed run's outputs are now
validated outside the timing interval, so a repetition that truncated its own output
contributes no measurement. Validation is deliberately outside the timing interval:
doing it inside would charge the harness's checking cost to the command.

## 7. Execution

### 7.1 Scheduling: adjacent matched blocks

Machine drift and thermal state move on the timescale of a single run, so pairing is
established by adjacency: repetition *i*'s two arms run back to back. Randomisation
is retained but moved **inside** each block, where it decides only which arm goes
first. That keeps the original intent, that no arm systematically occupies the
warmer slot, without destroying the pairing the intervals depend on.

Each run records its repetition ID, its position in the schedule, and a start
timestamp. The pairing is read from the schedule, never reconstructed by zipping two
independently ordered lists.

**Uncertainty** comes from a bootstrap that resamples whole blocks, preserving the
pairing, with the seed recorded. An interval on one shared-machine run does not
capture dataset or hardware generalisability, which is why §9 requires replication
across independently scheduled sessions.

**Precision rule, declared in advance:** begin with ten valid matched pairs per
command. If the 95% interval's relative width exceeds 0.25, add pairs in blocks of
five until it narrows or until 30 pairs, and record the stopping point. Repetitions
are not collected indefinitely, and the rule is fixed before the data are seen so
the stopping point cannot be chosen to flatter a result.

### 7.2 Timeouts and process groups

Every arm runs as a retained `Popen` in a new session, so its process group is
knowable from the object itself. On timeout the harness terminates **only** the group
it created: SIGTERM, then SIGKILL, then a reap.

The previous handler called `killpg(getpgid(0))`, and `subprocess.run`'s
`TimeoutExpired` carries no PID, so `getpgid(0)` returned the caller's group. A
timeout could therefore kill the harness and everything in its session. Executed
against version 1's code, the timeout test terminated the test runner itself with
SIGKILL rather than reporting a failure.

Output is captured to files rather than pipes. A pipe is inherited by every
descendant, so a grandchild surviving a timeout holds it open and the harness hangs
waiting for an EOF that never arrives. Files also preserve a killed run's output on
disk, which is the diagnostic a timeout investigation needs.

### 7.3 Hardware and environment

Dedicated or demonstrably controlled hardware, with recorded CPU model, core and
thread counts, frequency governor, load average, NUMA topology, storage device and
filesystem, kernel, and the full toolchain. Fixed CPU, thread and compression
budgets. Warm-cache primary study; a cold-cache or deployment-condition series is
run separately and labelled if it is claimed at all.

## 8. Provenance

Every run records, and every archived row carries:

- source revision and whether the tree was dirty;
- SHA256 of every candidate binary invoked, not a prefix;
- SHA256 of `compatibility/upstream.lock`, so a row cannot be attributed to an
  oracle other than the one recorded;
- SHA256 of every input file;
- the exact command line;
- comparator version;
- per-run raw measurements with repetition ID, schedule position and timestamp;
- toolchain, dependency and container versions.

Version 1 identified an older clean commit (`3a5b88c`) with no per-run binary hash,
and its scaling file held summaries with no per-run samples and no environment
manifest of its own.

## 9. Uncertainty and replication

Confidence intervals for time and memory as appropriate, computed from the corrected
paired blocks. All command rows appear, including losses.

Replication across **independently scheduled sessions** and representative datasets
is required. One bootstrap interval on one shared-machine run describes that run. A
claim about a production deployment needs the deployment's datasets, at least, and
the number of sessions is stated with the result.

## 10. Pipeline impact

Measure a real bulk-QC workflow, not just per-command speed: samples per hour,
CPU-hours per sample, memory at intended concurrency, output and I/O volume,
cold-start and install burden, and reproducibility of the QC decisions themselves.
Parallelism across samples may already be the right production design; measure that
deployment before optimising a single process.

If QC is a fraction of the complete RNA-seq pipeline, the overall benefit is
correspondingly bounded, and the report says by how much.

**Implemented and recorded** by `verification/measure_pipeline.py`
([`benchmarks/pipeline-impact-rat-8M.json`](../benchmarks/pipeline-impact-rat-8M.json)),
against the 8.2M-record rat alignment. The panel is the five commands a bulk RNA-seq
QC step consists of — `bam_stat`, `infer_experiment`, `read_distribution`,
`geneBody_coverage`, `read_duplication` — run in that order as one workflow.

One sample, serially: **48.8 s wall, 0.0135 CPU-hours, 1.06 GB peak, 0.1 MiB across 6
files.** The distribution is the finding: `geneBody_coverage` is 45% of the wall time,
`read_duplication` 33%, and the two streaming commands this project has optimised
(`bam_stat`, `infer_experiment`) are **1.5% between them and 6 MB of the peak**.

Whole panels in parallel, which is the deployment shape:

| Concurrent samples | Batch wall | Aggregate peak | Throughput | vs serial | CPU-hours/sample |
|---|---|---|---|---|---|
| 1 | 47.1 s | 2.0 GB | 76 samples/hour | 1.00x | 0.0133 |
| 2 | 50.1 s | 2.0 GB | 144 samples/hour | 1.88x | 0.0141 |
| 4 | 57.2 s | 4.5 GB | 252 samples/hour | 3.29x | 0.0160 |
| 8 | 77.8 s | 7.3 GB | 370 samples/hour | 4.84x | 0.0214 |

Throughput scales to 3.3x at four samples on 16 cores, then flattens: per-sample CPU
rises 61% from one to eight as the panels contend for memory bandwidth. **Running one
sample at a time leaves most of the machine idle.**

Each stage declares what its output must *contain*, and a failed stage stops the panel
rather than being timed and averaged in — a stage run on the output of a broken one
measures nothing, and including its time would understate a working panel's cost.

**Not yet measured, and named:** cold-start and install burden (the archive's
extraction-plus-first-invocation cost, which the archive smoke test exercises but does
not time), decision reproducibility as a *measured* property (the differential suite
establishes output equality, which is stronger, but is not framed as a workflow
property), and the ratio of QC to a complete RNA-seq pipeline. The single-sample,
single-organism, single-machine basis is stated with the figures rather than buried.

## 11. Memory measurement

GNU `time -v`'s child RSS is the **largest single child's** RSS, not the sum over a
concurrently running process tree (Linux documents
`RUSAGE_CHILDREN.ru_maxrss` this way). Reported as `peak_rss_mb` and labelled as
such. Helper-heavy workflows additionally get an explicitly labelled aggregate
process-tree measurement, and a separate per-process RSS measure where the two
differ materially.

## 12. Archive

Raw per-run samples with pair IDs, failures, command lines, full hashes, dependency,
container and toolchain versions, workload manifests, the analysis code, and the
publication figures. Raw scaling runs, not only summaries.

## 13. Stopping conditions

These are evidence gates, not calendar promises.

1. **Validation foundations:** converter truth cases, benchmark mutations, expected
   outputs, safe timeouts, correct pairing. Stop only when each deliberate defect is
   detected. *Status: met on 2026-10-01; see `benchmarks/test_bench_harness.py`,
   `verification/test_comparators.py`, `verification/test_run_diff.py`,
   `datasets/test_refgene_to_gtf.py`, `datasets/verify_refgene_frame.py`.*
2. **Release candidate:** pinned baseline, per-case provenance and mixed-input
   reporting, meaningful failure branches, artifact interoperability, executable CI.
   Stop when every advertised beta capability has a passing recorded profile.
3. **Whole-genome and production-size pilots:** correctly prepared independent
   inputs, constrained memory and concurrency, reviewed discrepancies. Stop when the
   declared deployment envelope has evidence.
4. **Collect the study once:** frozen candidate and protocol, matched blocks, equal
   work, raw samples, regenerated figures. Recollect affected rows only when code,
   methods or an unresolved discrepancy justify it.
5. **Package and review:** validate the actual archive and container plus an external
   installation; complete URLs, citation and support materials; release within the
   passed scope.
6. **Manuscript and submission:** human scientific and methods review, final archived
   study, clear scope, reproducible artifacts.

More testing is justified where it can reveal a meaningful undetected defect. More
benchmarking is justified only after correctness and measurement are qualified.
