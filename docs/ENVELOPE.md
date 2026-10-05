# Production envelope

What this project can be relied on to do at what scale, and what it cannot.

**Status: NOT a declared envelope.** This records what has been measured and what
remains unmeasured. A published limit needs a measured slope at a named scale *and*
the observed failure behaviour past it; neither is established for most commands yet.
Publishing the numbers below as a supported envelope would be the error the
2026-10-01 readiness audit found elsewhere in this project.

## Measured so far

Fixture: `datasets/aligned/tiny/tiny.bam` (621,783 records, 3,518,066 covered bases,
chr17 only), measured with `/usr/bin/time -v`, whose "Maximum resident set size" is
the largest **single** child's RSS — Linux documents `RUSAGE_CHILDREN.ru_maxrss` that
way, so it is not an aggregate process-tree figure. Reproduce with:

```bash
cargo build --workspace --release --locked
oracle/venv/bin/python3 verification/measure_memory.py --json memory.json
```

The default sweep (0.25 / 0.5 / 1.0 of the busiest contig) is too coarse to produce
a slope from this fixture, and the script says so rather than reporting one. The
measurement below therefore uses six windows:

```bash
oracle/venv/bin/python3 verification/measure_memory.py \
  --fractions 0.01 0.02 0.05 0.10 0.20 0.40 --min-span 2.0 \
  --json benchmarks/memory-sweep-tiny.json
```

Six nested windows of chr17, 21,387 → 621,783 records, a **29.07x** range (each the
peak of three repeats, with the spread recorded):

| Cost driver | Command | Retains | Peak RSS over the range | Fitted slope |
|---|---|---|---|---|
| streaming (control) | `bam_stat` | nothing beyond the current record | 3.0 → 2.9 MB | **−0.1 B/record (flat)** |
| duplicate counting | `read_duplication` | one entry per read | 56.5 → 85.1 MB | **43.1 B/record** |
| covered positions | `bam2wig` | a value per covered base | 7.5 → 151.7 MB | **250.4 B/record** |

Raw samples, per-window record and covered-base counts, the input's SHA256 and
per-command wall times are in
[`benchmarks/memory-sweep-tiny.json`](../benchmarks/memory-sweep-tiny.json).

**The bytes-per-record figures this document previously reported are withdrawn.**
They were 0.6, 135.6 and 270.4 B/record, fitted by least squares over the two
sizes above — a 1.22x range. A fit over a 1.22x range cannot distinguish a linear
cost from a step, a saturating curve, or a single allocation that happens to land
inside the window, and three of the measurements in the sweep do not even have
monotone denominators: restricting the fixture to 25%, 50% and 100% of chr17 gives
507,983, 621,783 and 621,783 records, because the fixture's reads are concentrated
and the last two windows contain the same reads. Two of the three points are
therefore the same measurement, and the slope was fitted through a repeated point.

`verification/measure_memory.py` now refuses to report a slope unless the measured
inputs span at least 2x (`--min-span`), records the input's own digest and record
count alongside the numbers, and says why it declined. A per-record cost is
re-published only from a real alignment sliced into windows that genuinely differ
in record count.

Two things the sweep does establish, and one it does not.

**The streaming claim now holds across a 29x range of record counts.**
`bam_stat` measures 3.0, 2.9, 2.9, 3.1, 2.9 and 2.9 MB at 21k, 43k, 84k, 246k,
359k and 622k records — a fitted slope of −0.1 B/record, i.e. flat within
measurement noise, over a range wide enough that a per-record cost of even
20 bytes would have been visible. That is the control that makes the other two
rows meaningful: the method can detect a slope where one exists.

**`read_duplication` costs about 43 bytes per record here, not the 135.6 previously
reported.** The earlier figure came from a 1.22x range whose two larger points held
*identical* record counts, so the fit passed through a repeated measurement. At 43
B/record the cost is consistent with a hash entry per read rather than a
per-base or per-quality-string structure.

**`bam2wig` costs about 250 bytes per record**, which on this fixture is 3,518,066
covered bases over 621,783 records — about 44 bytes per covered base. This is the
driver that bounds the envelope: a whole-genome human alignment at ordinary
coverage covers on the order of 10^9 positions, which at that rate is tens of
gigabytes.

**What a slope is not.** These windows are nested prefixes of one contig, so the
figures are a cost-per-record estimate on a single synthetic alignment. They say
nothing about behaviour that changes with alignment *shape* — read length, number of
contigs, read distribution — and nothing about the per-transcript drivers, whose cost
follows transcript count rather than record count. They are not a supported limit.
The next section repeats the sweep on a real 8.2M-record alignment, where
`read_duplication` costs 130.5 B/record against 43.1 here: **a per-record cost
measured on synthetic reads is a property of the generator.**

## On a real alignment: 8.2M records, 174x range

The same sweep on the rat alignment (`SRR1177982`, 8,200,844 records, 38,306,284
covered bases, 3 contigs) — the largest alignment measured anywhere in this project.
Six nested windows of chr1, 22,911 → 3,985,238 records, a **173.94x** range. Every
point is the **peak of three repeats**, and the spread across those repeats (0.1–0.5 MB
here) is recorded as this instrument's own noise floor for that command and input:

```bash
oracle/venv/bin/python3 verification/measure_memory.py \
  --bam datasets/heldout/aligned/SRR1177982/SRR1177982.bam \
  --bed datasets/heldout/reference/rn6.indexed.bed12 \
  --fractions 0.02 0.05 0.10 0.25 0.50 1.0 --min-span 2.0 \
  --json benchmarks/memory-sweep-rat-8M.json
```

| Cost driver | Command | Peak RSS over the range | Fitted slope |
|---|---|---|---|
| streaming (control) | `bam_stat` | 3.0 → 3.0 MB | flat within measurement repeatability; no slope published |
| duplicate counting | `read_duplication` | 56.4 → 530.5 MB | **130.5 B/record** |
| covered positions | `bam2wig` | 12.3 → 771.5 MB | **207.5 B/record** |

`bam_stat` moves 0.1 MB across a 174x range of record counts — less than its own
repeat-to-repeat spread of 0.2 MB, so the sweep withholds a slope for it entirely
rather than publishing a number indistinguishable from zero. The two rows below
detect their slopes easily, so the method demonstrably can; its failure to detect one
here is the finding.

**`read_duplication` costs about 130 bytes per record on real data**, against 43 on
the tiny fixture. The difference is what the fixture cannot show: a hash entry's
cost depends on read length, and this panel's reads are 101 bp against the fixture's
synthetic ones. A per-record cost measured on synthetic reads is a property of the
generator.

**`bam2wig` costs about 207 bytes per record here**, and 44 bytes per covered base
(18,317,951 covered bases over 3,985,238 records). It is the driver that bounds the
envelope. On the *full* 8.2M-record alignment it peaks at **1.69 GB** and takes 51 s.

### Re-measurement after cards B1/B2 (2026-10-04)

B1 re-keyed `read_duplication`'s maps by 128-bit hash and B2 replaced
`bam2wig`'s per-position map with sparse/dense chunks plus streamed
drain rendering. Single-point peaks on the same full 8.2M-record rat
alignment, same instrument (`/usr/bin/time -v`, `Maximum resident`):

| Command | Peak RSS before | Peak RSS after | Wall before | Wall after |
|---|---|---|---|---|
| `read_duplication` | 1,085,832 KB (1060 MB) | 543,472 KB (531 MB) | 15.1 s | 12.3 s |
| `bam2wig` | 1,725,856 KB (1685 MB) | 759,376 KB (741 MB) | 46.5 s | 42.6–43.1 s |

Outputs byte-identical before/after on both (both `.xls` files for
`read_duplication`, the 498 MB `.wig` for `bam2wig`, each `cmp`-clean over
three post-change runs for `bam2wig`). The per-record slopes above
(130.5, 207.5 B/record) were fitted on the old implementations and are
**stale until the sweep is re-run**; only the peaks in this section are
re-measured. The per-job maxima quoted elsewhere in this document
(1.69 GB for `bam2wig`, 1.06 GB for `read_duplication`, e.g. in the
concurrency section) are likewise pre-B1/B2 figures.

## What happens past the limit

A cost is half of a published limit. The other half is the behaviour when memory
runs out, which the audit requires to be either complete output or output that is
*identifiably* incomplete — a killed process leaving a plausible partial file fails
that, and a silent OOM fails it worse.

`verification/check_memory_failure.py` runs each non-flat command under `ulimit -v`
at 80%, 40% and 20% of its limit-free peak, in a fresh directory each time so a
previous run's output cannot be mistaken for a failed one's. Measured on the rat
alignment ([`benchmarks/memory-failure-rat-8M.json`](../benchmarks/memory-failure-rat-8M.json)):

| Command | Limit-free peak | At 80% / 40% / 20% | Output files written | Input changed |
|---|---|---|---|---|
| `bam2wig` | 1,685 MB | SIGABRT (134) at all three | **0** | no |
| `read_duplication` | 1,060 MB | SIGABRT (134) at all three | **0** | no |

**The contract property holds.** Neither command leaves a partial artifact, so no
consumer can read a truncated result as a complete one, and neither modifies its
input. Each failure prints something on stderr, so the run is not silent.

**The diagnostic is fixed for the allocator-abort path (B4, 2026-10-04).**
A `#[global_allocator]` in `crates/cli/src/lib.rs` (linked into all 33
binaries, verified via the message string in each) writes
`rseqc-rust: out of memory; see docs/ENVELOPE.md for per-command memory
costs` to stderr without allocating, then returns null so the normal
abort follows with the same exit status as before. Re-ran
`verification/check_memory_failure.py` on the rat alignment: zero output
files at every limit, exits unchanged (134 at 80/40/20% for
`read_duplication`; 134 at 80% for `bam2wig`), and the new message now
precedes the runtime's `memory allocation of N bytes failed` line. Two
caveats stay open: at 40% and 20% `bam2wig` still dies inside `zlib-rs`'s
inflate (`assertion left == right failed, left: MemError, right: Ok`,
exit 101) before any Rust allocation fails, so no actionable message can
appear there; and `bam_stat` wall time is unchanged (4.12/4.09/4.16 s
before, 4.12/4.18/4.12 s after, three runs each).

## Per-transcript drivers: a sweep that produced nothing, and why

Three commands retain per-transcript state — `geneBody_coverage` computes a percentile
list per transcript, `tin` an exon index, `junction_annotation` the reference's intron
sets — so their cost follows **transcript count**, not alignment-record count. The
per-record sweep cannot see that: restricting an alignment by coordinate does not change
how many transcripts a command builds state for. `verification/measure_memory.py
--transcript-drivers` therefore varies the BED12 while holding the 8.2M-record
alignment fixed ([`benchmarks/memory-transcripts-rat.json`](../benchmarks/memory-transcripts-rat.json)).

**The sweep's first run published 1234.5 bytes/transcript from fifteen runs that all
failed.** Every invocation exited 1 in 0.00s, because the alignment path was relative
and every command runs with its cwd set to a per-run scratch directory, so it resolved
against that directory. The slope function did not check exit codes, so it fitted a
line through fifteen instantaneous failures. That is precisely the failure mode the
script was written to stop — and it took a real measurement to hit it, not a review.

Two rules now make it structurally impossible: a slope is never reported over runs
with a non-zero exit, and input paths are resolved before use.

**With the runs actually succeeding, all three slopes came out negative**, and a
negative per-transcript cost is not a cost. The number is withheld, with the reason
recorded per command:

- `geneBody_coverage` and `tin` fall from 320 MB to 221 MB as the model grows 16x. Both
  hold a read-side window that cannot retire a read until every transcript it might
  belong to is scored, so a **smaller** model retains **more** reads. That effect runs
  opposite to the annotation cost and dominates it. The per-transcript cost is not
  measured here; the confound is.
- `junction_annotation` is essentially flat (75 → 74 MB) across the same 16x. It
  scores records in one streaming pass and holds no read-side window, so its annotation
  cost *should* dominate here — which says the per-transcript state is small next to
  the fixed cost of the input, not that no per-transcript cost exists.

**Holding the read stream small enough removes the confound, and the answer is mostly
"smaller than we can measure".** The same sweep against the 621,783-record tiny
fixture, where read retention cannot dominate
([`benchmarks/memory-transcripts-tiny.json`](../benchmarks/memory-transcripts-tiny.json)),
each point being the peak of three repeats:

| Command | 365 → 5,840 transcripts | Peak RSS | Repeat spread | Verdict |
|---|---|---|---|---|
| `tin` | 16x | 27.9 → 36.9 MB | 0.2–0.5 MB | **1673 B/transcript** (8.8 MB of movement, ~20x the noise) |
| `junction_annotation` | 16x | 73.8 → 74.9 MB | 0.1–0.6 MB | withheld: 1.8x the noise |
| `geneBody_coverage` | 16x | 56.0 → 56.0 MB | 0.2–0.4 MB | withheld: 0.5x the noise |

So `tin`'s per-transcript cost **is** measurable and is about 1.7 KB per transcript —
its exon block index, which is a `Vec` of blocks per transcript. The other two hold
per-transcript state small enough to be invisible at this scale: `junction_annotation`
at 1.1 MB across 16x, `geneBody_coverage` at 0.2 MB, both below the sweep's own
repeatability. That is a real answer, not an absence of one — but bounding those two
properly needs model sizes far larger than this panel's 5,359 transcripts, which means
a whole-genome annotation.

**The rule that produced these verdicts.** A slope is published only when the fitted
quantity moves at least 3x the repeat-to-repeat spread measured by running each point
three times. Every measurement in this document reports its spread. This replaced a
fixed 20 MB threshold, which was derived from `/proc/meminfo`'s ~90 MB system-wide
noise floor measured by `measure_concurrency.py` — the wrong instrument, applying the
wrong noise figure to a per-process measurement, and it discarded `tin`'s real 8.8 MB
signal.

## Concurrency: is per-job memory really additive?

`verification/measure_memory.py` and `check_memory_failure.py` both measure one
process at a time, and this project has stated throughout that "memory per concurrent
invocation is therefore additive: N concurrent jobs need N times the per-job figure".
That was an assumption. `verification/measure_concurrency.py` measures it, by sampling
`MemTotal - MemAvailable` at 20 Hz while N invocations run and subtracting a baseline
taken before the batch starts.

On the 8.2M-record rat alignment
([`benchmarks/concurrency-rat-8M.json`](../benchmarks/concurrency-rat-8M.json)):

| Command | ×1 | ×2 | ×4 | ×8 | Largest single job |
|---|---|---|---|---|---|
| `bam2wig` | 1.89 GB (1.12×) | 3.11 GB (1.85×) | 6.77 GB (4.02×) | 11.99 GB (7.11×) | 1.69 GB |
| `read_duplication` | 0.78 GB (0.73×) | 2.47 GB (2.33×) | 4.42 GB (4.17×) | 7.92 GB (7.47×) | 1.06 GB |

**The assumption holds.** Aggregate memory tracks concurrency linearly to within the
measurement's resolution at ×4 and ×8 (4.02× and 7.11× against per-job maxima of
1.69 GB), and wall time is nearly flat in concurrency — `bam2wig` takes 50.5 s at ×1
and 64.8 s at ×8, so 8 samples finish in 1.3× the time of one. That is the deployment
shape the audit asks about: parallelism across samples, not inside a command.

**Two caveats, both recorded rather than smoothed over.** The ×1 figures sit *below*
the single job's own peak RSS (0.73× and 1.12×) in one direction and above it in the
other, because `/proc/meminfo` has a noise floor of roughly 90 MB on an idle machine
here — page cache moves under `MemAvailable` — and the script measures that floor
rather than assuming a tolerance. And on a shared machine, memory in use *after* a
batch varies by several hundred MB between adjacent runs, so individual cells should
be read as an order of magnitude, not as a precise figure. The ×4/×8 agreement is
what makes the conclusion; any single cell is not.

**What this does not establish.** Both commands are single-process and allocate their
own working set, so additivity is the expected result. It is measured rather than
assumed now, and a command that later added shared memory (a mmap'd index, a
thread-pool arena) would break the pattern — which is why the measurement belongs in
the release's evidence rather than in a code comment.

## The whole QC panel, measured as a workflow

Per-command costs do not answer the deployment question, and the audit is explicit
that "if QC is only a fraction of the complete RNA-seq pipeline, the overall benefit is
correspondingly bounded". `verification/measure_pipeline.py` times the five-command
bulk QC panel as a workflow — `bam_stat`, `infer_experiment`, `read_distribution`,
`geneBody_coverage`, `read_duplication` — against the 8.2M-record rat alignment
([`benchmarks/pipeline-impact-rat-8M.json`](../benchmarks/pipeline-impact-rat-8M.json)).

**One sample, serially: 48.8 s wall, 0.0135 CPU-hours, 1.06 GB peak, 0.1 MiB of
output across 6 files.** The cost is not evenly spread:

| Stage | Wall | Share | Peak RSS |
|---|---|---|---|
| `geneBody_coverage` | 22.1 s | 45.2% | 221 MB |
| `read_duplication` | 15.9 s | 32.7% | 1060 MB |
| `read_distribution` | 6.3 s | 12.8% | 26 MB |
| `bam_stat` | 4.2 s | 8.6% | 2.9 MB |
| `infer_experiment` | 0.4 s | 0.7% | 3.3 MB |

**The two commands this project has optimised are 1.5% of the panel.** `bam_stat` and
`infer_experiment` are the streaming ones, and together they account for 15% of wall
time and 6 MB of peak memory. A speedup claim about either is bounded by the fact that
they are not where the time goes.

**Parallelism across samples is the deployment that works**, and it is measured rather
than assumed:

| Concurrent samples | Batch wall | Aggregate peak | Throughput | vs serial |
|---|---|---|---|---|
| 1 | 47.1 s | 2.0 GB | 76 samples/hour | 1.00x |
| 2 | 50.1 s | 2.0 GB | 144 samples/hour | 1.88x |
| 4 | 57.2 s | 4.5 GB | 252 samples/hour | 3.29x |
| 8 | 77.8 s | 7.3 GB | 370 samples/hour | 4.84x |

Throughput scales to 3.3x at four samples on 16 cores and then flattens — 4.84x at
eight, where per-sample CPU rises from 0.0133 to 0.0214 CPU-hours as the panels
contend for memory bandwidth. **A deployment that ran one sample at a time would leave
most of the machine idle**, which is the audit's own point that per-command optimisation
is the wrong lever here.

**What this is not.** One sample, one organism, one alignment, one machine, and a
five-stage panel chosen to represent a bulk QC step rather than to reproduce anyone's
specific pipeline. It is not a claim about total RNA-seq pipeline cost, and the ratio
between QC and everything else in a real pipeline is unmeasured.

## Known code-level limits, from source

These are read from the implementation, not inferred from timing, because the fixtures
available are far too small to demonstrate them.

| Limit | Where | Consequence |
|---|---|---|
| CRAM decodes container by container (B3) | `crates/formats/src/lib.rs`, `AlignmentRecords::CramStreaming` | One container's records buffered at a time, decoded via the public `read_container` / `slices` / `decode_blocks` / `slice.records` calls. Measured: `bam_stat` on an 8.2M-record embedded-reference CRAM (397 MB file, pysam `embed_ref=2` from the rat BAM) peaks at 51 MB vs 3 MB on the BAM, byte-identical output; a 1-record CRAM peaks at 3.7 MB. Slower than BAM (~79 s vs ~4 s, dominated by CRAM decode plus the per-container BAM round trip), but flat in file size. The external-reference error still surfaces at open with the historical message. |
| `bam2wig` retains covered positions | `crates/commands/src/bam2wig.rs` | Proportional to covered bases; see the measured slope above. |
| `read_duplication` retains one entry per read | `crates/commands/src/read_duplication.rs` | Proportional to alignment records; see the measured slope above. |
| SAM and CRAM round-trip through memory | `crates/formats/src/lib.rs`, `encode_decode` | Bounded by the chunk size for SAM (4,096 records) and by one container for CRAM (B3); whole-file buffering is gone on both paths. |
| Threading | whole workspace | The port exposes no `--threads` flag and has no first-party parallelism. Memory per concurrent invocation is therefore additive: N concurrent jobs need N times the per-job figure. |

## What is not established

Recorded as gaps rather than left to be assumed:

- **No production-size run at the audit's stated targets.** The rat alignment is
  4.1M read pairs (8,200,844 records) over a 3-contig index. The audit's production
  pilot calls for 10M, 50M and optionally 100M pairs on a full annotation; none has
  been run, and the contig subset is a hardware limit recorded in
  `datasets/manifest.yaml`, not a design choice.
- **No declared memory limit per command.** A cost per record and a measured peak are
  not a limit. A limit needs both, plus the behaviour past it — and past-it behaviour
  is now measured for two commands (above) but only against a *virtual* limit set by
  `ulimit -v`, not against a cgroup or an OOM kill, which behave differently.
- **The out-of-memory diagnostic is an internal abort**, not a message naming the
  command, its input and the reason. Established above; not fixed.
- **`junction_annotation`'s and `geneBody_coverage`'s per-transcript cost is bounded
  only from above** (see the section above): below this measurement's repeatability at
  a 16x model range. `tin`'s is measured at ~1.7 KB per transcript. Bounding the other
  two properly needs a whole-genome annotation, which is the same gap as item 1.
- **Windows are nested prefixes of one contig.** A cost per record is estimated well
  by that, but nested windows cannot reveal behaviour that depends on the *shape* of
  the input — read distribution across a contig, number of distinct chromosomes, or a
  per-reference structure.
- **Concurrent-invocation memory below ×4.** Measured at ×1, ×2, ×4 and ×8 on this
  machine for the two non-flat commands (above), and additive. Three samples of a
  shared machine is not a concurrency study: the x1 cells move by several hundred MB
  between adjacent runs, and levels above ×8 were not run.
- **CRAM at scale: measured once (B3).** Beyond the 10-record fixture, `bam_stat`
  on the 8.2M-record embedded-reference rat CRAM peaks at 51 MB with output
  byte-identical to the BAM run. One file, one command, one machine -- not a
  per-record CRAM slope.
- **The rat panel is one exposed sample.** It is the largest alignment measured, and
  it is also the stratum whose result drove development changes. Its cost
  characteristics are not a general claim about RNA-seq alignments.

## How to close this

1. Build or obtain an independently aligned whole-genome BAM at 10M and 50M pairs
   with a full annotation, and re-run `measure_memory.py` across at least four
   windows per driver.
2. Add a cgroup-based sweep alongside the `ulimit -v` one, so the OOM-kill path
   (exit 137, no chance to clean up) is measured too, not just the allocator-failure
   path.
3. ~~Repeat the whole sweep at 4 and 8 concurrent invocations.~~ Done for the two
   non-flat commands; the streaming ones need no aggregate measurement because their
   per-job figure is already flat.
4. ~~Re-run the per-transcript sweep against a small alignment.~~ Done: `tin` measures
   ~1.7 KB per transcript, and the other two fall below the measurement's
   repeatability at a 16x model range. Bounding them needs a whole-genome annotation,
   which is item 1.
5. Only then state per-command limits, with the scale, the annotation and the
   hardware they were measured on.

Until then, the honest position is: this project streams where it has been measured
to, and three commands have measured per-record costs that bound the workloads they
can be used on.
