# Post-freeze changes to the benchmark protocol and harness

`benchmarks/protocol.md` was frozen on 2026-09-29 before any measured run. Protocol
§1 states that any later change must be recorded here with its reason, and that any
result depending on a changed rule must be re-derived. This file is that record.

Two classes of change are distinguished:

- **HARNESS BUGS** — the harness was measuring the wrong thing. Every result produced
  by a buggy harness revision was discarded and the affected rows re-measured from
  scratch. The discarded numbers are named so they are not mistaken for results.
- **PROTOCOL CLARIFICATIONS** — the measurement rules were ambiguous and are now
  explicit. Applied uniformly to all rows.

---

## 1. Protocol clarifications (no result invalidated)

### 1.1 An emitted R plotting script is not a deliverable (2026-09-29)

The gate initially compared every file in both output trees byte-wise. The 14
`--skip-plot` commands emit an `.r` script, which embeds the absolute output paths of
the arm that produced it. Since each arm runs in its own directory, those paths
necessarily differ, so the gate failed on a path difference and discarded rows whose
actual data tables were identical (`read_quality`, `read_NVC`).

`.r` files are now excluded from content comparison. Their *presence* is still
enforced by the file-set comparison, and the data table each one plots is compared
normally. Result: unchanged timings, gate outcome corrected.

### 1.2 Non-deterministic commands are gated against a self-control, not equality (2026-09-29)

`RPKM_saturation` subsamples reads per percentile with no RNG seed. Verified
empirically: two consecutive **upstream** runs on identical input produce different
files. Byte-equality is therefore unsatisfiable for a correct port.

Rather than exempting the command, the harness now runs the reference twice and
compares it to itself. The port passes only if it agrees with upstream at least as
often as upstream agrees with itself. Measured: port-vs-upstream row agreement 59.3%
/ 59.4%, upstream-vs-upstream 58.4% / 58.5% — the port is marginally *more*
self-consistent with the reference than the reference is with itself.

### 1.3 `divide_bam` is gated on the partition union (2026-09-29)

DIV-0017 already documents that the two builds use different RNG algorithms, so
`--seed N` cannot assign the same query names to the same subset files. Per-file
membership is therefore not a valid equivalence criterion. The gate now checks the
union of the output subsets — that both arms perform the same partition — and the row
is labelled as a documented divergence. Per-file equality is explicitly not claimed.

### 1.4 `geneBody_coverage`'s `log.txt` is a log, not a result (2026-09-29)

Upstream writes `log.txt` (`geneBody_coverage.py:57`); the port logs to stdout. It is
excluded from the file-set check and reported on the row.

### 1.5 `sc_seqLogo` produces strictly more than upstream here (2026-09-29)

Upstream's plotting path crashes in the pinned oracle (logomaker/pandas incompat,
DIV-0016), so upstream emits only the count matrix while the port also renders SVG.
The two SVGs are recorded as extra candidate artifacts and the common deliverable
(count matrix) is compared. The port's end-to-end time therefore *includes* rendering
that upstream never performed, which makes its number conservative.

### 1.6 The stochastic self-control required strict dominance — 2026-09-29

The `RPKM_saturation` self-control originally passed only if the port's agreement with
upstream was *strictly greater* than upstream's agreement with itself. Both figures carry
sampling noise, so a 0.1-percentage-point difference decided the row: one run passed at
59.3% vs 58.4%, the next failed at 58.3% vs 58.4%, with no change in the port. A test
whose outcome depends on which side the noise falls on is not a test.

The gate now requires the port to be no more than 2 percentage points behind upstream's
own reproducibility (`STOCHASTIC_AGREEMENT_TOLERANCE`), and reports both figures on the
row either way. The measured situation is that the two are indistinguishable, which is
the accurate conclusion.

---

## 2. Harness bugs (affected results discarded and re-measured)

### 2.1 Output directories were shared across repetitions — INVALIDATED all rows, 2026-09-29

Each arm wrote into one directory reused by every repetition. Several upstream
commands refuse to overwrite an existing output and exit 2 (`split_bam`,
`split_paired_bam`, `divide_bam`), and the port exits 1, so **every repetition after
the first measured two immediate failures instead of the command**. The affected rows
looked like extreme wins (`split_bam` "0.79x") or losses for reasons that had nothing
to do with the code.

Fixed by giving each arm a fresh output directory per repetition
(`<cmd>/{py,rs}/repNN`), which is also what protocol §2/§5 intend by isolated
workdirs. All 29 rows were discarded and re-measured. `split_bam` went from a
meaningless 0.79x to 3.14x, and the gate failures it caused disappeared.

### 2.2 The point estimate and the confidence interval described different quantities — INVALIDATED all reported ratios, 2026-09-29

The ratio was computed from **raw** medians while the bootstrap interval was computed
from **floor-subtracted** (E1) pairs. For short commands these diverge badly: the
`bam_stat` 2k-read scaling point reported a 14.24x point estimate with a
[2.05, 4.02] interval — an interval that does not contain its own estimate, which is
the signature of the bug.

Fixed so the ratio and its interval are both computed from the E1 quantities, and the
raw end-to-end medians are reported separately and labelled. `analyze.py` recomputes
both from the stored raw per-run measurements, so the committed report is consistent
regardless of which harness revision produced the raw file.

### 2.3 The point estimate and the interval are now recomputed from raw data

`benchmarks/analyze.py` re-derives the E1 medians, both estimators and both intervals
from `runs` + `floor_s` in the raw JSON. No reported number is taken on trust from the
harness's own summary.

### 2.4 Workload completeness was not verified — INVALIDATED one scaling point, 2026-09-29

`scaling.py` skipped regeneration whenever `reads.bam` existed. The generator's own
reference validation had failed for the 150bp point (it rejects a mismatch rate above a
fixed threshold that did not scale with read length), leaving a partial workload with
no manifest and no FASTQ. Both arms then ran over missing input, agreed perfectly on
empty output, and the point was recorded as a clean pass with a nonsense 25.79x.

Fixed three ways: the generator's tolerance now scales with read length; `gen()`
validates that the full artefact set exists and is non-empty before reusing a
workload; and a missing workload now raises rather than silently measuring nothing.
The 150bp point was regenerated and re-measured (3.35x, not 25.79x).

### 2.5 A stray runaway process contaminated a full 10-rep run — INVALIDATED that run, 2026-09-29

An orphaned `cargo test --workspace --release` from 14:13 on 2026-09-29 was stuck in
`tests::open_alignments_decodes_a_real_no_reference_cram_fixture` and ran continuously
for 3+ hours at 100% of one core, with RSS growing past 20 GB (48 GB virtual, 6.3 GB
swapped). It was present for the entire duration of a full 10-repetition `main` run, so
every row in that run was measured on a machine with a persistent background load of
unknown interference.

The process was killed and the whole suite re-measured from scratch. The CRAM hang
itself was a transient mid-edit state that the later `crates/formats/src/lib.rs` change
already fixed: the test passes in 0.00 s now, three runs in a row.

**This is a second, independent instance of the same underlying class of error as §2.1
and §2.4 — a measurement that was never checked for whether it measured anything.** The
harness records the load average at start (`os.getloadavg()` in `environment()`), but
nothing compared load before against load after, so a machine that got busier *during* a
block, or was already busy before it started, produced a clean-looking run.

### 2.6 `tin` holds a whole-file read index where a sliding window suffices — results changed, re-measured, 2026-09-29

Not a harness bug, but recorded here because it moved a reported number by more than an
order of magnitude in the column that decides whether the port is usable.

`tin` builds its read index with `build_read_index`, which retains every read in the BAM
for the whole run. Each `IndexedRead` carries per-read heap buffers (query name,
qualities, sequence, CIGAR, two block lists), so this costs roughly 600 bytes per read:
600k reads measured 366 MB against upstream `pysam`'s 43 MB, and the port lost that row's
memory 8.5x even while winning on time 8.05x.

But `compute_tin` only ever consumes the reads whose *start* falls inside the current
transcript's span. The port now streams the BAM once in transcript coordinate order and
keeps only the reads that can still reach a transcript not yet scored, retiring a read
as soon as `end <= tx_start` (safe, because `start < end <= tx_start` puts it before
every later window and `reads_starting_in` would have excluded it anyway). The resident
set falls from 597,048 reads to a mean of 3,424 (max 16,438) on the benchmark workload.

This is exact, not an approximation: the per-transcript scoring arithmetic is a single
shared `score_sample` function, transcripts are scored in coordinate order but collected
by original index so row order and the summary's pairwise-summation order are unchanged,
and both output files are byte-identical to upstream on all 3000 transcripts on the
600k-read workload and all 26,590 on the 400k-read workload, with and without
`--subtract-background`. A non-coordinate-sorted input is detected and falls back to the
whole-file path, which is the only behavioural difference (pysam itself refuses to index
an unsorted BAM, so this is not reachable through the documented CLI).

---

## 3. Additional artefacts produced

- `benchmarks/generate_workload_real.py` — Tier A/B workload generator over real hg38
  sequence and real RefSeq annotation. Its post-write reference validation caught four
  independent generator bugs that all produced syntactically valid but
  scientifically wrong BAMs (see its own comments and RESULTS.md §6).
- `benchmarks/bench.py` — the harness.
- `benchmarks/scaling.py` — cost-driver sweeps.
- `benchmarks/analyze.py` — regenerates every table from raw JSON.
- `benchmarks/results-main.json`, `benchmarks/results-scaling.json` — raw measurements.
