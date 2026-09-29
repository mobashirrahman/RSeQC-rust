# Benchmark protocol — rseqc-rust vs upstream RSeQC

**Status: FROZEN.** Written before any measured run. Any change to this document
after measurement must be recorded in `benchmarks/CHANGES.md` with the reason and
the date, and any resulting claim must be re-derived from scratch.

- Protocol version: 1.0.0
- Frozen: 2026-09-29
- Reference implementation: upstream RSeQC 5.0.5, commit
  `59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24` (`compatibility/upstream.lock`),
  executed from the pinned oracle environment `oracle/venv`
- Candidate implementation: this repository, `cargo build --workspace --release --locked`
- Governing requirements: `testing.md` §12, `docs/PORTING_PLAN.md` Step 8 (gate G6a)

---

## 1. Purpose and claim scope

To quantify, per command, how much faster (and how much less memory-hungry) the
Rust port is than the Python implementation it replaces, **and to characterise the
shape of that advantage as a function of input size** rather than reporting a
single ratio.

The unit of analysis is a **(command, mode, workload)** triple. Per-command results
are primary. No cross-command aggregate is published as a headline number; if an
aggregate is shown it must name its command set, its weighting, and its missing
entries.

### 1.1 Claims this study may support

1. A per-command speedup ratio with a 95% confidence interval, on named workloads.
2. A per-command peak-memory ratio, reported separately from time.
3. A scaling exponent / cost curve along each cost driver in §6.
4. An explicit list of commands that cannot be benchmarked here, with the reason.

### 1.2 Claims this study may NOT support

1. **No publication-grade speedup claim.** See §2 on hardware isolation.
2. No claim about macOS or Windows performance.
3. No claim that a workload is *scientifically* equivalent merely because outputs
   are byte-identical; that is a T3/T4 question, not a T5 one. Equivalence here
   means "outputs match on this workload", nothing more.
4. No thread-scaling claim. The port is single-threaded (no `rayon` in first-party
   code; `rayon` appears only transitively via `noodles-bgzf`) and exposes no thread
   flag, so a scaling sweep is not possible. Recorded as a known gap in §8.

---

## 2. Hardware isolation — declared limitation

Measurements in this study are taken on a **shared, non-isolated** machine. They
are internal engineering evidence, useful for finding hotspots and regressions and
for comparing orders of magnitude. They are **not** publication-grade.

Recorded per run: CPU model, core/thread counts, governor, load average, kernel,
filesystem, `rustc` and Python versions, and executable hashes.

Mitigations applied even though hardware is shared:

- Pairing and randomisation (§5) so that machine drift cannot systematically favour
  one arm.
- Warm-up run discarded before measurement.
- Background load sampled before and after each block and recorded; blocks whose
  load varies by more than a stated threshold are flagged.
- The *primary* reported statistic is a paired ratio, which cancels slow-moving
  machine-level drift that a ratio of unpaired totals would not.

Reproducing these numbers on isolated hardware with locked frequency remains
outstanding work, and no README/CHANGELOG/release-note speedup claim may cite this
study without saying so.

---

## 3. Datasets

Three tiers. Tier A is the headline; Tier B is for cost-driver sweeps; Tier C is
the existing small synthetic generator, retained only for regression pinning.

### Tier A — real-model panel (headline)

Annotation is the **real** UCSC `hg38.ncbiRefSeq` BED12 gene model (real transcript
boundaries, real exon/intron structure, real contig naming). Read-bearing sequence is
drawn from the **real** UCSC hg38 human chromosomes.

Reads are **simulated from the real transcript sequences** with a realistic
Illumina-like quality profile, not random `ACGT` noise. This is a deliberate,
declared limitation: the *annotation and reference sequence are real*, but the
reads are not a real sequencer output. It is used in place of a real aligned BAM
because no suitably sized public aligned BAM with a matching BED12 was obtainable
in this environment; see §9.

| ID | Contigs | Transcript source | Role |
|---|---|---|---|
| `A-small` | chr17 | real RefSeq BED12 | small: exposes fixed-cost dominance |
| `A-mid` | chr1 + chr17 | real RefSeq BED12 | mid: typical QC run |
| `A-large` | chr1 + chr11 + chr17 | real RefSeq BED12 | large: asymptotic behaviour |

### Tier B — cost-driver sweeps

Same real model, varying exactly one driver at a time (§6).

### Tier C — existing synthetic generator

`benchmarks/generate_workload.py`. Retained for pinning historical numbers. **Not
used for any claim**: it emits one contig, 3 BED12 genes and fixed 100 bp
fragments, so it exercises none of the cost drivers in `testing.md` §12.3.

---

## 4. Experiment classes

Reported separately; never merged into one number.

| ID | Definition | Included |
|---|---|---|
| **E1 compute-only** | The command's own algorithmic work. Python **interpreter startup and import cost are excluded** by measuring and subtracting a per-invocation floor. | All commands |
| **E2 end-to-end** | Whole process tree as a user experiences it, **equal deliverables**: if upstream produces a plot or a BigWig, the candidate must too. | All commands |
| **E3 setup** | Index building / conversion cost, timed separately. | Commands with an index or setup step |

**E1 is the honest "language/algorithm" number. E2 is the honest "what the user
waits for" number.** Both are reported; E2 includes a fixed cost that has nothing
to do with the algorithm, and E1 excludes a fixed cost a real user always pays.
Reporting only one of them is a protocol violation.

### 4.1 The Python floor (measured, not assumed)

Upstream RSeQC pays a large per-invocation fixed cost before any input is read.
This study measures it explicitly rather than assuming it, because preliminary
probing showed it is not negligible and that it has a surprising property:

- Interpreter + import floor: ~0.10 s wall.
- **Additional ~1.0 s of CPU time per invocation is burned by OpenBLAS thread-pool
  initialisation** triggered by `import numpy` (measured: `import numpy` costs
  0.08 s wall but 0.90 s user on a 16-thread machine; with
  `OPENBLAS_NUM_THREADS=1` it costs 0.05 s user). Upstream RSeQC contains no
  `multiprocessing` and no explicit threading; this is dependency-internal.
- Consequence: an unpinned CPU-time column for upstream is ~90% BLAS
  initialisation noise, and wall-vs-CPU disagreement can be misread as parallelism
  that does not exist.

**Policy:** parity runs execute with `OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and `NUMEXPR_NUM_THREADS=1`, set
identically for both arms. A separate, clearly-labelled "upstream default env"
series may be run to characterise that default; it is never compared against a
pinned candidate.

---

## 5. Repetitions, ordering, and estimator

- **Pilot:** ≥3 untimed observations per (command, workload) to choose sizes and
  detect timeouts. Pilot data are not reported as results.
- **Measured:** **≥10 paired repetitions** for every primary deterministic
  comparison. Increased where pilot variance demands.
- **Ordering:** reference and candidate are **randomised and interleaved** within a
  block, one pair at a time, not run back-to-back. The schedule is generated from a
  recorded seed before collection and is not adjusted after seeing results.
- **Cache:** warm-cache only, and labelled as such. A fresh process is **not** a
  cold-cache run; no cache flushing is performed on a shared machine.
- **Estimator:** speedup = median(reference) / median(candidate) in wall time, with
  a 95% CI from a **paired block bootstrap** (10 000 resamples, resampling whole
  pairs, seed recorded). The exponentiated mean of paired log time ratios is also
  computed and reported as a secondary estimator; where the two disagree, both are
  reported.
- **A positive speedup claim requires the CI to lie entirely above 1.0.** An
  interval spanning 1.0 is reported as **inconclusive**, not as a win. Intervals
  below 1.0 (i.e. the port is slower) are reported as losses.
- No early stopping after a favourable measurement. The full schedule runs.
- **Multiplicity:** the per-command rows are the primary result and are not
  corrected against each other. Any summary claim spanning multiple commands is
  labelled **exploratory**.

---

## 6. Cost drivers swept (Tier B)

One driver varied at a time, all others held fixed, following `testing.md` §12.3:

1. Read count (total)
2. Transcript count (number of BED12 records)
3. Covered bases / mean depth
4. Mean transcript length and exon count
5. Contig count
6. Read length
7. Output volume (for the whole-output-buffering commands)

Each sweep records the **largest operating point measured** and the behaviour at or
beyond the resource limit. Scaling exponents are reported where a power law fits.

---

## 7. Equivalence gate (a speedup cannot be earned by doing less work)

A (command, workload) pair yields a reportable speedup **only if the equivalence
gate passes**. Failing the gate marks the row `unsupported` or `failed`; it never
yields a speedup.

The gate checks, in order:

1. Both arms exit 0.
2. **Structural output comparison**, not exit codes and not raw bytes alone. Text
   tables are compared line-wise with numeric-aware tolerance at the precision both
   tools actually emit. BAM outputs are compared as decoded record multisets
   (name, flag, reference, position, CIGAR, sequence, quality) — not as compressed
   bytes, which is meaningless across writers. BigWig/bedGraph compared as
   interval-value maps. FASTQ compared record-wise.
3. Any row where a *documented intentional divergence* applies
   (`compatibility/divergences.yaml`, e.g. RNG differences, timestamped log lines,
   BAM byte layout) is marked `diverged-documented` and reports timing but is
   excluded from the speedup table's headline statistics, with the divergence named.
4. Timeouts, OOMs, crashes and incomplete outputs are **retained in the table**.
   They are never replaced by an invented runtime, and the successful-only subset is
   never averaged without saying so.

This gate is deliberately stricter than the existing `benchmarks/run_benchmarks.py`,
which compares exit codes plus a raw file diff and covers only 5 commands.

---

## 8. Known exclusions, declared before measurement

| Command(s) | Reason | Handling |
|---|---|---|
| `FPKM-UQ.py` | Requires `htseq-count`, absent in this environment. Only a mock exists, which makes the dominant cost O(1) and the measurement meaningless. | Row reported as `unsupported-missing-dependency`. Not omitted. |
| `bam2wig.py` | Requires `wigToBigWig`, absent. Upstream also fails (silently, inside `subprocess.call(shell=True)` in a bare `except`). | Benchmarked in WIG-text mode only, with the absent BigWig deliverable stated on the row. No equal-deliverable BigWig claim. |
| `sc_editMatrix.py`, `sc_seqQual.py` | `pheatmap` R package absent. | Run with `--skip-heatmap`; the heatmap deliverable is excluded and stated. |
| 12 Rscript commands | R present, but end-to-end numbers would be dominated by Rscript startup, not RSeQC. | E1 (compute-only) is the primary reported figure; E2 measured separately and labelled as including Rscript. |
| All commands | Thread scaling impossible; port has no `--threads` and no first-party parallelism. | Recorded as a gap (§1.2). |
| All commands | No cold-cache series. | Warm-cache only, declared. |

---

## 9. Declared limitations of the dataset tier

- Reads are simulated from real transcript sequences, not a real sequencer run.
  Therefore: no platform-specific artefacts, no real PCR duplicates beyond
  coordinate collision, no real chimeric reads. Commands whose cost driver is
  *duplicate structure* (`read_duplication`) or *error/quality artefacts* are
  therefore measured on a favourable input and their ratios must be read as
  compute-only, not as a real-sample prediction.
- No real single-cell dataset is in the panel. The `sc_*` commands are run on
  synthetic single-cell inputs and their rows are labelled as such.
- The environment is shared (§2).

These limitations are part of the protocol, not an appendix, because they bound
which claims the resulting numbers can support.

---

## 10. Preregistered expectations

Stated before measurement so that a result contradicting them is a finding rather
than a redefinition.

- **E1 speedups will be far smaller than E2 speedups** for every command, because
  the candidate's win is diluted by a fixed Python floor that grows in relative
  importance as inputs shrink.
- **Speedup will decrease monotonically with input size** and approach a
  large-input asymptote. A single-point measurement is therefore expected to be
  uninformative, and the asymptotic ratio is expected to be the defensible number.
- **Peak memory will favour the port on the 12 streaming commands** and may
  **favour upstream on the 14 commands that decode the whole input eagerly**
  (`open_alignments` in `crates/formats/src/lib.rs` collects the entire file into a
  `Vec`; SAM/CRAM additionally round-trip through an in-memory BAM buffer). At
  least one such row is expected to be a loss or a wash.
- **The whole-output-buffering commands** (`bam2wig`, `normalize_bigwig`,
  `overlay_bigwig`, `inner_distance`, `tin`, `read_hexamer`) are expected to show
  output-volume-dependent memory growth.
- `read_hexamer` is expected to scale worse than linear in read length, because it
  allocates a `String` per base position (`read_hexamer.rs:83`).
- At least one command is expected to be **slower or equal** in the port, given
  that upstream delegates to `pysam`/`htslib` C and to `pyBigWig`, whereas the port
  uses pure-Rust `noodles`/`bigtools` and pays full BGZF decode cost in Rust. The
  literature (Rust vs C++ in bioinformatics) finds library maturity, not language,
  dominates this axis, so a loss is a live possibility and is not a failure of the
  study.

## 11. Reproducibility

Raw per-run measurements, environment manifests, the schedule seed, and the analysis
scripts that regenerate every table and figure are committed. Workload data and
large results are gitignored but regenerable from the recorded seed and manifest.
