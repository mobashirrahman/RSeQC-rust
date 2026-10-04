# RSeQC-rust

A from-scratch Rust reimplementation of [RSeQC](https://github.com/liguowang/RSeQC), the
RNA-seq quality-control toolkit, targeting a **drop-in replacement**: standalone binaries under
the original command names, no Python/R/helper-binary runtime dependency for the core CLI
workflow.

## Status

All 33 upstream commands are implemented, build cleanly (`cargo build --workspace`), and pass
`cargo test --workspace` / `cargo clippy --workspace --all-targets -- -D warnings`. Compatibility
is verified with a 90-case differential harness (`verification/run_diff.py`) that runs the real
upstream Python CLI and this port's compiled binary against the same fixture and diffs their actual
output — not just unit tests against this port's own expectations. The latest recorded run, on
2026-10-01 after the comparators were made failure-closed, is **90 of 90 passing**. The previous
invocation was 86 of 90 with four skipped for want of a fixture. Both combine substituted real
inputs with fixed synthetic/regression fixtures, so this is 90 compatibility checks and not 90
real-data validations. This is compatibility evidence for the exercised cases; see
[Limitations](#limitations) and the [readiness audit](docs/READINESS_AUDIT_2026-10-01.md).

**This is not yet a finished, published release.** See [Limitations](#limitations) below for what
still needs work before that's a fair claim, and `testing.md` / `docs/PORTING_PLAN.md` for the
full, authoritative punch list this summary is derived from.

## Upstream baseline

Ported against a frozen snapshot recorded in `compatibility/upstream.lock`:

- Source: `https://github.com/liguowang/RSeQC.git`
- Commit: `59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24`
- Declared package version: `5.0.5`

Compatibility claims in this repository are relative to that exact snapshot, not to "RSeQC" in
general — a later upstream release may behave differently. See `docs/PORTING_PLAN.md` Step 12
for how upstream-version drift is meant to be tracked going forward.

## Installation (build from source)

```sh
git clone <this-repo>
cd RSeQC-rust
cargo build --workspace --release
scripts/install-aliases.sh   # creates bam_stat.py, split_paired_bam.py, ... in target/release/
```

Add `target/release` to `PATH`, or copy its contents somewhere that already is. Running a command
by its original upstream name now runs this port's compiled binary — no Python interpreter or
`qcmodule` import needed for the core CLI.

A checksummed, self-contained distribution archive (binaries + aliases, no source tree required)
can be built with:

```sh
scripts/build-release-archive.sh
```

## Usage

Commands take the same flags as upstream. For example:

```sh
bam_stat.py -i sample.bam -q 30
read_distribution.py -i sample.bam -r reference.bed12
RPKM_saturation.py -i sample.bam -r reference.bed12 -o sample_saturation
```

Run any command with `--help` for its full flag list.

## Compatibility

| Command | Purpose | Status |
| --- | --- | --- |
| `bam_stat.py` | Summarize mapping statistics for a BAM or SAM alignment file | ✅ |
| `bam2fq.py` | Convert alignments in BAM or SAM format to FASTQ | ✅ |
| `divide_bam.py` | Randomly divide a BAM file into approximately equal subsets | ✅ |
| `split_bam.py` | Split a BAM file according to exon regions from a BED gene list | ✅ |
| `split_paired_bam.py` | Split a paired-end BAM into read-1, read-2, and unmapped BAM files | ✅ |
| `read_GC.py` | Calculate the GC-content distribution of aligned reads | ✅ |
| `read_NVC.py` | Calculate nucleotide frequency at each read cycle | ✅ |
| `read_quality.py` | Calculate per-cycle Phred quality-score distributions for aligned reads | ✅ |
| `read_duplication.py` | Calculate sequence-based and mapping-based read duplication rates | ✅ |
| `clipping_profile.py` | Estimate the clipping profile of RNA-seq reads from a BAM or SAM file | ✅ |
| `deletion_profile.py` | Calculate the distribution of deleted nucleotides across aligned reads | ✅ |
| `insertion_profile.py` | Calculate the distribution of inserted nucleotides across RNA-seq reads | ✅ |
| `mismatch_profile.py` | Calculate the distribution of mismatches across aligned reads | ✅ |
| `read_hexamer.py` | Calculate normalized hexamer frequencies from FASTA or FASTQ files | ✅ |
| `infer_experiment.py` | Infer RNA-seq library layout and strandedness from a SAM/BAM file | ✅ |
| `read_distribution.py` | Summarize read distribution across genomic annotation categories | ✅ |
| `inner_distance.py` | Estimate the inner distance between paired-end RNA-seq reads | ✅ |
| `RNA_fragment_size.py` | Calculate fragment-size statistics for each transcript or gene | ✅ |
| `junction_annotation.py` | Annotate splice junctions against a reference gene model | ✅ |
| `bam2wig.py` | Convert a sorted, indexed BAM file into WIG coverage files | ✅ |
| `geneBody_coverage.py` | Calculate RNA-seq read coverage across the gene body from BAM input | ✅ |
| `geneBody_coverage2.py` | Calculate RNA-seq coverage across the gene body from a BigWig file | ✅ |
| `normalize_bigwig.py` | Normalize a BigWig signal to a fixed total WIG sum | ✅ |
| `overlay_bigwig.py` | Apply an arithmetic operation to two BigWig signal tracks | ✅ |
| `FPKM_count.py` | Calculate fragment counts, FPM, and FPKM for BED12 transcript models | ✅ |
| `FPKM-UQ.py` | Calculate raw counts, FPKM, and upper-quartile normalized FPKM (TCGA-compatible) | ✅ |
| `RPKM_saturation.py` | Assess whether transcript RPKM estimates have reached sequencing saturation | ✅ |
| `junction_saturation.py` | Assess whether splice-junction discovery has reached sequencing saturation | ✅ |
| `tin.py` | Calculate transcript integrity number (TIN) per transcript or gene | ✅ |
| `sc_bamStat.py` | Report mapping statistics for single-cell RNA-seq BAM files | ✅ |
| `sc_editMatrix.py` | Visualize error-correction edits in cellular barcodes and UMIs | ✅ |
| `sc_seqQual.py` | Generate sequencing-quality matrices and a heatmap from a FASTQ file | ✅ |
| `sc_seqLogo.py` | Generate a DNA sequence logo from FASTA, FASTQ, or sequence-only input | ✅ |

✅ means the command's differential harness case(s) pass byte-identical (or equivalent, where a
format's own container isn't byte-comparable in principle — e.g. gzip/BAI headers) against real
upstream, for every case currently exercised. It does **not** mean every possible input/flag
combination has been tried — see [Limitations](#limitations).

`sc_seqLogo.py` is now fully implemented: its `.count_matrix.csv` output is computed and verified
byte-identical against real upstream, and `--oformat svg`/`png`/`pdf` all produce real, valid
sequence-logo images via native from-scratch renderers (`crates/render`). The renderers' known
gaps (approximate glyph placement, no `shade_below`/`fade_below` support) are documented in
`crates/render/src/seqlogo.rs`/`seqlogo_png.rs`/`pdf.rs`. No reference image is byte-comparable:
the pinned upstream's own `logomaker`+`pandas` combination crashes before rendering for every
format, so the differential case checks the byte-identical CSV plus candidate-side PDF artifacts.
See DIV-0016 in `compatibility/divergences.yaml`.

Every known, intentional behavioral difference from upstream — including ones that are permanent
by design (e.g. RNG algorithm differences, timestamped log lines) — is recorded in
`compatibility/divergences.yaml` with its own rationale and verification evidence. That file, not
this table, is the authoritative compatibility record.

## Limitations

- **`sc_seqLogo.py`'s sequence-logo image rendering is from-scratch and shares documented visual
  caveats** (DIV-0016) — `--oformat svg`/`png`/`pdf` are all implemented via native renderers
  (from-scratch; `pdf` wraps the PNG renderer's raster in a minimal single-page PDF image
  container). No two independent renderers can byte-match matplotlib's own image output regardless;
  a real reference image to visually compare against isn't currently available either (the upstream
  oracle used during development can't produce one in its own environment, due to an unrelated
  `logomaker`/`pandas` version incompatibility). Renderers use approximate glyphs and don't honor
  `shade_below`/`fade_below`.
- **`.cram` input is now supported** (all 14 commands that accept `.bam`/`.sam` also accept
  `.cram`), decoded with no external reference file — the common case, matching what `pysam`/
  `htslib` themselves fall back to when writing CRAM without one configured. A CRAM file that
  genuinely requires external reference resolution will fail to decode rather than silently
  producing wrong data; this port does not fetch references over the network the way `htslib` can,
  by design (see DIV-0002/0004).
  The CRAM path buffers the whole file (`AlignmentRecords::CramBuffered` in
  `crates/formats/src/lib.rs`), so CRAM input is bounded by available memory rather than
  by alignment size. This is a code fact, not a timing measurement: the CRAM fixture is
  10 records. See [`docs/ENVELOPE.md`](docs/ENVELOPE.md).
- **Scientific validation remains incomplete, and three verdicts have been withdrawn.** The
  differential suite is 90/90 (2026-10-01, failure-closed comparators), on the
  development panel with retained synthetic/regression fixtures. Separately, the held-out endpoint
  validator scores Rust binary output against pre-registered endpoints; it does not run
  upstream on held-out inputs, so it cannot compare the two implementations there.
  Three things were wrong with how its results were read, and are now corrected:
  - **E1's cross-lab PASS is retracted.** Its strand expectation came from ENA
    `library_selection`, which records how a library was amplified, not its strandedness.
    Without explicit protocol metadata E1 is now NOT_EVALUATED rather than guessed.
  - **A6 is withdrawn.** The rat junction "annotation density" explanation was written
    after the rat result was seen and rested on a mis-stated coverage figure (32% from
    transcript *spans*, which count introns as annotated; from merged exon bases it is
    1.8%).
  - **The rat annotation conversion was wrong and has been corrected**, independently
    confirmed against the genome sequence: UCSC refGene is half-open 0-based, not
    1-based inclusive. With the corrected BED12 and a rebuilt STAR index, the rat
    annotated-junction fraction is 0.585 against 0.026 before, clearing the
    pre-registered 0.50 bar.
  - **E3's estimand is narrower than a mechanistic reading.** It correlates one library's
    gene-body curve with the development panel's, which is a similarity between two
    samples and not evidence about any mechanism. `--estrand` now names the claim.
  - **Feature stratification now exists** (`--stratify`, by coverage, transcript length,
    exon count, GC and annotation ambiguity) and already shows that the aggregate E3 pass
    conceals stratum-level variation.
  All three held-out runs have been inspected and are consumed for confirmation; anything
  they prompt needs fresh independent data. See
  [the endpoint results](datasets/ENDPOINT_RESULTS.md) and `testing.md` sections 10-11.
- **Performance has been benchmarked, but not to publication standard.** A version-1
  preregistered run over 29 commands is in `benchmarks/RESULTS.generated.md`. It was
  measured on shared, non-isolated hardware, and the independent audit found defects in
  the harness that produced it: comparator false-pass paths, lost repetition pairing,
  and an output gate that accepted two empty directories. All are now repaired and
  covered by executable tests (`benchmarks/test_bench_harness.py`), and the replacement
  study is specified in [`benchmarks/protocol-v2.md`](benchmarks/protocol-v2.md). **No
  speedup figure from the version-1 run may be published**, because its gate passes do
  not establish equivalence and its intervals describe unpaired measurements. The
  version-1 relative orderings remain useful engineering evidence. The benchmark did
  find three commands that were *slower* than upstream (`infer_experiment`, `bam2fq`,
  `inner_distance`), all since root-caused and fixed, and a memory regression -- the
  alignment reader decoded whole files up front, extrapolating to ~23 GB for a
  50M-read-pair BAM -- since fixed for BAM/SAM by streaming, with `tin`'s and
  `geneBody_coverage`'s own indexes replaced by a sliding window.
- **`FPKM_count.py` no longer materialises the whole BAM, and this was the last
  command that used materially more memory than upstream.** Upstream splits its
  work: `count_total_fragments` streams the file once, but `count_transcript`
  re-queries the BAI with `samfile.fetch(chrom, tx_start, tx_end)` **once per
  transcript**, so pysam never holds more than one region's reads and pays for it
  in seeks. This port had loaded the file once and then scanned a start-sorted
  prefix per transcript -- the mirror image, 3.68x faster and memory-hungry.
  `compute_fpkm_rows_windowed` now walks transcripts in coordinate order and
  streams the BAM once, discarding a read at push time as soon as it can no
  longer reach the current or any later transcript. `count_transcript` is
  unchanged -- the driver hands it exactly the reads the whole-file prefix scan
  would have -- and five tests assert the windowed driver reproduces the
  whole-file rows byte for byte, including input row order, an absent chromosome
  (which still yields a zero row), and the out-of-order detection that routes to
  the whole-file fallback. On the 8.2M-record rat alignment: **153 MB -> 17 MB,
  now 2.9x LESS than pysam's 49 MB** and 4.75x faster, byte-identical. On the
  2.05M-read human alignment: 6.2 MB against pysam's 42.9 MB and 4.4x faster,
  byte-identical. The per-read record was also halved, from four `i64`s and five
  `bool`s (40 bytes padded) to `i32` coordinates and bit-packed flags (20 bytes);
  the `i32` is exact, not a narrowing, because POS/endpos/PNEXT are int32 on the
  wire in a BAM.
- **`tin`'s sliding window had the same defect `geneBody_coverage`'s did: it
  buffered the inter-transcript gap and trimmed afterwards, so peak RSS tracked
  the distance from the chromosome start rather than the local depth.** The
  window driver already existed and its own doc comment claimed a 174x resident
  reduction, but the reduction was measured on a workload whose first transcript
  sits near position 0. On the 2.05M-read human alignment, whose first chr1
  transcript starts at 114 Mb, the initial pull buffered 114 Mb of reads before
  the trim ran: **331 MB against pysam's 42.7 MB (7.9x worse)**, deteriorating to
  **11.2 MB (3.8x BETTER than pysam)** once reads are discarded at push time,
  with `tin.xls` byte-identical and 28x faster. The same fix was applied to
  `geneBody_coverage` (373 MB -> 14 MB). Neither change alters which reads the
  scoring functions see, only when unreachable ones are dropped.
- **What remains of `tin`'s memory is a real ultra-deep region, not a window
  bug.** On the rat alignment `tin` still peaks at 219 MB against pysam's 51 MB,
  and the cause is measured rather than assumed: one 4 kb transcript at
  chr1:80,612,893 is overlapped by **275,288 reads (68 reads per base)**, and the
  window holds that region in full. `pysam`'s `pileup` caps depth at 8000 per
  column during iteration, so it never sees the rest; this port must keep them
  because the depth-cap behaviour itself is emulated later (DIV-0024) and the
  reads it needs cannot be known in advance. A correction to an earlier draft of
  this note: it attributed several `tin`/`geneBody_coverage` memory numbers to a
  degenerate gene model with chromosome-scale transcripts. That was wrong -- it
  came from an `awk '$3-$1'` on a BED whose first column is a chromosome NAME,
  which awk read as 0 and which inflated every span by the transcript's start
  coordinate. The models are ordinary (mean spans 33 kb and 38 kb, max 2.1 Mb).
  The measured numbers were real; only the explanation was mistaken.
- **The rat endpoint figure is confirmed, and a third reference-preparation defect
  was found while confirming it.** The corrected annotation and a rebuilt index were
  both on disk and every digest matched — but the rat BAM had been aligned hours
  before the index was rebuilt, so the annotation was current and only the
  *alignment* was stale. Re-running the endpoint command reproduced an earlier figure
  while the recorded table said a different one, and nothing in the repository could
  distinguish them: the alignment's `align.json` named the index's parameters but not
  its annotation digest. `datasets/align_run.sh` now records that digest and refuses
  to align against an index whose stamp disagrees with the annotation beside it;
  `datasets/build_star_index.sh` gives each assembly its own index directory, so
  rebuilding the rat index cannot overwrite the human one; and
  `verification/check_rat_reference.py` asserts the whole chain (27 checks). Both
  implementations were then run on the identical corrected inputs and produced
  byte-identical junction tables, stdout and stderr — see
  [`datasets/ENDPOINTS.md`](datasets/ENDPOINTS.md) §6.3.
- **Capacity is measured per cost driver, and no envelope is declared.**
  [`docs/ENVELOPE.md`](docs/ENVELOPE.md) records what is measured. On a synthetic
  fixture over a 29x range of record counts: `bam_stat` flat at ~2.9 MB (the control,
  −0.1 bytes/record), `read_duplication` 43.1 bytes/record, `bam2wig` 250.4. An earlier
  version of this document reported 0.6, 135.6 and 270.4 bytes/record; those were
  fitted over a 1.22x range whose two larger points held *identical* record counts,
  and they are withdrawn.
  On the real 8.2M-record rat alignment over a 174x range: `bam_stat` flat,
  `read_duplication` **130.5** bytes/record, `bam2wig` 207.5 — a three-fold difference
  from the synthetic fixture, which is the point: a per-record cost measured on
  generated reads is a property of the generator.
  Behaviour past the limit is measured too. Under `ulimit -v` at 80%, 40% and 20% of
  their limit-free peaks, both non-flat commands abort with **no partial output and
  no input damage**, so a consumer cannot read a truncated result as a complete one —
  but the diagnostic is Rust's internal `memory allocation of N bytes failed`, not a
  message naming the command and its input. That is a known gap, recorded rather than
  closed. No whole-genome or production-size run at the audit's 10M/50M-pair targets
  has been performed, and no per-command memory limit is claimed.
- **Bounded property/fuzz suites exist.** Five `proptest` suites and a campaign report
  are in [verification/FUZZING.md](verification/FUZZING.md), covering BED parsing,
  CIGAR traversal, Python numeric formatting, FASTA/FASTQ parsing and the compressed
  input layer, with fixed regression seeds. They are bounded and deterministic; they do
  **not** establish coverage-guided fuzzing, resource-bounded robustness, or
  command-wide scientific correctness, and "no panic" is not evidence of correct
  biological meaning. Fixed seeds run in CI; varied-seed campaigns with recorded seeds
  and resource limits belong after shared-parser changes.
- **The validation machinery itself is tested.** 153 unit tests assert that the
  comparators and gates still reject what they must: corrupted BAM qualities, flags,
  tags, mates and headers; a truncated FASTQ record; a finite metric replaced by NaN;
  two empty output trees; a timeout that must kill only its own process group; a
  `CITATION.cff` placeholder that a TOML-only pattern had silently exempted. Every
  one asserts a *failing* result, so a comparator or check that stops rejecting
  corruption turns CI red rather than turning a regression into a reported speedup.
  On top of that, 120 command-contract checks (missing, empty, unreadable, malformed
  and truncated input; invalid flags; missing sidecars; existing output; spaces in
  paths; unwritable output; killed and concurrent runs) and 37 interoperability
  checks against htslib. See
  [`benchmarks/test_bench_harness.py`](benchmarks/test_bench_harness.py),
  [`verification/test_comparators.py`](verification/test_comparators.py),
  [`verification/test_run_diff.py`](verification/test_run_diff.py),
  [`scripts/test_release_metadata.py`](scripts/test_release_metadata.py),
  [`verification/check_command_contracts.py`](verification/check_command_contracts.py)
  and [`verification/check_interop.py`](verification/check_interop.py).
- **There is now a real measurement for 18 command rows, and it says which of them may
  be quoted.** [`benchmarks/results-rat-real-8M/results.json`](benchmarks/results-rat-real-8M/results.json)
  holds 18 rows on an 8.2M-record rat alignment, 5 matched blocks each, with per-run raw
  timings, failures, the interleaved schedule and a per-binary SHA256 — and that hash
  was checked against the binaries on disk *after* the run, so the rows provably measured
  the code being shipped. **16 of 18 gates pass.** The two that fail are `geneBody_coverage`
  and `tin`, which share the disclosed `max_depth` defect; their speedups are recorded
  and explicitly not claimed, because a gate that refuses is the gate working.

  | | ratio | 95% CI | gate |
  |---|---|---|---|
  | `inner_distance` | 14.79× | 14.45–15.36 | pass |
  | `read_NVC` | 11.19× | 11.12–11.54 | pass |
  | `read_distribution` | 4.58× | 4.50–4.68 | pass |
  | `bam_stat` | 4.46× | 4.35–4.49 | pass |
  | `clipping_profile` | 4.11× | 3.68–4.15 | pass |
  | `read_hexamer` … `RNA_fragment_size` | 3.52×–1.54× | see the JSON | pass |
  | `geneBody_coverage` | 12.36× | 11.83–12.81 | **fail — not claimed** |
  | `tin` | 7.79× | 7.78–7.92 | **fail — not claimed** |

  What this does **not** license is stated rather than implied: it is **one session on
  one machine**, which the protocol says a bootstrap interval cannot make generalisable;
  it is **one dataset**; and 5 blocks per row is short of the 10 the protocol starts
  from. These are indicative figures for a scoped beta, not publication figures. (Those
  rows are on the three-contig panel, which is retained because it is the same data the
  re-collected rows and the memory work used, so the two sets are comparable.)

- **There is a whole-genome rat alignment now, and it showed the panel above was 20% of
  the library.** A 58-contig rn6 index (150,217 annotated junctions) aligned the same
  17,168,681 input reads to **13,849,121 uniquely mapped (80.67%)** with 11,388,194
  splices, against **3,487,314 (20.31%)** and 2,255,750 on the three-contig panel — most
  rat reads fall outside `chr1`/`chr2`/`chr10`. The whole-genome BAM holds **32,626,178
  measured records = 8,584,340 read pairs**, recorded as the two separate quantities the
  protocol requires rather than as either alone. On it, `junction_annotation` is
  **byte-identical** between upstream and the port on all three data artifacts —
  `junction.xls` 5.9 MB, `junction.bed` 12.1 MB, `junction.Interact.bed` 29.0 MB — across
  **152,476 transcript rows**, and both arms report `total = 11388194`, exactly STAR's
  own spliced-junction count. See
  [`datasets/heldout/endpoint_results/SRR1177982_wholegenome/`](datasets/heldout/endpoint_results/SRR1177982_wholegenome/). The harness itself also earned three corrections during
  this study — three commands were gated on declarations that did not describe them, and
  a relative `--output-dir` silently invalidated the first attempt entirely.
- **Two commands' real-data equivalence is unproven, and the release says so.**
  `geneBody_coverage.py` and `tin.py` do **not** currently reproduce upstream's output on
  real data, from **one shared cause**: they both reach it through the same pileup
  primitive, and the port budgets `max_depth` per position over *visited* pileups where
  upstream budgets over the pileup *buffer* — so `is_del` pileups consume budget that
  upstream's does not. `geneBody_coverage`'s 100-bin curve has 76 of 100 bins differing
  on the 8.2M-record rat alignment, by up to 839 reads in both directions; `tin`'s TIN
  on transcript `NM_013162` is `93.75219481388487` upstream against `93.85965648789987`
  here, and the sample summary inherits it. Identified cause, not yet fixed
  (`compatibility/divergences.yaml` DIV-0024, open), and both are in the archive
  manifest's `known_limitations`. Every other pileup filter has been verified equivalent
  one mechanism at a time. The 90-case differential suite passed while both were wrong,
  because its fixtures are two orders of magnitude too shallow to bind a cap of 8000 —
  the benchmark harness's structural gate is what caught it. **No speedup or scientific
  claim is made for either command.** Two differential cases now assert this divergence
  *with a floor*, so a fix trips them rather than silently retiring the disclosure, and
  `verification/fixtures/check_gene_body_divergence.sh` re-checks on CI that the defect
  is still present.
- **Every shipped command has a measured capability record, not a prose claim.**
  [`verification/capability_matrix.py`](verification/capability_matrix.py) produces
  [`benchmarks/capability-matrix.json`](benchmarks/capability-matrix.json) by
  *executing* each command: feeding it a real file of each declared input format,
  re-running it with a `PATH` containing no `Rscript`, and running it twice to compare
  artifacts. That record is embedded in the archive manifest's `per_command` block, and
  `scripts/check_release_metadata.py` refuses a manifest whose record does not cover
  exactly the binaries the archive ships, **or whose recorded per-command binary hash
  is not the hash the manifest records for that binary**. Coverage alone cannot catch
  a stale hash: all 33 records once existed while all 33 hashes were stale after a
  rebuild, and every coverage check passed.
  What it establishes: **15 of 33 commands require `Rscript`** at runtime and report
  its absence, and a sixteenth needs a helper but not Rscript (`FPKM_UQ`, whose helper
  is `htseq-count`) — the record names *which* helper, because deriving "requires
  Rscript" from a single generic "needs a helper" flag mislabels exactly the command
  that does not need R; **24 accept BAM, 14 accept SAM, and 9 are BAM-only because upstream's own
  `validate_args` rejects any other extension**; `junction_saturation`,
  `RPKM_saturation` and `divide_bam` are non-deterministic by design (unseeded
  resampling, as upstream — `divide_bam` becomes reproducible with `--seed`); and
  `FPKM_UQ`, `sc_editMatrix` and `sc_seqQual` cannot be judged here because
  `htseq-count` and R's `pheatmap` are absent. No unexplained non-determinism remains.
- **The benchmark gate's own declarations are verified against live runs, not
  trusted.** [`verification/verify_stream_declarations.py`](verification/verify_stream_declarations.py)
  re-derives every `EXPECTED_STREAMS` entry — stdout presence, declared labels or
  substrings, produced artifacts — by running both arms once, and fails on any
  mismatch. This exists because three commands were gated on rules that did not
  describe them while their arms agreed byte-for-byte: `bam2wig` (a one-line stdout),
  `RNA_fragment_size` (a 369 KB table on stdout that a label-parsing probe had mistaken
  for no stdout at all), and `bam2wig` again via its artifact declaration. A gate built
  on a wrong declaration cannot report anything trustworthy about the commands it gates.
- **The QC panel is also delivered as a workflow, and the workflow is tested by running
  it.** [`workflows/qc_panel/Snakefile`](workflows/qc_panel/Snakefile) wires the
  five-command panel into Snakemake, and
  [`workflows/qc_panel/test/run_panel.sh`](workflows/qc_panel/test/run_panel.sh) runs
  the whole DAG on a real alignment, re-runs the same five commands directly, and diffs
  every artifact — so a wrapper bug fails a test instead of quietly changing somebody's
  QC numbers. Running it found five real errors in the Snakefile, each now documented at
  the mistake: a missing `-r`; `-o X` producing `X.NVC.xls` rather than `X`;
  `read_quality --skip-plot` producing *only* an R script and no data table; an awk
  `print` outside `BEGIN` never running and yielding an empty summary that the checker
  then passed, because a loop over zero lines checks nothing; and awk's
  `getline var < file` returning the line rather than a count. A corresponding Bioconda
  recipe is at [`recipes/rseqc-rust/meta.yaml`](recipes/rseqc-rust/meta.yaml), built
  `--locked` and installing the third-party notices, with
  `scripts/test_bioconda_recipe.py` rendering it so a recipe conda-build would reject
  cannot sit unnoticed.
- **The QC panel is measured as a workflow, and shows per-command speedups are the
  wrong lever.** [`verification/measure_pipeline.py`](verification/measure_pipeline.py)
  times the five-command bulk QC panel against the 8.2M-record rat alignment:
  **48.8 s wall, 0.0135 CPU-hours, 1.06 GB peak per sample.** *(The wall and
  CPU-hour figures stand. The aggregate-memory figure in the JSON is marked
  SUPERSEDED: it was measured with a baseline sampled after the panels started,
  which the harness has since fixed, and the panel included `geneBody_coverage`,
  which is not currently qualified. Recollect before citing it.)* The two streaming
  commands this project optimised — `bam_stat` and `infer_experiment` — are
  **1.5% of that**; `geneBody_coverage` and `read_duplication` are 78%. Running
  whole panels in parallel scales throughput to **3.3x at four samples** and 4.8x at
  eight on 16 cores, with per-sample CPU rising from 0.0133 to 0.0214 CPU-hours as they
  contend for memory bandwidth. The deployment that works is parallelism across
  samples; a per-command speedup is bounded by where the time actually is. This is one
  sample on one machine, not a claim about total pipeline cost.
- **Per-job memory is measured additive, not assumed.**
  [`verification/measure_concurrency.py`](verification/measure_concurrency.py) samples
  aggregate system memory while N invocations run. On the 8.2M-record rat alignment,
  `bam2wig` peaks at 1.89 GB at ×1 and 11.99 GB at ×8 (1.12× and 7.11× a 1.69 GB
  per-job maximum), and wall time is nearly flat in concurrency — 8 samples finish in
  1.3× the time of one. That is the parallelism-across-samples deployment the audit
  asks about, and it is measured rather than asserted; the caveats, including a ~90 MB
  noise floor in `/proc/meminfo` on this shared machine, are recorded in
  [`docs/ENVELOPE.md`](docs/ENVELOPE.md). The per-transcript cost of
  `geneBody_coverage` and `junction_annotation` is bounded only from above — below
  the measurement's own repeatability at a 16x model range. `tin`'s is measured at
  ~1.7 KB per transcript. A slope is published only when the fitted quantity moves at
  least 3x the repeat-to-repeat spread, and every point in the record is the peak of
  three repeats with its spread reported.
- **Release packaging covers one platform, and cannot publish yet.** CI builds and
  tests on macOS and Windows and builds the Docker image, but only a Linux x86_64
  archive is produced; there are no macOS/Windows release archives and no Python wheel.
  The archive builder now stages an explicit allowlist of 33 commands, bundles the
  LICENSE, README, CHANGELOG, citation metadata and third-party notices, writes a
  manifest with per-file digests and a declared glibc floor, and smoke-tests the
  *extracted* archive outside the source tree on a real workload and on paths
  containing spaces. **Publication is blocked**: `Cargo.toml` and `CITATION.cff` still
  carry `TBD` repository URLs, because the permanent URL and support channel are
  maintainer decisions, and `scripts/check_release_metadata.py` refuses a strict build
  while they are unresolved. `release-validation.yml` runs on every push, every pull
  request and every tag, and `release.yml` cannot publish without it.
- **The Python API compatibility layer is not implemented yet.** The `rseqc-python` crate is
  currently a stub; importing upstream's `qcmodule` API is not supported.
- **The repository license is settled.** The root `LICENSE` carries the canonical GPLv3 text and
  this project declares **GPL-3.0-or-later** (DIV-0003 resolved 2026-10-01). Upstream's own
  license metadata is internally inconsistent — its README says GPL-3.0-or-later while a packaging
  classifier says GPLv2 — and that inconsistency is recorded rather than reproduced.
- **Clean-room testing has been done once, on Linux x86_64 only** (2026-09-28): the built
  distribution archive was extracted and run with `PATH` limited to its own `bin/` directory (no
  Python/R/`wigToBigWig`). Core commands and `sc_seqLogo.py` work standalone; `bam2wig.py` falls
  back gracefully without `wigToBigWig`; the 15 `Rscript`-using commands now report upstream's own
  "Rscript executable not found" message (DIV-0021) — use `--skip-plot`/`--skip-heatmap` to avoid
  needing R. Not repeated on macOS/Windows.

## Development

- `docs/PORTING_PLAN.md` — the governing multi-step plan for this whole project.
- `testing.md` — the authoritative test/verification punch list.
- `compatibility/` — command inventory, API inventory, the frozen upstream baseline lock, and the
  full divergence log.
- `verification/run_diff.py` — the differential verification harness. Run with
  `oracle/venv/bin/python3 verification/run_diff.py` (requires the gitignored `oracle/` Python
  environment; see `docs/PORTING_PLAN.md` for how it's set up).
