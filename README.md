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

## Comparison with RustQC

Shared-hardware warning: the timing numbers below were measured on a shared machine (16 logical CPUs, 31 GB); every timing row records the load average before the run. Do not pool these numbers with any other study. Full evidence: `benchmarks/rustqc-comparison/RESULTS.md` (generated by `report.py` from `raw/` only).

## 2. Scope: tools covered by each implementation

| Implementation | RSeQC commands covered |
|---|---|
| Upstream RSeQC 5.0.5 (U) | 33 |
| This repository (P) | 33 |
| RustQC `rna` (R) | 8 (bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance, tin) |

## 3. Equivalence: tool x workload (P vs U, R vs U)

Verdicts: `byte-identical` > `numerically-identical` > `within-tolerance` > `not-comparable` / `differs` / `missing` (worst file per tool). Tolerance 0 unless noted. `not-comparable` only for junction_saturation (stochastic shuffle, no seed).

| Tool | W-rat-3c P | W-rat-3c R | W-hum-3c P | W-hum-3c R | W-rat-wg P | W-rat-wg R |
|---|---|---|---|---|---|---|
| bam_stat | byte-identical | byte-identical | byte-identical | byte-identical | byte-identical | byte-identical |
| infer_experiment | byte-identical | differs | byte-identical | differs | byte-identical | differs |
| read_duplication | byte-identical | byte-identical | byte-identical | byte-identical | byte-identical | byte-identical |
| read_distribution | byte-identical | differs | byte-identical | differs | byte-identical | differs |
| junction_annotation | byte-identical | differs | differs | differs | byte-identical | differs |
| junction_saturation | not-comparable | not-comparable | not-comparable | not-comparable | not-comparable | not-comparable |
| inner_distance | byte-identical | differs | byte-identical | byte-identical | byte-identical | differs |
| tin | differs | differs | missing | missing | differs | differs |

- W-hum-3c timeouts during equivalence: U:tin.py (recorded, not crashed).
- Detail (largest differences, first lines; full text in raw/equivalence.*.json):
  - W-rat-3c infer_experiment R=differs: SRR1177982.infer_experiment.txt: largest absolute 0.2618, largest relative 0.522764; first differences: / line count 5 vs 6 / line 2: 'This is PairEnd Data' vs '' / line 3: 'Fraction of reads failed to determine: 0.0009' vs 'This is PairEnd Data' / line 4: 'Fraction of reads explained by "1++,1--,2+-,2-+": 
  - W-rat-3c read_distribution R=differs: SRR1177982.read_distribution.txt: largest absolute 3.34404e+08, largest relative 0.659064; first differences: / line 3: 'Total Assigned Tags           8895639' vs 'Total Assigned Tags           9390784' / line 6: 'CDS_Exons           10963531            7425160             677.26            ' vs 'CDS_Ex
  - W-rat-3c junction_annotation R=differs: SRR1177982.junction.xls: differs: largest absolute 0, largest relative 0; first differences: / line 34989: 'chr2\t54832670\t54833250\t322\t annotated' vs 'chr2\t54832670\t54833250\t322\t complete_novel' / line 34990: 'chr2\t54833386\t5483 | SRR1177982.junction.bed: differs: largest absolute 0, largest r
  - W-rat-3c inner_distance R=differs: SRR1177982.inner_distance.txt: differs: largest absolute 1.01224e+06, largest relative 0.999831; first differences: / line 166: 'SRR1177982.10971043\t-52\tsameTranscript=No,dist=genomic' vs 'SRR1177982.10971043\t-52\treadPairOverlap' / line 167 | SRR1177982.inner_distance_freq.txt: differs: largest abso
  - W-rat-3c tin P=differs: SRR1177982.tin.xls=differs; SRR1177982.summary.txt=differs: largest absolute 3.69091, largest relative 0.0409556; first differences: / line 157: 'NM_013162\tchr1\t256806475\t256813678\t93.75219481388487' vs 'NM_013162\tchr1\t256806475\t256813678\t93.85965648789987' / line 654: 'NM_012559\tchr2\t18198723
  - W-rat-3c tin R=differs: SRR1177982.tin.xls: differs: largest absolute 2.7581e+08, largest relative 1; first differences: / line count 5360 vs 5181 / line 6: 'NM_013134\tchr2\t27480225\t27500654\t93.92694940851895' vs 'NM_013134\tchr2\t27480225\t27500654\t93 | SRR1177982.summary.txt: differs: largest absolute 1.52714, largest r
  - W-hum-3c infer_experiment R=differs: ERR10229623.infer_experiment.txt: largest absolute 0.1663, largest relative 0.401206; first differences: / line count 5 vs 6 / line 2: 'This is PairEnd Data' vs '' / line 3: 'Fraction of reads failed to determine: 0.1546' vs 'This is PairEnd Data' / line 4: 'Fraction of reads explained by "1++,1--,2+-,2-+":
  - W-hum-3c read_distribution R=differs: ERR10229623.read_distribution.txt: largest absolute 6.0416e+06, largest relative 1; first differences: / line 3: 'Total Assigned Tags           12288501' vs 'Total Assigned Tags           12282960' / line 6: 'CDS_Exons           23789748            9932742             417.52            ' vs 'CDS_Exons  
  - W-hum-3c junction_annotation P=differs: ja.junction.xls=byte-identical; ja.junction.bed=byte-identical; ja.junction.Interact.bed=differs
  - W-hum-3c junction_annotation R=differs: ERR10229623.junction.xls: byte-identical: 2569140 bytes identical | ERR10229623.junction.bed: byte-identical: 5198771 bytes identical | ERR10229623.junction.Interact.bed: differs: largest absolute 0, largest relative 0; first differences: / line 1: 'track type=interact name="Splice junctions" descript
  - W-hum-3c tin P=missing: U tin timed out; U outputs partial, no verdict
  - W-hum-3c tin R=missing: U tin timed out; U outputs partial, no verdict (R files listed, not judged)
  - W-rat-wg infer_experiment R=differs: SRR1177982.infer_experiment.txt: largest absolute 0.2604, largest relative 0.518622; first differences: / line count 5 vs 6 / line 2: 'This is PairEnd Data' vs '' / line 3: 'Fraction of reads failed to determine: 0.0004' vs 'This is PairEnd Data' / line 4: 'Fraction of reads explained by "1++,1--,2+-,2-+": 
  - W-rat-wg read_distribution R=differs: SRR1177982.read_distribution.txt: largest absolute 3.34404e+08, largest relative 0.659064; first differences: / line 3: 'Total Assigned Tags           8407135' vs 'Total Assigned Tags           8643988' / line 6: 'CDS_Exons           10963531            7310437             666.80            ' vs 'CDS_Ex
  - W-rat-wg junction_annotation R=differs: SRR1177982.junction.xls: differs: largest absolute 0, largest relative 0; first differences: / line 73351: 'chr2\t54832670\t54833250\t318\t annotated' vs 'chr2\t54832670\t54833250\t318\t complete_novel' / line 73352: 'chr2\t54833386\t5483 | SRR1177982.junction.bed: differs: largest absolute 0, largest r
  - W-rat-wg inner_distance R=differs: SRR1177982.inner_distance.txt: differs: largest absolute 667258, largest relative 0.999831; first differences: / line 7529: 'SRR1177982.10075785\t-16\tsameTranscript=No,dist=genomic' vs 'SRR1177982.10075785\t-16\treadPairOverlap' / line 7798: ' | SRR1177982.inner_distance_freq.txt: differs: largest abso
  - W-rat-wg tin P=differs: SRR1177982.tin.xls=differs; SRR1177982.summary.txt=differs: largest absolute 3.65601, largest relative 0.0405751; first differences: / line 157: 'NM_013162\tchr1\t256806475\t256813678\t93.71165044454774' vs 'NM_013162\tchr1\t256806475\t256813678\t93.81863111963806' / line 654: 'NM_012559\tchr2\t18198723
  - W-rat-wg tin R=differs: SRR1177982.tin.xls: differs: largest absolute 2.7581e+08, largest relative 1; first differences: / line count 5360 vs 5181 / line 6: 'NM_013134\tchr2\t27480225\t27500654\t93.93058017243231' vs 'NM_013134\tchr2\t27480225\t27500654\t93 | SRR1177982.summary.txt: differs: largest absolute 1.32015, largest r

## 4. Total cost for the eight tools (wall, CPU, peak memory)

Medians with (min-max) across repetitions in `raw/`; n shown. No confidence interval except where n>=5 (then range is still shown, not a CI). A single upstream repetition is labelled as a single run. P-par memory is the sampled aggregate (baseline-before-start method); U-seq/P-seq memory is the max single-process peak (sequential). R rows include non-RSeQC tools (see section 6).

| Workload | Config | n | Wall s median (min-max) | CPU s median (min-max) | Peak MB median (min-max) |
|---|---|---|---|---|---|
| W-rat-3c | U-seq | 3 | 320.3 (314.3-359.9) | 313.0 (310.0-345.9) | 1093.8 (1093.8-1094.1) |
| W-rat-3c | P-seq | 5 | 61.9 (59.5-72.4) | 61.1 (58.6-66.9) | 530.8 (530.2-1060.3) |
| W-rat-3c | P-par | 5 | 22.9 (22.1-47.7) | 61.3 (59.8-96.3) | 813.8 (348.0-4729.8) |
| W-rat-3c | R-1 | 5 | 59.2 (56.9-71.6) | 58.5 (56.9-66.9) | 747.4 (747.2-747.6) |
| W-rat-3c | R-8 | 5 | 36.3 (34.2-44.5) | 70.0 (66.4-70.5) | 1104.3 (1102.6-1105.3) |
| W-hum-3c | U-seq | 3 | 3771.1 (3770.6-3771.1) | 167.8 (167.3-167.9) | 1083.0 (1082.9-1083.4) |
| W-hum-3c | P-seq | 5 | 158.5 (157.0-159.4) | 158.1 (156.6-159.0) | 531.0 (530.8-531.1) |
| W-hum-3c | P-par | 5 | 113.5 (113.0-119.1) | 161.1 (160.4-168.3) | 815.0 (642.3-900.5) |
| W-hum-3c | R-1 | 5 | 177.7 (176.2-186.1) | 177.8 (176.2-186.1) | 1519.0 (1518.6-1519.3) |
| W-hum-3c | R-8 | 5 | 125.4 (124.7-125.6) | 192.9 (192.1-193.6) | 1784.3 (1775.4-1784.8) |
| W-rat-wg | U-seq | 1 (single run) | 645.9 (645.9-645.9) | 620.1 (620.1-620.1) | 3606.7 (3606.7-3606.7) |
| W-rat-wg | P-seq | 5 | 156.9 (156.4-166.9) | 156.5 (156.1-166.0) | 1322.6 (1322.3-1322.8) |
| W-rat-wg | P-par | 5 | 48.4 (48.2-50.0) | 159.8 (159.3-164.7) | 1353.2 (1066.7-1638.2) |
| W-rat-wg | R-1 | 5 | 156.4 (155.7-210.8) | 156.4 (155.7-194.0) | 1557.2 (1557.2-1557.5) |
| W-rat-wg | R-8 | 5 | 48.4 (48.1-48.6) | 193.0 (191.7-194.5) | 3099.6 (3052.4-3104.0) |

## Limitations

- **Sequence-logo images are from-scratch renders** with approximate glyphs and no `shade_below`/`fade_below` support; no upstream reference image is available ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **CRAM input is bounded by available memory** (whole-file decode), and CRAM files needing external reference resolution fail rather than fetching it ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **Held-out scientific validation is incomplete**, and three previously reported endpoint verdicts have been withdrawn ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **No publication-grade performance figures exist yet**; measured rows are single-session and indicative ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **`read_duplication` (~130 bytes/record) and `bam2wig` (~208 bytes/record, 1.69 GB on the full rat alignment) bound workload size**; no per-command memory limit is declared ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **Out-of-memory failures abort with Rust's internal allocation message**, not a diagnostic naming the command and input ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **`geneBody_coverage` and `tin` diverge from upstream past pysam's `max_depth` of 8000** (DIV-0024, open); no speedup or scientific claim is made for either ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **Property and fuzz suites are bounded and deterministic**, not coverage-guided fuzzing or command-wide correctness evidence ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **Real-data timing covers 18 commands on one alignment in one session** (16 of 18 gates pass); it is indicative, not generalisable ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **Fifteen commands require `Rscript` at runtime** (`FPKM_UQ` needs `htseq-count` instead); only a Linux x86_64 archive is produced and clean-room testing ran once on Linux ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **The Python `qcmodule` API layer is an unimplemented stub** ([details](docs/KNOWN_LIMITATIONS.md#limitations)).
- **Publication is blocked**: repository URLs are still `TBD`, so strict release metadata refuses to build ([details](docs/KNOWN_LIMITATIONS.md#limitations)).

## Development

- `docs/PORTING_PLAN.md` — the governing multi-step plan for this whole project.
- `testing.md` — the authoritative test/verification punch list.
- `compatibility/` — command inventory, API inventory, the frozen upstream baseline lock, and the
  full divergence log.
- `verification/run_diff.py` — the differential verification harness. Run with
  `oracle/venv/bin/python3 verification/run_diff.py` (requires the gitignored `oracle/` Python
  environment; see `docs/PORTING_PLAN.md` for how it's set up).
