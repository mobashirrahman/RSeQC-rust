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
