# T4 real-data panel

Real sequenced reads, aligned by a pinned third-party aligner. T4.1 compares
this port with upstream RSeQC on the development panel; T4.2 scores the Rust
port's output against scientific endpoints on held-out samples. The T4.2
validator does not run upstream on those held-out inputs.
Governed by `testing.md` §11 and `datasets/manifest.yaml`.

The panel has been aligned and evaluated. The T4.1 differential run is **90 of 90
passing** as of 2026-10-01, with failure-closed comparators; the earlier invocation
was 86 of 90 with four skipped for want of a fixture. Both retain fixed
synthetic/regression fixtures alongside real-panel substitutions, so the case count
is a count of compatibility checks, not of independently validated real libraries.
The separate T4.2 held-out
endpoint report records 3 raw validator passes, 3 failures, 3 inconclusive
outcomes, and 9 not-evaluated outcomes across endpoint/stratum rows. An independent technical
review found that the raw tally is not completed scientific validation:
E1's cross-lab strand expectation is unsupported, the required feature
stratification is absent, and E3/E5 need same-input port/upstream comparisons.
The subsequent independent coordinate audit found the rat refGene converter
still incorrect; its annotation-based results are invalid pending corrected
reference preparation and alignment. See the
[readiness audit](../docs/READINESS_AUDIT_2026-10-01.md).
See [the results and next steps](ENDPOINT_RESULTS.md) and
[the frozen protocol](ENDPOINTS.md).

All three held-out runs have been inspected. The rat E5 result was scored
before and after an annotation-converter correction. These samples are
consumed for confirmation: any follow-up change motivated by these outcomes
needs fresh independent data. The panel supports limited, workload-specific evidence; it
does not establish population-level, cross-aligner, whole-genome, or
single-cell claims.

## Why a third-party aligner

The obvious shortcut is to generate BAMs with this project's own workload
generator, the way the benchmark suite does. That is exactly the wrong move,
and this repo has already paid for learning why: the sliding-window `tin`
regression passed every benchmark workload because the generated reads are
cleanly block-separated by chromosome, and only the independently-authored
fixtures in `verification/` caught it. A BAM from our own generator would
encode our own assumptions about transcripts, strand and junctions, so
comparing the port against upstream on it tests the port against us.

## The one-way door

`manifest.yaml` splits runs into development and held-out at **study** level.
Development runs are for debugging and adding verification cases. Held-out
runs are inspected once, at final validation. If a held-out run is ever
opened to diagnose a defect, that stratum stops being held out and a fresh
run must be selected to replace it — `testing.md` §11.1 is explicit about
this, and `manifest.yaml` records `exposure:` per stratum so the history is
auditable rather than remembered.

## Pipeline

```bash
# 1. Reference (pinned by URL + SHA256; idempotent, resumable)
./datasets/fetch_reference.sh

# 2. STAR index (pinned params; ~30-60 min, ~16 GB RAM)
./datasets/build_star_index.sh

# 3. Reads (per-run MD5 from the ENA archive, not from this repo)
oracle/venv/bin/python3 datasets/fetch_fastq.py SRR1216016 SRR1216063

# 4. Align
./datasets/align_run.sh SRR1216016

# 5. Re-run the real-data differential matrix against the aligned development panel
RSEQC_REAL_DATA=datasets/aligned/real oracle/venv/bin/python3 verification/run_diff.py
```

Steps 1-3 are independent of step 4 and can run concurrently: downloads are
network-bound and the index build is memory-bound.

### Indices are per-assembly, and an alignment must match one

`build_star_index.sh` writes each index to its own directory under `star_index/`,
named after the genome, the annotation, the `sjdbOverhang` and the contig set. The
human and rat panels differ in all four, and a STAR index is only valid for the
inputs it was built from — the parts are overwritten in place, so sharing one
directory means "rebuild the rat index" silently replaces the human one with a
mixture of the two. That happened here once, and the result would have been a
human alignment against a corrupt index producing plausible-looking output.

`align_run.sh` derives the same directory name, and additionally refuses to proceed
unless the index directory's `.built-*` stamp records a `gtf_sha256` matching the
`annotation.gtf` actually sitting beside it. It writes that whole stamp into the
run's `align.json`, which is what makes "was this BAM aligned against that index?"
answerable afterwards. To align the rat stratum:

```bash
GENOME=datasets/heldout/reference/rn6.fa.gz \
GTF=datasets/heldout/reference/rn6.gtf \
CONTIGS="chr1 chr2 chr10" SJDB_OVERHANG=100 \
IDX_ROOT=datasets/heldout/star_index REF_DIR=datasets/heldout/reference \
RAW_DIR=datasets/heldout/raw OUT_ROOT=datasets/heldout/aligned \
  ./datasets/build_star_index.sh
```

`sjdbOverhang` is the single mate length minus one (101 bp reads → 100). Do not
read it off STAR's own progress output: its "Read length" column reports the summed
mate length, so a 2x101 library displays 202. A wrong overhang degrades
splice-junction detection silently, which is exactly what endpoint E5 measures.

## What is committed

The scripts, `manifest.yaml`, endpoint protocol and summary, and reference-curve
provenance are committed. The FASTQ/BAM data, indexes, per-run endpoint JSON,
and working logs are gitignored. They are re-derivable from the pinned digests
and recorded run accessions; a reviewer can reproduce the panel with the
commands above.
