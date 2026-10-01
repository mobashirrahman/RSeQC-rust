# T4 real-data panel

Real sequenced reads, aligned by a pinned third-party aligner, used to compare
this port with upstream RSeQC on data nobody in this project generated.
Governed by `testing.md` §11 and `datasets/manifest.yaml`.

The panel has been aligned and evaluated. The T4.1 differential run executed
86 of 90 cases and all 86 passed against upstream; four cases were skipped
because the panel lacks their required inputs. The separate T4.2 held-out
endpoint report records 3 passes, 3 failures, 3 inconclusive outcomes, and 9
not-evaluated outcomes across endpoint/stratum rows. These endpoints assess
scientific claims about upstream outputs, not port-versus-upstream agreement.
See [the results](ENDPOINT_RESULTS.md) and [the frozen protocol](ENDPOINTS.md).

All three held-out runs have been inspected once. They are consumed for
confirmation: any follow-up change motivated by these outcomes needs fresh
independent data. The panel supports limited, workload-specific evidence; it
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

## What is committed

The scripts, `manifest.yaml`, endpoint protocol and summary, and reference-curve
provenance are committed. The FASTQ/BAM data, indexes, per-run endpoint JSON,
and working logs are gitignored. They are re-derivable from the pinned digests
and recorded run accessions; a reviewer can reproduce the panel with the
commands above.
