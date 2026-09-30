# T4 real-data panel

Real sequenced reads, aligned by a pinned third-party aligner, used to check
that this port agrees with upstream RSeQC on data nobody in this project
generated. Governed by `testing.md` §11 and `datasets/manifest.yaml`.

**Nothing here supports a scientific claim yet.** No BAM has been produced
and no endpoint has been scored. See `manifest.yaml` for what is still
absent.

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

# 5. Full command matrix, port vs upstream, on the resulting BAM
```

Steps 1-3 are independent of step 4 and can run concurrently: downloads are
network-bound and the index build is memory-bound.

## What is committed

The scripts, `manifest.yaml`, and the logs. Not the data: `datasets/reference/`,
`raw/`, `aligned/` and `star_index/` are gitignored, because they are ~4 GB
and are fully re-derivable from the pinned digests plus the recorded run
accessions. A reviewer reproduces them with the commands above.
