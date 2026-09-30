# Changelog

All notable changes to this project are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/); this project has not yet made a
versioned release, so everything below is grouped under **Unreleased**.

For the full, granular history (including every individual bug fix and its verification
evidence), see `git log` and `compatibility/divergences.yaml` — this file summarizes at the level
of a feature or a closed divergence, not every commit.

## Unreleased

### Added

- A measured performance benchmark suite: a preregistered protocol
  (`benchmarks/protocol.md`), a workload generator over **real** hg38 sequence and real
  RefSeq BED12 annotation (`benchmarks/generate_workload_real.py`), a harness with
  process-tree resource measurement, a structural output-equivalence gate, randomised
  interleaved paired repetitions and paired block-bootstrap intervals
  (`benchmarks/bench.py`), cost-driver sweeps (`benchmarks/scaling.py`), and
  table regeneration from raw data (`benchmarks/analyze.py`). Results for 29 commands at
  10 repetitions are in `benchmarks/RESULTS.generated.md`: **all 29 are faster (1.51x to
  71.6x) and all 29 pass the output-equivalence gate**, and the port uses less memory
  than upstream on all 29. Those results trace to a clean commit (`3a5b88c`,
  `git_dirty: false`), satisfying PORTING_PLAN.md gate G6b. The run is not publication-grade (shared hardware), and no
  speedup figure from it may be published without repeating it on isolated hardware.
  The first run found three commands that were reproducibly *slower* (`infer_experiment`
  0.18x, `bam2fq` 0.35x, `inner_distance` 0.53x) and 14 that used more memory; all four
  issues have since been root-caused, fixed and re-measured (see **Fixed** below).
- Full Rust ports of all 33 upstream RSeQC commands, each with an original-name PATH alias
  (`scripts/install-aliases.sh`) so existing pipelines invoking commands by their upstream name
  work unmodified.
- A differential verification harness (`verification/run_diff.py`) that runs the real upstream
  Python CLI and this port's compiled binary against the same fixture and diffs their actual
  output — the primary correctness signal this project relies on, not unit tests alone.
- Shared alignment I/O (`rseqc_formats::open_alignments`) supporting BAM, plain-text SAM, and CRAM
  input uniformly across every command that advertises more than BAM upstream.
- Native `.bai` index writing (`rseqc_formats::write_bai_index`) for the commands that offer an
  `--index`/`--index-output` flag.
- Compressed (`.gz`/`.Z`/`.z`/`.bz`/`.bz2`/`.bzip2`) FASTA/FASTQ input support for the single-cell
  sequence-quality/logo commands.
- Native SVG and PNG sequence-logo renderers (`crates/render::seqlogo`/`seqlogo_png`) for
  `sc_seqLogo.py`, plus a native PDF renderer (`crates/render::pdf`) that wraps the same raster —
  so `--oformat svg`/`png`/`pdf` are all functional (DIV-0016 closed).
- GitHub Actions CI (build/test/clippy across Linux/macOS/Windows, a Docker image build, an
  advisory `cargo fmt` check) and a tag-triggered release workflow that publishes a checksummed
  distribution archive.
- A `Dockerfile` for a minimal container image.
- This project's own README, CONTRIBUTING guide, and this changelog.

### Fixed

- **`tin.py` no longer holds the whole BAM in memory.** It built a whole-file per-read
  index, retaining every read for the whole run, at roughly 600 bytes per read (query
  name, qualities, sequence, CIGAR and two block lists are all per-read heap buffers).
  A 600k-read BAM measured **366 MB** against upstream `pysam`'s 43 MB, because `pysam`
  answers each region query from a BAI index and never materialises the rest of the file.
  But `tin` only ever consumes reads whose *start* falls inside the current transcript's
  span, so the port now streams the BAM once in transcript coordinate order and keeps only
  the reads that can still reach a transcript not yet scored, retiring a read once
  `end <= tx_start`. Resident reads fall from 597,048 to a mean of 3,424 (max 16,438).
  Measured: peak RSS **366 MB -> 20.3 MB** and wall time 19.8 s -> 9.7 s, so `tin` goes
  from the port's worst memory loser to using less than half of upstream's. Exact, not an
  approximation: the per-transcript scoring is a single shared function, transcripts are
  scored in coordinate order but collected by original index so output row order and the
  summary's pairwise-summation order are unchanged, and both output files are
  byte-identical to upstream on all 3,000 transcripts (600k reads) and all 26,590
  (400k reads), with and without `--subtract-background`. A non-coordinate-sorted input is
  detected and falls back to the whole-file path.
  Re-measured after the fix: **16.32x** against upstream, 20.3 MB against upstream's
  42.3 MB.
- **The alignment reader now streams.** `open_alignments` returned
  `Vec<io::Result<Record>>`, decoding the entire input before any work, which the 14
  commands routed through it used as their working set. It now returns a streaming
  iterator: BAM decodes with `read_record` into one reusable buffer (O(1) memory), SAM
  converts in 4096-record batches, and the whole workspace rebuilt with no call-site
  changes because the `compute_*` functions were already generic over the iterator item
  type. Measured on `bam_stat` at 1.6M records: peak RSS 364 MB -> 2.9 MB (upstream:
  39 MB), and it got faster too (0.50 s -> 0.45 s). This was found by the benchmark suite, whose
  scaling sweep showed peak RSS growing linearly at ~226 bytes/record -- extrapolating
  to ~23 GB for a 50M-read-pair BAM, i.e. the port could not open a normal dataset. CRAM
  remains whole-file buffered because `noodles-cram` 0.99 exposes record iteration only
  as a single-use `records(&header)` that returns a spurious error when re-entered on a
  drained reader.
- **`geneBody_coverage.py` no longer holds the whole BAM in memory either**, via the
  same sliding window `tin` uses (`compute_coverage_windowed`). Exact rather than
  approximate, but for a different reason than `tin`'s: `geneBody_coverage` aggregates
  across *all* transcripts at once, so `tin`'s "a later transcript cannot need this read"
  argument does not apply -- but its accumulation is a sum plus an OR, both commutative,
  so visiting transcripts in coordinate order instead of BED order cannot change the
  result. Measured on the 400-transcript workload: peak RSS **75.1 MB -> 14.4 MB**
  (upstream 40.6 MB), and the command got faster too, **36.2x -> 46.6x**, with
  byte-identical output.
- **Fixed a regression the two windowed drivers introduced in the change above**: the pull
  loop treated *any* record not belonging to the chromosome being scored as "this block is
  over" and stopped, which stranded the record iterator and left every later transcript
  scoring 0.0. Triggered by ordinary input: an UNMAPPED record that still carries a
  reference id and a position, and a leftover from an earlier chromosome. Both are now
  consumed rather than treated as a stop signal.
  `verification/run_diff.py` caught this at 7 of 84 cases failing; the benchmark's own
  generated workloads could not, because they are cleanly block-separated by chromosome
  and never produce that record shape. See `benchmarks/CHANGES.md` 2.7.
- Together, the port now uses **less memory than upstream on all 29 benchmarked
  commands** (the lowest ratio is `tin` at 2.08x lighter; the highest absolute figures
  are `bam2wig` at 793 MB and `geneBody_coverage` at 14 MB).
- Three performance regressions found by `benchmarks/RESULTS.generated.md` and since
  fixed, each root-caused with a profile and re-measured through the same harness with
  outputs verified byte-identical against real upstream:
  - `infer_experiment.py` was **0.18x** the speed of upstream: it scanned every gene
    range on the chromosome once per sampled read (~1.8e9 comparisons at the 200k
    sample cap), where upstream uses `bx.intervals.Intersecter`, a bitset interval
    index. Replaced with an equivalent O(log n) index (sorted starts + running maximum
    of ends, per distinct strand), now **3.43x**.
  - `bam2fq.py` was **0.35x**: its output went to an unbuffered `File`, so a FASTQ
    record's 5 writes each became a syscall -- 5,600,008 `write` calls and 4.6 s of
    system time on 800k reads. Now wrapped in a 1 MiB `BufWriter` (5,600,008 -> 178
    syscalls), now **2.32x**.
  - `inner_distance.py` was **0.53x**: it built two `HashSet<String>` of transcript
    names per pair, twice per pair, by scanning all transcripts. Replaced with an
    O(log n) same-transcript query, now **11.18x**.
  The optimisation in each case is covered by a randomised test asserting it returns
  exactly what the code it replaced returned, and the full 84-case differential harness
  re-runs green.
- The 16 `Rscript`-invoking commands now pre-resolve `--rscript` (a `shutil.which` port,
  `rseqc_commands::exec_resolve`) and report upstream's `Rscript executable not found: <name>`
  instead of a raw OS spawn error (DIV-0021, found by clean-room testing).

Every entry in `compatibility/divergences.yaml` with `status: fixed-working-tree` began as a real
behavioral difference found via the differential harness against real upstream, then closed.
Notable categories (see the divergences file for the full, individually-verified account of each):

- SAM-text and CRAM input support for every command upstream advertises it for (was BAM-only).
- Real R-script generation and `Rscript` invocation for every plotting command (was previously
  producing no plot artifact, or an incomplete one, for several commands).
- `.bai` index writing for `split_paired_bam.py`, `split_bam.py`, and `divide_bam.py`.
- `bam2fq.py`'s `-c/--compress` gzip output.
- Several instances of a Rust CLI printing a raw Debug-formatted dump instead of upstream's real
  progress/report text, or omitting progress messages upstream always prints.
- A genuine scientific-correctness bug in `RPKM_saturation.py`'s percentile-resampling population
  (was rebuilt independently per iteration instead of accumulated cumulatively, producing an
  incorrect, non-monotonic saturation curve).
- Several Python duck-typing int-vs-float formatting quirks (`defaultdict(int)` promoted to float
  by a `+= 1.0` increment pattern) that this port previously rendered incorrectly.
- A genuine `htslib`/`noodles-cram` interoperability discrepancy in how an unmapped read's mapping
  quality round-trips through CRAM.

### Known limitations

See the README's own "Limitations" section for the current, up-to-date list — it is **not**
duplicated here to avoid the two going out of sync. As of this entry: no real biological/held-out
dataset validation or fuzzing has been performed; packaging covers Linux x86_64 distribution archives only (CI now
builds/tests on macOS and Windows too, and a Docker image build exists, but neither has been
verified to actually work yet — no Python-wheel publication); performance benchmarking **has**
been done, but only on shared, non-isolated hardware, so it is internal engineering evidence
and not a publishable speedup claim (`benchmarks/RESULTS.generated.md` §7); this project's own
release license is pending resolution of `compatibility/
divergences.yaml`'s DIV-0003 (upstream's own license metadata is internally inconsistent).
