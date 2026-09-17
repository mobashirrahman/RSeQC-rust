# Changelog

All notable changes to this project are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/); this project has not yet made a
versioned release, so everything below is grouped under **Unreleased**.

For the full, granular history (including every individual bug fix and its verification
evidence), see `git log` and `compatibility/divergences.yaml` — this file summarizes at the level
of a feature or a closed divergence, not every commit.

## Unreleased

### Added

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
  `sc_seqLogo.py`.
- GitHub Actions CI (build/test/clippy across Linux/macOS/Windows, a Docker image build, an
  advisory `cargo fmt` check) and a tag-triggered release workflow that publishes a checksummed
  distribution archive.
- A `Dockerfile` for a minimal container image.
- This project's own README, CONTRIBUTING guide, and this changelog.

### Fixed

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
duplicated here to avoid the two going out of sync. As of this entry: `sc_seqLogo.py`'s `pdf`
output format is not implemented (`svg`/`png` are); no real biological/held-out dataset validation
or fuzzing has been performed; packaging covers Linux x86_64 distribution archives only (CI now
builds/tests on macOS and Windows too, and a Docker image build exists, but neither has been
verified to actually work yet — no Python-wheel publication); no performance benchmarking has
been done; this project's own release license is pending resolution of `compatibility/
divergences.yaml`'s DIV-0003 (upstream's own license metadata is internally inconsistent).
