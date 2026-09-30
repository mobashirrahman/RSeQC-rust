# Publication gates plan — full vs restricted-first release

Prepared: 2026-09-17. Baseline: working tree at `6d088e5` (plus uncommitted
`read_hexamer` fix under review).

> **Status note (2026-09-30).** Sections 1 and 5 are the original 2026-09-17
> snapshot and are **stale in three specific ways**, recorded here so a reader does
> not act on them:
>
> 1. **`read_hexamer` is resolved.** The 2026-09-17 item 1 ("commit the fix or drop
>    it from R1") is done: fixed in `3b95212`, differential case added in `38760bb`.
> 2. **The 12 zero-coverage commands now have cases.** The matrix is **84 cases
>    spanning all 33 commands** (was 24 cases / 21 commands). All 84 pass.
> 3. **T5 is substantially done.** `benchmarks/protocol.md` is frozen and
>    preregistered; all 29 benchmarkable commands are measured at 10 paired
>    repetitions with a structural equivalence gate; gate G6a is met. What remains
>    is G6b — the committed results still trace to a commit with `git_dirty: true`,
>    so a final re-run against a clean tree is required before any performance
>    claim may cite them.
>
> Still open and unchanged: T4 (no `datasets/manifest.yaml`, no real-data panel, no
> held-out split), the LICENSE decision and DIV-0003, `crates/python` (still a
> 2-line stub), and the R1 matrix's missing **negative / option / overwrite**
> branches — 84 cases but only ~2 exercise a failure path, which is the gap that let
> a real `tin` regression ship (see §5.1).

This plan finishes the scoping decision required by the goal: either finish
all work scoped to **full 33-command + standalone + API publication gates**,
or ship a **restricted named-command/profile release first**. Evidence below
is from executed commands on this machine, not from plan documents alone.

Related inputs (checked, not re-argued here): `docs/PORTING_PLAN.md` (Steps
1–12, gates G0–G8), `testing.md` (gates T0–T6, AUD-01–16, Appendix A),
`docs/decisions/0001-target.md`, `compatibility/{commands,api,divergences,
upstream}.yaml|lock`, `verification/run_diff.py` + fixtures.

## 1. Executed evidence (2026-09-17)

- `cargo build --workspace --release --locked`: succeeds.
- `cargo test --workspace --locked`: succeeds (`rseqc-commands` 191 passed,
  `rseqc-formats` 37 passed, plus small CLI unit tests; `cargo fmt --check`
  and strict Clippy status were not re-run in this pass — re-run before any
  release candidate).
- `oracle/venv/bin/python3 verification/run_diff.py`: **24/24 cases PASSED**
  (tail output verified this session).
- `CASES` in `verification/run_diff.py` covers **21 distinct commands**:
  `bam_stat`, `read_NVC`, `read_GC`, `read_duplication`, `read_quality`,
  `FPKM_count`, `RNA_fragment_size`, `geneBody_coverage`, `tin`,
  `clipping_profile`, `insertion_profile`, `bam2fq`, `split_bam`,
  `infer_experiment`, `split_paired_bam`, `deletion_profile`,
  `mismatch_profile`, `junction_annotation`, `inner_distance`,
  `read_distribution`, `divide_bam`.
- **12/33 commands have zero differential cases**: `read_hexamer`
  (fix in working tree, uncommitted: `crates/cli/src/bin/read_hexamer.rs`,
  `crates/commands/src/read_hexamer.rs` + 2 untracked fixtures),
  `bam2wig`, `geneBody_coverage2`, `normalize_bigwig`, `overlay_bigwig`,
  `FPKM-UQ`, `RPKM_saturation`, `junction_saturation`, `sc_bamStat`,
  `sc_editMatrix`, `sc_seqLogo`, `sc_seqQual`.
- Standalone blockers (source-confirmed): `Rscript` subprocess paths remain
  in `genebody_coverage`, `genebody_coverage2`, `rpkm_saturation`,
  `sc_editmatrix`, `sc_seqqual`; `wigToBigWig` in `bam2wig.rs`;
  `htseq-count` in `fpkm_uq.rs`; `crates/render/src/lib.rs` (2 lines) and
  `crates/python/src/lib.rs` (2 lines) are doc-only stubs — see
  `testing.md` §2 and DIV-0005/DIV-0016.
- Oracle/contract gaps: `compatibility/upstream.lock` pins source commit
  `59a24c5` but explicitly leaves the executable environment unpinned (no
  pysam/htslib/R/wigToBigWig versions, container digest, locale/timezone)
  and no PyPI 5.0.5 cross-check. Open divergences include DIV-0002/0004
  (SAM-text input), DIV-0005 (native plotting), DIV-0006 (index output),
  DIV-0007 (bam2fq compress), DIV-0003 (license metadata).
- No real-data panel, no held-out validation, no benchmark protocol/runner,
  no clean-install/platform matrix exists yet (`testing.md` T4–T6 open).

## 2. Decision: restricted-first, full-track in parallel

**Ship a restricted named-command/profile release first (R1 below). Do not
claim the full 33-command + standalone + Python-API port until its gates
close.** Rationale:

- R1 can close on measured work (21 commands already have green
  differential cases + AUD-04/05/06/07 regressions fixed in-tree).
- Full scope is blocked by 12 missing command matrices, empty
  render/python crates, external-helper elimination, oracle environment
  lock, real-data/benchmark/packaging gates — weeks of implementation and
  independent review, not a paperwork exercise.
- `testing.md` §1 already permits this: “A publication using one validated
  command need not wait for every unrelated command” provided the release
  evidence explicitly identifies the restricted scope.

## 3. Restricted release R1 (proposed)

- **Name:** `rseqc-rust R1 — BAM-only data-table profile`.
- **Named commands (21):** the 21 listed in §1 with passing differential
  cases. `read_hexamer` joins R1 only after its dirty fix is committed plus
  a new `run_diff` case passes against the pinned oracle.
- **Enforced profile (precondition, not “rare in practice”):**
  BAM input only (SAM-text rejected, cf. DIV-0002/0004); R-script outputs
  are validated as **text artifacts only** — R execution / PDF rendering is
  out of scope and must be refused or skipped with a documented flag
  (no silent `Rscript` requirement); `divide_bam` seed gives
  same-implementation repeatability only (DIV-0017, cross-language seed
  equivalence explicitly excluded); BAM writers compared by decoded
  record-multiset + stdout report, not compressed bytes.
- **Explicitly excluded from R1 claims:** the 12 zero-coverage commands,
  native plots/PDFs, `wigToBigWig` conversion, `htseq-count`/`FPKM-UQ`
  subprocess behavior, `.bai` index writing, `bam2fq -c` compression,
  `import qcmodule` / `python *.py` invocation, and any performance claim.
- **Acceptance checklist (per command, plus profile):**
  1. Commit clean tree (commit the `read_hexamer` fix or drop it from R1).
  2. `cargo build --workspace --release --locked`,
     `cargo test --workspace --locked`, `cargo fmt --check`,
     `cargo clippy --workspace --all-targets --locked -- -D warnings` all green.
  3. `verification/run_diff.py` all R1 cases pass from a clean checkout with
     the oracle profile hash recorded (bind binary + fixture + oracle-env
     hashes; `RSEQC_KEEP_FAILURES=1` retention on failure).
  4. Close or scope every R1-relevant DIV entry; no “accepted” label covers
     unfinished R1 functionality (`testing.md` §14).
  5. Add negative/option-branch cases for each R1 command (empty input,
     invalid input, existing-output/overwrite, missing index where
     applicable) — currently the matrix is mostly happy-path.
  6. One small real-data spot check per R1 command family (bulk RNA-seq BAM
     + BED12) with provenance record; held-out split declared before tuning.
  7. Clean-install smoke test of the R1 artifact set outside the source
     tree (PATH aliases, spaces in paths, `--help`/exit codes).
  8. Independent reviewer sign-off on expectations, goldens, and the R1
     limitation statement; methods/commands/hashes archived per
     `testing.md` §14.1.
- **Work order for R1:** (a) commit/decide `read_hexamer`; (b) negative +
  option-branch cases; (c) oracle-env pin for the R1 profile only
  (Python/pysam/htslib versions, container digest); (d) real-data spot
  checks; (e) packaging smoke + report; (f) review + versioned release.

## 4. Full publication gates (retained scope, not closed by R1)

Full 33-command + standalone + API release requires all of `testing.md`
T0–T6, summarized as actionable gates:

- **T0 spec:** finish executable oracle env lock + PyPI cross-check;
  rebuild `commands.yaml`/`api.yaml` from pinned source + observed runs;
  resolve DIV-0003 license metadata; machine-readable capability map.
- **T1 runner:** provenance binding (source/bin/fixture/oracle/comparator
  hashes), per-case resource limits + process-group timeout, isolated
  workdirs + read-only inputs, output-tree allowlist, strict result states
  (`pass/fail/expected-upstream-failure/unsupported/infra/not-run`).
- **T2 foundation:** keep Appendix A + AUD-04–07 regressions; add
  flag/CIGAR/annotation/pileup truth tables (§6 fixture families).
- **T3 matrix:** differential + semantic BAM/BigWig/FASTQ comparators,
  full option/artifact/error branches, mutation + property tests for all 33
  rows (§7) — the 12 missing commands land here.
- **T4 science:** `datasets/manifest.yaml`, dev/held-out split, biological
  endpoints + downstream workflow pins (§11).
- **T5 perf:** `benchmarks/protocol.md` + paired, correctness-qualified
  measurements with uncertainty (§12); no speedup claim before T3–T4 pass.
- **T6 distribution:** standalone audit (no Python/R/`htseq-count`/
  `wigToBigWig` on PATH; child-process + `ldd`/bundled-component audit),
  native rendering in `crates/render`, Python bindings + `qcmodule` suite in
  `crates/python`, OS/arch matrix, CI layers, reproducibility bundle (§13).
- **Rule:** do not relabel missing work “accepted divergence” to close
  full gates (`testing.md` §14). R1-excluded behavior stays open on the
  full track.

## 5. Immediate next actions

### 5.0 As of 2026-09-30

Original list, annotated:

1. ~~Decide `read_hexamer`~~ **DONE** (`3b95212`, `38760bb`).
2. Extend the matrix with error/option/overwrite branches (§3 item 5).
   **STILL OPEN and now the highest-value item** — see §5.1.
3. Pin the R1 oracle profile. **PARTLY DONE**: `compatibility/upstream.lock`
   pins the source commit and a source-tree hash; the executable environment
   (pysam/htslib/R versions, container digest) is still unpinned.
4. ~~Start T3 for the 12 missing commands~~ **DONE**: all 33 commands have
   differential cases. `crates/render` remains a partial implementation
   (seqlogo SVG/PNG only).
5. ~~Freeze the benchmark protocol only after T3–T4 outputs qualify~~ **DONE**:
   the protocol was frozen and used; T4 is still open, so the *performance
   claims* remain internal engineering evidence rather than publishable.

Newly added, in priority order:

6. **Re-run the benchmark against a clean tree** to satisfy G6b. The committed
   `results-main.json` records `git_dirty: true`, which gate G6b forbids.
7. **Resolve DIV-0003 / ship a LICENSE file.** Upstream's own metadata is
   self-contradictory (GPL-3.0-or-later in README, GPLv3 in LICENSE, a GPLv2
   packaging classifier), and this project cannot choose its own release
   license until that is decided. This blocks packaging.
8. **T4**: create `datasets/manifest.yaml`, declare the dev/held-out split
   before any tuning, and add biological endpoint checks.

### 5.1 Why the branch-case gap is now the top item

On 2026-09-30 a sliding-window rewrite of `tin`'s and `geneBody_coverage`'s
read index was committed (`78776e8`, then fixed in `6cd93e6`). It broke `tin`
on **7 of 84** differential cases: the pull loop treated any record not
belonging to the chromosome being scored as end-of-block, stranding the
iterator so every later transcript scored 0.0.

Two things about how that happened are worth recording as process, not history:

- **The benchmark could not have caught it.** The generated workloads are
  cleanly block-separated by chromosome, so they never produce the record
  shape that triggers the bug. A benchmark whose inputs the project generates
  tests the project against its own assumptions.
- **The differential harness caught it only by luck.** It passed because
  `verification/fixtures/synthetic/pe.bam` happens to interleave an
  UNMAPPED-but-positioned record into each chromosome's block. No case was
  *designed* to exercise that transition; the fixture set just contains it.
  The same class of bug in a command whose fixtures lack that shape would have
  shipped.

Hence item 2. A negative case per command family (missing input, missing
`.bai`, pre-existing output, malformed BED, empty input) converts "the fixture
happened to catch it" into "a test was written to catch it".

## 6. What “done” means

- **R1 done:** §3 checklist green on a clean revision, reviewer-signed,
  with a versioned release whose notes name exactly the 21 (or 22)
  commands + BAM-only/R-text profile + exclusions. No full-port language.
- **Full done:** T0–T6 evidence + G7/G8 packaging/publication gates, zero
  unexplained mismatches, per-command benchmarks with raw data, clean
  installs on every advertised platform, standalone audit without helper
  executables, Python API suite green, reproducible report regeneration.
