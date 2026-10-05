# Agent execution plan

Prepared 2026-10-04. This turns the limitations review into task cards that a
small LLM agent can execute one at a time. An orchestrator (the main session or
the maintainer) hands out cards, verifies each result, and merges.

Each card is self-contained: an agent needs this file's "Rules for every agent"
section plus its own card, nothing else.

## Tiers

| Tier | Who | Use for |
|---|---|---|
| S | small model (Haiku class, free swarm backends) | docs, config, mechanical edits following an existing pattern |
| M | mid model (Sonnet class or stronger) | any change to Rust under `crates/`, any parity debugging |
| L | orchestrator or maintainer | design decisions, long measurements, anything outward-facing |

Small models have, in this repository, stopped half-way, committed TODOs,
narrowed test inputs to hide failures, and reported checks clean when they were
not. Do not give an S agent a tier-M card.

## Rules for every agent

1. Work only in the worktree and branch you were given (`agent/<card-id>`).
   Never commit to `main`, never push, never tag.
2. Edit only the files listed under **Files** on your card. If the task cannot
   be done without touching another file, stop and report which file and why.
3. Start every shell with `export PATH="$HOME/.cargo/bin:$PATH"`.
4. Never change expected output, fixtures, comparators, tolerances or test
   inputs to make a check pass. Never add `#[ignore]`, `todo!()` or a `TODO`.
5. Output compatibility is the priority: upstream typos and quirks are
   reproduced, not corrected (see `CONTRIBUTING.md`).
6. One repair attempt. If the acceptance checks still fail after one fix
   attempt, stop and report the failing output. Do not keep trying.
7. Do not edit `benchmarks/results-*`, `benchmarks/*.json`,
   `datasets/heldout/`, `oracle/` or `compatibility/upstream.lock`.
8. Commit with specific `git add <file>` paths, never `git add -A`.

### Standard checks

Run all of these unless the card is docs-only:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
cargo build --workspace --release --locked
cargo test --workspace --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
oracle/venv/bin/python3 verification/run_diff.py <case names from the card>
oracle/venv/bin/python3 verification/synthetic_sweep.py --only '<regex from the card>'
```

`run_diff.py` takes case names as positional arguments and uses the binaries in
`target/release` of the current worktree.

### Report format

End with exactly this, and paste real output, not a summary of it:

```
CARD: <id>
STATUS: done | blocked
COMMITS: <git log --oneline main..HEAD>
GIT STATUS: <git status --short>
CHECK OUTPUT: <last 5 lines of each standard check that was run>
CARD-SPECIFIC EVIDENCE: <what the card's Acceptance section asks for>
NOT DONE / DOUBTS: <anything skipped, assumed or unverified>
```

## Rules for the orchestrator

- Create worktrees by hand (the built-in worktree isolation has failed here):
  `git worktree add <scratchpad>/<card-id> -b agent/<card-id> main`, then
  symlink the gitignored inputs into it: `oracle/upstream-src`, `oracle/venv`,
  `oracle/python313`, `datasets/heldout`, `datasets/aligned`.
- Before merging a card, re-run the standard checks yourself, read the diff
  (`git diff main...agent/<card-id> --stat` first: files outside the card's
  list are a rejection), and look for `*.proptest-regressions` files.
- Run cards marked **serial** one at a time. Cards in the same phase that are
  not marked serial have disjoint files and can run in parallel.
- After each phase, run the full `run_diff.py` and `synthetic_sweep.py` once on
  `main`.

## Decisions the maintainer must make first

These block the cards named; no agent may decide them.

| # | Decision | Blocks |
|---|---|---|
| D1 | Retire the Python API (`crates/python`) and native plots from the v0.1 scope, or keep them as goals | A3 |
| D3 | Whether the single-pass driver (phase C) is in the first paper or a later release | phase C ordering |

The only `open` entry in `compatibility/divergences.yaml` is DIV-0024, which was
scoped rather than fixed on 2026-10-04; no card touches it.

## Phase A: documentation and scope (tier S, parallel)

### A1. Make the README a user document

- **Goal:** a reader can tell in one screen what is supported, how to install
  it, and what the known limits are.
- **Files:** `README.md`, `docs/KNOWN_LIMITATIONS.md` (new).
- **Steps:**
  1. Move every bullet of the README's `## Limitations` section, verbatim, into
     `docs/KNOWN_LIMITATIONS.md` under the same headings. Delete nothing.
  2. Replace the README section with at most 12 short bullets, one per limit,
     each ending with a link to its anchor in `docs/KNOWN_LIMITATIONS.md`.
     State the current limit only; no history of how it was found.
  3. Leave every other README section unchanged.
- **Acceptance:** `wc -l README.md` is under 260. Every sentence removed from
  the README appears in the new file (show `git diff --stat` and confirm the new
  file's line count is at least the number of lines removed).

### A2. Add RustQC to related work

- **Goal:** the manuscript outline names the closest existing tool.
- **Files:** `docs/MANUSCRIPT_OUTLINE.md`.
- **Steps:** under "Statement of need and related work", add a paragraph and a
  comparison table using only these facts: RustQC (Seqera, Bioconda package
  `rustqc`, nf-core module `rustqc`) reimplements 15 RNA-seq QC tools in one
  single-pass binary, including eight RSeQC tools: `bam_stat`,
  `infer_experiment`, `read_duplication`, `read_distribution`,
  `junction_annotation`, `junction_saturation`, `inner_distance`, TIN. It
  describes its outputs as format-compatible with upstream and MultiQC. This
  project covers all 33 RSeQC commands under their original names with
  byte-identical output against a pinned upstream. Mark every speed comparison
  cell "to be measured (card E1)".
- **Acceptance:** no number about RustQC's speed appears. Source link
  `https://seqera.io/blog/rustqc/` is cited.

### A3. Draft ADR 0002 (blocked on D1)

- **Goal:** record the v0.1 scope decision next to `docs/decisions/0001-target.md`.
- **Files:** `docs/decisions/0002-v0.1-scope.md` (new).
- **Steps:** same headings as ADR 0001. Status `proposed`. Content is whatever
  the maintainer decided in D1; do not argue a position.
- **Acceptance:** file exists, under 40 lines, ADR 0001 unmodified.

### A4. Mark the 3.7x figure provisional

- **Goal:** remove the contradiction between the README ("no speedup figure
  from the version-1 run may be published") and the outline's "median 3.7x".
- **Files:** `docs/MANUSCRIPT_OUTLINE.md`.
- **Steps:** in the "Scope decisions" section, state which measurement the 3.7x
  comes from (find it with `grep -rn "3.7x" CHANGELOG.md README.md benchmarks/`),
  and add that the paper will cite only the protocol-v2 study (card E2). If you
  cannot find the source measurement, say so in the report and change nothing.
- **Acceptance:** report quotes the source line found.
- **Serial with A2** (same file).

## Phase B: memory envelope (tier M, parallel, one command each)

Baseline numbers are in `docs/ENVELOPE.md`. Each card must report peak RSS
before and after on the rat alignment with:

```bash
/usr/bin/time -v target/release/<binary> <args> 2>&1 | grep 'Maximum resident'
```

using `datasets/heldout/aligned/SRR1177982/SRR1177982.bam`. Do not edit
`docs/ENVELOPE.md`; the orchestrator updates it from the reports.

### B1. `read_duplication`: store hashes, not strings

- **Why:** about 130 bytes per record on real data; a 50M-pair library needs
  over 10 GB.
- **Files:** `crates/commands/src/read_duplication.rs`.
- **Design:** `compute_duplication` keeps `seq_dup` keyed by the full sequence
  string and `pos_dup` keyed by a `"chrom:start:exon_boundary"` string. Only the
  occurrence histograms are ever output, never the keys. Key both maps by a
  128-bit hash of the same bytes instead. Use a deterministic hash with fixed
  keys; use a crate already in `Cargo.lock` if one fits, otherwise two
  `std::collections::hash_map::DefaultHasher` passes with distinct prefixes.
  Do not add a dependency.
- **Before coding:** confirm by reading the file and
  `crates/cli/src/bin/read_duplication.rs` that no output contains a key. If
  one does, stop and report.
- **Acceptance:** standard checks with cases `read_duplication_basic
  read_duplication_sam_text` and sweep regex `read_duplication`. Output on the
  rat BAM byte-identical before and after (`cmp` on both `.xls` files). A new
  unit test showing that two different keys produce different hashes and equal
  keys equal hashes. RSS before and after.

### B2. `bam2wig`: chunked dense coverage storage

- **Why:** about 44 bytes per covered base; whole-genome human is tens of GB.
- **Files:** `crates/commands/src/bam2wig.rs`.
- **Design:** `ChromWig.forward` and `.reverse` are `BTreeMap<i64, f64>`, one
  node per covered position. Replace the storage with fixed-size chunks (for
  example 4096 positions of `f64` per chunk, allocated on first touch, held in a
  `BTreeMap<i64, Box<[f64; N]>>` keyed by chunk index) plus a parallel
  "touched" bitmap if a position with value 0.0 must be distinguishable from an
  untouched one. Keep `f64`. Iteration order and every rendered value must be
  unchanged. Keep the public function signatures callers use, or update the one
  caller in `crates/cli/src/bin/bam2wig.rs` and say so in the report (that file
  is then allowed).
- **Before coding:** read `render_chrom_block` and determine whether a touched
  position holding 0.0 is rendered differently from an untouched one. State the
  answer in the report.
- **Acceptance:** standard checks with cases `bam2wig_basic
  bam2wig_synthetic_stdout bam2wig_with_wigsum
  bam2wig_missing_output_dir_is_upstream_inconsistent` and sweep regex
  `bam2wig`. `.wig` output on the rat BAM byte-identical before and after. RSS
  before and after; target is under half of the 1.69 GB baseline.

### B3. CRAM: stream instead of buffering the file

- **Why:** `AlignmentRecords::CramBuffered` in `crates/formats/src/lib.rs`
  (around line 398) decodes the whole file before the first record is used.
- **Files:** `crates/formats/src/lib.rs`.
- **Design:** yield records container by container (or slice by slice) through
  the existing `AlignmentRecords` iterator, keeping the existing
  missing-reference error message exactly. No command code changes.
- **Acceptance:** standard checks with every `run_diff.py` case whose name
  contains `cram` (list them with `grep -o 'name="[^"]*cram[^"]*"'
  verification/run_diff.py`), plus the `crates/formats` fuzz tests. Because the
  committed CRAM fixture is 10 records, also build a larger CRAM from the rat
  BAM with `oracle/venv/bin/python3` and pysam in the scratch directory, and
  show `bam_stat` output on it is identical to `bam_stat` on the BAM, with RSS
  for both.
- **If streaming is not achievable in one attempt:** stop, report why, change
  nothing.

### B4. Actionable out-of-memory message

- **Why:** past the memory limit every command dies with an allocator abort and
  no explanation (`docs/ENVELOPE.md`, "What happens past the limit").
- **Files:** `crates/cli/src/lib.rs`.
- **Design:** a `#[global_allocator]` that wraps `std::alloc::System`; when the
  inner allocation returns null, write a fixed message to file descriptor 2
  without allocating (for example `rseqc-rust: out of memory; see
  docs/ENVELOPE.md for per-command memory costs`), then return null so the
  normal abort follows. Exit status must stay as it is today. Confirm all
  binaries under `crates/cli/src/bin/` link this crate; if any does not, stop
  and report.
- **Acceptance:** standard checks (full `run_diff.py`, since every binary
  changes).
  `oracle/venv/bin/python3 verification/check_memory_failure.py --json <scratch>/oom.json`
  still shows zero output files written and the same exit status, and the new
  message appears on stderr. Wall time of `bam_stat` on the rat BAM, three runs
  before and after, to show no slowdown.

## Phase C: single-pass driver (blocked on D3)

A pipeline runs many of these commands on the same BAM and each one re-reads
it. One new binary reads the file once and feeds several commands.

### C1. Design and two-command pilot (tier L)

- **Files:** new `crates/cli/src/bin/rseqc_multi.rs`, new
  `crates/commands/src/multi.rs`, `crates/cli/Cargo.toml`.
- **Design constraint:** the existing `compute_*` functions take a record
  iterator and must not change. One reader thread decodes batches of records and
  sends each batch (shared, not copied per consumer) over bounded channels to
  one worker thread per selected command; each worker runs the unchanged
  `compute_*` over a channel-backed iterator. Bounded channels give
  backpressure so memory stays flat.
- **Contract:** `rseqc_multi -i x.bam [-r x.bed] -o <prefix> --run bam_stat,read_GC`
  writes exactly the files each command would write on its own, and each
  command's stdout and stderr to `<prefix>.<command>.stdout` / `.stderr`.
- **Pilot:** `bam_stat` and `read_GC`. Each command's CLI `run` logic must be
  callable as a function taking a record iterator; document the pattern used at
  the top of `multi.rs`, because the C2 cards copy it.
- **Acceptance:** a new script `verification/check_multi.py` that runs each
  selected command alone and through `rseqc_multi` on the same input and
  compares every file and stream byte for byte. It must fail when a file is
  missing on either side.

### C2. One command per card (tier M, serial merge)

One card each, in this order (streaming commands first): `read_NVC`,
`read_quality`, `clipping_profile`, `insertion_profile`, `deletion_profile`,
`mismatch_profile`, `infer_experiment`, `read_distribution`, `inner_distance`,
`junction_annotation`, `junction_saturation`, `read_duplication`.

- **Files:** `crates/cli/src/bin/<command>.rs` and
  `crates/commands/src/<command>.rs` for that command only. The agent does not
  edit `multi.rs`; it reports the one registration line and the orchestrator
  adds it.
- **Steps:** follow the pattern documented at the top of `multi.rs`.
- **Acceptance:** standard checks for that command's `run_diff.py` cases, and
  `verification/check_multi.py` passing for that command once the orchestrator
  has registered it.

Windowed commands (`tin`, `geneBody_coverage`, `FPKM_count`,
`RNA_fragment_size`) are out of scope until C2 is complete and measured.

## Phase D: upstream drift (tier S)

### D1. Scheduled drift check

- **Files:** `.github/workflows/upstream-drift.yml` (new).
- **Steps:** weekly `schedule` plus `workflow_dispatch`. Read the pinned commit
  from `compatibility/upstream.lock`, run
  `git ls-remote https://github.com/liguowang/RSeQC.git HEAD`, and fail with a
  message naming both commits if they differ. No other action.
- **Acceptance:** the YAML parses (`oracle/venv/bin/python3 -c "import yaml,sys;
  yaml.safe_load(open(sys.argv[1]))" <file>`; if `yaml` is missing, say so). The
  report must state that the workflow has not been run on GitHub.

## Phase E: evidence (tier L, not for delegation)

### E1. RustQC comparison

Approved 2026-10-04. The full card set is in
[`benchmarks/RUSTQC_COMPARISON_PLAN.md`](../benchmarks/RUSTQC_COMPARISON_PLAN.md)
(cards R1 to R8); most of it is delegable.

### E2. Benchmark study v2

Run `benchmarks/protocol-v2.md` once, at a tagged release candidate, after
phases B and (if D3 says so) C are merged. A whole-genome rat run at 32M already
exists in `benchmarks/results-rat-wholegenome-32M`; the missing workload is a
human whole-genome alignment at 10M and 50M pairs (`docs/ENVELOPE.md`, open
item 1). Needs a quiet machine; an agent cannot provide that.

## Phase F: release (maintainer only)

### F1. Beta tag

After A1, B1 to B4 and D1 are merged: tag
`v0.1.0-beta.1`, Linux x86_64 only, and let `release.yml` build it. Pushing and
tagging are the maintainer's actions.

## Order

1. Decisions D1 and D3.
2. Phase A and phase D in parallel (docs and config only).
3. Phase B in parallel, one worktree per card.
4. F1 beta tag.
5. Phase C if D3 puts it in scope, then E1 and E2.
