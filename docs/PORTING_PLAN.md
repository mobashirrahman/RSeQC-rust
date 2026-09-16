# RSeQC full-port execution plan

Status: planning only; no implementation or GitHub publication has been performed.
Prepared: 2026-09-16.

## Objective and decisions

Produce an independently buildable replacement for RSeQC, with documented CLI and Python API compatibility, complete command coverage, reproducible verification, per-command benchmarks, and a public GitHub release suitable for scientific citation. A full port means implementing the behavior, not invoking the original Python implementation as a hidden fallback.

The selected target language is **Rust**, prioritizing performance and standalone binaries. The primary distribution must run all 33 commands without Python, R, or separately installed helper executables. Native format libraries may be linked or bundled where their licensing and target-platform support permit; inspect the final runtime dependencies before describing an artifact as standalone. A standalone executable is not necessarily fully statically linked.

Retain a separate optional Python compatibility distribution with thin bindings and launchers: a native binary alone cannot replace `import qcmodule` or `python command.py`. These adapters must call the Rust implementation, not reproduce the algorithms in Python. Python/R and upstream helper programs remain permitted in the isolated oracle and verification environments. Native plotting and replacement of external conversions are required implementation work. If a renderer cannot meet the agreed compatibility contract, record that as a release blocker or explicit limitation rather than silently requiring R.

Deliver in two independently tested forms: standalone Rust CLI artifacts first, then Python-package compatibility artifacts before claiming the full package port is complete. Original command names, including `.py` names executed through PATH, must remain available in the standalone distribution. Invocation through `python command.py` belongs to the optional Python distribution.

Default scope: all 33 installed commands, their option branches and artifacts, plus an explicitly inventoried Python API. Do not silently reduce “full” to a few popular commands. Start with Linux; audit upstream platform behavior and select the supported OS/architecture matrix before release claims. Do not assume identical behavior across upstream GitHub and PyPI releases.

## Inspected baseline

Upstream: https://github.com/liguowang/RSeQC

Inspected commit: `59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24`.
Package metadata version: `5.0.5`. This is a candidate reference, not proof that the published 5.0.5 package has identical contents.

Evidence from this checkout:

- `setup.py` installs 33 Python scripts.
- `src/qcmodule` contains 20 Python modules, including its empty initializer. Together, command scripts and these modules contain approximately 21,300 lines.
- `SAM.py` and `BED.py` contain substantial shared implementation; splitting work solely by source filename would create bottlenecks.
- Declared dependencies include pysam, bx-python, NumPy, pyBigWig, and logomaker. Imports and subprocesses reveal further dependencies, including pandas, matplotlib, Rscript, and wigToBigWig; execution must establish the complete dependency graph.
- Several routines use random shuffling/sampling. `divide_bam.py` uses a seeded RNG; other routines need individual inspection.
- The checked-in CI performs installation/import/help checks and distribution builds. No test-named files were found in this checkout; there is no demonstrated functional compatibility suite to inherit.
- README/license metadata are inconsistent: README describes GPL-3.0-or-later, LICENSE contains GPLv3, and a project classifier says GPLv2. Record and resolve applicable notices before choosing release metadata; do not copy the inconsistent classifier blindly.

Source anchors: [installed commands](https://github.com/liguowang/RSeQC/blob/59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24/setup.py), [metadata](https://github.com/liguowang/RSeQC/blob/59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24/pyproject.toml), [shared modules](https://github.com/liguowang/RSeQC/tree/59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24/src/qcmodule), [CI](https://github.com/liguowang/RSeQC/blob/59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24/.github/workflows/ci.yml).

## Step 1 — Freeze the compatibility contract

Owner: coordinator, assisted by one bounded inventory agent.

1. Record the selected Rust target and performance/standalone priorities. Freeze the allowed build, linked, bundled, and runtime dependencies; use `docs/decisions/0001-target.md` as the architecture decision.
2. Pin the oracle source commit, source archive hash, Python version, exact dependency versions, R packages, external executables, locale, timezone, and container digest.
3. Independently inspect the matching PyPI release. Choose GitHub or PyPI as the primary baseline; any additional versions become separately tested profiles.
4. Inventory every installed command, option/default, input format, output filename/schema, stdout/stderr behavior, exit code, sorting requirement, index requirement, and external process.
5. Inventory importable modules, public functions/classes, signatures, return types, exceptions, mutation behavior, generators, and observable object state. Include unreferenced public helpers; command coverage does not establish API coverage. Audit `psyco_full.py`, `update.pl`, and packaged non-code files explicitly rather than silently treating all repository files as runtime interfaces.
6. Define exact compatibility levels: invocation compatibility; deterministic textual compatibility; numeric equivalence; binary-file semantic equivalence; visual equivalence; Python API equivalence.
7. Establish a known-divergence ledger. Upstream failures and questionable scientific behavior require explicit decisions. Do not silently “correct” an algorithm in compatibility mode.

Deliverables: `compatibility/upstream.lock`, `compatibility/commands.yaml`, `compatibility/api.yaml`, `compatibility/divergences.yaml`, `docs/compatibility.md`, `docs/decisions/0001-target.md`.

Gate G0: every installed command and candidate API surface has an owner and a disposition; compatibility claims and baseline are unambiguous.

## Step 2 — Build the executable reference environment

Owner: oracle agent. Can run alongside Step 3 once G0 is available.

1. Build the locked upstream environment and run all 33 help entry points and one valid workload for each command.
2. Record undeclared dependencies and runtime failures. Preserve the unmodified reference; any necessary patch creates a separately identified oracle variant with an explicit diff.
3. Execute each case in a fresh working directory with bounded CPU, memory, disk, and time. This also isolates upstream shell commands that use wildcard filenames and cleanup.
4. Capture argv, environment, input hashes, stdout, stderr, exit status/signal, output tree, executable versions, and elapsed resources.
5. Repeat selected runs to identify timestamps, ordering differences, RNG behavior, and other nondeterminism before generating goldens.

Deliverables: `oracle/`, environment locks, reference runner, immutable oracle artifacts with provenance.

Gate G1a: each command has a successful reference run or a reproducible, explicitly classified upstream failure. Such failures are not counted as successful functional coverage.

## Step 3 — Build fixtures and independent scientific expectations

Owner: fixture agent; does not implement production algorithms.

1. Generate tiny, human-auditable SAM/BAM, FASTQ, FASTA, BED12, WIG/BigWig, chromosome-size, and single-cell inputs as required by the inventory.
2. Cover single/paired reads; all relevant flags; mapped/unmapped mates; secondary/supplementary alignments; duplicates; QC failures; MAPQ thresholds; missing tags; mixed lengths; missing qualities; and strand orientation.
3. Cover CIGAR operations M/I/D/N/S/H/P/=/X, splice boundaries, overlapping mates, chromosome boundaries, BED half-open coordinates, negative-strand transcripts, multiple isoforms, overlapping annotations, and contig-name mismatches.
4. Cover empty inputs, zero denominators, malformed/truncated files, absent indexes, invalid sorting, unwritable output paths, and paths with spaces. Specify applicability per command.
5. Add hand-calculated expected counts, distributions, normalization values, and transcript metrics. Implement independent checks of formulas; never derive every expected value from the port itself.
6. Select redistributable public bulk and single-cell datasets, plus degraded RNA and spliced-read cases for TIN/junction coverage. Record accession, URL, license/access conditions, checksums, preprocessing, and reference genome/annotation versions.
7. Reserve held-out real datasets for final verification. Synthetic downsampling is useful for scaling, but does not replace biological diversity.

Deliverables: fixture generators and seeds, tiny committed fixtures, `datasets/manifest.yaml`, download/preprocessing scripts, mathematical expectations.

Gate G1b: boundary cases can reveal coordinate/filtering errors independently of upstream output; fixtures can be regenerated and validated.

## Step 4 — Build verification before production porting

Owner: verification agent, independent of the command implementer.

Create one runner that executes reference and candidate against the same case and reports actionable differences. Store expectations outside production code. Comparator tests must deliberately introduce bad outputs and confirm that mismatches are detected.

| Output/behavior | Required comparison |
|---|---|
| CLI | Original executable names, `.py` launch semantics where promised, flags/defaults, help, streams, exit codes, output locations |
| Deterministic text | Byte comparison where meaningful; explicitly named normalization only for paths/timestamps or documented unstable fields |
| Counts/categories | Exact integer equality, including zero/missing categories and ordering when observable |
| Floating point | Metric-specific absolute/relative tolerance justified before testing; explicit NaN/Inf, zero, rounding, percentile and interpolation rules |
| BAM | Decoded headers, records, sequence/qualities/tags, order where promised, valid indexing; compression bytes may differ |
| BigWig | Chromosomes, intervals, values, missing regions versus zero, and summaries used by consumers; assess zoom/aggregation behavior |
| FASTQ | Headers, sequences, qualities, pairing, record order, gzip-decoded contents |
| Plots | Underlying data and generated scripts where part of the interface, plus rendered labels/layout/curves in a locked renderer; ignore only declared metadata |
| Python API | Imports, signatures, types, values, exceptions, mutability/state and iterator behavior |
| Invalid inputs | Exit status/error class, streams and relevant messages, partial artifacts, cleanup |

For randomized algorithms, first determine whether the reference has a usable seed. Identical seeds across languages do not imply identical samples: exact mode may require the same PRNG, sampling algorithm, and traversal order. When reference instrumentation is necessary, label those runs separately and retain unmodified black-box runs. Predefine repeated-seed statistical checks and uncertainty for stochastic equivalence; do not call that byte-exact compatibility.

Add invariants and metamorphic checks where their assumptions hold: partition counts sum correctly; coverage mass is conserved under the specified filters; reference renaming preserves values; duplicating reads doubles additive counts; order-invariant metrics survive permutation. Document exceptions for sampling, duplicate-sensitive and paired-order algorithms.

Deliverables: `verification/`, comparator unit tests, case manifests, golden provenance, machine-readable and human-readable reports.

Gate G2: a deliberately broken candidate fails each comparator class; no broad normalization or tolerance can conceal scientific discrepancies.

## Step 5 — Establish shared interfaces and a vertical slice

Owner: core agent; coordinator owns shared architecture and build metadata.

1. Define shared alignment records/streams, optional tags, reference dictionary, CIGAR blocks, interval queries, annotation models, numeric primitives, seeded sampling, output writers, and structured errors.
2. Document reference/query coordinate conventions and read-versus-fragment counting. Keep command-specific filtering configurable; do not force all tools to share one filter policy.
3. Evaluate mature alignment/BigWig libraries with small decoding/writing probes against reference behavior before adoption. Pin versions and inspect redistribution requirements.
4. Separate computation from CLI parsing and plot rendering so the CLI and Python bindings use the same implementation.
5. Implement `bam_stat.py` end to end, including the relevant `qcmodule` binding if required, packaging, one valid and one invalid case, and differential verification.
6. Freeze interface version 1 only after this slice passes. Create stubs so command agents can work against stable contracts.
7. Use a Cargo workspace with separate core, format I/O, command, rendering, CLI, and Python binding crates. Keep the core independent of Python and CLI parsing. Use explicit integer widths, checked coordinate conversions, and documented floating-point accumulation order. Preserve observable ordering and seeded sampling behavior when adding concurrency.
8. Prototype alignment/BigWig I/O and native PDF/image rendering early: validate the hardest input/output formats and plotting features before distributing all command tasks. Do not promise pure Rust dependencies or universal static linkage before these probes pass.
9. Assign one owner to workspace manifests, dependency locks, toolchain pinning, and shared traits. Command workers own task-specific modules. Maintain a small number of coherent crates rather than requiring one crate per command.

Deliverables: core interfaces, build system, launchers/binding skeleton, first fully verified command.

Gate G3: the first command passes in a clean installed environment, and independent command modules can compile/test without editing shared internals.

## Step 6 — Port every command in bounded parallel work packages

All names below retain their original `.py` suffix. Groups are scheduling queues, not single oversized agent tasks. Assign one command or one coherent helper per task; split complex commands into algorithm and serialization tasks after defining the contract.

| Queue | Commands (33 total) | Prerequisites |
|---|---|---|
| A: alignment and file operations | `bam_stat.py`, `bam2fq.py`, `divide_bam.py`, `split_bam.py`, `split_paired_bam.py` | Alignment I/O, flags, writers; intervals for split; RNG for divide |
| B: read-level profiles | `read_GC.py`, `read_NVC.py`, `read_quality.py`, `read_duplication.py`, `clipping_profile.py`, `deletion_profile.py`, `insertion_profile.py`, `mismatch_profile.py`, `read_hexamer.py` | Alignment/sequence readers, CIGAR/tag semantics, histograms, plotting |
| C: annotation and fragment metrics | `infer_experiment.py`, `read_distribution.py`, `inner_distance.py`, `RNA_fragment_size.py`, `junction_annotation.py` | BED12, interval queries, strand/pair logic, junction extraction |
| D: coverage tracks | `bam2wig.py`, `geneBody_coverage.py`, `geneBody_coverage2.py`, `normalize_bigwig.py`, `overlay_bigwig.py` | Alignment/BigWig I/O, coverage/percentile rules, normalization, plotting |
| E: expression and sampling | `FPKM_count.py`, `FPKM-UQ.py`, `RPKM_saturation.py`, `junction_saturation.py`, `tin.py` | Annotation, coverage, quantiles, sampling, formula verification |
| F: single-cell tools | `sc_bamStat.py`, `sc_editMatrix.py`, `sc_seqLogo.py`, `sc_seqQual.py` | BAM tags, FASTQ, cell/barcode handling as actually used, matrices, logo/heatmap rendering |

Run a separate API queue for public `qcmodule` behavior not exercised by these commands: BED/SAM helper methods, FASTA/FASTQ utilities, CIGAR variants, annotation, statistics/quantiles, k-mers, ORFs, wiggle handling, and other inventoried helpers. Similar-looking legacy helper implementations may behave differently; merge them only after proving equivalence.

Each task must deliver source, relevant CLI/API adapters, targeted differential cases, independent edge checks, and a short discrepancy report. Plotting and optional argument branches are part of completion. No “TODO: plotting” command is counted as finished.

Gate G4: 33/33 commands implemented; every inventoried API item implemented or explicitly classified; all required option and artifact families covered. Unimplemented public API is a disclosed limitation and blocks an unqualified full-package claim.

## Step 7 — Integrate and resolve discrepancies

Owner: coordinator/integrator; separate reviewer for each completed task.

1. Run the full differential matrix on the integrated commit, not only individual branches.
2. Classify every mismatch: candidate bug, upstream bug, reference-environment failure, permitted serialization difference, stochastic variation, or explicit compatibility exception.
3. Fix candidate bugs. For upstream scientific errors, preserve documented compatibility behavior where appropriate and provide corrected behavior only as an explicit separate mode/version. Record unavoidable incompatibilities.
4. Run held-out datasets and workflow tests that use original command names and consume generated tables/BAM/BigWig/plots. Select real consumers during discovery and pin their versions.
5. Test Python imports/calls independently of CLI tests, installation into a clean environment, missing optional dependencies, and supported platform builds.
6. Verify that production execution never imports or shells out to the upstream implementation. Python bindings may delegate to the new core; the upstream package is confined to oracle environments.

Gate G5: zero unexplained scientific mismatches; full report published; every exception visible in the compatibility statement. Matching a broken reference alone is not evidence of scientific validity.

## Step 8 — Establish fair baseline benchmarks

Owner: benchmark agent. Harness work can start earlier; performance conclusions require G5.

1. Write the benchmark protocol and resource budget before measuring final speedups.
2. Include every command, its principal execution modes, and small/medium/large workloads appropriate to its inputs. Select sizes by pilot measurements to avoid an unbounded cross-product. Explicitly label unsupported combinations and failed runs.
3. Measure end-to-end wall time, user/system CPU, peak memory for the full process tree, I/O volume where available, output size, failures, and throughput. Also report compute versus plotting/serialization separately where measurable.
4. Use the same hardware, filesystem, inputs, required outputs, compression settings, and equivalent threading/resource limits. Record CPU model, RAM, OS/kernel, toolchains, dependency versions, build flags, environment, container digest, and commit.
5. Compare one-thread parity first; add scaling runs such as 1/2/4/8 threads up to hardware limits. Report upstream defaults as well as any tuned comparison.
6. Separate warm-cache and genuinely controlled cold-cache runs. Do not call a fresh process a cold-cache run. If caches cannot be controlled, disclose that limitation.
7. Randomize/interleave candidate/reference order; use a pilot and at least five measured repetitions for primary cases, increasing repetitions when variance warrants it. Repeat stochastic workloads over recorded seeds. Store all measurements, including failures.
8. Report median, spread, and confidence intervals with a stated method. Speedup = reference time / candidate time; report memory changes separately. Validate outputs for benchmark cases so a fast incomplete run cannot count as a win.
9. Use isolated hardware for publication results; shared GitHub runners can provide smoke checks, not stable performance evidence.

Deliverables: `benchmarks/protocol.md`, runners, workload manifest, raw JSON/CSV, environment manifest, scripts generating all result tables/figures.

Gate G6a: every command has a correctness-qualified benchmark result or an explicit failure entry. No selective omission of slowdowns. Any aggregate identifies its command/workload set and weighting; per-command data remains primary.

## Step 9 — Optimize only verified bottlenecks

Owner: one optimization agent per measured hotspot; reviewer retains comparator control.

1. Profile workloads and rank CPU, allocation, memory, and I/O bottlenecks.
2. Optimize one bounded component at a time; retain a simple verified reference path if useful for testing.
3. Rerun affected semantic tests and benchmarks on each change; rerun the integrated full suite for the release candidate.
4. Reject speedups caused by skipped reads, changed precision, output suppression, altered sampling, or different filtering unless released as an explicitly different mode.
5. Set performance acceptance targets after the baseline pilot. Scientific equivalence is mandatory; a universal speedup promise is not justified in advance.

Gate G6b: all final performance claims refer to the final verified commit and reproducible artifacts.

## Step 10 — Package and automate release verification

Owner: release agent; coordinator controls central workflow/build files.

Planned Rust workspace organization:

```text
Cargo.toml, Cargo.lock, rust-toolchain.toml
crates/core/              # scientific models, algorithms, numeric primitives
crates/formats/           # alignment, annotation, sequence, coverage I/O
crates/commands/          # independently owned command modules
crates/render/            # native plots and legacy script serialization
crates/cli/               # dispatch, flags, original executable aliases
crates/python/            # optional native Python bindings
bindings/python/qcmodule/ # Python namespace adapters and script launchers
compatibility/            # inventories, baseline lock, divergences
oracle/                   # upstream environment and runner only
verification/             # comparators, differential and API tests
fixtures/                 # small fixtures and generators
datasets/                # manifests and reproducible fetch scripts
benchmarks/               # protocol, runners, result analysis
reports/                  # generated summaries, release evidence
docs/                    # architecture, compatibility, user docs
.github/workflows/        # verification, packaging, scheduled suites
LICENSE, NOTICE, CITATION.cff, README.md, CONTRIBUTING.md
```

Implement PR checks for Rust formatting/lints, locked workspace builds, targeted unit/property tests, tiny differential cases, API imports, CLI smoke checks, and distribution installation. Run medium/full matrices on a schedule and release candidates. Keep expensive reference results cached only under keys including source, environment, case, input, and runner hashes.

Produce versioned source archives, supported binary distributions, Python packages/bindings where promised, and container images as appropriate. Test artifacts after installation outside the repository. Use a distinct project/distribution name and document how its compatible executables and import namespace replace upstream in an isolated environment; do not assume both packages can safely coexist.

For standalone artifacts, run all command smoke workloads and representative plot/conversion workloads in a minimal environment with Python, R, wigToBigWig, and other oracle helpers absent. Audit dynamic-library dependencies and bundled components. Test archive extraction, executable aliases, CPU compatibility, and paths containing spaces. Specify the minimum OS/ABI for each target; test static-link targets separately if offered. Test the optional Python wheel in a distinct environment.

Gate G7: clean-install end-to-end tests pass for every advertised platform/artifact; commands and bindings work without source-tree paths or the original package. Standalone CLI workloads, including plotting and file conversion, pass without separately installed language runtimes or helper executables.

## Step 11 — Prepare GitHub and scientific publication assets

Owner: documentation/release agent; coordinator approves factual claims.

1. Create a clear README with project status, supported upstream baseline, installation, original-command examples, compatibility table, performance table, and limitations.
2. Include source attribution and resolved license metadata; inventory third-party licenses and dataset redistribution conditions. Preserve applicable notices.
3. Add architecture/developer documentation, reproducibility instructions, issue/PR templates, contribution guidance, changelog, citation metadata, and a release checklist.
4. Store small evidence in Git; publish large datasets/results through versioned release artifacts or an appropriate archive with checksums and persistent identifiers. Do not commit huge BAMs into ordinary Git history.
5. Generate a methods/results report from raw evidence: baseline, hypotheses, fixtures, dataset selection, formulas, comparator rules, hardware, uncertainty, all results and limitations. A publishable repository supports evaluation; it does not guarantee journal acceptance.
6. Create the public GitHub repository under the selected owner/name, configure CI and branch rules available to that account, push a reviewable release candidate, and publish the versioned release after gates pass. Actual publication is a later execution phase, not part of this planning task.
7. Archive the release and reproducibility bundle for citation, then keep upstream-change monitoring separate from the frozen baseline.

Gate G8: a fresh checkout can rebuild, verify, and regenerate the report using documented commands and accessible inputs; every reported number traces to a raw result and exact source commit.

## Step 12 — Maintain compatibility after publication

Track upstream versions as distinct targets. New upstream behavior gets an inventory diff, new fixtures, an explicit compatibility profile, and re-verification before any claim is extended. Re-run release benchmarks after algorithm/dependency changes. Maintain regression fixtures for every resolved user-reported mismatch.

## Agent orchestration protocol

The coordinator owns architecture, task assignment, dependency ordering, shared-file changes, discrepancy decisions, integration, and release claims. Workers receive bounded evidence and code context rather than the entire repository.

Use this task ticket for each assignment:

```yaml
id: PORT-read_GC
objective: Implement read_GC.py and its required compatibility surface
baseline: 59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24
depends_on: [G3, profile-contract-v1, read_GC-oracle-cases]
read_first: [command-spec, relevant-upstream-functions, shared-interface-spec]
owned_paths: [command-module, task-specific-tests]
forbidden_paths: [shared-core, dependency-locks, golden-results, central-CI]
required_cases: [default, option-boundaries, reverse-strand, empty, invalid-input]
acceptance: [differential-pass, independent-check-pass, no-upstream-fallback]
handoff: [commit, changed-files, exact-checks-and-results, open-discrepancies]
```

Rules:

- Workers are not alone in the codebase: never revert another worker's edits; accommodate concurrent changes.
- Give each task exclusive file ownership and a separate branch/worktree. Do not run simultaneous cherry-picks or merges into the integration tree.
- Shared changes are requested through an interface-change ticket, implemented once by the owner, and rebased into dependent work.
- Implementers cannot loosen tolerances, rewrite goldens, or waive their own failures. Comparator/oracle changes need independent review and provenance.
- Bound tasks by an acceptance check that fits a small agent's context. Split a task when it spans multiple unrelated algorithms or shared modules.
- Review actual diff and executed results; an agent's “done” message is not acceptance evidence.
- Keep a machine-readable status board: blocked, ready, active, review, verified, merged. Only the coordinator marks merged acceptance.
- A blocked worker reports its minimal reproducer and interface need, then yields; do not repeatedly retry an architectural blocker.

### Scheduling with three available agent slots

The present environment permits a coordinator plus two workers. The work queues remain useful at larger concurrency, but do not assume unlimited workers.

| Wave | Worker 1 | Worker 2 | Coordinator |
|---|---|---|---|
| 0 | CLI/dependency inventory | API/algorithm inventory | Baseline/scope decisions and contract |
| 1 | Locked oracle runner | Fixtures and independent expectations | Approve case schemas and discrepancy policy |
| 2 | Differential harness | Core interfaces and bam_stat slice | Integrate/review interfaces and comparators |
| 3+ | Next ready command/helper task | Independent command/helper or review task | Integrate completed tasks; unblock shared work |
| Coverage completion | API gaps and bindings | Plots/remaining command branches | Audit all 33 commands and coverage matrix |
| Validation | Held-out/workflow verification | Benchmark harness and pilot | Triage discrepancies and freeze protocol |
| Optimization | One measured hotspot | Independent validation/benchmarking | Accept only verified improvements |
| Release | Packaging/install matrix | Reproducibility docs/report | Final audit and publication |

Stagger implementation and review so a worker never reviews its own patch. If only one worker slot is free, prioritize the critical path rather than creating idle agents. Larger teams can run more command tasks after G3; shared-core ownership remains singular.

Critical path: scope/baseline → executable oracle + fixtures → verified harness/core slice → all commands/API/plots → integrated verification → final benchmarks → installable artifacts → reproducible release.

## Definition of done

- All 33 installed commands and all agreed public API surfaces accounted for, with no hidden upstream fallback.
- CLI defaults, options, output artifacts, error behavior, plotting, and dependency behavior covered by the declared contract.
- Deterministic/numeric/semantic/stochastic compatibility separately and honestly reported; all exceptions visible.
- Independent scientific checks and held-out real datasets pass; zero unexplained mismatches.
- Benchmarks cover every command, retain raw data and failures, and reproduce all published tables/figures.
- Advertised release artifacts pass clean installation and actual workloads. Standalone Rust binaries require no Python/R/helper executable installation; optional Python bindings use the same Rust algorithms.
- Attribution, citation, accessible source, reproducibility bundle, and versioned GitHub release are complete.

No test suite proves equivalence for every possible input. Publish the tested scope and evidence, rather than claiming universal proof. Estimate calendar time only after the inventory, executable oracle, and first vertical slice reveal the actual compatibility burden.
