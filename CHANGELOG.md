# Changelog

All notable changes to this project are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/); this project has not yet made a
versioned release, so everything below is grouped under **Unreleased**.

For the full, granular history (including every individual bug fix and its verification
evidence), see `git log` and `compatibility/divergences.yaml` — this file summarizes at the level
of a feature or a closed divergence, not every commit.

## Unreleased

### Added

- A **real-data validation panel** (`datasets/`) implementing `testing.md` section 11.1:
  real RNA-seq reads from public archives, fetched with per-run MD5 verification from the
  archive, aligned by a pinned third-party aligner (STAR 2.7.11b) rather than by this
  project's own tooling -- a BAM from our own generator would encode our own assumptions
  about transcripts, strand and junctions, so testing the port against upstream on it
  would test the port against us. The 90-case differential matrix ran against the real
  panel via `RSEQC_REAL_DATA`: **86 ran and passed, 4 were skipped** because the panel has
  no equivalent fixture, with zero failures. **Re-run 2026-10-01 with the
  failure-closed comparators and the corrected rat annotation: 90 of 90 pass.** Both
  invocations include fixed synthetic/regression cases alongside real-panel
  substitutions, so the count is compatibility checks, not independent real-data
  validations. A study-level development/held-out split is
  recorded in `datasets/manifest.yaml`. The held-out runs have since been evaluated
  against the scientific endpoints; raw outcomes are mixed and are reported separately in
  `datasets/ENDPOINT_RESULTS.md`, which also records the three corrections made to how
  those outcomes were read (E1's strand expectation withdrawn, A6 withdrawn, E3's
  estimand narrowed to curve-to-curve similarity). **Differential agreement does not
  establish scientific validity**: the held-out panel exists to test the science on
  data neither arm was tuned against, not to re-confirm parity.
- A measured performance benchmark suite: a preregistered protocol, a workload
  generator over **real** hg38 sequence and real RefSeq BED12 annotation
  (`benchmarks/generate_workload_real.py`), a harness with process-tree resource
  measurement, a structural output-equivalence gate, adjacent matched-block
  repetitions and block-bootstrap intervals (`benchmarks/bench.py`), cost-driver
  sweeps (`benchmarks/scaling.py`), and table regeneration from raw data
  (`benchmarks/analyze.py`).
- **A replacement benchmark protocol** (`benchmarks/protocol-v2.md`), frozen 2026-10-01.
  The version-1 harness was found to have comparator false-pass paths, lost
  repetition pairing and an unsafe timeout handler, which makes its recorded gate
  passes and intervals uninterpretable; version 1's results are retained as
  historical development evidence and **no speedup figure from it may be
  published**. The replacement study has not yet been run at scale.
- **Tests for the validation machinery itself**, because a comparator that cannot
  detect a defect converts a regression into a reported speedup:
  `benchmarks/test_bench_harness.py` (49 tests),
  `verification/test_comparators.py` (46),
  `verification/test_run_diff.py` (21),
  `scripts/test_release_metadata.py` (28),
  `datasets/test_refgene_to_gtf.py` (9),
  `verification/check_command_contracts.py` (120 contract checks over ten commands),
  `verification/check_interop.py` (41 checks against htslib, including a clean
  refusal on a reference-based CRAM). Every comparator test asserts a **failing**
  result for a deliberate defect.
- **A per-command capability matrix derived by executing the binaries**
  (`verification/capability_matrix.py`, `benchmarks/capability-matrix.json`), embedded
  in the archive manifest and required to cover every shipped binary. The audit's
  Stage A asks for per-command records of accepted inputs, helper-dependent outputs,
  stochastic behaviour and known quirks; `compatibility/` held those in prose across
  three files, which no test verified. This feeds each command a real file of each
  declared format, re-runs it with a `PATH` containing no `Rscript` at all, and runs it
  twice to compare artifacts.

  Writing it surfaced **eleven probe defects**, each of which would have made the
  matrix a statement about the harness rather than about the commands:

  - a bare `-i <bam>` argv that failed five commands for reasons unrelated to the
    format under test, three of which do not read alignments at all (`bam2wig` needs
    `-s chrom.sizes`; `overlay_bigwig`, `normalize_bigwig` and `geneBody_coverage2`
    read BigWig; `FPKM_UQ`'s flags are entirely `--bam/--gtf/--info/--output`);
  - a "minimal PATH" assembled by editing a PATH *string*, which still found `Rscript`
    because `/usr/bin` was on it; replaced by a directory of symlinks to the tools a
    shell needs and nothing else;
  - a determinism probe that gave each run a different output prefix, so every
    generated `.r` script differed by construction -- every command with a plot contract
    looked stochastic;
  - `--sequencing` passed as two values (`pa 100`) when the flag takes a single
    `SE`/`PE` layout, so `clipping_profile` and `insertion_profile` never ran;
  - a per-format scratch directory keyed only by format, not by command, so a file one
    command wrote named `out` broke the directory `tin` needs there;
  - an output path built as `work / "out"` while the probe's cwd was the run directory,
    producing `out/out`, whose parent did not exist -- `RPKM_saturation` did all its
    work and then failed with a bare ENOENT;
  - `overlay_bigwig --action mean`, a plausible synonym for upstream's `Average`;
  - a BigWig/BED pairing whose chromosome names did not match, failing with "no valid
    BED12 records matched chromosomes in the BigWig";
  - PDF comparisons that reported every R-backed command as non-deterministic because
    R stamps `/CreationDate` into every plot;
  - a hand-written probe list that had drifted from the spec table and from the
    archive's allowlist, so 32 of 33 shipped commands were probed and one
    (`RPKM_saturation`) was silently absent while another appeared twice. The order and
    the specs are now cross-checked against the allowlist at run time.

  What the matrix found, having none of those in the way: **16 of 33 commands require
  `Rscript`** at runtime and say so when it is missing; **24 accept BAM, 14 accept SAM,
  and 9 are BAM-only because upstream's own `validate_args` rejects any other
  extension** (recorded per command with the mechanism, since a port accepting more
  than upstream would disagree on an input the user was told was invalid); three
  commands are **non-deterministic by design** -- `junction_saturation` and
  `RPKM_saturation` resample with an unseeded `random.shuffle`, `divide_bam` is
  unseeded unless `--seed` is given (DIV-0017) -- and three could not be judged at all
  because the helper they need is absent (`FPKM_UQ` needs `htseq-count`,
  `sc_editMatrix` and `sc_seqQual` need R's `pheatmap`). No unexplained
  non-determinism remains.
- **Per-job memory is now measured additive rather than assumed.** The documentation
  has long said "memory per concurrent invocation is additive: N concurrent jobs need
  N times the per-job figure". `verification/measure_concurrency.py` measures it by
  sampling `MemTotal - MemAvailable` while N invocations run: on the 8.2M-record rat
  alignment `bam2wig` peaks at 1.89 GB alone and 11.99 GB at eight concurrent, and
  `read_duplication` at 0.78 GB and 7.92 GB -- additive to within the instrument's
  resolution, with wall time nearly flat in concurrency (eight samples in 1.3x the
  time of one). That is the parallelism-across-samples deployment the audit asks about.
  Building the measurement took two corrections to itself. The baseline was first
  sampled *after* the processes were spawned, on the reasoning that their allocations
  should be "already included" -- which subtracts the commands' entire working set from
  their own peak and produced an aggregate (895 MB) below the single job's own peak RSS
  (1060 MB), a physically impossible result. And an "aggregate below per-job" flag
  against a fixed tolerance is not a real check: `/proc/meminfo` has a noise floor here
  of roughly 90 MB with nothing running, so the script measures that floor on the same
  machine before comparing. Post-run memory is recorded signed, because memory left in
  use after a batch and memory correctly returned to the page cache are different
  findings.
- **A memory sweep published a cost per transcript that was fitted through fifteen
  failed runs.** `measure_memory.py --transcript-drivers` varies the gene model while
  holding the alignment fixed, which is the only way to see the per-transcript state
  that `geneBody_coverage`, `tin` and `junction_annotation` retain -- a cost the
  per-record sweep cannot see, because restricting an alignment by coordinate does not
  change how many transcripts a command builds state for. Its first run reported
  **1234.5 bytes/transcript**: every invocation exited 1 in 0.00s, because the
  alignment path was relative and each command runs with its cwd set to a per-run
  scratch directory, so it resolved there. The slope function did not check exit
  codes, so it fitted a line through fifteen instantaneous failures -- the exact class
  of figure this script was written to stop publishing, reached by running it rather
  than by reading it. A slope is now never reported over runs with a non-zero exit, and
  input paths are resolved before use.

  With the runs actually succeeding, **all three slopes came out negative**, which is
  not a cost. `geneBody_coverage` and `tin` fall from 320 MB to 221 MB as the model
  grows 16x, because the read-side window they hold cannot retire a read until every
  transcript it might belong to is scored: a *smaller* model retains *more* reads, and
  that effect dominates the annotation cost the sweep is trying to isolate.
  `junction_annotation` is flat (75 -> 74 MB) and holds no such window, so there the
  reading is that its per-transcript state is small next to the fixed cost of the
  input. Both slopes are withheld with the reason recorded per command; a negative
  per-unit cost is now refused outright rather than reported.

  Holding the read stream small enough to remove the confound gives a real answer for
  one of the three: **`tin` costs about 1.7 KB per transcript** (27.9 -> 36.9 MB
  across a 16x model range, ~20x the measurement's own noise), which is its exon
  block index. `junction_annotation` (1.1 MB of movement) and `geneBody_coverage`
  (0.2 MB) fall below that noise and are bounded from above only. Getting there also
  replaced a bad rule: a fixed 20 MB minimum signal, derived from the ~90 MB
  system-wide noise floor `/proc/meminfo` shows on an idle machine -- the wrong
  instrument, applying its noise figure to a per-process measurement -- which
  discarded `tin`'s real 8.8 MB signal. Every sweep point is now the peak of three
  repeats with its spread recorded, and a slope is published only when the fitted
  quantity moves at least 3x that spread.
- **The bulk QC panel is now measured as a workflow, which shows per-command speedups
  are the wrong lever here.** `verification/measure_pipeline.py` runs the five commands
  a bulk RNA-seq QC step actually consists of -- `bam_stat`, `infer_experiment`,
  `read_distribution`, `geneBody_coverage`, `read_duplication` -- as one workflow
  against the 8.2M-record rat alignment: **48.8 s wall, 0.0135 CPU-hours, 1.06 GB peak,
  0.1 MiB of output per sample.** The two streaming commands this project has
  optimised, `bam_stat` and `infer_experiment`, are **1.5% of the panel's wall time and
  6 MB of its peak**; `geneBody_coverage` and `read_duplication` are 78%. Running whole
  panels in parallel -- the deployment the audit suggests may already be correct --
  scales throughput to **3.3x at four samples and 4.84x at eight** on 16 cores, with
  per-sample CPU rising from 0.0133 to 0.0214 CPU-hours as the panels contend for
  memory bandwidth.

  Each stage declares what its output must *contain*, not merely that a file exists,
  and a failed stage stops the panel: a stage run on the output of a broken one
  measures nothing, and including its time would understate what a working panel costs.
  Two of the five stages needed their output predicates written from what the commands
  actually emit -- `read_distribution` takes no output flag and prints to stdout, as
  upstream does, so the first version's `--out-prefix` was rejected outright.
- **The published workflow-impact memory figure was wrong, and is now marked
  superseded rather than quietly carried forward.** `verification/measure_pipeline.py`
  sampled its memory baseline 1.0 s *after* the panels started, so their allocations
  landed inside the baseline and were then subtracted from the peak — which measures
  how much the panels had **not** yet allocated, and is not conservative either, since
  a panel that has not peaked by 1 s inflates it. The symptom was 2,026 MB of "aggregate
  peak" for a single panel whose own peak is 1.06 GB. The concurrency harness had
  already been corrected for exactly this and the two disagreed, which is how the flaw
  was found; this file kept the old ordering. Both now sample the baseline before
  anything starts. The JSON file carries an explicit `SUPERSEDED` status block naming
  both this and the second problem — its panel includes `geneBody_coverage`, whose cost
  belongs to an unqualified command — so neither figure can be cited before recollection.
- **Stage C now has both artifacts, and the workflow one is executed rather than
  plausible.** [`workflows/qc_panel/Snakefile`](workflows/qc_panel/Snakefile) wires the
  measured five-command QC panel into Snakemake, and
  `workflows/qc_panel/test/run_panel.sh` runs the whole DAG on a real alignment, re-runs
  the same five commands directly, and diffs every artifact — so a wrapper bug fails the
  test instead of quietly changing somebody's QC numbers. Running it found five real
  errors in the Snakefile itself, each of which is now a comment where the mistake was:
  `infer_experiment` and `read_distribution` need `-r`; `-o X` produces `X.NVC.xls`,
  not `X`; with `--skip-plot`, `read_quality`'s only deliverable is `X.qual.r` and it
  writes no data table at all; awk runs its main block once per **input line**, so a
  `print` outside `BEGIN` with no input files never executes and the summary marker came
  out empty — and the checker then reported success over a zero-line file, because a
  loop over no lines checks nothing; and awk's `getline var < file` returns the *line*,
  not a count, so a `wc()` helper stored file content in every "count" field.
  [`recipes/rseqc-rust/meta.yaml`](recipes/rseqc-rust/meta.yaml) is the Bioconda recipe,
  built `--locked` and installing `THIRD_PARTY_NOTICES.txt` per Bioconda's Rust guidance,
  with `scripts/test_bioconda_recipe.py` rendering it and failing on anything conda-build
  would reject. Both carry the unresolved repository URL, so both are reported by the
  same release-metadata gate as `Cargo.toml` and `CITATION.cff`.
- **A third command was gated on a wrong declaration, and the declarations are now
  re-derived from live runs instead of trusted.** `RNA_fragment_size`'s entire
  per-transcript table is its stdout — upstream writes it to `sys.stdout` when no `-o`
  is given — and `EXPECTED_STREAMS` said `stdout: False`, so the gate rejected two arms
  whose stdout was byte-identical across 369,382 bytes. That declaration had been
  written from an observation, and the observation was wrong in a way worth recording:
  a probe looked for `label value` pairs, found none in a table whose *header row*
  carries the column names, and concluded there was no stdout. "No labelled stdout" was
  read as "no stdout". With `bam2wig` that is two commands gated on rules that do not
  describe them, and it is the audit's P0 gate finding recurring in the declarations.
  `verification/verify_stream_declarations.py` now re-derives every declaration from a
  live run of both arms — stdout presence, declared labels or substrings, and produced
  artifacts — and fails on any mismatch, so the class cannot recur silently. A stream
  failure also now retains the compared stdout: `RNA_fragment_size`'s mystery was a
  ten-second diagnosis once it did, and an unfindable one before.
- **The build invocation is part of the release artifact.**
  `cargo build -p rseqc-cli --release` and `cargo build --workspace --release` produce
  **different bytes for the same binary on the same source** — `bam_stat` hashes
  `d9841652…` under `--workspace` and `c497e29f…` under `-p`. So a capability matrix
  probed after a per-package build describes different bytes than the archive ships,
  even though the revision is identical, and a study run against one cannot be quoted
  alongside the other. The archive builder now spells out `--workspace --release
  --locked` with the reason, and `scripts/check_release_metadata.py` compares every
  command's recorded hash against the binary actually staged in the archive — a check
  that already earned its place by catching exactly this.
- **The whole-genome index built successfully, and then the script failed writing the
  record of what built it — for the same reason it had failed before, one line over.**
  `build_star_index.sh` put the contig set into both the index directory name *and* the
  `.built-…` marker name inside it. Bounding only the directory let STAR finish a 29 GB
  whole-genome rat index and then fail on `mkdir` of the marker, which left an index on
  disk with no record of its parameters — precisely the state the marker exists to
  prevent, and the same class of damage as the human rat index this project already
  suffered. Both names now take the bounded slug, whose contents still record the
  complete contig list. The test covers both names and is a real negative control:
  reverting just the marker line fails exactly the marker test. Separately,
  `--print-index-key` first sat *above* the variables it printed, so `set -u` made it
  exit 1 — caught only because the test invokes it rather than trusting it.

- **The whole-genome rat pilot now exists, and it shows the old panel was a heavily
  biased 20% of the library.** A 58-contig rn6 STAR index (150,217 annotated junctions,
  `genomeSAindexNbases` 11) aligned the same 17,168,681 input reads to
  **13,849,121 uniquely mapped (80.67%)** with **11,388,194 splices** — against
  **3,487,314 (20.31%)** and **2,255,750** on the three-contig panel. Roughly 80% of rat
  reads fall outside `chr1`/`chr2`/`chr10`, so every endpoint computed on that panel was
  a subset selected by where the reads happened to land, and "8.2M records" was never a
  record count for this library: it was the part that mapped at all. The whole-genome BAM
  holds **32,626,178 measured records = 8,584,340 read pairs**, reported as the two
  separate quantities protocol-v2 section 3 requires rather than as either alone.

  On that substrate `junction_annotation` is **byte-identical** between upstream and the
  port on all three data artifacts — `junction.xls` (5.9 MB), `junction.bed` (12.1 MB)
  and `junction.Interact.bed` (29.0 MB) — across 152,476 transcript rows. Both arms
  report `total = 11388194`, exactly STAR's own spliced-junction count, so neither is
  inventing or dropping junctions relative to the aligner. Only the generated R script
  differs, and only in the output directory it embeds in its own `pdf()` calls.
  Evidence and provenance: `datasets/heldout/endpoint_results/SRR1177982_wholegenome/`.

- **A scientific command silently produced a wrong answer where upstream refuses.**
  `verification/verify_stream_declarations.py` runs each declared command once and
  compares its outcome against upstream's, so it catches the *reverse* of what it was
  built for: `geneBody_coverage2` **succeeded** where upstream **refused**. Upstream calls
  pyBigWig's `values(chrom, pos-1, pos)`, which raises `Invalid interval bounds!` and
  exits 1 for a gene model extending past a chromosome's last base. The port's BigWig
  reader returns NaN there instead of failing, and that NaN was mapped to **0.0** — so
  the uncovered tail contributed nothing and the command emitted a confident coverage
  curve. That is the audit's "no silent metric loss and no successful corrupt output"
  failure, in a scientific command. It now refuses, naming the sampled position, the
  chromosome and its length. The exit status matches upstream; the message is clearer
  than pyBigWig's, and that difference is recorded (DIV-0027) rather than papered over
  by reproducing a worse diagnostic for byte parity. The new differential case fails on
  the old binary and passes on the fixed one, so it is a real negative control.

- **My own tool silently ignored five of six `--workload-for` overrides**, which is why
  the declaration sweep sat at 24/29 while every override worked when passed alone.
  `argparse` with `nargs="*"` and the default `store` action *replaces* the list on each
  occurrence rather than appending, so only the last flag survived. The consequence was
  misleading rather than loud: the affected commands were reported as "UNVERIFIED
  because the workload does not supply their input shape", which reads as a property of
  the commands. `action="append"` fixes it, and the sweep now verifies **28 of 29**
  declarations against live runs.

- **The whole-genome claim is now checked by the reference checker, not just recorded.**
  `verification/check_rat_reference.py` verifies the 58-contig index's build record and
  both of its digests, that the alignment is bound to *that exact* index build (compared
  as parsed key/value sets, because `align_run.sh` flattens the marker to one line and a
  raw string compare fails on formatting alone), and — the point of the exercise — that
  the recorded byte-identity is true of the files on disk, plus that both arms' junction
  total equals STAR's own spliced-junction count, so neither is inventing or dropping
  junctions relative to the aligner. 39/39 checks pass. Both new check families were
  mutation-tested: falsifying the recorded hash fails the hash check, and appending one
  byte to one arm fails the identity check. Two of the first-draft checks failed for
  reasons of my own making — the stamp records the digest of the UNPACKED `genome.fa`
  STAR read, not of the re-fetchable `rn6.fa.gz`, so comparing it to the `.gz` was
  simply the wrong comparison.

- **The study is now 18 rows on real data, and the README says which of them may be
  quoted.** 16 of 18 gates pass; the two failures are `geneBody_coverage` and `tin`,
  sharing the disclosed DIV-0024 defect, and their speedups are recorded and explicitly
  not claimed. Three things bound the claim and are stated rather than implied: the
  alignment carries reads on **3 of 58 annotated contigs**, so it is a diagnostic panel
  rather than the whole-genome pilot protocol-v2 section 2 asks for; it is **one session
  on one machine**, which section 6 says a bootstrap interval cannot make generalisable;
  and 5 blocks per row is short of the 10 the protocol starts from. These are indicative
  figures for a scoped beta, not publication figures — which is a narrower and more
  useful statement than the previous blanket "no performance figure may be published",
  and it names the three measurements that would change it.

- **A whole-genome rat index could not be built at all, because the index directory
  name embedded the contig set and overflowed the 255-byte filename limit.** That is
  the scope the production pilot needs, and the failure was `mkdir: File name too long`
  from inside the build script — after it had located STAR and begun unpacking the
  genome, so 90 seconds into a 45-minute job rather than at argument-parsing time. A
  long contig set is now summarised rather than dropped: the first two contigs stay
  readable, the count of the remainder follows, and a digest of the full set keeps the
  name unique, so rn6 and rn7 whole-genome indexes still get different directories. The
  complete contig list is recorded in the index's marker file either way, so nothing is
  lost by abbreviating. `build_star_index.sh --print-index-key` exposes the rule without
  building anything, and `datasets/test_star_index_naming.py` checks it: a short set is
  unchanged (otherwise every existing index directory would become unreachable), a
  whole-genome set fits in a filename, the summary says what it summarises, and two
  different contig sets can never collide — which is the same silent-damage mode that
  destroyed the human rat index earlier.

- **A partial re-collection silently shrank the study summary from 18 rows to 3.**
  Recollecting a subset of rows is a normal protocol operation — protocol-v2 says to
  recollect affected rows only when code, methods or an unresolved discrepancy justify
  it — and `results.json` was written from the current invocation alone, so re-measuring
  three rows left a summary that read as if only three commands had ever been
  benchmarked. The per-command records beside it survived, so no measurement was lost,
  but the aggregate is what a reader opens first. `merge_summary` now folds the new
  rows into any existing summary, carrying un-remeasured rows forward **with the
  revision that produced them** so a merged summary never presents old and new numbers
  as one run, and `--summarise-only` rebuilds the aggregate from the per-command records
  without re-measuring anything. Seven tests cover the merge and the rebuild, including
  that a corrupt summary is not fatal and that a re-measured row replaces its earlier
  version.

- **The gene-body divergence is two commands, not one, and the disclosure was wrong
  because of it.** The 18-row real-data study also failed `tin` — on the *same*
  transcript, `chr1:256806475-256813678`, `NM_013162`, TIN `93.75219481388487`
  upstream against `93.85965648789987` in the port, with the sample summary inheriting
  the difference. Both commands reach that behaviour through the same shared primitive,
  `tin::genebody_coverage_with_visited`, so it is one cause and two commands; confirmed
  on the synthetic `depth_cap` fixture with no real-data confound at all (TIN
  `39.570534359417486` against `39.2610946473031`). DIV-0024, the archive manifest's
  `known_limitations`, the README and this changelog all named only
  `geneBody_coverage`, and all four are corrected. Neither command is claimed. This is
  also the lesson worth keeping: a defect recorded against the first command it was
  observed in will silently under-report the second, so the differential suite now pins
  both with floor-asserted cases rather than one.

- **`geneBody_coverage.py` does not reproduce upstream's curve on real data, and that
  was previously unknown.** Running the repaired benchmark harness against the
  8.2M-record rat alignment produced a gate failure the 90-case differential suite could
  not: **76 of 100 bins differ, by up to 839 reads in one bin, in both directions.**
  Bisected to a single 904-base minus-strand transcript
  (chr1:256806475-256813678, NM_013162), committed as
  `verification/fixtures/genebody_divergence_minimal.bed12` so it cannot change
  silently.

  Established by construction, not by inference. A 100-base transcript was built in
  which each 10-base range exercises exactly one pileup knob -- plain reads,
  sub-threshold base quality, duplicate flag, an overlapping proper pair with identical
  bases in the overlap, an overlapping proper pair with mismatched bases, a
  paired-but-not-proper orphan, a `20M20D20M` deletion, secondary and supplementary
  flags, and a plain tail. BED12 base *i* maps to output bin *i*, so a disagreement
  localises to one mechanism: **all nine knobs agree exactly.**

  The cap was then isolated the same way. 8,000 / 8,001 / 12,000 uniform reads both
  arms cap at exactly 8000, so the cap's *value* is not the defect. Three versus 9,000
  copies of a single `20M20D20M` read are byte-identical at three and diverge at
  9,000, once the cap binds: **50 of 100 bins differ**, the port low at every one of
  them by 1..20 reads, reporting a flat 0 where clean reads demonstrably cover.

  So the cause is pysam/htslib `max_depth` semantics. The port reimplements the cap as
  a per-position budget over the **visited** set (`crates/commands/src/tin.rs`), so
  `is_del` pileups consume budget that upstream's **buffer** budget does not; upstream
  counts 8001..8010, 11..20, 8020 where the port reports a flat 8000, a flat 0, 8000.
  Which of htslib's effects dominates is not yet pinned down. **The mechanism is
  identified; the fix is not yet written.**

  Recorded as **DIV-0024, open**, and disclosed rather than filed as an accepted
  difference: the archive manifest's `known_limitations` states it, no speedup is
  claimed for the command, and `verification/fixtures/check_gene_body_divergence.sh`
  re-checks on every CI run that the defect is *still present* -- it goes green when the
  divergence is fixed, at which point the ledger entry and the CI step must change
  together. The differential suite needs a deep-coverage fixture; until it has one, its
  90/90 does not speak to this command on real data.
- **The differential suite could not have caught that, and now can.** Its
  `geneBody_coverage` fixture has depth around 40, which cannot reach a cap of 8000,
  so its 90/90 said nothing about this command on real data. Two changes close that:
  `verification/synthetic_data.py` now emits a `depth_cap` fixture of 8,100 `20M20D20M`
  reads that does bind the cap, and `run_diff.py` gained a `divergent_files` mechanism
  that asserts a known divergence **with a floor** instead of tolerating it. The new
  case `genebody_coverage_depth_cap_divergence` passes only while the two outputs still
  differ by at least the recorded amount; when the cap semantics are fixed the files
  become identical, the case fails, and it names the ledger entry to retire. A
  divergence that silently changes shape fails it too. Five unit tests cover that in
  both directions, because a check that only ever passes is decoration.
- **19 of the binaries read the whole alignment, computed every metric, and then
  discarded it.** Upstream's `validate_args` refuses an output prefix whose parent
  directory does not exist, before reading any input, via `parser.error()` and exit 2.
  Nineteen port binaries had no such check: they produced
  `No such file or directory (os error 2)` and exit 1 *after* computing everything --
  an error naming neither the directory nor the flag, indistinguishable from a missing
  input file. That is the audit's "no silent metric loss" failure in its purest form.
  Fixed by `rseqc_cli::require_existing_output_parent_or_exit`, one shared helper
  rather than twenty-two copies of a rule a new binary would forget. `bam2wig` and
  `bam2fq` are the **only** two upstream scripts without the check, so the helper is
  deliberately *not* applied to them (DIV-0026): matching upstream's inconsistency is
  required, because otherwise the port would reject an input upstream accepts.
  Five new differential cases cover the contract, including that counter-example.
- **A relative `--output-dir` silently invalidated an entire real-data study.** Each
  benchmark arm runs with `cwd` set to its own run directory, so a relative
  `--output-dir` was re-resolved *inside* that directory: every command taking
  `--out-prefix` was handed a path whose parent did not exist and failed in both arms.
  It read as two broken commands rather than one bad invocation -- and because the
  failure record kept only the tail of stderr, the retained diagnosis of a ten-rep
  failure was the truncated text `imum resident set size`. Both are fixed: the path is
  resolved once at parse time, and failures now retain `stderr_head` as well as
  `stderr_tail`, because an argparse error is printed at the *start* of stderr and GNU
  time's report at the end.
- **A failed benchmark row destroyed the report for every other row.** A row whose
  every repetition failed has `median_python = None`, and `{None:.3f}` raises
  `TypeError` -- so one bad command aborted the report and lost the numbers for the
  commands that passed, which is the opposite of the protocol's "retain failures,
  include all command rows". A missing measurement now prints as `MISSING`.
- **The benchmark gate failed `bam2wig` for a reason that was the gate's own
  declaration.** Both arms print the identical one-line `wigToBigWig` command to stdout;
  `EXPECTED_STREAMS` declared no stdout for it, so the gate rejected a byte-identical
  pair with "unexpected stdout content for a command that declares none" -- the audit's
  P0 gate finding recurring in the declaration rather than in the comparison, where a
  wrong declaration gates a command on a rule that does not describe it. Fixed by
  declaring the stream, and by a `stdout_mode: text` path that requires declared
  substrings *and* stream equality rather than parsed `label value` metrics. (Declaring
  it as a metric label instead -- the obvious fix, tried first -- makes the gate demand
  parseable metrics from a sentence and fails both arms.) `bam2wig` now passes its real
  data gate at 2.66x.
- **Independent truth cases** that a differential suite cannot supply:
  `crates/formats/tests/semantic_truth.rs` (CIGAR M/=/X/N/D rules, interval overlap
  boundaries, coordinate-frame invariants) and
  `crates/commands/tests/scientific_semantics.rs` (strandedness truth table, paired
  exon membership). Both are hand-derived rather than recordings of current
  output, and both name the upstream divergences they deliberately preserve.
- **A pinned, drift-checked oracle environment.** `compatibility/upstream.lock` now
  records the interpreter, the full `pip freeze`, the importable versions and each
  external helper's status; `verification/check_oracle_env.py` exits non-zero when the
  live environment differs. A recorded differential result can now be attributed to
  a specific reference environment.
- **Mandatory CI for this project's own validation.**
  `.github/workflows/release-validation.yml` runs on every push, every pull request
  and every tag; `release.yml` cannot publish without it. Previously no CI job ran
  the differential suite, the comparators, or the truth cases.
- **Release artefacts that are validated, not merely built.** The archive builder
  stages an explicit allowlist of 33 commands (it previously staged every executable
  it found in `target/release`), bundles `LICENSE`, `README`, `CHANGELOG`, citation
  metadata and generated third-party notices, writes a manifest with per-file SHA256
  digests and a declared glibc floor, and smoke-tests the *extracted* archive
  outside the source tree on a real workload and on paths containing spaces.
  `scripts/check_release_metadata.py` refuses a strict build while the repository
  URL is still a `TBD` placeholder, so a release cannot publish metadata that
  resolves to nobody.
- **A measured capacity record** (`docs/ENVELOPE.md`,
  `benchmarks/memory-sweep-tiny.json`, `benchmarks/memory-sweep-rat-8M.json`,
  `benchmarks/memory-failure-rat-8M.json`, `verification/measure_memory.py`,
  `verification/check_memory_failure.py`). Over a 174x range of record counts on the
  real 8.2M-record rat alignment: `bam_stat` flat at ~3 MB -- the control that makes
  the other rows meaningful -- `read_duplication` 130.5 bytes/record, `bam2wig` 207.5,
  with `bam2wig` peaking at 1.69 GB on the full alignment.
  `measure_memory.py` now refuses to report a slope over a range narrower than 2x and
  drops input windows that turn out to hold identical record counts: the first version
  of it fitted a confident bytes/record figure over a 1.22x "range" whose two larger
  points were the same measurement, and those figures (0.6 / 135.6 / 270.4) are
  withdrawn. A per-record cost measured on the synthetic fixture is 43.1 against 130.5
  on real data, which is why the fixture's numbers are not the headline.
  **What happens past the limit is now measured too**: under `ulimit -v` at 80%, 40%
  and 20% of their limit-free peaks, `bam2wig` and `read_duplication` both abort with
  no partial output and no input damage -- satisfying the contract's "complete output
  or identifiably incomplete output" -- but with Rust's internal
  `memory allocation of N bytes failed` rather than an actionable message. Recorded as
  a known gap. **No envelope is declared**.
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

- **The UCSC refGene -> BED12/GTF conversion was wrong in both outputs, and is now
  corrected against the genome rather than against documentation.** The converter
  treated refGene as 1-based inclusive; it is half-open **0-based**, with `txStart` and
  `exonStarts` as 0-based offsets and `txEnd`/`exonEnds` exclusive. Two successive
  revisions were wrong in opposite directions -- one shifted the transcript start by
  -1 while lengthening exons by +1, the next lengthened exons by +1 while leaving GTF
  starts unshifted -- which is why the errors did not cancel: BED12 and GTF disagree
  about inclusivity, so a length computed with an inclusive end is wrong even when its
  start is right. The audit's probe row (`exons [100,200), [300,400)`) produced BED12
  start 99 with sizes 101,101 and GTF starts 100,300; it now produces start 100 with
  sizes 100,100 and GTF starts 101,301. `datasets/verify_refgene_frame.py` confirms
  the frame from sequence: 60/60 transcripts whose CDS refGene marks complete begin
  with ATG at their 5' end and end with a stop codon under this reading, and 0/60
  under a 1-based reading. `datasets/test_refgene_to_gtf.py` pins the arithmetic
  against hand-derived cases including minus strand and the zero-coordinate
  boundary.
- **The rat index builder silently reused a stale contig subset.** `build_star_index.sh`
  regenerated `annotation.subset.gtf` only when the file was *absent*, so after the
  annotation above was corrected the corrected GTF was written and then ignored: STAR
  built a junction database from the old coordinates and reported success. Found while
  propagating the converter fix, and it is a worse failure than the original because
  the rebuild was indistinguishable from a correct one. The script now records a digest
  of its inputs and rebuilds when they change, and `verification/check_rat_reference.py`
  verifies the annotation on disk against the digests in `datasets/manifest.yaml`.
- **The rat alignment was never bound to the index it was produced from, and the
  recorded 0.617 figure was not reproducible.** The STAR index was rebuilt from the
  corrected annotation at 21:23; the rat BAM on disk had been aligned at 11:30
  against the previous index. Every annotation digest matched, because the annotation
  on disk really was the corrected one -- only the *alignment* was stale. Re-running
  `junction_annotation` on that BAM reproduced 25,851 of 44,176 annotated junctions,
  the intermediate figure, while the recorded table said 26,905 of 43,594, and
  nothing in the repository could say which was current. The alignment's `align.json`
  named the index's parameters but not its annotation digest, so "was this BAM
  aligned against this index?" had no checkable answer. `align_run.sh` now records
  the digest (it was copying a stamp that only began recording digests after the
  index-rebuild fix), and `verification/check_rat_reference.py` asserts the digest in
  `align.json` equals the stamp of the index on disk and that the BAM is not older
  than that stamp. The index's exon records are additionally required to be exactly
  the exon records of the corrected `rn6.gtf` restricted to the stamp's contigs,
  rather than a digest of a copy sitting beside it. The rat BAM was re-aligned, and
  both implementations were then run on the identical corrected inputs: the junction
  table, stdout and stderr are byte-identical, and a third invocation through the
  validator reproduces the same table digest.
- **`build_star_index.sh` gave the human and rat panels one index directory.** "Build
  the rat index" overwrote the human index's `Genome`, `SA`, `SA_*` and exon tables in
  place, and `align_run.sh` only checks that `SA` exists, so the next human alignment
  would have run against a corrupt index and produced plausible-looking output. This
  was not hypothetical: it happened while re-deriving the rat junction database, and
  the human index had to be rebuilt from hg38 afterwards. Each assembly now gets its
  own directory, named after the genome, annotation, overhang and contig set, because
  the two rat runs already differed by `sjdbOverhang` alone and a wrong overhang
  degrades splice-junction detection silently -- which is what E5 measures.
- **The rat assembly metadata conflated two assemblies.** `manifest.yaml` read
  `rn6 / mRatBN7.2`. UCSC serves rn6 as Rnor_6.0 and mRatBN7.2 as rn7, a different
  assembly; the field now names rn6 = Rnor_6.0 and says so explicitly.
- **A CRAM needing an external reference crashed with a Rust panic instead of an
  error.** `noodles-cram` 0.99's slice reader does
  `repository.get(name).transpose()?.expect("invalid slice reference sequence name")`,
  and this build resolves no external reference by design, so the `expect` always
  fired for such a file. The user saw
  `panicked at .../noodles-cram-0.99.0/src/io/reader/container/slice.rs:355` and
  exit 101, naming neither the input file nor the fact that a reference was the
  problem -- while the format's own doc comment promised "surfaces as a decode error".
  This is the shape most real CRAMs have: `samtools view -C -T ref.fa` records the
  reference URI and M5 in the header and does not embed the bases, so the defect only
  reached users who had done the ordinary thing. The decode now runs inside
  `catch_unwind` with the default panic hook silenced, and that one known panic is
  converted to an `io::Error` naming the file and the cause; any other panic is
  re-raised with the hook restored so a genuine bug keeps its location. Both branches
  are tested against htslib-written fixtures (`cram_no_reference.cram` must decode,
  `cram_external_reference.cram` must fail cleanly), and the refusal is also asserted
  in `verification/check_interop.py`. The exit status is what the interop check
  asserts, not the absence of panic text -- the hook is silenced, so a text-based
  check passed while the defect was still present.
- **`scripts/check_release_metadata.py` silently exempted `CITATION.cff`.** Its
  assignment pattern matched only TOML's `key = "value"`, not YAML's
  `key: "value"`, so the two placeholder-bearing files in the repository are written
  in two syntaxes and the check read as though it covered both while exempting the one
  citation tooling actually reads. The defect was invisible to the obvious test --
  run the check, see it fail, stop -- because `Cargo.toml` uses `=`, so the check
  still failed. It now matches both syntaxes, and `scripts/test_release_metadata.py`
  asserts a *failing* result for a `CITATION.cff`-only placeholder. The archive
  builder also runs the check against the staged `doc/` copy rather than the source
  tree, so what is verified is the file a user unpacks, and it bundles `Cargo.toml`
  alongside `CITATION.cff`.
- **The benchmark comparators accepted corrupted output.** The BAM signature read only
  `flag & 0xC0` and omitted base qualities, mate coordinates and tags, so a changed
  duplicate flag, an NM tag, a quality string or a mate position compared equal. The
  FASTQ reader stepped by four lines over `len - 3`, so a truncated final record was
  discarded and compared equal to a complete file. Text comparison normalised `nan`
  and `inf` like any other token, so a finite metric replaced by NaN compared equal --
  including symmetrically, where both arms had collapsed. Comparators now live in
  `verification/comparators.py`. The differential runner now imports the same module
  rather than keeping its own numeric-cell and numeric-table comparison, and
  `verification/test_comparators.py` asserts that it does -- a "shared" claim that no
  test checks is how the two copies drifted apart in the first place.
- **The benchmark output gate could pass two empty directories, and never compared
  stdout.** Only the last repetition's trees were gated, so nine of ten repetitions per
  row were unvalidated, and a command whose entire deliverable is on stdout received a
  gate pass with no metric checking whatsoever. Each command now declares its expected
  stdout, required labels and artifacts; a command with no declaration fails; every
  timed run's artifacts are validated outside the timing interval; and two empty files
  no longer count as agreement.
- **The benchmark timeout handler could kill the harness itself.** `subprocess.run`'s
  `TimeoutExpired` carries no PID, so the fallback `killpg(getpgid(0))` signalled the
  *caller's* process group. Executed against the previous code, the timeout test
  terminated the test runner with SIGKILL rather than reporting a failure. The runner
  now retains its `Popen`, terminates only the group it created, reaps children, and
  captures output to files so a surviving grandchild cannot hold a pipe open and hang
  the harness.
- **The benchmark's repetitions were not paired.** The schedule globally shuffled
  individual arms and the results were zipped by completion order; at the recorded seed
  20260929 all ten "pairs" had different repetition IDs, so the intervals described
  unpaired measurements. Repetition *i*'s arms now run back to back, with the
  within-block order randomised from the seed, and the pairing is read from the
  schedule rather than reconstructed from two lists.
- **`tin.py` no longer scores an alignment before discovering its output directory does
  not exist.** Upstream validates `-o` up front and exits 2; the port scored the whole
  BAM and then failed with a bare ENOENT minutes later, naming neither the directory
  nor the fact that the check was missed.
- **Two command-contract checks were added and one was measuring nothing.** A
  truncated alignment is now rejected rather than partially processed (a command that
  reads to the first decode error and reports what it got exits 0 with metrics
  covering an unknown fraction of the input), and an unwritable output directory must
  produce a diagnostic rather than a successful exit with the metric on stdout. The
  truncated-alignment check first passed against the 347-byte
  `bam_stat_basic.bam` fixture because its cut point (`max(1024, 70%)`) fell past
  the end of the file, so the "truncated" copy was byte-identical to the original; it
  now reports a fixture too small to truncate as a skip. 120 checks pass.
- **`geneBody_coverage.py` no longer aborts the process on a sample with no coverage.**
  Found by the T4 real-data panel, not by the synthetic matrix. When a sample has no
  coverage over the gene model, every coverage value is 0, so `max - min` is 0 and
  upstream's `(value - min) / (max - min)` raises `ZeroDivisionError`, which it reports
  as `geneBody_coverage.py: error: float division by zero` and exit 1. Rust's float
  division yields NaN instead of raising, so the NaN skewness then reached
  `partial_cmp().unwrap()` in the sample sort and **panicked**, aborting with SIGABRT
  (exit 101) where upstream printed one clean line and still wrote its coverage table.
  The port now reproduces upstream's error and exit code exactly, and the sort comparator
  is total (NaN last, tie-broken by name) so a NaN can never abort the process again.
  The coverage `.txt` it writes before failing is byte-identical to upstream's.
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

See the README's "Limitations" section for the current list, and
[`docs/ENVELOPE.md`](docs/ENVELOPE.md) for what has and has not been measured about
capacity.

**Closed since the 2026-10-01 [readiness audit](docs/READINESS_AUDIT_2026-10-01.md):**
the rat coordinate conversion and its frame are corrected and confirmed against the
genome; the benchmark comparators, output gate, timeout handler and pairing are
repaired and covered by executable tests; the oracle environment is pinned and
drift-checked; CI runs this project's own validation on every push and every tag; and
the release archive is allowlisted, bundled with its licence and citation metadata, and
smoke-tested outside the source tree. The differential matrix is **90/90**.

**Still open, and each is a real limit rather than a pending nicety:**

- **The replacement study has been run, and what it does not yet license is narrower
  than "no figure may be published".** Version-1 numbers remain void (defective
  comparators, unpaired repetitions). Under `benchmarks/protocol-v2.md` there are now 18
  command rows on real data, 5 matched blocks each, with per-run raw timings, failures,
  the interleaved schedule, and a per-binary SHA256 that was verified against the
  binaries on disk afterwards. **16 of 18 gates PASS; the two failures are
  `geneBody_coverage` and `tin`, which share the disclosed DIV-0024 cap defect, so no
  speedup is claimed for either.** Three things still bound the claim, and each is a
  specific missing measurement rather than a general reservation: the alignment covers
  **3 of 58 annotated contigs**, so it is a diagnostic panel rather than the
  whole-genome pilot protocol-v2 section 2 asks for; it is **one session on one
  machine**, and section 6 says a bootstrap interval on a single shared-machine run does
  not capture dataset or hardware generalisability; and there is **one dataset**. The
  protocol requires 10 valid matched pairs as a starting point and a predeclared
  precision rule; 5 blocks per row with intervals reported is short of that, so these
  are indicative figures for a scoped beta, not publication figures.
- **Three endpoint verdicts are withdrawn or narrowed.** E1's cross-lab strand
  expectation came from the wrong archive field and is now NOT_EVALUATED without
  explicit protocol metadata; A6's annotation-density explanation is withdrawn; E3's
  estimand is curve-to-curve similarity, not a mechanism. Feature stratification now
  exists and shows the aggregate E3 pass concealing stratum-level variation. The rat E5
  figure (0.617) is computed on corrected, index-bound reference preparation and
  reproduced by three independent invocations, but the stratum has been inspected five
  times across three preparation states, so it is diagnostic rather than confirmatory.
- **No capacity envelope is declared.** A cost per record and a measured peak are not a
  limit, and no run at the audit's 10M/50M-pair production targets has been performed.
  The whole-genome rat index needed for a realistic pilot is being built and has already
  cost one real obstacle: `build_star_index.sh` was killed mid-sort when the invoking
  shell exited, so the build must be detached with `setsid`. That is recorded here
  because it is the kind of thing that looks like an intermittent STAR failure and is
  not.
- **Publication is blocked** on the permanent repository URL and support channel, which
  are maintainer decisions; `scripts/check_release_metadata.py` enforces that rather
  than letting a `TBD` URL ship.
- **DIV-0003 (licence variant) is open.** Upstream's own metadata is self-contradictory;
  `doc/LICENSE` is the repository's existing file and the archive records the ambiguity
  rather than resolving it.
- **The release archive covers Linux x86_64 only.** No macOS, Windows or ARM artifact
  has been validated on its platform, and the Python API remains a stub.
- **Bounded property/fuzz suites exist; coverage-guided fuzzing and mutation-coverage
  measurement remain open.** "No panic" is not evidence of correct biological meaning.
- **The fixture-panel environment** (STAR 2.7.11b in a micromamba prefix outside the
  repository) cannot be rebuilt from a fresh checkout, so the real-data panel is not
  reproducible without that step.
