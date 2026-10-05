# Manuscript preparation outline

Prepared 2026-10-01. This is an outline for a future qualified software release; numerical results, author details and release identifiers remain pending. Use [the readiness audit](READINESS_AUDIT_2026-10-01.md) as the evidence and work-order reference.

## Working title and claim

**RSeQC-rust: a native implementation for RNA-sequencing quality control**

The proposed contribution is preservation of specified RSeQC CLI quantities with measured runtime/memory and deployment benefits on realistic workflows. Define the command/mode/input profile before writing the abstract. The original full native-plot/Python target remains a later milestone; scope the first paper to actually qualified capabilities.

## Scope decisions (2026-10-04, maintainer)

- **The claim is parity plus performance, not scientific validity on held-out
  data.** What is proven and therefore claimable: the port produces the same
  numbers as RSeQC 5.0.5 (byte-identical on 24/24 runnable commands against a
  held-out human library, 90/90 and 100/100 differential cases) and runs faster
  while using less memory (median 3.7x wall speedup, no command slower, no
  command more memory-hungry after the 2026-10-04 fixes). The held-out
  scientific-endpoint validation the manifest calls for is a different claim
  about biological interpretation, and its incompleteness does not block a paper
  whose contribution is preserved quantities plus measured cost.
- **DIV-0024 is scoped, not fixed.** `geneBody_coverage` and `tin` are claimed
  identical to upstream **except in regions deeper than pysam's `max_depth`
  of 8000 reads per position**, where the two implementations' pileup-cap
  semantics differ by construction (DIV-0024, open). The divergence is bounded
  (1 of 5359 rat transcripts, TIN off by 0.107; synthetic off by 0.309),
  disclosed, and pinned by a dedicated fixture so a fix would trip the test. No
  code change; the scope statement above is what makes the identity claim
  complete and honest.

## Statement of need and related work

Explain why existing RNA-seq QC is useful, and which installation, runtime or memory constraints a native implementation addresses. Identify target users and a concrete bulk-QC workflow. Compare with upstream RSeQC as the direct compatibility baseline; discuss other QC tools as related capabilities rather than treating different metrics as interchangeable speed comparisons.

Cite Wang, Wang and Li (2012), **RSeQC: quality control of RNA-seq experiments**, Bioinformatics 28(16):2184–2185, DOI [10.1093/bioinformatics/bts356](https://doi.org/10.1093/bioinformatics/bts356). Attribute biological algorithms to their original sources.

### Related tool: RustQC

RustQC (Seqera, Bioconda package `rustqc`, nf-core module `rustqc`) reimplements 15 RNA-seq QC tools in one single-pass binary, including eight RSeQC tools: `bam_stat`, `infer_experiment`, `read_duplication`, `read_distribution`, `junction_annotation`, `junction_saturation`, `inner_distance`, and TIN. It describes its outputs as format-compatible with upstream and MultiQC ([announcement](https://seqera.io/blog/rustqc/)). This project covers all 33 RSeQC commands under their original names with byte-identical output against the pinned upstream (see `compatibility/upstream.lock`).

| Tool | Scope | Output compatibility | Speed comparison |
|---|---|---|---|
| Upstream RSeQC 5.0.5 | 33 commands | baseline | sequential baseline; see `benchmarks/rustqc-comparison/RESULTS.md` section 4 (shared hardware) |
| RustQC v0.2.1 | 15 tools in one binary, including 8 RSeQC equivalents | equivalence verdicts in `benchmarks/rustqc-comparison/RESULTS.md` section 3 | single-pass `rna` at 1 and 8 threads; see `benchmarks/rustqc-comparison/RESULTS.md` section 4 (shared hardware; R time includes dupRadar/featureCounts/preseq/Qualimap/samtools) |
| This project | all 33 RSeQC commands under their original names | equivalence verdicts in `benchmarks/rustqc-comparison/RESULTS.md` section 3 | sequential and 8-concurrent medians; see `benchmarks/rustqc-comparison/RESULTS.md` section 4 (shared hardware) |

## Design and implementation

Describe the shared Rust command implementation, binary/alias surface, streaming BAM/SAM paths, interval indexing and sliding-window changes. Explain trade-offs, retained upstream quirks, stochastic policies, native/helper-dependent outputs, and the whole-file CRAM limitation if it remains in the release.

State versions and targets accurately. Pin the released source and binaries, the reference environment and comparison profiles. Provide a build/capability manifest.

## Validation methods

Separate output compatibility from independent scientific expectations. Describe test inputs, study-level split, complete artifact/stream checks, tolerance rules, command errors and discrepancy classification.

Include independently specified annotation-coordinate, paired-exon, pileup and CIGAR expectations. Archive the rat converter defect and corrected preparation history; existing rat scores do not qualify a biological result.

**The three-defect chain is a method result, not a footnote.** In this project the
same failure mode -- trusting a file rather than checking the link between two
artifacts -- produced three separate defects on the rat stratum: a refGene conversion
in the wrong coordinate frame; a STAR index that silently reused a stale contig subset
so a "successful" rebuild was built from the old coordinates; and a BAM that outlived
the index it was aligned against, so every annotation digest matched while the
endpoint result was not reproducible. Each was found by the same check (digests and
derivations, not trust) and each was invisible to the check before it. That is
reportable as a validation-method contribution, and the checks themselves
(`datasets/test_refgene_to_gtf.py`, `datasets/verify_refgene_frame.py`,
`verification/check_rat_reference.py`, `scripts/test_release_metadata.py`) are the
evidence.

The rat stratum must still be described as exposed: it was inspected five times
across three preparation states, so its 0.617 junction figure is diagnostic. A fresh
independent rat sample is required before it can carry a biological claim.

Describe real whole-genome data with actual library-preparation evidence. Use identical alignment/annotation/options for reference and candidate comparisons. Use fresh confirmation for conclusions changed after held-out exposure. Keep confounded chemistry/lab observations exploratory.

## Performance methods

Use the repaired, frozen benchmark protocol. Define actual matched blocks, primary elapsed-time estimator, resource limits, equal delivered work, output validation, hardware/storage/cache settings and raw sampling.

Include representative production sizes, memory-sensitive cost drivers and a full workflow at intended concurrency. Label pairs versus alignment records, data-only versus plotting/conversion modes, and failures/unsupported modes. Explain any startup-subtraction sensitivity analysis.

**The measurement harness needed its own validation, and that is reportable.** Three
separate defects in this project produced measurements that looked valid and were not,
and none was caught by looking harder at the numbers:

- The gate rejected *correct* command pairs. `bam2wig` prints a one-line stdout and
  `RNA_fragment_size`'s whole per-transcript table is stdout; both were declared as
  writing none, so byte-identical arms failed on "unexpected stdout content for a
  command that declares none". The second declaration had been written from a probe
  that looked for `label value` pairs, found none in a table whose header row carries
  the column names, and concluded there was no stdout. Declarations are now re-derived
  from live runs (`verification/verify_stream_declarations.py`).
- A relative `--output-dir` made every command taking `--out-prefix` fail in both arms,
  because each arm runs with its working directory set to its own run directory. It
  presented as two broken commands rather than one bad invocation, and an entire
  real-data study was launched before it was caught.
- A workflow memory baseline sampled 1.0 s *after* the panels started, so their
  allocations sat inside the baseline and were then subtracted from the peak. The
  concurrency harness had already been corrected for exactly this, and the two
  disagreed, which is how it was found.

Each is a measurement-validity result, and each is the kind that a paper claiming
performance numbers has an obligation to report: they are properties of how the
measurements were taken, they change reported magnitudes, and in two cases they changed
a verdict from fail to pass.

**Two scientific commands are not qualified, and the paper must say so in the methods,
not a footnote.** `geneBody_coverage` and `tin` do not reproduce upstream's output on
real data, from one shared cause: they reach it through the same pileup primitive, and
the port budgets `max_depth` per position over visited pileups where upstream budgets
over the pileup buffer. `tin` is the more consequential of the two to disclose, because
transcript integrity number is a widely used QC metric, so a reader who trusts a `tin`
figure from this port is trusting a number that is wrong in the fourth significant
figure on a deep transcript. The cause was isolated by construction rather than inference: a
transcript was built in which each base range exercises exactly one pileup mechanism,
and all nine -- quality threshold, duplicate flag, overlap rewriting, orphan handling,
deletion skipping, secondary and supplementary flags -- agree exactly. The divergence is
in pysam/htslib `max_depth` semantics, where the port's per-position budget over
*visited* pileups differs from upstream's budget over the pileup *buffer*. The
differential suite passed while this was wrong, because its fixture was two orders of
magnitude too shallow to bind the cap; a dedicated deep fixture now pins the divergence
numerically, so a fix trips the test instead of quietly retiring the disclosure. The
command is excluded from every speedup claim.

## Results placeholders

- Compatibility matrix: [pending qualified release run]. The 90-case differential
  suite passes, but it is a compatibility result, not a scientific one.
- Independent truth cases: available now — 17 format/semantic truth cases
  (`crates/formats/tests/semantic_truth.rs`), 15 scientific-semantics cases
  (`crates/commands/tests/scientific_semantics.rs`), 9 refGene converter cases, and a
  60/60-vs-0/60 coordinate-frame confirmation against the rat genome sequence.
- Real-data concordance on corrected, index-bound inputs: **available for one command on
  a whole-genome rat alignment** — `junction_annotation` byte-identical between upstream
  and port on all three data artifacts (`junction.xls` 5.9 MB, `junction.bed` 12.1 MB,
  `junction.Interact.bed` 29.0 MB) across 152,474→152,476 transcript rows, with both
  arms reporting `total = 11388194`, exactly STAR's own spliced-junction count
  (`datasets/heldout/endpoint_results/SRR1177982_wholegenome/`).
  It is not a whole-stratum validation, and the sample is **exposed**: it was inspected
  across several preparation states while the depth-cap defect was found on it.

  The substrate itself is a finding worth stating rather than burying. The earlier
  three-contig rat panel retained only **20.31%** of uniquely mapped reads, because most
  rat reads fall outside `chr1`/`chr2`/`chr10`; the whole-genome alignment of the same
  library retains **80.67%**. So endpoints computed on a contig subset are not a smaller
  sample of the same thing — they are a sample selected by where reads happened to land,
  and "8.2M records" was the part that mapped rather than a record count for the library.
  A paper reporting a per-transcript or per-junction statistic must say which substrate
  produced it, and a contig subset is a bias, not a sample.
- Primary end-to-end/runtime and memory results: **available for 18 command rows on
  real data** ([`benchmarks/results-rat-real-8M/results.json`](../benchmarks/results-rat-real-8M/results.json)),
  5 matched blocks each, 16 of 18 gates PASS. The two failures are `geneBody_coverage`
  and `tin`, which share the disclosed `max_depth` defect, so no speedup is claimed for
  either. The gate caught three separate declaration defects during this study that
  would otherwise have appeared as command failures. What it is NOT: whole-genome
  alignment (3 of 58 annotated contigs carry reads), a second session, or a second
  dataset — protocol-v2 section 6 requires those before a generalisable claim.
- Memory cost drivers: measured on the 8.2M-record rat alignment over a 174x record
  range, with repeats and a rule that refuses to report a slope from failed runs,
  negative costs, or signal below 3x the repeat noise
  ([`benchmarks/memory-sweep-rat-8M.json`](../benchmarks/memory-sweep-rat-8M.json),
  `docs/ENVELOPE.md`) — `bam_stat` flat over the range; `read_duplication`
  **130.5 B/record**; `bam2wig` **207.5 B/record** with a full-rat peak near
  **1.69 GB**. Per transcript, `tin` is about **1.7 KB**, while `geneBody_coverage`
  and `junction_annotation` fall below measurement repeatability at this model scale,
  so no per-transcript cost is claimed for them. Report as a cost-driver
  characterisation with its stated range, not as an envelope.
  An earlier tiny-fixture sweep (29x range, 621,783 records) is retained as
  development history and is superseded; its per-record figures must not be quoted.
- Scaling and maximum measured operating point: [pending production-size pilot].
- Pipeline impact: **the wall-time and CPU-hour figures stand; the memory figure does
  not.** A five-command bulk QC panel costs **48.8 s wall and 0.0135 CPU-hours** on the
  8.2M-record rat alignment, with `geneBody_coverage` at 45% of the wall time and
  `read_duplication` at 33%; the two streaming commands this project optimised are
  1.5% of the panel. Running whole panels in parallel scales throughput 3.3x at four
  samples and 4.84x at eight on 16 cores. That bounds any per-command speedup claim:
  the optimised commands are not where the time is. The **aggregate memory column in
  [`benchmarks/pipeline-impact-rat-8M.json`](../benchmarks/pipeline-impact-rat-8M.json)
  is marked SUPERSEDED** and must not be quoted: its baseline was sampled after the
  panels started, so it subtracted their allocations from their own peak, and the
  panel's composition includes `geneBody_coverage`, which is not qualified. Recollect
  before any memory figure from this file is used. Report the rest as one sample, one
  machine, one panel — not as a total-pipeline figure.
- Failures, disclosed divergences, and limits: [populate from the exact study].

Do not transfer the version-1 benchmark confidence intervals into these
placeholders: its comparators, gate and repetition pairing were defective, so its
intervals are uninterpretable. Do not transfer the rat annotation-density
conclusion (A6, withdrawn) either.

## Main figure and supplement

For an application note, plan one composite main figure: concordance, runtime/memory scaling and workflow impact. Put the full command/mode matrix, methods, negative controls, raw paired measurements, failure records and additional stress curves in the supplement/archive.

For a longer paper, separate compatibility, cost-driver scaling and deployment impact into figures only if the additional research contribution justifies the format.

## Discussion and limits

Explain the extent to which an unchanged QC quantity yields unchanged downstream interpretation. Describe the supported operational envelope, untested modes and retained upstream behavior. Avoid extrapolating one hardware platform or study into universal industry-scale claims.

## Availability and authorship

Provide permanent public source, release tag/commit, archive DOI, container digest, license, example data, test/benchmark rerun instructions, and support/contribution routes. Ensure the source/test data supporting figures are accessible.

**Blocked on maintainer decisions, not on analysis.** The permanent repository URL and
support channel are unset (`Cargo.toml` and `CITATION.cff` still carry `TBD`, and
`scripts/check_release_metadata.py` refuses a strict build while they do). No DOI or
public history exists yet. JOSS additionally requires more than six months of public
development and demonstrated research use, neither of which this project has: the
local history runs 2026-09-16 to 2026-10-01 with no tags. Venue route stays
Bioinformatics Advances, with Bioinformatics Application Notes as the alternative if
the demonstrated advance fits that scope. The technical reviews in this repository
are not human domain-expert or journal-reviewer sign-off, and must not be described
as such.

**AI assistance disclosure.** This work has been developed with AI assistance
throughout. Both target venues require it to be disclosed, and JOSS additionally
requires human review of assisted work. The role of assistance across the code,
tests, analysis and manuscript drafting should be stated plainly rather than in a
generic acknowledgement.

Human authors should confirm affiliations, contributor roles, funding, conflicts and their own scientific review. Document AI tools/versions where known, their assistance in code/docs/manuscript, and human verification; do not claim verification that has not occurred. [Bioinformatics author guidance](https://academic.oup.com/bioinformatics/pages/author-guidelines), [JOSS disclosure policy](https://joss.readthedocs.io/en/latest/submitting.html).

## Submission decision

Prefer an application-note venue once the realistic study, artifacts and external pilot are complete. Evaluate JOSS only after its current public-history and research-use conditions are met. Confirm venue rules and costs again at submission. A draft narrative can be written now; measured results and publication-ready claims require the evidence gates.

