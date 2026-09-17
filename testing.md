# Scientific validation, compatibility, and performance testing plan

Prepared: 2026-09-17. Audited implementation: `b52114315e98e84b201dc769ddf750125570b1a0`.

This is an implementation and acceptance plan, **not a statement that validation has been completed**. It expands Steps 1–10 of [PORTING_PLAN.md](docs/PORTING_PLAN.md) into testable obligations. The [standalone architecture decision](docs/decisions/0001-target.md), [command inventory](compatibility/commands.yaml), [API inventory](compatibility/api.yaml), and [divergence ledger](compatibility/divergences.yaml) remain inputs to the work. Their assertions must be checked against executable evidence.

The objective is to make a precisely scoped claim defensible: a named command, operating mode, input domain, baseline, and released artifact produces scientifically understood results, with measured compatibility and performance. Testing cannot establish correctness for every possible input. Unexamined domains and unresolved discrepancies remain visible limitations.

Navigation: [audit and shared-code impact](#2-evidence-at-the-audited-revision), [runner requirements](#5-make-the-verification-harness-capable-of-failing-correctly), [33-command matrix](#7-all-command-acceptance-matrix), [independent mathematics](#8-independent-mathematics-properties-and-comparators), [biological validation](#11-real-data-and-biological-validation), [performance](#12-fair-performance-and-scalability-measurements), [work order and gates](#14-ordered-implementation-work-and-release-gates), [audit recipes](#appendix-a-reproducible-audit-cases-to-promote-into-permanent-tests).

## 1. Separate the claims and their acceptance criteria

| Claim | Evidence required | Evidence that is insufficient |
|---|---|---|
| The implementation builds | Locked clean build and installed-artifact smoke tests on each supported platform | A successful build on the developer's machine |
| Invocation compatibility | Original names, options, defaults, streams, exits, and output locations exercised by real processes | Help text or aliases alone |
| Output compatibility | Per-artifact comparison against an identified upstream environment, including error and optional-output branches | Matching selected columns or matching two failures |
| Scientific correctness within a stated domain | Independent expectations, boundary tests, biological controls, and downstream impact assessment | Agreement with an upstream implementation that shares the same bug |
| Stochastic equivalence | Defined sampling unit, reproducibility policy, invariants, and adequately powered equivalence tests | Same integer seed, similar-looking curves, or a nonsignificant difference |
| Performance improvement | Correctness-qualified, paired end-to-end measurements with uncertainty and matched work | A faster run that omits plots, conversion, indexing, or reads |
| Standalone distribution | Installed commands complete their advertised workloads without separately installed Python, R, or helper executables | Rust source code or absence of dynamic Python linkage |
| Full package compatibility | All 33 command contracts and the inventoried Python API pass their own gates | All binary names exist |

A publication using one validated command need not wait for every unrelated command. Its methods and release evidence must explicitly identify that restricted scope. An unqualified full-port claim requires the whole inventory. No single aggregate “percent compatible” replaces these distinctions.

A scope restriction must be identifiable from recorded inputs and parameters, with an enforced precondition or a reviewed way to detect violations. Calling an affected case “rare in practice” is not a restriction. Each study run needs a provenance record connecting its actual BAM/annotation/input hashes, options, artifact, and validation profile; qualification of one build does not establish that a later study used it correctly.

Maintain two independent result fields: **compatibility status** and **scientific validity status**. A preserved upstream bug can pass the first and fail the second. A deliberate correction can pass the second and fail the first. If a corrected behavior is introduced, give it an explicit, versioned profile and separate expectations; never change the meaning of an existing profile silently. Such a profile is proposed policy, not an existing implementation feature.

## 2. Evidence at the audited revision

The repository review established:

- A locked release build succeeds, and 224 Rust unit tests pass. This does not measure semantic branch coverage or biological validity.
- The five existing differential cases pass: `bam_stat_basic`, `read_NVC_basic`, `read_GC_basic`, `read_duplication_basic`, and `read_quality_basic`.
- Those cases reuse one 10-record BAM. Every record has a `20M` CIGAR, a 20-base all-A sequence, and Phred-40 qualities. The fixture declares coordinate sorting but places an unpositioned unmapped read before subsequent mapped reads. It is useful for sequential flag-counting tests, but must not serve as evidence for sorted/indexed workloads.
- The differential suite has no semantic BAM/BigWig comparator, independent biological dataset panel, or comparator mutation tests. Four cases skip stream comparison; only selected output files are checked.
- `cargo fmt --all -- --check` fails. Strict Clippy reports 14 diagnostics before stopping; that is not proof that only 14 issues exist.
- Rendering and Python binding crates contain only documentation. Several commands still invoke `Rscript`, `wigToBigWig`, or `htseq-count`. `sc_seqLogo` always fails after producing its count matrix.
- The oracle lock identifies a source snapshot but does not lock the complete execution environment. The local diagnostic probes used pysam `0.24.1` and its reported samtools version `1.24`; these observations are not an environment lock.

The original five passes are **smoke evidence**, not command-wide qualifications. The audit below found failures outside their input space without changing production code.

### 2.1 Confirmed loopholes and defects

`P0` blocks reliance on affected scientific results or on a verification pass. `P1` blocks an affected compatibility, performance, or release claim. “Reproduced” refers to the audited revision and local oracle, not all upstream releases.

| ID | Priority / evidence | Finding and consequence | Required closure |
|---|---|---|---|
| AUD-01 | P0, reproduced | [run_diff.py](verification/run_diff.py) accepts equal exit-1 failures with empty stdout: it reports a pass with zero matched labels. | Declare expected exit and nonempty schema per case; missing metrics and failed positive cases must fail. |
| AUD-02 | P0, reproduced | Its label parser treats `1e-3` and `1e+3` as the same value `1`; duplicate labels overwrite earlier values. | Strict complete-token parsing, duplicate rejection, explicit field set, and comparator negative controls. |
| AUD-03 | P0, reproduced | Selecting a valid case plus an unknown case name returns success after silently omitting the unknown case. | Reject every unknown selector; report requested, selected, executed, skipped, and passed counts separately. |
| AUD-04 | P0, reproduced | [FPKM_count](crates/commands/src/fpkm_count.rs), `count_total_fragments`, omits the mapped mate's exon-overlap condition. With `--only-exonic`, the audit fixture has an upstream denominator of 1 and a Rust denominator of 2; FPKM changes from 20,000,000 to 10,000,000. | Add the paired exon-membership truth table and preserve the Appendix A reproducer before fixing the implementation. |
| AUD-05 | P0, reproduced | [reference_span](crates/formats/src/cigar.rs) ignores `=`/`X` but is used to emulate indexed fetch in [RNA_fragment_size](crates/commands/src/rna_fragment_size.rs). A `20=` read starting exactly at the transcript boundary is dropped: upstream reports one fragment of length 40; Rust reports zero. | Distinguish a legacy CIGAR helper from standards-based reference span; verify every caller against its actual upstream path. |
| AUD-06 | P0, reproduced | [gene-body coverage](crates/commands/src/genebody_coverage.rs) uses the TIN coverage helper without upstream's default depth cap. At one site, upstream reports 8000 and Rust 8001. | Pin and reproduce the applicable pileup contract, or explicitly exclude this compatibility domain. Test both consumers and depth/quality/order interactions. |
| AUD-07 | P0, documented and source-confirmed | [TIN](crates/commands/src/tin.rs) lacks overlapping-mate handling; gene-body coverage shares it. DIV-0011 understates the scope by focusing on TIN. | Test overlap agreement/disagreement, adjusted qualities, pair flags, missing mates, and both commands. Quantify the effect; do not assume a universal direction of TIN bias. |
| AUD-08 | P0, source-confirmed risk; output reproducer pending | FPKM's index stores `end = start + sequence_length` for region retrieval. Upstream's indexed fetch uses alignment span before the command's separate query-length-based counting rules. Splicing/deletions can change candidate selection. | Add reads starting before a transcript whose reference span reaches an exon or mate inside it. Keep fetch semantics separate from legacy fragment formulas. |
| AUD-09 | P1, source-confirmed | `divide_bam` uses Rust `StdRng`, while upstream uses Python `random.Random`; the seed alone does not establish identical assignment. Saturation commands also use different shuffling implementations. | Decide exact seeded compatibility versus a disclosed stochastic profile; implement the appropriate tests in Section 9. |
| AUD-10 | P1, source-confirmed | The runner executes both processes with `cwd=REPO_ROOT`, has no timeout/resource limits, inherits environment, checks selected artifacts, deletes scratch evidence even on failure, and does not bind binaries to the tested source revision. | Implement the runner contract in Section 5 before treating green results as release evidence. |
| AUD-11 | P1, inspected | One unrepresentative BAM and mostly implementation-local expectations leave nucleotide composition, splice/CIGAR, quality, index, and real-data behavior underconstrained. | Independent fixtures and the 33-command matrix; validate the fixture generators themselves. |
| AUD-12 | P1, source-confirmed | `split_bam` ignores `--overwrite`; BAM writers rely on drop in several CLIs; some unimplemented requested outputs produce only warnings. | Existing-output and input-alias tests; propagate finalization errors; assert every promised artifact exists and is readable. Test unsupported branches explicitly. |
| AUD-13 | P1, source-confirmed | The standalone and full-package contracts are incomplete: external helpers, empty renderer/bindings, BAM-only paths, missing compression/index branches. | Maintain capability records per artifact/option and complete clean-install gates. A skipped feature is not a pass. |
| AUD-14 | P1, source-confirmed risk | Whole-BAM indexes, per-base `BTreeMap` coverage, repeated chromosome scans, and whole-output `String` buffers can dominate time/memory. FASTQ sequences/qualities are loaded before sampling limits apply. | Scale reads, covered bases, transcripts, overlap density, output size, and sequence count independently; measure peak process-tree resources. |
| AUD-15 | P1, source-confirmed risk | `python_str_float` explicitly omits scientific-notation compatibility. NaN-sensitive sorts, zero-variance normalization, unchecked coordinate casts, and differing summation orders need audit. | Numeric boundary and empty/constant-input suites, including magnitude extremes and integer overflow. |
| AUD-16 | P1, inspected | Inventories and documentation have drifted. For example, the `bam2wig` inventory omits conversion dependencies and does not match stranded output casing. Some ledger classifications are outside its stated vocabulary; “accepted” includes unfinished functionality. | Rebuild contracts from pinned source plus observed executions; validate schema and separate decision acceptance from implementation/test completion. |

Preserving a bug in an upstream helper does **not** justify applying that helper to a different upstream API. AUD-05 is a concrete example. Treat comments such as “exactly equivalent,” “same defaults,” and “lossless simplification” as hypotheses requiring tests.

### 2.2 Shared-code impact map

Every shared change must trigger the listed consumer tests, in addition to local unit tests. Keep this map machine-readable once the harness is extended; recompute it when imports or call paths change.

| Shared behavior | Consumers requiring regression coverage |
|---|---|
| BAM reader, records, flags, tags, coordinates, headers | Every BAM-consuming command, including the single-cell commands and all BAM writers |
| `fetch_exon_blocks` and its legacy soft-clip/`=`/`X` behavior | `read_duplication`, `read_distribution`, `inner_distance`, `bam2wig`, `RPKM_saturation` |
| `fetch_intron_blocks` | `inner_distance`, `junction_annotation`, `junction_saturation` |
| `reference_span` | `RNA_fragment_size`; audit future callers before reuse |
| BED extraction and interval set operations | Annotation/fragment commands, split classification, FPKM, TIN background, exonic BigWig normalization |
| TIN read index and position coverage | `tin`, `geneBody_coverage` |
| Percentile helpers | `geneBody_coverage`, `geneBody_coverage2`; separately test TIN and saturation's different percentile definitions |
| Strand-rule parser in FPKM module | `FPKM_count`, `bam2wig`; independently test the separate saturation parser |
| Python numeric formatting | FPKM variants, TIN, gene-body variants, saturation, BigWig normalization, single-cell matrices, and every new caller |
| BigWig reader | `geneBody_coverage2`, `normalize_bigwig`, `overlay_bigwig` |
| Heatmap script generator | `sc_editMatrix`, `sc_seqQual` |

Do not infer cross-command equivalence merely because two functions share a name. Query coverage, fragment coverage, pileup entries, aligned bases, exon blocks, and read-start counts are distinct quantities.

### 2.3 Checks available today

From the repository root, the current baseline can be inspected with:

```bash
cargo build --workspace --release --locked
cargo test --workspace --locked
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
PYTHONDONTWRITEBYTECODE=1 oracle/venv/bin/python3 verification/run_diff.py
```

Run the commands separately so a formatting failure does not suppress other diagnostics. The last command requires the separately provisioned local oracle. These commands do not execute the future matrix, property tests, biological validation, or publication benchmarks described below. The historical audit table above records the runner's former AUD-01/02/03/10 weaknesses; the current working-tree runner has executable checks for those false-pass paths, but provenance binding and resource isolation still need to be completed for a publication release.

#### Working-tree first executable gate (2026-09-17)

The first regression gate is now executable and reproducible from the repository:

- `verification/run_diff.py` runs 11 cases: the original five smoke cases plus six focused scientific regressions for paired FPKM membership, indexed reference-span selection, `=` CIGAR spans, gene-body depth capping, gene-body overlapping mates, and TIN overlapping mates.
- `verification/test_run_diff.py` contains eight harness tests covering strict numeric parsing, duplicate labels, unknown/failed/empty cases, required labels, semantic numeric comparison, and timeout handling.
- The gate is fail-closed: expected-success cases reject non-zero exits, empty labeled output, missing required metrics, unknown selectors, and timed-out processes. Failed-run scratch evidence can be retained with `RSEQC_KEEP_FAILURES=1`.
- The corresponding implementation regressions are fixed in the working tree and each has a committed generator plus deterministic BAM/BED fixture. The regression command outputs match the local pinned oracle byte-for-byte except where the case explicitly uses the declared numeric-table comparator.
- Locked workspace tests, strict Clippy, and a release build are required after every gate update. These checks are evidence for the restricted regression profile only; they do not close the full all-command, API, standalone, real-data, or publication gates below.

The AUD table remains an historical record of findings at the audited revision. Closure status must be migrated into the versioned divergence ledger after independent review; a passing local regression is not by itself a compatibility or scientific-validity claim.

## 3. Validation structure and independence

Use complementary sources of evidence for every scientific metric:

1. **Analytical expectations:** tiny hand-auditable input, an explicit formula/filter table, and manually or independently calculated values.
2. **Pinned upstream execution:** unmodified real CLI and, where useful, public API calls, with complete process and artifact capture.
3. **Independent implementations or readers:** a small declarative enumerator, independently implemented formula, or a suitable third-party tool with its differing semantics accounted for.
4. **Properties and biological controls:** metamorphic relations with stated assumptions, controlled perturbations, and held-out real data.

Expectations must not call production helpers to calculate the answer. A round trip through the same writer and reader is useful but is not independent format validation. Pysam and samtools share htslib, so count them as a common implementation lineage, not two independent votes.

Separate implementation review from approval of scientific expectations, changed goldens, and tolerances. One person can prepare both when resources require it, but release acceptance needs another reviewer to inspect the actual evidence. Reviewer independence is a release requirement, not an instruction to launch agents during this planning task.

Track each requirement as:

`requirement -> command/profile/options -> case IDs -> independent expectation -> comparator -> run evidence -> reviewer -> allowed claim`

Track semantic branches explicitly: each relevant filtering decision, coordinate convention, denominator, sampling unit, error mode, and artifact family needs a test. Line/branch coverage is a diagnostic for omissions; neither a high percentage nor a large test count proves scientific validity. Document unreachable branches and remaining uncovered behavior rather than hiding them with exclusions.

## 4. Freeze a reproducible oracle and contract

Before expanding goldens:

- Obtain the complete upstream revision specified in [upstream.lock](compatibility/upstream.lock), including package metadata and auxiliary files. Verify the checkout against that revision, not just an existing directory name.
- Record a canonical file-content manifest and hashes. The existing hash of a host-generated tar stream is sensitive to archive ordering and metadata; specify an exact reproducible archive recipe or prefer sorted relative-path/content hashes. Store source archives with persistent access.
- Cross-check the published release separately. A GitHub commit declaring version `5.0.5` is not evidence of PyPI `5.0.5` compatibility. Keep distinct oracle profiles if contents differ.
- Lock Python, pysam **and bundled htslib/samtools**, bx-python, NumPy, pandas, pyBigWig/libBigWig, logomaker, matplotlib, R, R packages, fonts, `htseq-count`, and UCSC helpers as applicable. Capture build options, resolved executable paths, package hashes, OS/architecture, locale, timezone, and image digest.
- Provide a scripted environment rebuild from a clean checkout and exercise all 33 commands with real inputs, including plotting/conversion. Import/help checks do not establish a functioning oracle.
- Run a reference case repeatedly to establish deterministic fields and nondeterminism before creating a golden. Record random state and environment such as `PYTHONHASHSEED` where relevant.
- Hash before and after execution to detect oracle/fixture mutation. Never patch the primary reference to make it pass. A necessary patch or seeded instrumentation creates a separately named profile with a diff and justification.
- Freeze the actual API and CLI surfaces from source and execution: defaults, aliases, negative values, optional branches, sorting/index requirements, files, stream placement, error exits, warning behavior, and overwriting. Inspect package assets and unreferenced public helpers as well as CLI call paths.

Do not extend compatibility claims when dependencies or upstream behavior change without rerunning the corresponding matrix. A container provides an environment boundary, not independent assurance that the source or assumptions are correct.

## 5. Make the verification harness capable of failing correctly

### 5.1 Required execution contract

Each case must declare input hashes, oracle/candidate profiles, expected outcomes for each side, exact argument vectors, environment, expected artifacts, comparators, and resource bounds.

- Build the candidate from the claimed revision with `--locked`. Capture toolchain, target, build flags, lockfile hash, source-tree hash, and binary hashes. Record dirty-tree content if deliberately tested; release qualification uses a clean, immutable revision.
- Run each side in its own fresh working directory, with isolated temporary/config/cache directories and read-only fixture/oracle inputs where feasible. Relative outputs and side files such as `log.txt` must be captured there.
- Use argument vectors, controlled environment, bounded CPU/memory/disk/wall time, and bounded or streamed log capture. A timeout must terminate the whole process group. Record signal, timeout, OOM, dependency failure, and ordinary exit separately.
- A positive case requires its declared successful exit **and** all required artifacts/fields. Both programs failing is never a positive pass. A negative case requires its specified failure behavior; it does not count toward valid-workload coverage.
- Capture stdout and stderr separately as bytes. Decode under the declared contract only for semantic comparison. Do not concatenate streams and accidentally hide wrong routing or lose framing.
- Compare the complete output tree to an allowlist: required, optional under an explicit condition, and forbidden artifacts. Check file types, sizes, nonempty schemas where required, and side effects outside the sandbox. Do not accept a stale artifact from an earlier run.
- Validate fixture hashes even when files already exist. Missing fixture, golden, schema, comparator, dependency, or selected case is an infrastructure failure, never success or automatic regeneration of an approved expectation.
- Enforce strict result states: `pass`, `fail`, `expected-upstream-failure`, `unsupported`, `infrastructure-error`, `not-run`. Only actual successful tests contribute to the relevant pass count. Required `not-run` cases block release.
- Every expected failure links to a specific issue/ledger entry and a narrow expected error. An unexpected pass forces review; broad expected-failure decorators cannot conceal unrelated failures.
- Preserve failed runs with inputs or content-addressed references, commands, raw streams, output tree, diagnostics, and comparator versions. Preserve concise successful evidence and immutable input/output hashes. Apply an explicit storage-retention policy after archiving.
- Cache only by source, binary, oracle environment, inputs, case definition, comparator, tolerance profile, and runner hashes. A stale cache cannot satisfy a changed requirement.

Different arguments on the two sides are allowed only as an explicit capability-scoped case. For example, a reference `--skip-plot` versus a candidate with no plotting option tests selected data outputs, not invocation or complete output compatibility.

### 5.2 Case-manifest design

Implement a validated schema and version it. The following is an illustrative **future** record, not a file understood by today's runner:

```yaml
schema_version: 1
id: fpkm_paired_exon_denominator
command: FPKM_count.py
profile: github-59a24c5
kind: positive
requirements: [AUD-04, paired-exon-membership, fpkm-denominator]
inputs: [fixture:paired-exon-membership-v1, fixture:single-exon-bed-v1]
# Each fixture reference resolves to a reviewed manifest with real hashes.
oracle_argv: ["-i", "{bam}", "-r", "{bed}", "-o", "{work}/out", "-e"]
candidate_argv: ["-i", "{bam}", "-r", "{bed}", "-o", "{work}/out", "-e"]
expected_exit: {oracle: 0, candidate: 0}
artifacts:
  - path: out.FPKM.xls
    required: true
    comparator: fpkm-table-v1
    expected_rows: 1
    unique_key: [chrom, st, end, accession, gene_strand]
    expectation: analytical:paired-exon-membership-v1
    numeric_policy: exact-for-this-fixture
streams: {stdout: empty, stderr: fpkm-progress-v1}
normalizations: []
limits: {wall_seconds: 30, memory_mib: 1024, output_mib: 20}
```

Resource limits above are a starting allocation for that tiny case, not defaults for whole-genome jobs. A valid manifest must resolve every ID; incomplete examples and placeholders must fail validation.

### 5.3 Test the tests

Create comparator/runner tests before trusting expanded coverage. Deliberately introduce each error below and require detection:

- Equal nonzero exits; two empty outputs; header-only results; missing mandatory categories; duplicated labels/rows; unknown or misspelled case selectors; zero selected/executed cases.
- `1e-3` versus `1e+3`; trailing junk after a number; changed sign; one count changed by one; NaN versus finite/zero; positive versus negative infinity; signed zero where observable; near-zero errors that a relative tolerance would hide.
- Missing or extra rows, swapped samples/columns, repeated identifiers, changed order when ordered, truncated last row, missing newline, wrong delimiter, wrong stdout/stderr channel.
- A one-base coordinate shift, dropped alignment/tag/mate, changed flag, damaged header, wrong reference dictionary, missing BGZF footer, valid-looking but unqueryable index, malformed compressed body.
- BigWig missingness replaced by zero; lost interval; boundary shift; altered summary/zoom result; all-negative signal incorrectly clamped; WIG/bedGraph coordinate conversion errors.
- Missing plot/conversion output despite exit 0; blank rendered page; swapped axes or sample labels; omitted series; a legend inconsistent with the underlying values.
- Stale binary, wrong source revision, changed fixture content, mismatched oracle dependency, stale cached pass, failure outside the selected output tree, child process surviving timeout.

Add algorithm mutation tests for off-by-one comparisons, reversed strand routing, removed flag filters, changed denominators, merged duplicate points, read-versus-fragment substitutions, and incorrect sampling populations. Every preselected scientifically meaningful mutation must be caught or receive an independently reviewed equivalence explanation. Archive surviving mutations; do not improve a score by deleting difficult mutants.

## 6. Fixture system and boundary coverage

### 6.1 Fixture construction

Commit small human-readable SAM/BED/FASTA/FASTQ/track descriptions with generation seeds and hashes of generated binaries. Keep the biological event model separate from its binary encoding. Generate BAMs using an independent tool and validate them by a second implementation where practical. Check sorting, index queries, QNAME grouping, mate consistency, CIGAR/query lengths, and contig lengths; a header declaration is not validation.

Keep valid edge inputs separate from deliberately malformed inputs. Cover valid unusual BAM states, such as positioned unmapped mates, without assuming every surprising flag combination is invalid. Do not replace existing fixtures silently; version corrections and retain a reproducer if a prior fixture exposed a defect.

Use a layered design: single-factor truth tables, targeted higher-order interactions, combinatorial pairwise coverage for the remaining option space, and randomized generated inputs with shrinking. Pairwise coverage alone is insufficient for known interactions such as mate overlap × base quality × pair flags.

### 6.2 Shared fixture families

| Family | Required cases |
|---|---|
| Alignment flags | Single/paired, each orientation, proper/improper pairs, first/last/neither segment, unmapped/positioned-unmapped/mate-unmapped, secondary, supplementary, QC failure, duplicate; meaningful combinations and supplementary records sharing QNAME |
| MAPQ and quality | MAPQ 0, cutoff−1/cutoff/cutoff+1, 254/255; distinguish unavailable MAPQ from an ordinary score. Base qualities 0, 12/13/14, valid maxima, missing quality, and malformed sentinel mixtures; qualities varying by cycle and strand |
| CIGAR | Each of `M I D N S H P = X`, compositions at both ends, consecutive indels, multiple junctions, long skips, no aligned bases, long CIGAR representation, consistency with sequence length; exact exon/transcript/chromosome boundaries |
| Pairing | Overlapping and nonoverlapping mates, equal/unequal qualities, agreeing/disagreeing overlap bases, unequal lengths, mates on different chromosomes, orphaned/missing mates, distant mates, negative/zero template lengths, QNAME collisions across read groups |
| Coordinates | Position 0 in internal/annotation coordinates; SAM position 1; half-open end exclusion; touching/nonoverlapping intervals; one-base intervals; chromosome end; large legal coordinates; negative/reversed/overflowing values as invalid cases |
| BAM access | Sorted, query-name sorted, unsorted, false sort declaration, equal-coordinate ties; `.bam.bai` and `.bai`, CSI if advertised; missing/stale/truncated/wrong-BAM indexes; absent contig; truncated file, corrupt blocks, EOF/footer errors; SAM where promised |
| Tags | Absent/present/wrong-type `MD`, `NM`, `NH`, `IH`, `H0/H1/H2`, read-group and cell/UMI tags; arrays and integer subtypes for format preservation; duplicated malformed tags; optional strings with delimiters |
| BED/annotation | BED6 versus BED12 where applicable; both strands; one/many exons; CDS/UTR boundaries; overlapping isoforms and genes; duplicated names/models; noncoding/zero CDS; trailing commas; inconsistent counts; blank/comment/track/browser lines; contig naming/case mismatches |
| Sequence text | A/C/G/T/N and IUPAC ambiguity; lowercase; unequal lengths; multiline FASTA; FASTQ qualities made only of A/C/G/T/N; empty/truncated records; CRLF; final newline absent; `.gz`/`.bz2` if advertised; wrong extensions/magic; sequence and quality length mismatch |
| Tracks | BigWig empty/missing/zero/negative/positive regions, chromosome subsets, unequal dictionaries, adjacent equal runs, floating-point extremes, zoom levels, exact and approximate summaries; WIG fixed/variable steps and spans where accepted; bedGraph chunk boundaries |
| Scale and sparsity | Empty, one observation, all filtered, constant coverage, entirely missing contig, sparse islands, high depth at 7999/8000/8001 and well beyond, long reads/transcripts, many isoforms/barcodes, high-output-volume cases |
| Parameters and files | Defaults and every option, zero/negative/huge/nonfinite values, incompatible options, invalid enums, paths with spaces/quotes/Unicode, duplicate basenames, unwritable outputs, existing outputs, output equal to input or an alias, interrupted writes |

Pin coordinate and encoding interpretations to the [SAM/BAM specifications](https://samtools.github.io/hts-specs/) and [UCSC format definitions](https://genome.ucsc.edu/FAQ/FAQformat.html). Upstream deviations from those definitions need separate compatibility expectations.

The [pysam pileup documentation](https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignmentFile.pileup) exposes depth, filtering, overlap, quality, and reference-dependent behavior. Verify which settings actually apply to the pinned library, stepper, inputs, and call site; a documentation default or source comment alone is insufficient. In particular, exercise depth caps, orphans, overlap quality adjustments, deletion/refskip entries, and BAQ with/without a supplied reference where relevant.

## 7. All-command acceptance matrix

Every row additionally requires CLI/default/error tests, independent expectations, fixture validation, repeated deterministic execution, and all advertised artifact families. Expand each row into named manifest cases; these are minimum obligations, not assertions of complete option coverage. Compare generated scripts and rendered images separately. The existing inventory is a discovery aid, not authoritative golden data.

### 7.1 Alignment and file operations

| Command | Scientific/semantic obligations | Artifact and branch obligations |
|---|---|---|
| `bam_stat.py` | Hand-count every reported category and denominator; flag precedence; supplementary versus secondary; spliced reads; cross-contig mates; exact MAPQ threshold and 255 | Full stdout/stderr contract, empty/all-filtered files, BAM and promised SAM support |
| `bam2fq.py` | Reverse complement and reversed qualities; ambiguity characters; missing sequence/quality; suffixes; orphan records; secondary/supplementary inclusion | Single/paired output order and record multiplicity, exact decompressed FASTQ, gzip mode, error on failed final write |
| `divide_bam.py` | All records with a QNAME remain together; every eligible record appears exactly once; unmapped filtering; seeded assignment versus stochastic profile | Subset counts 1/2/non-power-of-two and invalid 0; decoded BAMs, empty subsets, index option, repeatability, write/finalization failures |
| `split_bam.py` | Read-start/mate-start classification; exon boundaries; no-overlap and junk categories; contig case and mate chromosome behavior | Three BAMs conserve the eligible multiset; report consistency; overwrite protection and index output, including preexisting index files |
| `split_paired_bam.py` | Read-1/read-2/unmapped partition; exact historical flag transformation and mate fields; unpaired/abnormal flags | Header, tags, sequence, quality, ordering, counts, index/overwrite/verbose branches, independently decodable complete BAM files |

### 7.2 Read-level profiles

| Command | Scientific/semantic obligations | Artifact and branch obligations |
|---|---|---|
| `read_GC.py` | 0/25/50/75/100% GC, ambiguity denominator, mixed lengths, filters and rounding boundary cases | GC table, script, rendered plot, empty/zero-length sequence handling, skip-plot behavior |
| `read_NVC.py` | Cycle-by-base counts, reverse orientation, ambiguous bases, last-record-length behavior and filtered last record | Full table schema, `--nx` behavior, script/image artifacts, variable lengths and order sensitivity |
| `read_quality.py` | Per-cycle histograms with nonuniform qualities; actual MAPQ-only filter; missing qualities; reverse orientation; all-zero qualities | Full R vectors and rendered summaries, reduction factor including invalid zero, min/max score axes, last-read-length behavior |
| `read_duplication.py` | Sequence versus positional duplicate definitions; spliced/clipped CIGARs; QNAME is not sequence identity; reverse reads; histograms conserve multiplicities | Both tables, cutoff's actual scope, script/plot, paired and multimapping cases |
| `clipping_profile.py` | Per-cycle soft-clipping counts, both ends, reverse reads, single/paired routing, hard clipping distinctions and mixed lengths | Sequencing-mode enum, table and denominators, script/image, no-clipping and all-filtered cases |
| `insertion_profile.py` | Insertions at different cycles, multiple runs, reverse reads, read-1/read-2 denominators, combinations with clipping/deletion | Single/paired tables, script/image, zero insertions, variable read lengths |
| `deletion_profile.py` | Deletion location versus length; multiple events; read/query/reference coordinates; read-length filter and sampled qualifying-read count | Tables, script/image, cap boundary, `=`/`X` and clipped reads, empty qualifying set |
| `mismatch_profile.py` | `MD`/`NM` missing and wrong type, adjacent mismatches, MD deletions, zero NM, ambiguity, reverse complement, CIGAR consistency | Every substitution category and cycle, sample cap and natural exhaustion output, script/image, invalid tags |
| `read_hexamer.py` | Hand-count all windows on short sequences; 4096-bin completeness, N exclusion and denominator, record boundaries; FASTQ quality strings resembling DNA | FASTA/FASTQ and compressed claims, reference genome/gene branches, input lists, duplicate display names, missing-input policy, stdout/file identity |

### 7.3 Annotation and fragment metrics

| Command | Scientific/semantic obligations | Artifact and branch obligations |
|---|---|---|
| `infer_experiment.py` | Known forward/reverse/unstranded mixtures; all pairing orientations; overlapping opposite-strand genes; unassigned denominator; early sampling bias | Complete reported fractions/layout, sample-size boundaries, empty or mixed-layout input, annotation mismatch |
| `read_distribution.py` | CDS/UTR/intron/intergenic priority and ambiguity; exon-block midpoint versus read counts; overlapping annotations and region lengths | Every category including zeros, total reads/tags/assigned counts, tags-per-kb denominator, boundary fixtures |
| `inner_distance.py` | Negative/zero/positive distance, intron removal, mate chromosome and orientation, first-transcript-per-chromosome legacy omission | Raw distance rows, histogram edge inclusion, bounds/step/sample cap, unknown chromosome text, script/image |
| `RNA_fragment_size.py` | Transcript-boundary `=`/`X` reproducer; CIGAR reference span for fetch; query length in fragment formula; mean/median/population SD; fragment cutoff | Header and numeric types in stdout/file, invalid BED hard failures, index requirements, cross-contig and unequal-length pairs |
| `junction_annotation.py` | Known/both-end-novel/one-end-novel categories; event versus unique counts; multiple junctions/read; intron cutoff; annotation contig filtering | Table, R script, plot, BED12 and Interact tracks; `--skip-bed`/`--skip-interact`, coordinate/strand/color validity and rounding |

### 7.4 Coverage and track operations

| Command | Scientific/semantic obligations | Artifact and branch obligations |
|---|---|---|
| `bam2wig.py` | Per-base and integrated coverage; duplicate/secondary filters; tag-based normalization filter versus MAPQ-based signal filter; strand rules and reverse sign | Exact WIG coordinates/order/casing; both stranded files; normalized sums with rounding budget; BigWig conversion and semantic output, missing converter behavior |
| `geneBody_coverage.py` | Pileup semantics including AUD-06/07; percentile interpolation/ties; negative-strand reversal; transcript identity; empty/constant coverage; raw skewness | BAM/list/directory input and index variants; table, scripts, curves/heatmap for 1/2/3+ samples; sample-name collisions; no NaN-sort panic |
| `geneBody_coverage2.py` | BigWig missingness, float precision, strand, percentile positions; first-exon-shorter-than-100 legacy exclusion | Data, script, every advertised plot format, absent contigs, all transcripts excluded; cross-check BAM/BigWig siblings only on a proven common domain |
| `normalize_bigwig.py` | Genome versus merged-exon denominator; target sum; zeros/missing/negative values; contig mismatch and out-of-range exons | WIG/bedGraph, chunk size 1 and boundary-spanning equal runs, rounding, overwrite handling; semantic chunk invariance separately from text segmentation |
| `overlay_bigwig.py` | Each action and operands, missing chromosome versus missing interval, negatives/NaN/Inf, mismatched sizes; legacy Division failure | Complete WIG, chunk effects, overwrite, per-base independent arithmetic; Division classified as upstream failure, not successful functional coverage |

### 7.5 Expression, sampling, and integrity

| Command | Scientific/semantic obligations | Artifact and branch obligations |
|---|---|---|
| `FPKM_count.py` | AUD-04 and AUD-08; all four mate/exon-membership states, weights 0/0.5/1, unpaired and strand rules, total/exonic denominators, duplicate and spliced reads | Full FPKM table and row identity, both normalization modes, zero transcript length/denominator, rejected invalid weights/rules, case-sensitive fetch |
| `FPKM-UQ.py` | Hand-calculated counts/lengths/FPKM/UQ/log2; protein-coding denominator; zeros, ties, quantile method, missing/duplicate gene IDs | Actual counting from BAM/GTF, count table plus normalized table, fixed helper arguments, print-command branch, helper failure; standalone implementation is a separate unmet gate |
| `RPKM_saturation.py` | Exon-block midpoint sampling unit, duplicates, strand populations versus denominator, percentile interval overlap/gaps, exact 100% endpoint | Raw-count/RPKM tables and row order, quartile/error-summary logic, seed instrumentation, script/image, zero sample/length and invalid percentile combinations |
| `junction_saturation.py` | Per-read intron population, cumulative sampling, known recurrence threshold versus novel behavior, annotation chromosome filtering, terminal 100% | All series/axes in script and plot, deterministic supplied permutations, stochastic equivalence, start/step/ceiling interactions and empty population |
| `tin.py` | Entropy formula, strict distinct-start cutoff, duplicate filtering asymmetry, mate overlap/depth cap, sampled duplicates, background truncation and subtraction | Per-transcript and mean/median/SD summaries, directory/list/single BAM forms, repeated basenames, output-directory behavior, negative/zero coverage controls |

### 7.6 Single-cell tools

| Command | Scientific/semantic obligations | Artifact and branch obligations |
|---|---|---|
| `sc_bamStat.py` | Actual configurable tag/RE/confident-mapping rules; mitochondrial naming; splicing categories; read, barcode, and UMI denominators | Complete report including zero categories, absent/wrong-type tags, empty denominators, index validity versus mere sidecar presence |
| `sc_editMatrix.py` | Raw/corrected CB and UMI pairs; identical/missing/different states; unequal lengths; position/base substitution counts; suffixes | Both frequency tables and edit matrices, tie ordering, CSV types, sample limits, heatmap script/image/log2/style branches and dependency failures |
| `sc_seqQual.py` | Per-cycle Phred counts and fractions; variable lengths, low-quality tails, incomplete FASTQ; limit applies at the correct stage | Count and percentage matrices, column sums, compression, all heatmap branches and image formats, memory under small limits on large files |
| `sc_seqLogo.py` | Per-cycle base matrix, first-seen column order, exclude-N/limit interaction, lowercase/ambiguity, multiline FASTA legacy behavior | Count matrix **and** logo output, highlights, stack/style options and supported formats, compression, invalid/empty sequences; current always-failing command remains incomplete |

The matrix has 33 distinct commands. CI must verify that its machine-readable successor remains in one-to-one correspondence with the declared command inventory and installed aliases. No omitted row is implicitly out of scope.

## 8. Independent mathematics, properties, and comparators

### 8.1 Analytical controls and metamorphic tests

For each property, record its applicability and exclusions in the case. Do not force an attractive identity onto an algorithm with different semantics.

- **Counts and partitions:** eligible input record multisets equal the union of output partitions, with exact multiplicities and declared flag transformations. QNAME grouping alone cannot detect a dropped supplementary record. Conservation applies to the eligible set after specified exclusions.
- **Intervals:** on a small artificial chromosome, enumerate covered positions directly and compare union/intersection/subtraction and overlap lengths. Check commutativity/idempotence and `|A union B| + |A intersection B| = |A| + |B|`. Test raw-multiplicity queries separately from set queries.
- **Coverage:** for a controlled uncapped dataset without mate-overlap adjustment, coverage mass equals the sum of eligible aligned-base contributions. Duplication doubles raw counts under additive filters; it need not double normalized metrics, capped pileups, unique counts, or TIN.
- **FPKM:** independently calculate `C * 10^9 / (N * L)` using exact counts and rational lengths where applicable. Verify `C`, `N`, and `L` separately. Duplicating all eligible records preserves FPKM when numerator and denominator scale equally; changing the exon-only population changes the denominator by its own contract.
- **TIN:** for `n` distinct sampled positions, uniform positive coverage at all positions gives 100; one positive position gives `100/n`; uniform support at `k` positions gives `100*k/n`. Use `H = -sum(p*ln(p))`, `TIN = 100*exp(H)/n`, excluding zero terms. Keep duplicate sampling positions, background adjustment, and eligibility outside these simplified controls. Multiplying positive coverage uniformly preserves TIN when no filter/cap/background changes.
- **GC/NVC/quality:** hand-enumerate short reads with asymmetric bases and qualities; reverse-complementing/reorienting inputs transforms cycles exactly as specified. Histograms sum to the eligible observations for each position, accounting for variable lengths.
- **Track operations:** `Add(x, 0) = x`, `Subtract(x, x) = 0`, and multiplication/scaling identities on finite supported values. Missingness and sparse serialization are compared separately. Normalized mass must meet a derived serialization-error bound, not necessarily equal the target after decimal rounding.
- **Saturation:** test a fixed population and explicit permutation so each cumulative subset is independently enumerable. Monotonicity is only required for metrics/schedules that actually accumulate unique discoveries; RPKM values need not increase monotonically. A true full-population endpoint must agree independently of shuffle where the algorithm covers that population exactly.
- **Renaming/permutation:** bijective contig renaming across all relevant inputs preserves metrics on the supported case-sensitive domain. Input permutation is tested only for order-invariant algorithms; last-record-length bugs, first-transcript omission, finite sampling, and depth caps need explicit order-sensitive expectations.
- **Paired cross-tool checks:** compute BAM-to-track-to-gene-body relationships only after aligning filters, mate handling, annotations, sampling, missingness, and rounding. Otherwise disagreement between sibling commands may be expected and must not be erased.

For independent TIN method context, use the [original TIN paper](https://link.springer.com/article/10.1186/s12859-016-0922-z). Its biological results do not validate this port or establish that every preserved legacy behavior is scientifically desirable.

### 8.2 Numeric and textual acceptance

Define a versioned policy **per metric and artifact before examining candidate differences**. No global `1e-6` or correlation-only pass rule.

| Value/artifact | Comparison policy |
|---|---|
| Integer counts, coordinates, flags, identifiers, category membership | Exact, including zeros, multiplicity, row keys, and denominators. Never convert large integer counters to floating point for comparison. |
| Deterministic text | Byte comparison where claimed; only named, anchored path/timestamp/version fields may be normalized. Do not sort all rows, strip arbitrary whitespace, drop warnings, or discard headers globally. |
| Floating-point metric | Predeclared absolute/relative budget with units and derivation: `abs(actual-expected) <= atol + rtol*abs(expected)`. Exact known analytical values stay exact when representable. Near-zero and large-magnitude controls must prevent permissive tolerances. |
| Exceptional numbers | Declare when NaN, missing, infinity, and signed zero are legitimate. Match missing/NaN masks and infinity signs explicitly; neither NaN comparisons nor replacement by zero may create a pass. Unexpected nonfinite output fails. |
| Decimal serialization | For a specified rounding-to-nearest format with `d` decimal places, account for at most `0.5*10^-d` per rounded value, with the specified tie rule. Propagate by interval width when comparing integrated track mass. Exact text compatibility still requires the exact serialization. |
| BigWig numeric values | Compare decoded stored precision and missing intervals independently from downstream f64 arithmetic. Do not demand impossible extra precision or tolerate changed coordinates because values are close. |
| Quantiles/statistics | Explicit interpolation method, ties, sample versus population SD, small-N behavior, sort order, and accumulation order. Test each distinct helper rather than assuming one “percentile” definition. |

Derive bounds from independent high-precision/rational calculations, accumulation error, stored precision, and scientific impact. A tolerance must be small enough to catch the minimum scientifically relevant error, including a one-read change at important boundaries. If numerical noise and the intended decision margin cannot be separated, mark the result inconclusive or narrow the supported domain.

Report mismatch count, missing/extra identities, maximum absolute/relative discrepancy with its location, distribution of discrepancies, and decision changes. High correlation, low mean error, or agreement after aggregation cannot hide an incorrect subset of transcripts/samples. Golden or tolerance changes require a reason and independent review, never “update snapshots until green.”

### 8.3 Artifact-specific comparators

- **BAM:** fully decode headers and all records using an independent implementation. Compare reference dictionaries/order, read names, coordinates, mate fields, CIGAR, sequence, quality, flags, and typed tags with an explicit policy for representation-only differences. Preserve sequence order when promised; otherwise compare full record multisets including duplicates. Validate EOF, actual sorting, index construction, and randomized boundary region queries. Merely running a quick header/footer check is insufficient.
- **BigWig/WIG/bedGraph:** compare chromosome lengths, interval boundaries, span-weighted signal, zeros versus absent data, strand sign, and values. For large files use a complete interval-union comparison without expanding the whole genome. Validate consumer summary queries and zoom behavior separately. [pyBigWig documents](https://github.com/deeptools/pyBigWig#compute-summary-information-on-a-range) that summary queries may use approximations; exact values and approximate summaries are different contracts.
- **FASTQ/FASTA:** validate structure, read names/suffixes, sequence/quality length, record order, pairing, and decompressed content. Do not mistake gzip-byte differences for sequence differences. Test both encoding and parser errors.
- **Tables/matrices:** schema, key uniqueness, row/column order, types, all category cells, index names, missingness, and denominators. Check cell-level values before derived totals. CSV quoting and identifiers containing delimiters need dedicated cases.
- **Plots:** validate source data, axes, units, strand direction, sample identity, legend, limits, normalization, and all series before rendering. Render in a pinned environment; compare structure and images with justified renderer-specific tolerances, and inspect representative publication figures manually. A visually similar plot with wrong labels/data fails. New native rendering may justify visual differences but not scientific differences; scripts remain separate artifacts where promised.
- **Errors:** expected exit/error class, streams, partial artifacts, and cleanup. Error wording may have a scoped normalization, but the meaning, failing condition, and output safety remain testable.

## 9. Randomness and statistical equivalence

Inventory randomness in `divide_bam`, `junction_saturation`, `RPKM_saturation`, and any public API helper. First establish what is sampled: query-name groups, alignments, fragments, exon blocks, or junction observations. Their probability distributions differ.

For an exact seeded-compatibility claim, reproduce the PRNG, seed interpretation, integer sampling/rejection method, shuffle algorithm, traversal order, and draw schedule. An implementation that uses `StdRng` cannot infer Python seed compatibility from repeated agreement with itself. Pin the RNG implementation/version for reproducibility across dependency updates.

For unseeded upstream behavior:

1. Test counting/subsampling logic deterministically with explicit permutations, including adversarial ordering and awkward percentile schedules.
2. Use a separately identified instrumented oracle for matched-seed or matched-permutation investigations; preserve unmodified black-box runs as additional evidence.
3. Define distributional equivalence endpoints, scientifically meaningful margins, error control, and sample-size/power analysis before comparing results. Pilot at least 30 independent seeds, then select a fixed validation seed count to achieve the preregistered power target, normally at least 90%. Thirty seeds alone is not an acceptance criterion.
4. Evaluate means, variance, tails, allocation bias, rare features, and whole-curve behavior where relevant. Account for dependence between percentile points and multiplicity across endpoints; per-point nonsignificance does not establish curve equivalence.
5. Use an equivalence procedure with predeclared margins and matching confidence levels. Failure to detect a difference is inconclusive, not equivalent. Report seeds and all runs, including failures; do not keep sampling until a favorable p-value appears.

Exact invariants still apply on every seed: no dropped/duplicated eligible records, QNAME groups not split, valid sampling counts, and independently checked endpoints where applicable. Seed replicates estimate algorithmic randomness, not biological replication.

## 10. Robustness, fuzzing, and failure injection

Prioritize functions whose failure can silently alter scientific data:

- Fuzz BED parsing, CIGAR traversal, optional tags, MD parsing, FASTA/FASTQ parsing, strand rules, matrix serialization, numeric formatting, and CLI parameter validation. Use structured generators for valid inputs and byte mutations for malformed inputs.
- Assert termination, bounded resources, no panic/overflow, valid output or an explicit error, and invariants where their premises hold. Use the locked upstream only within resource limits; an upstream hang is not a behavior the candidate must reproduce.
- Test debug and release builds: unchecked casts and integer operations can behave differently. Exercise counts beyond 32-bit ranges with generated summaries, not necessarily billions of input reads. Include very large coordinate additions, signed-to-unsigned conversion, zero step/chunk size, and nonfinite numeric arguments.
- Inject short writes, write errors on final flush/BGZF finalization, full disk, permission changes, truncated input midstream, broken pipes, and child-process failures. Do not rely on `Drop` to communicate write failures. A success exit must not leave an unreadable or incomplete required artifact.
- Protect existing outputs and inputs, including symlink/hardlink aliases. Check multi-output partial failure, cancellation, concurrent output-prefix collisions, and rerun behavior. Define recoverable partial-output handling explicitly.
- Test path and sample-label escaping in generated R, CSV, and shell/subprocess arguments. A legitimate quote in a path must not change execution or invalidate output.
- Test large files with small sample limits: establish whether parsing should stop early and whether memory is proportional to the declared workload. Evaluate decompression/corrupt-record handling after the sampling stop according to the contract.

Run a bounded fuzz campaign before each release, record duration/seeds/corpus/coverage, minimize every failure, and commit a deterministic regression case. The release criterion is no unresolved relevant failures, not merely “fuzzing ran for N hours.”

## 11. Real data and biological validation

### 11.1 Dataset selection before measuring results

Create `datasets/manifest.yaml` with accession and exact file IDs, resolvable URLs, source/release dates, checksums, usage/redistribution terms, biological metadata, library protocol, read layout/length, genome/annotation versions, aligner/options, preprocessing, and reproducible subset recipes. Archive accessible evidence; controlled-access data cannot be the sole support for a publicly reproducible claim.

Minimum coverage domains, with multiple independent libraries where a biological population claim is intended:

| Domain | What it challenges |
|---|---|
| Bulk paired-end stranded RNA-seq | Orientation, fragments, junctions, overlapping mates, exon membership |
| Bulk unstranded and single-end RNA-seq | Layout/strand inference and different denominators |
| Different preparation methods and read lengths | Annotation sensitivity, variable-length handling, 5′/3′ coverage patterns |
| Experimentally degraded RNA or a measured integrity series | TIN and gene-body response against known sample conditions |
| Reference mixtures/spike-ins with known controls | Count/normalization response and independent expression controls |
| Droplet single-cell BAM and barcode/UMI FASTQ | Real tag conventions, correction matrices, huge barcode cardinality |
| A second organism/reference structure | Contig assumptions, chromosome naming, transcript structure |
| Deep/sparse/high-duplication and multimapping samples | Pileup caps, sampling, performance and memory extremes |

Candidate sources to select and pin include the [SEQC reference-RNA study](https://www.nature.com/articles/nbt.2957), datasets described by the [TIN study](https://link.springer.com/article/10.1186/s12859-016-0922-z), and the [10x Genomics public dataset collection](https://www.10xgenomics.com/datasets). These are selection leads, **not downloaded or validated project datasets**. Verify actual availability, protocol, reference, file identity, and permitted use before admitting a dataset.

Split development and held-out validation at the original sample/study level before tuning algorithms or margins. Crops, downsampled libraries, and differently aligned copies of the same reads belong to the same partition. Once a held-out case is inspected to guide a fix, record that exposure and use fresh independent confirmation for claims requiring untouched validation.

Aligner/reference dependence is part of the result. Include representative different alignment encodings and tag policies; do not claim robustness across aligners from one aligner's BAMs. Do not fabricate metadata or treat a known library label as perfect truth when a mixed library is possible.

### 11.2 Scientific endpoints and downstream impact

Before running held-out validation, define which conclusions the outputs will support and what change would matter scientifically. Examples:

- Known library orientation: correct dominant strand/layout and bounded unassigned fraction under the specified annotation; assess sampling uncertainty.
- Controlled synthetic degradation or coverage skew: response of TIN and gene-body summaries under analytically defined perturbations. For real degradation series, estimate association with measured quality while considering preparation, depth, and transcript-length effects. Do not require perfect monotonicity in heterogeneous biological samples.
- Mixture/spike-in controls: expected changes in counts and normalization where the command's estimand makes that expectation valid. FPKM is not automatically a validated substitute for every expression-analysis method.
- Junction/count results: exact differences on matched BAM/annotation inputs, then impact on a preregistered downstream decision or report. Validate metadata, skipped transcripts, and sample matching as carefully as the values.
- Single-cell outputs: manual tagged-read controls plus real-data summary consistency; confirm that different chemistry/tag versions are inside the advertised domain.

Report per-sample and per-feature errors, stratified by coverage, transcript length, GC, splice complexity, overlap, annotation ambiguity, and relevant quality flags. Investigate outliers, not only global correlations. When confidence intervals describe biological generalization, use sample/donor/study as the independent unit or a justified hierarchical model; millions of transcripts from one sample are not millions of independent biological replicates.

Pin at least one actual downstream workflow per advertised output family: consumption of tables by a QC report, BAM indexing/querying and downstream alignment handling, tracks queried by a viewer/analysis tool, and plots reviewed for publication. A parser smoke test alone is insufficient when the claim concerns unchanged scientific decisions.

## 12. Fair performance and scalability measurements

Performance claims start only after the exact workload passes the relevant scientific and compatibility gates. Missing functionality or an unresolved numerical discrepancy makes that workload ineligible for an equivalence-qualified speedup claim; still report its failure or unsupported status.

### 12.1 Protocol to freeze before benchmarking

- Define each command/mode, small/medium/large inputs, all required outputs, setup/indexing costs, thread limits, compression settings, and success criteria. Include every advertised command in the coverage table, even when its result is unsupported or failed.
- Measure data-only, compute-only, and complete end-to-end work as separate named experiments. Compare equal deliverables. If upstream draws a plot or writes BigWig, its end-to-end competitor must do the same.
- Match preprocessing and input files exactly, and record prebuilt indexes. Do not omit the port's in-memory indexing cost while including upstream setup. Report installation/build time separately from runtime.
- Record CPU, core count, affinity, RAM, NUMA, storage/filesystem, OS/kernel, power/frequency policy, resource limits, compiler/target features, optimization flags, library versions, executable hashes, environment, and storage/cache state.
- Measure wall time, total user/system CPU, peak concurrent memory for the **whole process tree**, I/O, output size, and throughput with a precise denominator. Summing child peak RSS is not peak concurrent memory. Account for temporary files and external R/counting/conversion processes.
- Begin with one-thread parity, explicitly controlling implicit threading in dependencies and external helpers. Add scaling only for supported modes with identical correctness. Keep default-versus-default and tuned-versus-tuned results separate.
- Use isolated hardware for publication numbers. Shared CI measurements detect gross regressions but cannot support a precise published speedup on their own.

### 12.2 Repetitions and analysis

Use separate pilot runs to choose runtime/resource budgets and repetition counts. As an initial protocol, use at least 3 untimed/pilot observations where relevant and at least 10 measured paired repetitions for a primary deterministic comparison; increase the fixed sample size when pilot variation requires it. Freeze the schedule before collecting the final measurements. Expensive cases can use a different justified preregistered design; do not silently drop them.

Randomize/interleave reference and candidate within comparable blocks. Measure warm-cache workloads separately from genuinely controlled cold-cache workloads, documenting how cache state was established. A new process is not a cold cache. Never flush shared machine caches casually.

For each command/workload, publish raw paired measurements, median times, spread, and a 95% confidence interval for a defined speedup estimator. One defensible estimator is the exponentiated mean of paired log time ratios, with resampling at the independent run/block level; if using a ratio of medians, name and implement that estimator consistently. State bootstrap/interval methodology and sample size. Do not pool unrelated workloads as independent replicates.

A positive speedup claim requires the preregistered uncertainty criterion, normally a confidence interval entirely above 1, and scientific equivalence of the measured workload. An inconclusive interval is reported as inconclusive. Avoid stopping after a favorable measurement. Apply multiplicity control or label exploratory claims when testing many speedup hypotheses.

Keep timeouts, OOMs, crashes, and incomplete outputs in the report. Do not replace timeouts with an invented exact runtime, omit slowdowns, or average only the successful subset without identifying it. A geometric aggregate must name its included workloads, weighting, missing entries, and uncertainty; per-command results remain primary.

### 12.3 Scale along the implementation's actual cost drivers

Measure independent sweeps in read count, covered genomic bases, maximum depth, transcript count, exon count/overlap density, number of contigs, distinct QNAMEs/barcodes, input length, output volume, and annotation size. Preserve expected results under each scaling operation.

In particular, test `bam2wig`'s per-base tree, TIN/gene-body/FPKM whole-BAM indexes, per-transcript scans, interval intersection/subtraction on dense models, saturation population copies, and full-text output buffers in BigWig operations. Synthetic duplication can expose a cost curve but is not a substitute for realistic biological inputs. Measure the largest supported operating point and the behavior beyond its resource limit.

After optimization, rerun affected semantic consumers, adversarial cases, and the benchmark workload outputs. Speed gained by dropped reads, changed sampling, reduced precision, weaker filters, or skipped artifacts is a changed computation and must not be presented as an equivalent optimization.

## 13. Packaging, Python API, and automation

Test installed release artifacts outside the source tree, under a fresh user configuration and without the oracle environment on `PATH`/`PYTHONPATH`.

- Run all claimed `.py` PATH aliases and native names, with correct case/hyphenation, from directories containing spaces. Record resolution to the intended executable; accidentally invoking an installed upstream script invalidates the test.
- For the standalone profile, exercise all required workloads in an image with Python, R, `htseq-count`, `wigToBigWig`, and other helper executables absent. Trace child processes and audit dynamic libraries/bundled components. Help-only checks and `ldd` alone do not establish standalone behavior.
- Test every advertised OS/architecture/ABI, archive installation, executable permissions, aliases, and CPU baseline. Unsupported platforms are explicit. An artifact built for the host CPU cannot imply portability without testing.
- For optional Python compatibility, cover all inventoried modules and public functions/classes: imports, signatures/defaults, returns/types, exceptions, state/mutation, iteration/laziness, file handles, side effects, and packaged helpers/assets. Test wheels in clean environments and supported Python versions. A `.py` executable alias does not make `python command.py` or `import qcmodule` work.
- Validate release metadata, license/attribution files, third-party dependency notices, citation metadata, source archives, and reproducibility bundle. Record dataset redistribution decisions rather than bundling unreviewed source material.

Automate the following layers. Time allocations are engineering targets, not correctness exceptions:

| Layer | Required work | Failure policy |
|---|---|---|
| Every change | Locked build/test, formatting/lint, schema/selector validation, comparator negative controls, impacted deterministic tiny differential cases | Fail closed; publish failed artifacts |
| Scheduled | Full tiny matrix, property/mutation/fuzz campaigns, medium data, nondeterminism and resource probes | File and retain failures; do not silently quarantine scientific defects |
| Release candidate | Full command/API/capability matrix, held-out validation, installed-artifact/platform checks, complete figures/conversions, benchmark qualification | Required missing/unsupported cases block the corresponding advertised claim |
| Publication snapshot | Frozen protocol, clean rebuild by a reviewer, raw results and report regeneration, final claim-to-evidence audit | Unreproducible or unexplained results excluded from claims and explicitly reported |

Do not make network availability of an upstream service a hidden source of skipped release tests. Cache verified immutable source/data artifacts under appropriate distribution terms. CI summaries must distinguish tests attempted from requirements satisfied.

## 14. Ordered implementation work and release gates

The files/directories below are proposed deliverables unless already present. These gates refine the original plan; they do not mark its historical gates as achieved.

| Gate | Work and deliverables | Acceptance evidence |
|---|---|---|
| T0 — trustworthy specification | Audit inventories, versioned oracle environment/build recipe, fixture/result schemas, machine-readable capability/requirement map, repaired status documentation | Clean oracle rebuild; every command/API item has scope, owner, and disposition; no unresolved profile identity |
| T1 — trustworthy runner | Extend `verification/run_diff.py` or replace it with the declared runner; comparator modules/tests; output-tree capture; provenance and failure preservation | All Section 5 negative controls fail for the intended reason; empty/failed/unknown cases cannot pass; deliberate output corruption is detected |
| T2 — scientific regression foundation | Versioned Appendix A fixtures plus flag/CIGAR/annotation/pileup truth tables; independent math checks; review the whole divergence ledger | AUD-04/05/06 reproduce before their fixes and pass afterward; documented overlap and other relevant gaps are fixed or explicitly excluded from restricted claims |
| T3 — complete command validation | Manifest cases for all 33 rows, full option/artifact/error coverage, semantic BAM/BigWig/FASTQ comparators, mutations and properties | No unexplained numerical mismatches; every required branch mapped; all requested outputs actually produced for a full-command claim |
| T4 — scientific qualification | Dataset manifest, development/held-out split, biological-control protocol, downstream impact report | Predeclared scientific endpoints and margins met within the advertised domain; reviewer signs off on limitations and exposed holdouts |
| T5 — equivalent performance evidence | `benchmarks/protocol.md`, runner, raw paired measurements, resource reports, reproducible analysis | Each claimed speedup is tied to passing output evidence and the measured released binary; failures/slowdowns retained |
| T6 — distributable result | Clean-install/platform tests, standalone runtime audit, optional Python API suite, CI and publication archive | Every advertised artifact/profile reproducible outside the repository; final source/build/data/report identifiers agree |

Immediate execution order:

1. Fix harness false-pass paths and lock the oracle; commit regression fixtures for the three reproduced scientific mismatches.
2. Resolve the paired denominator, fetch/CIGAR, and pileup semantics; audit their downstream consumers before broad optimization.
3. Expand the tiny fixture matrix across every command, including error/option/artifact paths; add independent controls and comparator mutation tests as each comparator lands.
4. Complete missing claimed capabilities and renderers, while building a frozen real-data panel and installed-artifact tests.
5. Run held-out scientific validation, then correctness-qualified performance measurements and release reproducibility review.

A gate can close for a clearly named restricted command/profile while other work remains, but the full-package and standalone milestones retain their original scope. Do not relabel missing work “accepted divergence” to close those gates.

Each work ticket must specify owned paths, requirement/case IDs, upstream call sites, independent expectations, affected consumers, exact acceptance commands, and evidence location. Production changes, goldens, and tolerance decisions need visible review boundaries. Every fixed defect receives a regression case and a release-note entry if it affected previous outputs.

### 14.1 Per-command publication acceptance checklist

- [ ] Baseline, dependencies, input domain, options, artifacts, and profile are frozen.
- [ ] Scientific definitions and denominators are independently reviewed.
- [ ] Deterministic fixture, boundary, negative, property, and mutation obligations pass.
- [ ] All expected artifacts and output identities are checked; no vacuous passes or stale binaries.
- [ ] Relevant shared primitives and their other consumers have been retested.
- [ ] Known divergences are quantified and either closed or excluded explicitly from the claim.
- [ ] Real-data/held-out endpoints and downstream impact meet predeclared criteria.
- [ ] Installed release artifact reproduces the result on supported environments.
- [ ] Performance claims, if any, have paired raw measurements and passing equivalent-output evidence.
- [ ] Methods, exact commands, versions/hashes, figures, limitations, and executable report-generation scripts are archived.
- [ ] An independent reviewer can rebuild and regenerate the relevant report without the developer's working directories.

No checklist item is satisfied merely because this document exists.

## 15. Loopholes that must remain under active review

Even after these gates, retain explicit residual-risk reporting for:

- Shared bugs between upstream, fixtures, and an allegedly independent oracle; reviewed handwritten controls reduce but do not eliminate this risk.
- Correct implementation of an unsuitable scientific definition, annotation, or library assumption. Compatibility cannot repair study design.
- Unseen combinations, new aligner/tag conventions, long-read inputs outside the tested domain, platform/library drift, and numerical extremes.
- Benchmark selection bias, nonrepresentative hardware, cache effects, hidden preprocessing, and unsupported command branches.
- Pseudoreplication, repeated inspection of held-out data, post hoc tolerance changes, and unreported failed analyses.
- Misidentified samples/references, duplicate gene identifiers, annotation mismatch, and wrong plot labels despite correct arithmetic.
- Differences that are small on average but change a boundary decision or a rare biologically relevant subset.

For each residual risk, record the domain affected, detection method, evidence available, and claim restriction. Add newly observed production failures to a versioned regression corpus. Revalidate affected claims after changes to scientific logic, I/O libraries, compiler options, randomness, rendering, or the oracle. A published release report remains attached to its exact artifact; later fixes do not retroactively validate earlier outputs.

## Appendix A. Reproducible audit cases to promote into permanent tests

These recipes produced the AUD-04/05/06 results at the audited revision. They are deliberately small and do not require a public dataset download. Preserve the event descriptions and expected values independently of the generator implementation. The current working tree promotes them, plus the paired-overlap and indexed-fetch cases, through `verification/fixtures/make_regression_fixtures.py`; the generated BAM/BED assets are permanent first-gate fixtures, while the oracle environment and broader matrix remain open.

Common header: SAM `VN:1.6`, `SO:coordinate`, one reference `chr1` of length 1000. Write records in ascending reference-start order and generate a valid BAM index with the pinned oracle tool. Coordinates below are **0-based starts**, not SAM text POS. All records have MAPQ 40, 20 A bases, and 20 Phred-40 qualities. Paired records set the mate reference to `chr1`; template length is the outer span with the appropriate sign.

Common BED12 line, using tabs:

```text
chr1	100	200	tx1	0	+	100	200	0	1	100,	0,
```

### A.1 Paired exon denominator

| QNAME | Start | Mate start | FLAG | CIGAR |
|---|---:|---:|---:|---|
| both | 110 | 140 | 99 | 20M |
| one | 120 | 500 | 99 | 20M |
| both | 140 | 110 | 147 | 20M |
| one | 500 | 120 | 147 | 20M |

Run each `FPKM_count` implementation with `-i fixture.bam -r model.bed -o out -e` in separate fresh directories. Both exit 0. Both count two transcript fragments and mRNA length 100. Upstream counts only the pair with both ends exonic in the normalization denominator: `N=1`, FPM `2000000.0`, FPKM `20000000.0`. The audited Rust code counts both pairs: `N=2`, FPM `1000000.0`, FPKM `10000000.0`.

This is an upstream-compatibility expectation under that command's specific counting definitions, not a general recommendation to mix different fragment populations in expression normalization. Independent scientific review must retain that distinction.

### A.2 Reference span at a transcript boundary

| QNAME | Start | Mate start | FLAG | CIGAR |
|---|---:|---:|---:|---|
| eq | 100 | 120 | 99 | 20= |
| eq | 120 | 100 | 147 | 20= |

Run `RNA_fragment_size` with `-i fixture.bam -r model.bed -n 1`. Both exit 0. Upstream's row has fragment count `1`, mean `40.0`, median `40.0`, population SD `0.0`. The audited Rust row has count `0` and zero summaries. The all-`=` reference span is incorrectly zero in the reused legacy helper, so the first read fails the half-open overlap check at the transcript start.

Also test the same reads with `20M`, `20X`, mixed `M/=/X`, and transcript boundaries shifted by one base. A first attempt with the read starting before the transcript produced zero on both sides because the command's later containment filter rejected it; this demonstrates why a reproducer must exercise the complete call path.

### A.3 Pileup depth cap

Create 8001 distinct unpaired records named `depth0` through `depth8000`, all at start 110 with FLAG 0 and `20M`. Add one unpaired `tail` record at start 150 with the same CIGAR/qualities. Run `geneBody_coverage` with `-i fixture.bam -r model.bed -o out --skip-plot`.

Both exit 0 and produce `out.geneBodyCoverage.txt`. Its maximum coverage is `8000` upstream and `8001` in the audited Rust implementation. The sparse tail avoids an all-constant profile. Repeat at 7999, 8000, 8001, and substantially higher depth, and vary read start/order and quality so a fix cannot simply clamp the final count without reproducing pileup admission/filtering semantics. Test the shared TIN consumer separately; a coverage discrepancy does not imply a fixed-direction TIN discrepancy.

### A.4 Harness negative controls

Using an isolated import and mocked process results, without changing production files:

| Probe | Audited behavior | Required behavior |
|---|---|---|
| Both sides return exit 1, empty stdout, error stderr in a positive stream case | `stream comparison PASS (0 labeled values matched)` and success | Fail the positive case |
| Compare `value: 1e-3` against `value: 1e+3` | Both parse as `{'value': '1'}` | Parse the complete numbers and fail |
| Two `value:` lines with different numbers | Last silently overwrites first | Reject unexpected duplicate field |
| Select `bam_stat_basic` and `typo_missing_case` together | Only known case runs; process returns success | Reject unknown case selector |

Promote these into runner unit tests that invoke the real selection/comparison paths. Mocked diagnostics establish a loophole, but permanent acceptance must also prove the repaired executable harness rejects deliberately broken candidate programs.
