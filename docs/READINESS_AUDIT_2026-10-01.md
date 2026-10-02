# RSeQC-rust readiness audit and delivery plan

Audit date: 2026-10-01. Source revision: `b39ef0a0545b7e32f0249c5256e9bb0f8298bc41`, with existing documentation edits in the working tree. This audit read implementation code, test definitions, raw benchmark JSON, recorded differential logs, CI and release configuration, and publisher policies. It executed targeted comparator mutations, a mocked timeout probe, and a hand-specified annotation-conversion example. It did not re-run the full Cargo, differential, real-data, or performance suites. Existing pass records are historical evidence, not fresh qualification of a release.

## Assessment

This is a substantial functional implementation suitable for controlled pilot evaluation. Production reliability and publishable performance remain unqualified. The most urgent work is fixing the validation and benchmark machinery and the rat annotation converter; repeating the current benchmarks would preserve their weaknesses.

A useful first release and software paper can describe a validated, explicit CLI profile. Completing Python bindings, every native plot, every organism, and every biological endpoint is required only for claims covering those features. The accepted full standalone target in [ADR 0001](decisions/0001-target.md) remains unfinished: many plots still invoke R, optional helper programs remain, and the Python crate is a stub. A scoped beta is an intermediate milestone, not completion of that target.

## What the evidence currently establishes

| Area | Observed evidence | Limit |
|---|---|---|
| Implementation | All 33 command names have differential cases; 90 cases in the current matrix | Command presence and selected cases do not qualify every option, artifact, or input format |
| Rust testing | 369 `#[test]` declarations and five bounded property/fuzz suites are present | Declaration count is not an executed pass count or a semantic-coverage measure |
| Differential results | Latest recorded real-panel-enabled invocation ran 86 cases, with four skipped and no reported failures | Many cases retain fixed synthetic/regression inputs. This is a mixed-input invocation, not 86 real-data validations or independent samples |
| Error coverage | Three cases declare a shared nonzero exit; two additional cases explicitly expect upstream failure and Rust success | This is limited failure-path coverage; successful empty/malformed cases also exist |
| Scientific endpoints | 18 raw rows: 3 PASS, 3 FAIL, 3 INCONCLUSIVE, 9 NOT_EVALUATED | E1 metadata is unsupported, E3 is development-curve similarity, stratification is missing, and rat annotation preparation is incorrect |
| Performance | 29 command rows, ten runs per arm, raw measurements, older clean source revision `3a5b88c` | All rows are labelled `compute`; machine is shared, pairing is broken, comparators are incomplete, and ten upstream logo failures are recorded |
| Scale | Main workload has 800,000 alignment records; smaller TIN/gene-body workloads use restricted transcript sets; read-count sweep reaches 800,000 requested pairs | No demonstrated whole-genome, tens-of-millions-of-pairs operating envelope |
| Distribution | Linux archive builder, container Dockerfile, three-OS build/test CI configuration, root LICENSE and CITATION.cff | No local version tags; no verified hosted CI history reviewed; archive omits documentation/license; public URLs are placeholders |

The CLI regression suite is materially stronger than early project plans suggest. Conversely, some current summaries overstate real-data coverage and say fuzzing has not happened even though bounded suites and a campaign report exist. [PUBLICATION_PLAN.md](PUBLICATION_PLAN.md), [T4_FINDINGS.md](T4_FINDINGS.md), and parts of the benchmark README contain historical or contradictory status assertions.

## Concrete blockers and their closure

| Priority | Finding and evidence | Required closure |
|---|---|---|
| P0: scientific evidence | `datasets/refgene_to_gtf.py` treats UCSC refGene coordinates as one-based inclusive. Probe input with exons [100,200), [300,400) yields BED start 99 and sizes 101,101, and GTF starts 100,300. Correct BED start is 100, sizes 100,100; GTF starts are 101,301 | Correct and independently test conversion, including minus strand and zero-coordinate boundaries; regenerate BED and GTF, rebuild the affected STAR index, re-align, and compare both implementations on identical corrected inputs. Current rat results cannot support scientific claims |
| P0: benchmark gate | File-tree gate can accept two empty output directories. Captured stdout is not compared, so stdout-only metrics can receive a gate pass without metric checking | Require per-command expected streams, labels, schemas and artifacts. Validate outputs for each timed run, outside the timing interval; equal missing outputs must fail |
| P0: benchmark comparators | Executed probes accept changed BAM quality scores, duplicate flags and NM tags; a FASTQ with a trailing partial record; and finite text replaced by NaN | Compare required BAM fields, headers, mate fields, tags, record order/sorting where required; validate FASTQ framing and lengths; explicitly reject unexpected nonfinite values. Mutation tests must reject every relevant corruption |
| P0: benchmark timeout | Mocked `TimeoutExpired` has no PID; handler calls `getpgid(0)`, then `killpg` on that group, risking termination of the harness's own group | Use a retained `Popen` object, terminate only its newly created process group, reap children, preserve diagnostics, and test that the runner survives |
| P1: statistical pairing | Schedule globally shuffles individual arms, appends results by completion order and then zips them. At seed 20260929, all ten zipped pairs have different repetition IDs | Randomize order inside adjacent matched blocks; retain pair IDs, order and timestamps; bootstrap actual blocks. Recollect measurements for interval claims |
| P1: benchmark provenance | Main evidence identifies an older clean commit, no exact candidate binary hash per run; scaling file contains summaries without per-run samples or its own environment manifest | Bind each run to source revision, binary SHA256, oracle lock, exact inputs/options, comparator version and full raw measurements; rerun at the release candidate |
| P1: scientific interpretation | E1 uses PCR selection as strand expectation; E3/E5 causal interpretation is confounded; rat A6 was added after exposure | Resolve protocol metadata, restrict E3 to its actual estimand, perform required stratification, and retire the rat density explanation pending correct-coordinate evaluation. Fresh independent samples are required for changes prompted by exposed data |
| P1: baseline/CI | `upstream.lock` leaves executable dependencies unpinned; current CI/release jobs do not run Python differential or comparator tests | Freeze a runnable reference environment and add a small mandatory differential/comparator CI layer, plus an archived release validation job |
| P1: production envelope | CRAM branch collects the entire file before processing. Duplicate counting retains unique sequences/positions; WIG generation retains covered positions | Measure memory at representative scale for each cost driver; stream CRAM or explicitly restrict it; define supported memory limits and test failure behavior. Streaming BAM alone does not bound every command's memory |
| P1: artifacts | Archive builder stages arbitrary top-level executables and includes no LICENSE, README or citation file. CI uses `ubuntu-latest` without an explicit ABI support floor; CFF/Cargo use `TBD` URLs | Package an allowlisted command set plus metadata/license/dependency notices, define target/ABI support, smoke-test the actual archive and container, and replace public metadata placeholders before publication |

The coordinate finding is grounded in [UCSC's documented database conventions](https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/), not in fitting a held-out result. The source converter's claim that refGene is one-based is incorrect. Its earlier “+1 correction” repaired neither conversion completely. Rat assembly metadata also conflates rn6 with mRatBN7.2: UCSC identifies rn6 as Rnor_6.0 and mRatBN7.2 as rn7. [UCSC assembly records](https://genome.ucsc.edu/goldenPath/releaseLog.html), [rn7 announcement](https://genome.ucsc.edu/goldenPath/newsarch.html).

These harness findings establish false-pass paths and invalid pairing, not that every native command or historical point estimate is wrong. The measured improvements remain useful engineering leads; their current qualification is insufficient for a publication headline. Raw probes are in [readiness-audit-evidence.json](readiness-audit-evidence.json).

## Tests that are necessary, and tests that can wait

Tests and benchmarks answer different questions. Additional benchmark repetitions cannot repair untested semantics or a comparator that accepts corruption. The next tests should target independent truth and operational failures, not increase the unit-test count indiscriminately.

| Layer | Minimum useful next work | Completion condition |
|---|---|---|
| Harness credibility | The mutations and timeout scenario above; missing artifacts/streams; altered headers/order/mates; wrong index; unexpected exits; incorrect selector; wrong pairing | Each intentional defect produces a failing result; runner survives and records the reason |
| Shared scientific semantics | Independent coordinate examples, pileup overlap/quality/depth boundaries, CIGAR M/= /X/N/D rules, paired exon-membership, strandedness truth table | Hand-calculated or independently generated expectations pass in debug and release; upstream comparison is a separate field |
| Command contracts | Risk-based empty/malformed input, missing index, invalid flags, existing output, paths with spaces, unreadable input and output failure, killed jobs and concurrent invocations | Defined exit/status, complete outputs or explicitly identifiable incomplete outputs; no silent metric loss or successful corrupt output |
| Format interoperability | Read generated BAM/BAI with htslib; validate sorted order and indexed fetch; FASTQ framing/qualities; BigWig interval values and chromosome dictionary; representative SAM/CRAM modes | Independent consumers read every advertised artifact and recover the expected content |
| Real-data qualification | Whole-genome, explicitly documented stranded/unstranded bulk libraries across independent studies; fresh confirmation where exposed results drove changes | Same BAM/BED/options in both arms; per-command artifact checks and clinically irrelevant but operationally important outlier review; no unexplained discrepancies |
| Production pilot | Representative 10M, 50M and, if in scope, 100M pairs; full annotation; low/high duplication and coverage; multiple jobs under fixed resource budgets; real workflow integration | Declared envelope completes or fails predictably, preserves input data, stays within its published budget, and produces reproducible usable outputs |

The proposed sizes are design targets to match intended deployments, not a universal industry certification rule. Start with a small qualification panel; expand where a failure or a specific advertised domain warrants it. Use public, independently aligned BAMs with verified provenance where building whole-genome STAR indexes locally is impractical.

For an initial bulk CLI paper, a new single-cell study, Python API tests, complete native plot parity, every aligner, and an ERCC/degradation experiment are conditional on those claims. A reimplementation does not have to rediscover all of RSeQC's biological methods. It must demonstrate preservation of the relevant quantities, independent boundary correctness, and unchanged interpretation within its declared scope. New accuracy, chemistry, degradation or cross-organism claims need their own appropriate controls and replication.

Bounded proptest suites already exist. Keep fixed regression seeds in normal CI; add varied-seed campaigns with recorded seeds and resource limits after shared-parser changes. Coverage-guided fuzzing and mutation-coverage measurement are useful subsequent work, but are not prerequisites for every low-risk CLI change. “No panic” alone is not evidence of correct biological meaning.

## Benchmark study version 2

Freeze a new dated protocol after repairing the runner; retain the old results as historical development evidence.

1. **Primary estimator:** actual command invocation elapsed time, CPU and memory with equal delivered work. Keep compute/data-only/plotting modes distinct. Treat subtraction of a `--help` floor as a sensitivity analysis: it is an imperfect proxy for startup and adds uncertainty.
2. **Data:** real whole-genome RNA-seq and controlled stress workloads. Report alignment records and read pairs separately. Retain full annotation for primary realistic workloads; reduced transcript panels are diagnostic sweeps.
3. **Scale:** small/intermediate/production sizes and selected high-risk drivers: unique sequences, covered bases, local depth, transcripts/exons, introns, output volume, CRAM buffers. The existing “20,000 transcript” point actually contains 9,179; label measured counts, not requested counts.
4. **Execution:** dedicated or otherwise demonstrably controlled hardware, recorded load/frequency/NUMA/storage, fixed CPU/thread/compression budgets, warm-cache primary study and separately defined cold-cache or deployment conditions if claimed. Start with ten valid matched pairs; use a predeclared precision rule to increase repetitions rather than collecting indefinitely.
5. **Validity:** validate expected streams/artifacts for every run. Preserve timeout/OOM/failure records. Upstream logo failures and absent helper programs require explicitly scoped modes; the present 29 rows are all compute-class runs, not full plotting/conversion pipeline measurements.
6. **Uncertainty:** correct paired blocks, confidence intervals for time and memory as appropriate, failures retained, all command rows including losses. Use replication across independently scheduled sessions and representative datasets; a bootstrap interval on one shared-machine run does not capture dataset or hardware generalizability.
7. **Pipeline impact:** measure a real bulk-QC workflow, samples/hour, CPU-hours/sample, memory at intended concurrency, output/I/O volume, cold-start/install burden, and reproducibility of decisions. If QC is only a fraction of the complete RNA-seq pipeline, overall benefit is correspondingly bounded.
8. **Resources:** GNU time's child RSS is not aggregate simultaneous process-tree memory. Use an explicitly labelled aggregate measurement for helper-heavy workflows, with a separate process RSS measure if needed. Linux documents `RUSAGE_CHILDREN.ru_maxrss` as the largest child's RSS. [Linux resource-usage documentation](https://man7.org/linux/man-pages/man2/getrusage.2.html).
9. **Archive:** raw samples with pair IDs, failures, command lines, full hashes, dependency/container/toolchain versions, workload manifests, analysis code and publication figures. Store raw scaling runs as well as summaries.

Additional performance optimization, internal multithreading and many repetitions on the current small workloads can wait. Parallelism across samples may already be the appropriate production design; first measure that deployment.

## Release sequence

**Stage A — release preparation.** Resolve the P0 evidence defects and the affected converter. Pin reference environment and capability manifest. Record per-command support for inputs, helper-dependent outputs, stochastic behavior and known upstream quirks. Bind a clean commit to the final binaries and evidence. Decide the permanent repository URL, support channel and maintainer policy.

**Stage B — scoped beta, proposed `v0.1.0-beta.1`.** Start with a declared Linux x86_64 ABI and a tested container. Ship versioned source, checksummed allowlisted binaries/aliases, LICENSE, third-party notices, README, citation metadata, machine-readable build/capability manifest and small example inputs/outputs. Validate the extracted archive outside the source tree using real workloads, not only `--help`. Gate tag publication on tests, comparators and the appropriate differential profile. Mark other included capabilities experimental where qualification is absent. The beta does not meet the full native-plot/Python milestone.

**Stage C — workflow pilot and `v0.1.0`.** Obtain external users in a research group and a representative workflow environment; collect reproducible issue reports and fixes. Add a tested Nextflow or Snakemake example and a Bioconda recipe once artifact behavior stabilizes. Linux-only is acceptable initially; ship macOS/Windows/ARM artifacts only after actual installed-artifact validation on each platform. Bioconda's Rust guidance includes locked builds and bundled dependency licenses. [Bioconda recipe guidance](https://bioconda.github.io/contributor/guidelines.html).

**Stage D — stable release and paper snapshot.** Freeze the claimed compatibility profile and document versioning/behavior changes, support expectations, upgrade guidance and rollback. Archive the exact paper source/release/data evidence, assign a release DOI and record container digests. Zenodo can archive GitHub releases and provides guidance for citation metadata. [Zenodo release archiving](https://help.zenodo.org/docs/github/archive-software/github-upload/).

Avoid announcing full-package compatibility while the Python API is a stub, or full runtime-independent plotting while R is required. Those are legitimate later milestones. Public publishing/tagging remains a separate action after the concrete release candidate has passed its gates.

## Publication strategy

The defensible contribution is a scientifically compatible native implementation with demonstrated runtime/memory benefits, simpler deployment, and rigorous interoperability evidence. Explain the architectural work—streaming, interval indexes and sliding windows—and show why it improves a real research workflow. A language change and developer-machine timing alone leave a weak paper.

My recommended application-note route is **Bioinformatics Advances**, with **Bioinformatics Application Notes** considered if the demonstrated advance and study fit its scope. This is an assessment of fit, not an acceptance prediction. The former explicitly values useful software, use cases, quality and documentation. [Bioinformatics Advances author guidelines](https://academic.oup.com/bioinformaticsadvances/pages/author-guidelines). Bioinformatics requires accessible software/test data, archived versions and reproducible results; its Application Notes are short software descriptions. [Bioinformatics author guidelines](https://academic.oup.com/bioinformatics/pages/author-guidelines).

**JOSS is a later option, not an immediate submission shortcut.** Current submission rules require more than six months of public development, active iterative history and demonstrated research use. This local history contains 197 commits from 2026-09-16 to 2026-10-01, no tags, and placeholder public URLs; it does not establish eligibility. Verify actual public history. Its AI policy requires disclosure and human review of assisted work, and limits AI use in author-editor/reviewer conversations. [JOSS submission rules](https://joss.readthedocs.io/en/latest/submitting.html). Its review also evaluates need, comparison to alternatives, design, impact, tests, packaging and support pathways. [JOSS review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html). The technical reviews in this repository are not human domain-expert or journal reviewer sign-off.

Prepare the narrative and methods now; write numerical results only after the study is qualified. The paper should contain:

- A declared command/mode/input scope and the original RSeQC citation.
- A compatibility and scientific-expectation table with discrepancies, stochastic modes and failures.
- Real-data validation plus independent truth cases; retire the current rat conclusion until reference preparation is fixed.
- Runtime/memory scaling and at least one realistic workflow impact example.
- Limitations, supported operating envelope, reproducibility/availability statements, contributor roles, funding and AI assistance disclosure.

An application note can use one composite primary figure and supplementary matrices, raw measurements and methods. A longer methods paper needs a substantive additional methodological contribution and a larger evaluation. A DOI-backed release or preprint can precede journal submission; it does not replace qualification or peer review. [Manuscript outline](MANUSCRIPT_OUTLINE.md).

## Ordered next work and stopping criteria

1. **Repair validation foundations:** converter truth cases, benchmark mutations, expected outputs, safe timeouts and correct pairing. Stop only when each deliberate defect is detected.
2. **Qualify a release candidate:** pinned baseline, per-case provenance and mixed-input reporting, meaningful failure branches, artifact interoperability and executable CI. Stop when every advertised beta capability has a passing recorded profile.
3. **Run whole-genome and production-size pilots:** reuse correctly prepared independent inputs, constrain memory and concurrency, review discrepancies. Stop when the declared deployment envelope has evidence.
4. **Collect the publication study once:** frozen candidate/protocol, matched blocks, equal work, raw samples and regenerated figures. Recollect affected rows only when code, methods or unresolved discrepancies justify it.
5. **Package and review:** validate actual archive/container and an external installation; complete URLs, citation/DOI and support materials; release beta/stable within the passed scope.
6. **Complete manuscript and submit to an eligible venue:** human scientific/methods review, final archived study, clear scope and reproducible artifacts.

These are evidence gates rather than calendar promises. More testing is justified where it can reveal a meaningful undetected defect; more benchmarking is justified after correctness and measurement qualification. Completing every feature in the long-term roadmap is not a prerequisite for a carefully scoped, useful software release.

