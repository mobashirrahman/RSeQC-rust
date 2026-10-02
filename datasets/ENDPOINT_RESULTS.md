# Endpoint results

T4.2 held-out endpoint outcomes. The pre-registered specification is
[`ENDPOINTS.md`](ENDPOINTS.md); the validator is
`verification/validate_endpoints.py`. Every number here is traceable to a command
output kept under `datasets/heldout/endpoint_results/`.

**Read this first.** The raw verdict tally is not completed scientific validation,
and what a verdict supports is narrower than its name suggests. Three of these rows
were corrected after an independent technical audit on 2026-10-01, and the
corrections are stated rather than quietly applied:

- **E1's cross-lab PASS is withdrawn.** Its strand expectation came from ENA
  `library_selection`, which records how a library was amplified, not its strandedness.
- **A6 is withdrawn.** The rat junction density explanation was written after the rat
  result was seen and rested on a mis-stated coverage figure.
- **The rat annotation conversion was wrong** and is now corrected, independently
  confirmed against the genome sequence rather than against documentation.
- **E3's estimand is narrower** than its mechanistic reading: a correlation between one
  library's curve and the development panel's is a similarity between two samples.

The validator now enforces each of these: E1 is NOT_EVALUATED without explicit
protocol metadata, E3 records its estimand unless `--estrand` is given, and feature
stratification exists and already shows the aggregate E3 pass concealing
stratum-level variation.

## Summary

| Stratum | Run | E1 | E3 | E5 | Raw tally |
|---|---|---|---|---|---|
| cross_lab | ERR10015758 (Charité, HiSeq 2500) | **PASS — claim retracted, see below** | FAIL | PASS | 2 pass, 1 fail |
| cross_chemistry | ERR10229623 (BGISEQ-500) | PASS | PASS | PASS | 3 pass |
| cross_organism | SRR1177982 (rat, rn6) | NOT_EVALUATED | FAIL | **PASS (0.617, corrected coordinates)** | 1 pass, 1 fail, 1 not evaluated |

The rat E5 figure was produced after the coordinate correction, the index rebuild and
the re-alignment described below. It is the only endpoint row here computed on
corrected reference preparation.

E2, E4 and E6 are NOT_EVALUATED on all three strata, with reasons.

## What each verdict actually supports

### E1 — cross_lab is retracted

**The PASS is withdrawn.** It was obtained by evaluating `library_selection=PCR` as
an unstranded-library expectation. That field records *how the library was
amplified*; it does not record the library's strandedness, which comes from the
library-preparation protocol. The audit found this inference and it does not hold.

`validate_endpoints.py` no longer infers a strand expectation from library
selection. Without explicit `--protocol-strandedness`, E1 is **NOT_EVALUATED** and
the measured fractions are recorded as observations only. A forward or reverse
expectation combined with `library_selection=PCR` is now rejected as
self-contradictory rather than evaluated.

Consequence: cross_lab's E1 is **NOT_EVALUATED** pending protocol metadata. The
dominant-fraction figure (0.369/0.365 on the development panel, well inside the
0.65 balance bound) remains in the raw results and is consistent with an unstranded
library, but that is an observation consistent with the hypothesis, not a test of
it.

E1 on the rat stratum remains NOT_EVALUATED for the same reason: `library_selection`
is `cDNA`, which does not determine strandedness either.

### E3 — the passing strata support a narrower claim than "chemistry matters"

The statistic is the Pearson correlation between a held-out library's gene-body
coverage curve and the **development panel's** curve, smoothed at a width frozen
before the held-out runs were scored.

That measures similarity between two samples. It is **not** evidence about
degradation, reverse transcription, or library-preparation chemistry: a mechanism
and a non-mechanism produce the same correlation. The BGISEQ-500 versus HiSeq-2500
contrast (r = 0.948 versus 0.627) is additionally confounded — those two libraries
differ in many respects besides chemistry, and the contrast cannot attribute the
difference to chemistry.

`validate_endpoints.py --estrand` now names the estimand being claimed. Without it,
the endpoint records the similarity estimand explicitly and refuses the mechanistic
reading. The supported claim is: *these libraries' gene-body coverage profiles
resemble the development panel's by more or less*. Nothing wider.

The HiSeq-2500 FAIL (r = 0.627) and the rat FAIL stand as measurements; what they
do not license is a mechanistic attribution.

### E5 — the rat result is superseded, and the density explanation is retired

The rat annotation was regenerated with the corrected coordinate conversion, the STAR
index was rebuilt from it, and the reads were re-aligned. The annotated fraction is
**0.617** of junctions (26,905 / 43,594) and **0.948** of splice events
(2,138,286 / 2,255,750), against 0.026 and 0.003 for the superseded run.

Both clear the pre-registered 0.50 bar on their own terms, so the density-based
INCONCLUSIVE is unnecessary. The port is not implicated: the change is in reference
preparation, and both implementations read the same BAM.

**The figure was confirmed by running both implementations on the identical corrected
inputs**, which is the condition the audit set before any rat number could be quoted.
`junction_annotation.py` and `target/release/junction_annotation` were run on the
same re-aligned BAM and BED12 with the same options and produced byte-identical
`junction.xls`, stdout and stderr
([`heldout/endpoint_results/SRR1177982_same_input/`](heldout/endpoint_results/SRR1177982_same_input/)).
The `junction.xls` is also byte-identical to the one recorded above, so 0.617 is
reproduced rather than restated.

Confirming it required finding a third defect in the same family. The BAM on disk had
been aligned at 11:30 against an index rebuilt at 21:23, and every digest check passed
because the annotation on disk really was the corrected one — only the alignment was
stale. Re-running the endpoint on that BAM reproduced the *intermediate* 0.585 while
the recorded table said 0.617, and nothing recorded could distinguish them: the
alignment's own `align.json` named the index's parameters but not its annotation
digest. `verification/check_rat_reference.py` now asserts that binding, and
`datasets/build_star_index.sh` gives each assembly its own index directory so
"rebuild the rat index" cannot overwrite the human one. See
[`ENDPOINTS.md` §6.3](ENDPOINTS.md).

The three-way comparison is in [`ENDPOINTS.md` §6.2](ENDPOINTS.md). It records an
intermediate state in which the BED12 was corrected but the index still carried the
old coordinates, because the index builder silently reused a stale contig subset.
That is worth keeping visible: it shows the annotation and the junction database
needed separate correction, and that a "successful" rebuild had produced an
indistinguishable-from-the-original index.

**A6 is withdrawn.** The explanation added on 2026-10-01 — that the rat annotation
covers too little of the indexed sequence for E5 to be evaluable — was formulated
*after* the rat result was seen, and the coordinate conversion underneath it was
independently found incorrect. It is neither confirmed nor refuted on corrected
data. The 32% figure it rested on also needs restating: it was computed from
transcript **spans** (chromStart..chromEnd), which counts introns as annotated.
Computed from merged **exon bases**, the corrected rat annotation covers **1.8%**
of the indexed sequence (chr1 2.0%, chr2 1.1%, chr10 3.0%).

The correction is large but not the whole story. The port is not implicated: the
change is in the reference preparation, and both implementations read the same
BAM.

### Cross-lab coverage rests on one stratum

The rat arm's submitting centre is Fudan, the same centre as the development panel,
so `cross_organism` is a cross-**organism** stratum and not a cross-lab one.
Cross-lab coverage genuinely rests on the Charité run alone — whose E1 claim is the
one retracted above.

### Aligner dependence is unexamined

All three runs were aligned by STAR 2.7.11b with the same options. `testing.md`
§11.1 forbids claiming robustness across aligners from one aligner's output, so no
such claim is made.

### Contig-subset indexes confound one figure

The indexes cover `chr1 chr17 chrM` (human) and `chr1 chr2 chr10` (rat) because
whole-genome GENCODE indexing does not fit in this machine's 31 GB and swap cannot
be enabled. E1's absolute unassigned fraction is confounded by off-index mates on
two of three strata; the threshold is preserved unchanged and re-specifying it would
be a separate pre-registration.

## Coverage domains still absent

| Domain | Status |
|---|---|
| Bulk single-end | development only (`se.bam`), not held-out |
| Experimentally degraded RNA / integrity series | **not covered** — E2's re-alignment was not run for the held-out strata |
| Reference mixtures / spike-ins | **not covered held-out** — zero ERCC in all three runs |
| Droplet single-cell | **not covered** — no dataset exists; no claim made |
| Deep / sparse / high-duplication / multimapping | **not covered** |
| Whole-genome index | **not covered** — contig subsets only |

## Feature stratification

`validate_endpoints.py --stratify` now implements the breakdown §3 of `ENDPOINTS.md`
requires: coverage, transcript length, exon count, GC content and annotation
ambiguity. A global pass with a stratum failing is a **PARTIAL PASS**.

On the development panel this is already informative: the aggregate scores
r = 0.81 while six of thirteen strata fall below 0.80, including exon-count strata
at r = 0.42 and 0.76 and the shared-annotation stratum at r = 0.57. So the aggregate
E3 PASS conceals substantial stratum-level variation even on the panel where the
aggregate passes. That is a finding about the endpoint's aggregate statistic, not
about the port: an aggregate over heterogeneous features is not a summary of them.

The GC stratum needs a genome FASTA (`--genome`); it is reported as unavailable
rather than dropped, since a silently missing stratum is the failure §3 exists to
prevent.

## Retired dispositions

| ID | Disposition |
|---|---|
| A6 (E5 INCONCLUSIVE on annotation density) | **WITHDRAWN.** Post-exposure and resting on a mis-stated coverage figure. Superseded by the corrected-coordinates run. |
| E1 strand expectation from `library_selection` | **WITHDRAWN.** Wrong field; see above. |
| rat E3 FAIL as a chemistry/biology signal | **RETIRED.** Confounded and post-exposure; the stratum is not independent confirmation of anything. |

## What would be needed to accept scientific claims

1. Library-preparation protocol metadata for each stratum, so E1 has an expectation
   that is not inferred.
2. Fresh independent samples for anything prompted by these results. All three
   held-out runs are consumed: each was inspected during validation.
3. Whole-genome indexes, or an explicit acceptance that contig-subset indexes bound
   which claims are available.
4. A second aligner, if any cross-aligner robustness claim is wanted.
5. Port-versus-upstream comparison on the same held-out inputs and options for E3
   and E5. The current held-out rows score the port's output against an endpoint;
   they do not compare the two implementations on held-out data.
