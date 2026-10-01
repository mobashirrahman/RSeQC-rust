# T4 findings: scientific qualification

Status: T4 **partially complete**. The held-out panel is built, aligned and
evaluated; four of the eight coverage domains in `testing.md` §11.1 remain
uncovered. This document records what was found, what it means, and what is
still owed. It is the honest state, not a summary that flatters the result.

Artifacts: `datasets/manifest.yaml` (split and provenance),
`datasets/ENDPOINTS.md` (pre-registered spec, signed 2026-10-01, amendments
6.1 A1–A6), `datasets/ENDPOINT_RESULTS.md` (verdict table),
`verification/validate_endpoints.py` (harness),
`verification/run_diff.py` (byte-level differential).

## 1. The result that qualifies the port

Port vs upstream RSeQC (`oracle/upstream-src`), same held-out BAM
(`ERR10015758`, 350 MB, ~19 M reads, Charité / Illumina HiSeq 2500):

| metric | port | upstream |
|---|---|---|
| Fraction of reads failed to determine | 0.2024 | 0.2024 |
| Explained by `1++,1--,2+-,2-+` | 0.3971 | 0.3971 |
| Explained by `1+-,1-+,2++,2--` | 0.4005 | 0.4005 |

Identical to four decimals on held-out data. `geneBody_coverage` and
`junction_annotation` differentials on the same BAM were still running at the
time of writing (upstream is pure Python; ~85 min for 50,725 transcripts) and
are **not yet claimed**.

Separately, T4.1 holds: 86/90 differential cases run on real reads with 86
passing, 4 explicit skips, 0 failures, 118 byte-identical assertions.

## 2. Held-out endpoint results

| endpoint | cross_lab | cross_chemistry | cross_organism |
|---|---|---|---|
| E1 layout/orientation | PASS 0.4005 | INCONCLUSIVE | INCONCLUSIVE |
| E3 gene-body skew shape | FAIL r=0.627 | PASS r=0.948 | FAIL r=0.670 |
| E5 junction classification | PASS 0.706 | FAIL 0.411 | INCONCLUSIVE |
| E2 TIN under degradation | NOT_EVALUATED | NOT_EVALUATED | NOT_EVALUATED |
| E4 FPKM vs spike-in | NOT_EVALUATED | NOT_EVALUATED | NOT_EVALUATED |
| E6 single-cell (non-endpoint) | NOT_EVALUATED | — | — |

**E3 and E5 are not port tests.** Given §1, a divergence from the development
reference is a property of the library, not a defect. Read as port defects they
would be false accusations; read as library characterisation they are correct
and useful:

- **E3 splits on library prep.** BGISEQ-500 reproduces the dev 3'-skew shape
  closely (r=0.948, peak at the 80th percentile vs the reference's 71st);
  HiSeq-2500 and rat do not (0.627, 0.670). A single threshold across
  chemistries measures the chemistry, not the software.
- **E5 splits on library prep too, and this one is a real failure.** 0.706
  (Charité) vs 0.411 (Tartu) against the *same* GENCODE annotation on the
  *same* three contigs, so Tartu's 0.411 is not an annotation artifact.

## 3. Bugs found, and what each one would have looked like

Three, all caught by held-out work rather than by tests.

**(a) `datasets/refgene_to_gtf.py` — BED12 off-by-one, mine.** refGene exons
are 1-based inclusive, so a BED12 exon size is `end - start + 1`. The
converter emitted `end - start`, making every rat exon one base short and every
exon end one base early. Junction matching collapsed to 7 annotated of 44,176
(0.02%). Fixing it moved the value to 1,134 (2.57%) — a 162x improvement,
which is what confirmed the off-by-one rather than leaving it a plausible
guess. The converter now agrees with `make_bed12.py` base-for-base, and all
19,160 rat transcripts are asserted to reconstruct their span exactly. The
STAR-facing GTF was never affected, so no re-index or re-alignment was needed.
*Would have appeared as:* a finding that rat is terribly annotated.

**(b) Rat STAR index — `sjdbOverhang` 149 on 101 bp reads.** The index was
built with the human default. A too-large overhang sizes STAR's junction
database for reads longer than it will see, degrading splice detection at the
long end — precisely what E5 measures. *Would have appeared as:* a finding
about rat annotation quality, when it was an alignment misconfiguration. Caught
by STAR's progress line, not by any check that existed; the value is now an
explicit override and is part of the index stamp key.

**(c) `datasets/build_star_index.sh` — `gunzip -c in > out` leaves a
zero-byte `out` when gunzip fails, and `set -e` does not remove it.** One
failed unpack therefore leaves a file every later run trusts as "already
unpacked", and an index is silently built from it. *Would have appeared as:* a
0-byte junction database, i.e. a port that cannot see any junction.

A fourth trap is recorded rather than fixed, because it is a documentation
hazard: **STAR's progress `Read length` column reports the summed mate
length** (202 for a 2x101 library, 300 for 2x150). The overhang must come from
the FASTQ, not that column.

## 4. Methodological findings: three endpoints were unevaluable as written

All three were found on **development** data, with the held-out BAMs not yet
existing, and fixed as dated amendments (`ENDPOINTS.md` §6.1) rather than
edited into the pre-registered text. A pre-registered document whose
thresholds quietly change is not pre-registered.

**E3's threshold sat above its own metric's noise floor.** Splitting the dev
panel's reads into two independent halves and correlating their gene-body
curves gives **r=0.744** — below E3's own 0.80 threshold — and the peak
percentile moves from 63 to 75 between two halves of the same library. So the
endpoint would have failed for reasons unrelated to biology or the software.

| smoothing | split-half r |
|---|---|
| raw 100 points | 0.744 |
| moving average, width 9 | 0.801 |
| **moving average, width 15** | **0.842** |
| 25 bins instead of 100 | 0.770 |

Width 15 is frozen; width 9 clears by 0.001, which is not a margin. Bin
reduction does not work, so smoothing is the mechanism. **Circularity stated:**
the width was chosen because it is the smallest putting the metric's *measured*
reliability below the threshold.

**E5's absolute junction window was depth-dependent.** The pre-registered
1e5–1e7 total assumed full depth; the 324k-read dev panel has 2,038 junctions
and failed it while its substantive claim passed comfortably. E5 is now the
annotated fraction ≥ 0.50, with the count kept as a diagnostic. The floor sits
below the dev value (0.697) deliberately, because rat's annotation is UCSC
refGene and a lower annotated fraction there is a property of the annotation.

**E5 is not evaluable when the annotation is too sparse to annotate the
index.** The rat index annotates **32%** of its own sequence (chr1 32%, chr2
27%, chr10 43%); the human index annotates **503%**, because GENCODE places
many overlapping transcripts per locus and rat refGene is a thin set. STAR
invents introns across the unannotated remainder — the rat BAM genuinely
contains `55M498576N46M` — and the port reports them as novel, correctly. E5
returns INCONCLUSIVE below 50% coverage, with the figure attached. The 0.5 bar
is a priori and leaves the human strata untouched.

**A hypothesis I had to discard:** giant introns looked like the rat-specific
cause until the human BAM turned out to have *more* of them (2,352,546 vs
1,378,374) while scoring 0.706 against rat's 0.026. Intron size was a red
herring; annotation density is the explanation.

**A trap the harness caught in itself:** the first E3 implementation returned
`r=1.0` on the dev panel, because the reference curve was derived from the
sample being scored. The harness now detects that case from recorded
provenance and returns NOT_EVALUATED. The dev PE/SE cross-check gives
r=0.999, but those two panels are views of the same reads, so that agreement is
nearly guaranteed and is **not** evidence E3 survives a change of lab,
sequencer or species.

## 5. Coverage against `testing.md` §11.1 — what is still owed

| §11.1 domain | status |
|---|---|
| Bulk paired-end stranded / unstranded | covered (SEQC dev; Charité `library_selection=PCR`; Tartu) |
| Different preparation methods and read lengths | covered (101 bp rat, 150 bp human; PCR vs cDNA vs SEQC) |
| A second organism / reference structure | **partial** — transcript structure and coordinate conventions yes; **chromosome naming is not covered at all**, because rn6 uses the same `chr1`-style naming as hg38 |
| Bulk single-end | dev only (`se.bam`), not held-out |
| Experimentally degraded RNA / integrity series | **not covered** — E2's re-alignment was not run |
| Reference mixtures / spike-ins | **not covered held-out** — zero ERCC in all three runs; only dev has spike-ins |
| Droplet single-cell | **not covered** — no dataset exists; no claim made |
| Deep / sparse / high-duplication / multimapping | **not covered** |

Four domains outstanding. T4 cannot be called passed on this evidence.

## 6. Known limitations of the held-out panel

- **Aligner dependence is unexamined.** All three runs were aligned by the
  same STAR 2.7.11b with the same options. §11.1 forbids claiming robustness
  across aligners from one aligner's BAMs, so no such claim is made — but the
  cross-lab and cross-chemistry strata differ from development *only* in
  library, not in how the BAM was produced.
- **Contig-subset indexes** (`chr1 chr17 chrM`; `chr1 chr2 chr10`) because
  whole-genome GENCODE indexing OOMs at 31 GB and swap cannot be enabled. This
  confounds E1's absolute unassigned fraction on two of three strata, and it
  is the direct cause of the rat junction situation in §4.
- **E1's expectation is metadata-dependent.** `library_selection=PCR` implies
  unstranded, so cross-lab is evaluated as a strand-*balance* check
  (dominant ≤ 0.65, fixed from the dev reference's 0.369/0.365). `cDNA`
  selection for rat does not determine strandedness, so rat's strong claim is
  not evaluable at all. Unstated protocol is recorded as *not evaluated*,
  never as *passed*.
- **Lab verification came from `center_name`, not accession inference**, which
  changed a conclusion: the rat arm's submitter is Fudan, the same centre as
  development, so cross-lab coverage genuinely rests on the Charité run rather
  than on the rat arm providing it.

## 7. Open items unrelated to T4

- **`DIV-0003` is open.** GPL-3.0-or-later was chosen and is declared in
  `Cargo.toml` and `CITATION.cff`, but no root `LICENSE` file has been added.
- **Benchmark artifacts trace to commit `3a5b88c`** (clean tree), while the
  branch head is later. The intervening changes are T4 data work plus a
  `geneBody_coverage` error-path fix that does not touch the benchmarked path,
  but a release claiming performance at the final tag needs a rerun.
