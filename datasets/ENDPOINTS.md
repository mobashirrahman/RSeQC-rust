# T4.2 — scientific endpoint specification (PRE-REGISTERED)

Prepared: 2026-09-30. Status: **pre-registered; reviewed and signed off
2026-10-01. Validation not yet run.**

**This file exists because `testing.md` §11.2 requires it to exist first:**

> Before running held-out validation, define which conclusions the outputs
> will support and what change would matter scientifically.

It is written *before* any RSeQC command has been run against the held-out
strata. `datasets/manifest.yaml` records those strata as
`exposure: NOT YET INSPECTED`, and that is deliberate: an endpoint defined
after seeing which samples pass it is not a held-out endpoint.

If any endpoint below turns out to be wrong once the held-out data is
opened, that is a finding to be recorded, not an invitation to rewrite the
threshold.

---

## 0. What T4.1 already established, and what this file is not about

T4.1 ran the 90-case differential matrix against real reads: **86 run, 86
pass, 4 skipped for want of a fixture, 0 failures, 118 byte-identical
assertions.** That answers "does the port reproduce upstream's outputs?".

This file is about the different question §11.2 poses: **are those outputs
scientifically right?** A port can be byte-perfect against an upstream that
is itself wrong on a given input, and on this project that is a live
possibility — `testing.md` §14 explicitly forbids relabelling missing work as
an accepted divergence, and several upstream behaviours in `tin.py` are
documented as probable bugs preserved verbatim (e.g. `read_NVC.py`'s table
sized by the last record processed, `mystat.percentile_list`'s round-half-even
interpolation). Preserving a bug is correct for a port and wrong for a
science claim, and this file is where the two are separated.

Consequence: **no endpoint below is a test of the port.** The port is
assumed byte-identical to upstream (T4.1). These are tests of whether
upstream's own answers survive contact with real data. Where upstream is
wrong, the finding is recorded against upstream, and whether the port should
diverges is a separate decision requiring its own evidence.

---

## 1. Panel, and the split that governs it

| stratum | organism | run | role |
|---|---|---|---|
| development | human | SRR1216016, SRR1216063 | used for all debugging and threshold-setting to date |
| held_out.cross_lab | human | ERR10015758 | validated once, at the end |
| held_out.cross_chemistry | human | ERR10229623 | validated once, at the end |
| held_out.cross_organism | rat | SRR1177982 | validated once, at the end |

Independent unit: **the sample**. §11.2 is explicit that "millions of
transcripts from one sample are not millions of independent biological
replicates", so every per-transcript statistic below is reported with the
sample as the unit of replication, and no interval is computed over
transcripts.

---

## 2. Endpoints

### E1 — Library layout and orientation inference is correct

**Command:** `infer_experiment.py`
**Claim tested:** the port (and upstream) recover the dominant strandedness
of a real library, and bound the unassigned fraction.

- **Statistic:** the fraction of reads assigned to the 1+/1−/2+/2−/unassigned
  buckets.
- **Endpoint:** the dominant bucket is the one the library's protocol
  predicts, **and** the unassigned fraction is below 20%.
- **What would falsify it:** a run where the dominant bucket is not the
  predicted one, or where unassigned exceeds 20%.
- **Known confound to state, not to resolve:** SEQC libraries are stranded
  and the exact protocol differs between the development arm and
  `cross_lab`. The expected bucket is taken **from the archive's
  `library_selection`/protocol metadata**, not assumed to be the same for
  both strata. If the archive does not state strand for a stratum, E1 is
  **not evaluated** on it rather than guessed.
- **Marginal case:** `cross_chemistry` (BGISEQ-500). If its protocol is not
  stated as stranded, E1 reduces to "unassigned fraction bounded" only.

### E2 — TIN responds monotonically to controlled 3′ depletion

**Command:** `tin.py`
**Claim tested:** the transcript-integrity signal tracks degradation in the
direction the biology predicts, so the port's TIN is measuring what it is
supposed to measure.

- **Design:** synthetic 3′-degradation applied to *real reads* before
  alignment. Not a degraded real sample — a controlled perturbation with an
  analytically known effect on the 3′ end.
- **Statistic:** per-transcript TIN.
- **Endpoint:** TIN is **non-increasing** as degradation increases, for the
  large majority of transcripts. No monotonicity threshold is asserted here
  because degradation is applied to already-aligned coordinates and the exact
  perturbation semantics are settled in the implementation task, not here.
- **What would falsify it:** TIN flat or increasing with degradation, which
  would mean TIN is not sensitive to 3′ loss at all.
- **Explicit non-requirement, from §11.2:** "Do not require perfect
  monotonicity in heterogeneous biological samples." This endpoint is on a
  *controlled* perturbation precisely so that a failure means the metric is
  broken, not that the biology is noisy.

### E3 — Gene-body coverage recovers a known bias

**Command:** `geneBody_coverage.py`
**Claim tested:** the percentile-binned gene-body profile recovers an imposed
  coverage skew.

- **Design:** coverage skew imposed by downsampling reads from a
  deterministic coordinate subset before alignment.
- **Statistic:** the 100-percentile aggregated coverage vector.
- **Endpoint:** the recovered skew's direction matches the imposed skew, with
  Pearson r ≥ 0.8 on the normalised vector.
- **Why r and not a stricter bound:** `geneBody_coverage`'s own aggregation
  is bin-indexed, and its `pearson_moment_coefficient` "skewness"
  centres on the middle-indexed value rather than the mean — both documented
  upstream quirks. A strict equality test would test the quirk rather than
  the science.
- **What would falsify it:** r < 0.8, or the wrong direction.

### E4 — Count normalisation responds to spike-in dilution

**Command:** `FPKM_count.py`
**Claim tested:** FPKM tracks a known fold-change.
**This endpoint has an explicitly limited claim, per §11.2:** "FPKM is not
automatically a validated substitute for every expression-analysis method."
The endpoint tests only that FPKM *responds monotonically* to an imposed
dilution, not that its absolute values are correct or that it is an adequate
expression measure.

- **Design:** the SEQC samples are ERCC-spike-in containing. Where an ERCC
  transcript set is present in the model, impose a known dilution series.
- **Statistic:** FPKM per transcript across the series.
- **Endpoint:** rank correlation between imposed log-fold-change and
  reported log-FPKM is ≥ 0.9 over the transcripts present in both.
- **What would falsify it:** correlation below 0.9, or a non-monotone
  response to a monotone perturbation.

### E5 — Junction calls agree with annotation on real data

**Command:** `junction_annotation.py`
**Claim tested:** junction classifications are recovered at the rates the
  method supports on real spliced reads.

- **Statistic:** fraction of called junctions annotated / total / non-canonical.
- **Endpoint:** the annotated fraction exceeds the non-canonical fraction,
  and total junctions land in the range 10⁵–10⁷ for this read depth.
- **Why this range and not an exact count:** junction counts are aligner- and
  depth-dependent, and pinning an exact number here would be asserting
  something about STAR, not about RSeQC.
- **Cross-stratum interest:** `cross_chemistry` is where an implausible
  junction profile would first show, since a DNBSEQ-specific quality or
  splice artefact would appear here and nowhere else.

### E6 — Single-cell outputs are NOT covered

No endpoint. The panel has no single-cell sample (§11.1 domain "droplet
single-cell"), and the 4 skipped differential cases include the only case
needing `sc.bam`. **No single-cell claim may be made**, and this is recorded
as a known gap rather than deferred silently.

---

## 3. Stratification

§11.2 requires errors reported "stratified by coverage, transcript length,
GC, splice complexity, overlap, annotation ambiguity, and relevant quality
flags", and "investigate outliers, not only global correlations". For each
endpoint above, per-sample results are therefore binned by:

- mean coverage over the transcript (low / medium / high)
- transcript length quartile
- GC content quartile
- exon count (as a splice-complexity proxy)
- ambiguity flag: whether the transcript's exons overlap another
  transcript's in the same model

**Rule:** a global pass with a stratum failing is recorded as a **partial
pass**, not a pass. Any stratum failing on more than one held-out stratum
blocks the corresponding claim outright.

---

## 4. What a failure means

| outcome | consequence |
|---|---|
| endpoint passes on all three held-out strata | claim supported, scope limited to §5 |
| passes on development, fails on a held-out stratum | the port still reproduces upstream (T4.1); the finding is against upstream's algorithm. Recorded, and whether the port should diverge becomes its own decision with its own evidence. |
| fails on a controlled perturbation (E2, E3) | the metric is insensitive to what it claims to measure. This is a finding against upstream regardless of the port, and the strongest kind available here. |
| cannot be evaluated (metadata absent, e.g. strand unstated) | recorded as **not evaluated**. Never as a pass. |

---

## 5. Scope limits on any resulting claim

Whatever the outcomes, a claim from this panel is limited to:

- **three contigs of one assembly** (chr1, chr17, chrM of hg38 for the human
  strata) — a hardware limit, recorded in `datasets/manifest.yaml`
- **one aligner** (STAR 2.7.11b). §11.1 forbids claiming cross-aligner
  robustness from one aligner's BAMs.
- **one library per stratum.** §11.1 requires multiple independent
  libraries for a population claim, so no population-level claim is available
  from this panel at all.
- **simulated degradation/coverage skew**, not naturally degraded samples.
  E2 and E3 say so on their face.
- **no single-cell claim** (E6).

`testing.md` §14: none of the above may be closed by relabelling missing
work "accepted divergence".

---

## 6. Sign-off

This specification is pre-registered. It required review before the held-out
strata were opened, because after they are opened the thresholds above stop
being predictions.

- [x] Endpoint definitions reviewed and accepted
- [x] Thresholds and falsification conditions reviewed (especially the two
      that are qualitative: E1's strand precondition and E2's monotonicity)
- [x] Stratification bins reviewed
- [x] Confirmed that no held-out run has been inspected as of sign-off
- [x] Signed off 2026-10-01; held-out data may now be prepared and validated

At the moment of sign-off: the three held-out runs were fetched and
MD5-verified, the rat reference was on disk, and **no RSeQC command had been
run against any held-out run**. The human held-out alignment was in progress
and is data preparation, not inspection.

### 6.1 Amendment after sign-off (2026-10-01)

Two items were flagged as under-specified at pre-registration and are
resolved here. This is an **amendment to a pre-registered document**, made
after sign-off and before any held-out output was seen, so it is recorded as
a dated change rather than folded silently into the text above. No threshold
was loosened; one threshold was made concrete and one definition was made
computable.

**A1 — E2's monotonicity threshold, now concrete.** At pre-registration E2
asserted only that TIN is non-increasing as degradation increases, with no
threshold, because the perturbation's exact semantics were to be settled at
implementation time. Settled as: for each degradation level, the fraction of
transcripts whose TIN is ≤ its level-0 value must be **≥ 0.90**, and the
mean TIN across transcripts must be non-increasing between consecutive
levels. The 0.90 allows for transcripts whose sampled positions do not span
the degraded tail; the mean-level check is what catches a metric that is
merely noisy rather than non-responsive. A metric that is flat (mean change
< 1% across the full degradation series) is a distinct failure from a metric
that responds in the wrong direction, and both are recorded separately.

**A2 — "annotation ambiguity", now defined.** Stratification needs a bin that
two people compute the same way. Defined as: a transcript is **ambiguous** if
any of its exons overlaps by ≥ 1 bp an exon of a different transcript in the
same model. Computable from the BED12 alone with no external tool, and it is
the coarsest defensible definition — it will not catch isoforms that differ
only in UTR structure, which is a known false-negative and is recorded here
rather than hidden.

**A3 — one substantive limit stated more plainly.** E1's strand expectation
is taken from archive metadata. For `cross_lab` that metadata must actually
be checked before the endpoint is evaluated rather than assumed, because
"unassigned fraction bounded" is a much weaker claim than "correct strand
recovered" and the difference is exactly the kind that gets lost in
summarising.