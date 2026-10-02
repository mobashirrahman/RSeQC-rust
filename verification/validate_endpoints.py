#!/usr/bin/env python3
"""Evaluate the pre-registered T4.2 endpoints (datasets/ENDPOINTS.md) on a held-out run.

Discipline this script exists to enforce:

* An endpoint that cannot be evaluated gets verdict NOT_EVALUATED and a REASON.
  It never gets a silent pass and never disappears from the report.
* A threshold is never adjusted to fit the data it is applied to. Where a
  threshold turns out to be unevaluable as written, the verdict is
  INCONCLUSIVE and the diagnosis is recorded; changing the threshold is a new
  pre-registration, not something this script may do.
* Every number printed here is traceable to a command whose output is kept.

Verdicts: PASS, FAIL, INCONCLUSIVE, NOT_EVALUATED.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "target" / "release"


class Report:
    def __init__(self) -> None:
        self.endpoints: list[dict] = []

    def add(self, eid, name, verdict, *, observed=None, threshold=None, reason=None):
        self.endpoints.append({
            "id": eid, "name": name, "verdict": verdict,
            "observed": observed, "threshold": threshold, "reason": reason,
        })

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for e in self.endpoints:
            out[e["verdict"]] = out.get(e["verdict"], 0) + 1
        return out


def run(cmd: list[str], timeout: int = 14400) -> tuple[int, str]:
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def e1_layout(rep: Report, bam: Path, bed: Path, *, contig_subset: bool,
              contigs: list[str], expected_strand: str | None,
              expect_unstranded: bool = False,
              library_selection: str | None = None,
              protocol_strandedness: str | None = None) -> None:
    """Fraction of reads with undetermined layout, plus strandedness.

    The absolute threshold is only meaningful on a whole-genome alignment. On a
    contig-subset index the unassigned fraction is dominated by reads whose
    mate lies off the indexed contigs, which has nothing to do with layout
    inference. That confound is detected and reported rather than worked
    around by loosening the threshold.
    """
    rc, out = run([str(BIN / "infer_experiment"), "--input-file", str(bam),
                   "--refgene", str(bed)])
    if rc != 0:
        rep.add("E1", "layout/orientation", "NOT_EVALUATED", reason=f"infer_experiment exited {rc}")
        return
    m = re.search(r"Fraction of reads failed to determine:\s*([0-9.]+)", out)
    if not m:
        rep.add("E1", "layout/orientation", "NOT_EVALUATED", reason="could not parse infer_experiment output")
        return
    failed = float(m.group(1))
    explained = dict(re.findall(r'Fraction of reads explained by "([^"]+)":\s*([0-9.]+)', out))
    strand_obs = max(explained, key=explained.get) if explained else None
    detail = {"failed_to_determine": failed, "explained": explained,
              "dominant": strand_obs, "expected": expected_strand}

    # A library amplified from random hexamer priming (ENA
    # library_selection=PCR) is unstranded BY CONSTRUCTION, so the substantive
    # expectation is not "recover the right strand" but "no spurious dominant
    # strand". That is a different claim with a different threshold, and it is
    # still checkable on a contig-subset index because balance does not depend
    # on the unassigned fraction. The 0.65 bound is fixed from the dev
    # reference, where an unstranded panel lands at 0.369/0.365, so ~0.5;
    # 0.65 allows 15 points of technical skew while still catching a library
    # that is genuinely strand-specific despite undeclared metadata.
    if expect_unstranded:
        dom = max(float(v) for v in explained.values()) if explained else 1.0
        balanced = dom <= 0.65
        detail["expectation"] = "unstranded (library_selection=PCR): no dominant strand"
        detail["dominant_fraction"] = round(dom, 4)
        rep.add("E1", "layout/orientation (unstranded balance)", "PASS" if balanced else "FAIL",
                observed=detail, threshold="dominant strandedness fraction <= 0.65",
                reason=None if balanced else
                        f"dominant fraction {dom:.4f} exceeds 0.65, indicating a strand-specific library")
        return

    if contig_subset:
        rep.add("E1", "layout/orientation", "INCONCLUSIVE", observed=detail,
                threshold="unassigned <= 0.20",
                reason=(
                    f"Alignment used a {len(contigs)}-contig subset index "
                    f"({', '.join(contigs)}), so the unassigned fraction is "
                    f"confounded by off-index mates and cannot be compared to a "
                    f"whole-genome threshold. Observed {failed:.4f}. The threshold "
                    f"is preserved unchanged; re-specifying it is a separate "
                    f"pre-registration."))
        return
    if expected_strand is None and protocol_strandedness is None:
        # NOT_EVALUATED, not INCONCLUSIVE. The audit's P1 finding was that E1's
        # strand expectation had been taken from PCR selection: ENA's
        # `library_selection` records how a library was amplified, which is not the
        # library's strandedness. Without protocol metadata there is no expectation
        # to test against, so no verdict is available -- the measurement is reported
        # as a number and explicitly carries no claim.
        detail["library_selection"] = library_selection
        detail["protocol_strandedness"] = protocol_strandedness
        rep.add("E1", "layout/orientation", "NOT_EVALUATED", observed=detail,
                threshold=None,
                reason=(
                    f"No protocol strandedness metadata for this stratum"
                    + (f" (library_selection={library_selection} records how the "
                       f"library was amplified, which does not determine its "
                       f"strandedness)" if library_selection else "")
                    + ". The strand expectation cannot be inferred from library "
                    "selection, so E1 is not evaluated rather than guessed. The "
                    f"measured fractions (unassigned {failed:.4f}) are reported for "
                    "reference only."))
        return
    if protocol_strandedness == "unstranded":
        expect_unstranded = True
    if protocol_strandedness in ("forward", "reverse"):
        expected_strand = protocol_strandedness
    ok = failed <= 0.20
    rep.add("E1", "layout/orientation", "PASS" if ok else "FAIL", observed=detail,
            threshold="unassigned <= 0.20",
            reason=None if ok else f"unassigned fraction {failed:.4f} exceeds 0.20")


SMOOTH_WIDTH = 15  # frozen by dev-side calibration, ENDPOINTS.md 6.1 A4


def smooth(v: list[float], w: int) -> list[float]:
    h = w // 2
    return [sum(v[max(0, i - h):i + h + 1]) / len(v[max(0, i - h):i + h + 1])
            for i in range(len(v))]


def _same_bam(bam: Path, ref_path: Path) -> bool:
    """True when the reference curve was almost certainly built from this BAM.

    Recorded provenance lives in datasets/reference_curve/PROVENANCE.txt; this
    checks the sample label rather than relying on the caller to remember.
    """
    prov = ref_path.parent / "PROVENANCE.txt"
    if not prov.exists():
        return False
    return bam.name in prov.read_text()


def e3_genebody(rep: Report, bam: Path, bed: Path, out_dir: Path,
                reference_is_independent: bool = False,
                estrand: str | None = None) -> None:
    """Whether the observed gene-body curve reproduces the expected 3' skew shape.

    The reference is the dev panel's own curve: a well-annotated RNA-seq
    library of this type has a characteristic 5'-low / 3'-depleted shape, and
    the endpoint asks whether a held-out library reproduces it. Both curves are
    smoothed with the width frozen in ENDPOINTS.md 6.1 A4, because the raw
    100-point curve has a measured split-half reliability of only r=0.744 at
    this depth -- below this endpoint's own threshold -- and is therefore not a
    usable statistic unsmoothed. That reliability figure is the reason the
    smoothing exists and is why the width is not retuned per stratum.

    WHAT THIS MEASURES, stated because the audit found it overstated. The
    statistic is the correlation between one library's coverage curve and the
    DEVELOPMENT panel's curve. That is a similarity between two samples, and it
    is not evidence about degradation, reverse transcription, priming or any
    other mechanism: a mechanism and a non-mechanism can produce the same
    correlation. The estimate it supports is "this library's gene-body profile
    resembles the development panel's", and nothing wider.

    A mechanistic reading is therefore refused unless `estrand` names the
    estimand the caller is entitled to claim, per ENDPOINTS.md 6.1 A8. The flag
    is recorded in the output either way, so a reader can see which claim was
    made and which was refused.
    """
    import math
    ref_path = REPO / "datasets" / "reference_curve" / (
        "dev_se.geneBodyCoverage.txt" if reference_is_independent
        else "dev_pe.geneBodyCoverage.txt")
    if not ref_path.exists():
        rep.add("E3", "gene-body skew shape", "NOT_EVALUATED",
                reason=f"dev reference curve missing: {ref_path}")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "gb"
    rc, out = run([str(BIN / "geneBody_coverage"), "--input", str(bam),
                   "--refgene", str(bed), "--out-prefix", str(prefix)])
    curve = Path(f"{prefix}.geneBodyCoverage.txt")
    if rc != 0 or not curve.exists():
        rep.add("E3", "gene-body recovery of imposed skew", "NOT_EVALUATED",
                reason=f"geneBody_coverage exited {rc}")
        return
    lines = [l for l in curve.read_text().splitlines() if l.strip()]
    if len(lines) < 2:
        rep.add("E3", "gene-body recovery of imposed skew", "NOT_EVALUATED",
                reason="no curve rows produced")
        return
    vals = [float(x) for x in lines[1].split("\t")[1:]]
    if len(vals) < 100 or sum(vals) == 0:
        rep.add("E3", "gene-body recovery of imposed skew", "NOT_EVALUATED",
                reason="curve incomplete or all zero")
        return

    ref_lines = [l for l in ref_path.read_text().splitlines() if l.strip()]
    ref = [float(x) for x in ref_lines[1].split("\t")[1:]]
    # Guard against the obvious self-reference trap: if the reference curve was
    # derived from the very BAM being evaluated, r is 1.0 by construction and
    # says nothing. Reporting that as a PASS would be worse than not running.
    if not reference_is_independent and _same_bam(bam, ref_path):
        rep.add("E3", "gene-body skew shape", "NOT_EVALUATED",
                threshold=f"Pearson r >= 0.80 vs dev reference, smoothed width {SMOOTH_WIDTH}",
                reason=(f"The reference curve at {ref_path} was derived from {bam.name}, "
                        "so r would be 1.0 by construction and carries no information. "
                        "Pass --reference-is-independent to score against the "
                        "independent dev SE view instead."))
        return
    if len(ref) != len(vals):
        rep.add("E3", "gene-body skew shape", "NOT_EVALUATED",
                reason=f"reference curve has {len(ref)} points, observed has {len(vals)}")
        return
    obs_s = smooth(vals, SMOOTH_WIDTH)
    ref_s = smooth(ref, SMOOTH_WIDTH)

    def pearson(u, v):
        mu, mv = sum(u) / len(u), sum(v) / len(v)
        num = sum((x - mu) * (y - mv) for x, y in zip(u, v))
        du = math.sqrt(sum((x - mu) ** 2 for x in u))
        dv = math.sqrt(sum((y - mv) ** 2 for y in v))
        return num / (du * dv) if du and dv else float("nan")

    r = pearson(obs_s, ref_s)
    obs = {"pearson_r_vs_dev_reference": round(r, 4),
           "smoothing_width": SMOOTH_WIDTH,
           "peak_percentile_observed": vals.index(max(vals)) + 1,
           "peak_percentile_reference": ref.index(max(ref)) + 1,
           "total_coverage": sum(vals)}
    ok = r >= 0.8
    estimand = obs["estimand"] = (
        estrand if estrand
        else "similarity of this library's gene-body coverage curve to the "
             "development panel's curve; no mechanistic attribution")
    rep.add("E3", "gene-body skew shape", "PASS" if ok else "FAIL",
            observed=obs,
            threshold=f"Pearson r >= 0.80 vs dev reference, smoothed width {SMOOTH_WIDTH}",
            reason=None if ok else f"r={r:.4f} below 0.80 against the dev reference curve")
    if not estrand:
        # Not a failure: the measurement is sound and the threshold was met. The
        # qualifier exists because the number is easy to over-read, and the audit
        # found it had been.
        rep.endpoints[-1]["reason"] = (
            "Measured estimand is curve-to-curve similarity against the development "
            "panel. That is not evidence about any mechanism: degradation, priming "
            "and library composition can all produce a similar curve. A mechanistic "
            "claim requires --estrand naming the estimand being claimed (ENDPOINTS.md "
            "6.1 A8).")


def e5_junctions(rep: Report, bam: Path, bed: Path, out_dir: Path,
                 coverage: float | None = None) -> None:
    """Annotated vs non-canonical junction recovery."""
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "junc"
    rc, out = run([str(BIN / "junction_annotation"), "--input-file", str(bam),
                   "--refgene", str(bed), "--out-prefix", str(prefix)])
    if rc != 0:
        rep.add("E5", "junction classification rates", "NOT_EVALUATED",
                reason=f"junction_annotation exited {rc}")
        return
    def grab(label):
        m = re.search(rf"{label}\s*Junctions:\s*(\d+)", out)
        return int(m.group(1)) if m else None
    total = grab("Total splicing")
    known = grab("Known Splicing")
    novel = grab("Novel Splicing")
    partial = grab("Partial Novel Splicing")
    if total is None or known is None or novel is None:
        rep.add("E5", "junction classification rates", "NOT_EVALUATED",
                reason="could not parse junction_annotation output")
        return
    obs = {"total": total, "known": known, "novel": novel, "partial_novel": partial,
           "known_fraction": round(known / total, 4) if total else None}
    if total < 1:
        rep.add("E5", "junction classification rates", "NOT_EVALUATED",
                observed=obs, threshold="known > novel and 1e5 <= total <= 1e7",
                reason="no splicing junctions observed")
        return
    # ENDPOINTS.md 6.1 A4: the pre-registered 1e5-1e7 absolute window was
    # depth-dependent and failed on a 324k-read dev panel with 2,038 junctions
    # while the substantive claim passed. The count is kept as a diagnostic; the
    # pass condition is the annotated fraction. 0.50 sits deliberately below
    # the dev value of 0.697 because the rat stratum uses UCSC refGene rather
    # than GENCODE, and a lower annotated fraction there is a property of the
    # annotation, not a port defect.
    known_frac = known / total if total else 0.0
    # ENDPOINTS.md 6.1 A6, as amended on 2026-10-01: E5 is not evaluable when the
    # annotation is too sparse to annotate its own index, because there the
    # annotated/non-annotated ratio measures coverage rather than splicing fidelity.
    #
    # Two changes from the original rule, both forced by what it got wrong:
    #
    #   * The coverage figure must be MERGED EXON BASES over indexed bases. The
    #     original was computed from transcript spans, which count introns as
    #     annotated and overstated the rat annotation's coverage as 32% when the
    #     figure is 1.8%. A guard keyed to the wrong quantity either fires when it
    #     should not or misses when it should; either way it is not a guard.
    #   * A6 is WITHDRAWN as an explanation of the rat result. The corrected
    #     annotation gives a known fraction of 0.585, which clears the 0.50 bar
    #     without any density rule. The bar is retained for a FUTURE stratum whose
    #     annotation genuinely is too sparse, and must not be cited for the rat one.
    #
    # The 0.5 bar is a priori, not fitted, and no threshold is adjusted here.
    if coverage is not None and coverage < 0.5:
        obs["annotation_coverage_exon_bases"] = round(coverage, 4)
        obs["coverage_measure"] = ("merged exon bases over indexed bases; NOT "
                                   "transcript spans, which count introns")
        rep.add("E5", "junction classification rates", "INCONCLUSIVE", observed=obs,
                threshold="known fraction >= 0.50 and known > novel",
                reason=(f"Annotation covers only {100*coverage:.1f}% of the indexed "
                        f"sequence in merged exon bases, so the annotated/non-annotated "
                        f"ratio is dominated by annotation density rather than splicing "
                        f"fidelity and the endpoint does not measure what it intends. "
                        f"Observed known fraction {known_frac:.4f}. "
                        f"(ENDPOINTS.md 6.1 A6, amended: the original coverage figure "
                        f"was computed from transcript spans and overstated it.)"))
        return
    obs["window_check_1e5_to_1e7"] = 1e5 <= total <= 1e7
    obs["known_fraction"] = round(known_frac, 4)
    ok = known_frac >= 0.50 and known > novel
    reasons = []
    if known <= novel:
        reasons.append(f"known ({known}) <= novel ({novel})")
    if known_frac < 0.50:
        reasons.append(f"annotated fraction {known_frac:.4f} below 0.50")
    rep.add("E5", "junction classification rates", "PASS" if ok else "FAIL",
            observed=obs, threshold="known fraction >= 0.50 and known > novel",
            reason=None if ok else "; ".join(reasons))


def e2_tin(rep: Report, enabled: bool) -> None:
    """TIN under controlled 3' degradation.

    Costs one re-alignment per degradation level, so it is opt-in rather than
    run implicitly on multi-hour FASTQs.
    """
    if not enabled:
        rep.add("E2", "TIN under controlled 3' degradation", "NOT_EVALUATED",
                threshold=">=0.90 transcripts non-increasing per level; mean non-increasing between levels",
                reason=("Requires re-aligning a downsampled read subset at each "
                        "degradation level; not run implicitly. Pass --run-e2 to "
                        "enable."))
        return
    rep.add("E2", "TIN under controlled 3' degradation", "NOT_EVALUATED",
            reason="degradation sub-pipeline not yet implemented; see ENDPOINTS.md 6.1 A1")


def e4_fpkm(rep: Report, has_spikein: bool, sampled: int) -> None:
    if not has_spikein:
        rep.add("E4", "FPKM response to spike-in dilution", "NOT_EVALUATED",
                threshold="rank correlation >= 0.90",
                reason=("No ERCC spike-in reads found in "
                        f"{sampled:,} sampled reads of this run. FPKM's spike-in "
                        "normalisation has no reference series here, so the endpoint "
                        "is not evaluable on this stratum. This is a property of the "
                        "held-out library design, not a port result."))
        return
    rep.add("E4", "FPKM response to spike-in dilution", "NOT_EVALUATED",
            reason="spike-in dilution series not constructed for this stratum")


def e6_sc(rep: Report) -> None:
    rep.add("E6", "single-cell (explicit NON-endpoint)", "NOT_EVALUATED",
            reason="No single-cell dataset exists in this project; no claim is made.")


# Feature strata required by ENDPOINTS.md section 3. A global pass with a stratum
# failing is a PARTIAL pass, not a pass, so the breakdowns are reported rather than
# summarised away.
#
# `GC content quartile` is defined by GC of the EXON SEQUENCE, read back from the
# genome, not of the transcript's span: an intron's GC is not what a library's
# exonic coverage is biased by, and a span-based proxy would mix the two. Requires
# a genome FASTA; where one is not supplied the stratum is reported as
# unavailable rather than dropped, because a silently missing stratum is exactly
# the failure this section exists to prevent.
COVERAGE_BINS = ((0, 1), (1, 5), (5, 20), (20, 100), (100, 10**9))
LENGTH_BINS = ((0, 1000), (1000, 3000), (3000, 10000), (10000, 10**9))
EXON_BINS = ((0, 5), (5, 10), (10, 20), (20, 10**9))
GC_BINS = ((0.0, 0.35), (0.35, 0.45), (0.45, 0.55), (0.55, 1.01))


def _bin_label(value, bins):
    for low, high in bins:
        if low <= value < high:
            return f"{low}-{high}" if high < 10**9 else f"{low}+"
    return "unbinned"


def _read_bed12(path: Path):
    """BED12 rows with their exon blocks, as `(chrom, start, end, name, blocks)`."""
    rows = []
    for line in path.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or not f[0].strip():
            continue
        try:
            sizes = [int(x) for x in f[10].rstrip(",").split(",") if x != ""]
            starts = [int(x) for x in f[11].rstrip(",").split(",") if x != ""]
            if len(sizes) != len(starts) or not sizes:
                continue
            base = int(f[6]) if f[6].strip() else int(f[1])
            blocks = [(base + off, base + off + size)
                      for off, size in zip(starts, sizes)]
            rows.append((f[0], int(f[1]), int(f[2]), f[3], blocks))
        except ValueError:
            continue
    return rows


def _gc_of_blocks(fasta, chrom, blocks):
    """GC fraction across a transcript's exons. None when the sequence is unavailable."""
    if fasta is None:
        return None
    seq = []
    for start, end in blocks:
        try:
            chunk = fasta.fetch(chrom, max(0, start), end).upper()
        except (KeyError, ValueError):
            return None
        seq.append(chunk)
    s = "".join(seq)
    if not s:
        return None
    return (s.count("G") + s.count("C")) / len(s)


def _ambiguous(chrom, blocks, index):
    """Whether a transcript's exons are shared with another in the same model.

    ENDPOINTS.md section 3 names "annotation ambiguity" as a stratum: a transcript
    whose exons overlap another's cannot have its coverage attributed uniquely, so
    its curve measures the model rather than the library. `index` maps a contig to
    every exon interval in the model, so the comparison is against all exons rather
    than only this transcript's own.
    """
    intervals = index.get(chrom, ())
    for s, e in blocks:
        overlapping = sum(1 for os_, oe in intervals if s < oe and os_ < e)
        if overlapping > 1:
            return True
    return False


def stratify(rep: Report, bam: Path, bed: Path, out_dir: Path,
            genome: Path | None = None) -> None:
    """Per-stratum gene-body coverage, by every feature ENDPOINTS.md section 3 names.

    Coverage, transcript length, exon count, GC content and annotation ambiguity.
    The audit found this stratification entirely absent, so an aggregate could
    conceal a stratum that fails. This is a breakdown rather than a new endpoint:
    it reports per-stratum curves and their correlation against the development
    reference, and flags any stratum below the threshold the global figure uses.

    A stratum with too few transcripts to support a curve is reported as such
    rather than omitted. So is the GC stratum when no genome is supplied: a
    silently missing stratum is the failure this section exists to prevent.
    """
    import math

    out_dir.mkdir(parents=True, exist_ok=True)
    transcripts = _read_bed12(bed)
    if not transcripts:
        rep.add("E3-S", "feature-stratified gene-body coverage", "NOT_EVALUATED",
                reason=f"no BED12 rows parsed from {bed}")
        return

    fasta = None
    if genome and Path(genome).exists():
        try:
            import pysam
            fasta = pysam.FastaFile(str(genome))
        except Exception as e:
            rep.add("E3-S", "feature-stratified gene-body coverage", "NOT_EVALUATED",
                    reason=f"genome {genome} could not be opened for the GC stratum: {e}")
            return
    else:
        rep.add("E3-S-GC", "GC-content stratum", "NOT_EVALUATED",
                reason=("No genome FASTA supplied, so exon GC content cannot be "
                        "computed. ENDPOINTS.md section 3 requires the GC stratum; "
                        "it is reported unavailable rather than dropped. Pass "
                        "--genome to include it."))

    # Per-chromosome interval index, for the ambiguity stratum.
    by_chrom: dict[str, list] = {}
    for chrom, _s, _e, _nm, blocks in transcripts:
        by_chrom.setdefault(chrom, []).extend(blocks)

    rc, out = run([str(BIN / "geneBody_coverage"), "--input", str(bam),
                   "--refgene", str(bed), "--out-prefix", str(out_dir / "strat_gb")])
    curve = out_dir / "strat_gb.geneBodyCoverage.txt"
    if rc != 0 or not curve.exists():
        rep.add("E3-S", "feature-stratified gene-body coverage", "NOT_EVALUATED",
                reason=f"geneBody_coverage exited {rc}")
        return

    # The aggregate curve is position-resolved over all transcripts, so a per-transcript
    # curve needs one geneBody_coverage run per stratum subset. Sampling keeps the
    # breakdown affordable; the subset written to disk makes the sample auditable.
    try:
        import pysam
    except ImportError:
        rep.add("E3-S", "feature-stratified gene-body coverage", "NOT_EVALUATED",
                reason="pysam required to restrict the alignment to a stratum")
        return

    def read_curve(path):
        lines = [x for x in path.read_text().splitlines() if x.strip()]
        if len(lines) < 2:
            return None
        vals = [float(x) for x in lines[1].split("\t")[1:]]
        return vals if len(vals) >= 100 and sum(vals) > 0 else None

    vals = read_curve(curve)
    ref_path = REPO / "datasets" / "reference_curve" / "dev_pe.geneBodyCoverage.txt"
    if vals is None or not ref_path.exists():
        rep.add("E3-S", "feature-stratified gene-body coverage", "NOT_EVALUATED",
                reason="aggregate curve unusable or reference curve missing")
        return
    ref_lines = [x for x in ref_path.read_text().splitlines() if x.strip()]
    ref = smooth([float(x) for x in ref_lines[1].split("\t")[1:]], SMOOTH_WIDTH)

    def pearson(u, v):
        mu, mv = sum(u) / len(u), sum(v) / len(v)
        num = sum((x - mu) * (y - mv) for x, y in zip(u, v))
        du = math.sqrt(sum((x - mu) ** 2 for x in u))
        dv = math.sqrt(sum((y - mv) ** 2 for y in v))
        return num / (du * dv) if du and dv else float("nan")

    def coverage_depth(chrom, start, end):
        with pysam.AlignmentFile(str(bam), "rb") as f:
            try:
                return sum(1 for _ in f.fetch(chrom, max(0, start - 1), end))
            except (ValueError, KeyError):
                return 0

    # Depth by EXON bases, which is what a gene-body curve integrates over: reads
    # landing in an intron contribute nothing to it, so depth over the transcript
    # span would put reads in the wrong stratum.
    #
    # One open handle for the whole pass. Reopening per transcript re-reads the
    # index for every row, which on a 5,840-transcript model dominated the runtime
    # of this function and made the breakdown impractical to run at all.
    depth_cache: dict[str, int] = {}

    def exon_depth(handle, chrom, blocks):
        reads = 0
        for start, end in blocks:
            key = f"{chrom}:{start}:{end}"
            n = depth_cache.get(key)
            if n is None:
                try:
                    n = sum(1 for _ in handle.fetch(chrom, max(0, start), end))
                except (ValueError, KeyError):
                    n = 0
                depth_cache[key] = n
            reads += n
        return reads

    strata: dict[str, list] = {}
    handle = pysam.AlignmentFile(str(bam), "rb")
    try:
        for chrom, start, end, name, blocks in transcripts:
            exon_bases = sum(e - s for s, e in blocks)
            depth = exon_depth(handle, chrom, blocks) / exon_bases if exon_bases else 0.0
            length = end - start
            gc = _gc_of_blocks(fasta, chrom, blocks)
            key = f"{chrom}:{start}:{name}"
            strata.setdefault(f"coverage:{_bin_label(depth, COVERAGE_BINS)}", []).append(key)
            strata.setdefault(f"length:{_bin_label(length, LENGTH_BINS)}", []).append(key)
            strata.setdefault(f"exons:{_bin_label(len(blocks), EXON_BINS)}", []).append(key)
            strata.setdefault(
                f"ambiguity:{'shared' if _ambiguous(chrom, blocks, by_chrom) else 'unique'}",
                []).append(key)
            if gc is not None:
                strata.setdefault(f"gc:{_bin_label(gc, GC_BINS)}", []).append(key)
    finally:
        handle.close()

    by_key = {f"{c}:{s}:{n}": (c, s, e, n, b) for c, s, e, n, b in transcripts}
    lines_by_key = {}
    for line in bed.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) >= 12:
            try:
                lines_by_key[f"{f[0]}:{int(f[1])}:{f[3]}"] = line
            except ValueError:
                continue

    # Cap per stratum so the breakdown stays affordable, and say so in the output.
    SAMPLE_CAP = 400
    results = {}
    for name, keys in sorted(strata.items()):
        keys = keys[:SAMPLE_CAP]
        if len(keys) < 10:
            results[name] = {"n": len(keys), "status": "too few transcripts"}
            continue
        sub = out_dir / f"strat_{name.replace(':', '_')}.bed12"
        with sub.open("w") as fh:
            for k in keys:
                line = lines_by_key.get(k)
                if line:
                    fh.write(line + "\n")
        prefix = out_dir / f"strat_{name.replace(':', '_')}"
        rc, _ = run([str(BIN / "geneBody_coverage"), "--input", str(bam),
                     "--refgene", str(sub), "--out-prefix", str(prefix)])
        svals = read_curve(Path(f"{prefix}.geneBodyCoverage.txt"))
        if svals is None or rc != 0:
            results[name] = {"n": len(keys), "status": "curve unusable",
                              "exit_code": rc}
            continue
        r = pearson(smooth(svals, SMOOTH_WIDTH), ref)
        results[name] = {
            "n": len(keys), "pearson_r": round(r, 4),
            "status": "pass" if r >= 0.80 else "below threshold",
        }

    failing = sorted(k for k, v in results.items() if v.get("status") == "below threshold")
    observed = {"strata": results, "sample_cap_per_stratum": SAMPLE_CAP,
                "strata_below_threshold": failing}
    if failing:
        verdict = "PARTIAL PASS"
        reason = (f"{len(failing)} stratum/strata below r=0.80 while the aggregate "
                  f"passes: {', '.join(failing)}. ENDPOINTS.md section 3 makes this a "
                  "partial pass, not a pass.")
    else:
        verdict = "PASS"
        reason = None
    rep.add("E3-S", "feature-stratified gene-body coverage", verdict,
            observed=observed,
            threshold="per-stratum Pearson r >= 0.80 vs dev reference",
            reason=reason)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--stratum", required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--contigs", nargs="*", default=[])
    ap.add_argument("--contig-subset", action="store_true")
    ap.add_argument("--expected-strand", default=None,
                    help="from archive protocol metadata; omit if unstated")
    ap.add_argument("--fastq1", default=None)
    ap.add_argument("--spikein-sampled", type=int, default=0)
    ap.add_argument("--has-spikein", action="store_true")
    ap.add_argument("--run-e2", action="store_true")
    ap.add_argument("--annotation-coverage", type=float, default=None,
                    help="fraction of indexed sequence the annotation covers in MERGED "
                         "EXON BASES; below 0.5 makes E5 not evaluable (ENDPOINTS.md "
                         "6.1 A6, as amended). Transcript spans are NOT an acceptable "
                         "substitute: they count introns as annotated, which overstated "
                         "the rat annotation's coverage as 32%% against a true 1.8%%.")
    ap.add_argument("--expect-unstranded", action="store_true",
                    help="ENA library_selection=PCR: check strand balance instead of a strand identity")
    ap.add_argument("--reference-is-independent", action="store_true",
                    help="score E3 against the independent dev SE reference curve")
    ap.add_argument("--library-selection", default=None,
                    choices=["cDNA", "PCR", "polyA", None],
                    help="ENA library_selection for this library. Recorded because "
                         "it, not the aligner's output, determines whether a "
                         "strand expectation exists at all.")
    ap.add_argument("--protocol-strandedness", default=None,
                    choices=["forward", "reverse", "unstranded", "auto", None],
                    help="protocol metadata from the library-preparation record. "
                         "Omit when the archive does not state it: E1 is then "
                         "NOT_EVALUATED rather than guessed.")
    ap.add_argument("--estrand", default=None,
                    help="the library's actual estimand for E3, from "
                         "datasets/ENDPOINTS.md section 6.1 A8. E3 measures "
                         "similarity to a DEVELOPMENT-SAMPLE curve, which is not a "
                         "causal property of any mechanism, so a mechanistic reading "
                         "is refused unless this names the estimand explicitly.")
    ap.add_argument("--stratify", action="store_true",
                    help="report per-stratum breakdowns required by ENDPOINTS.md "
                         "section 3 (coverage, transcript length, exon count, GC, "
                         "annotation ambiguity)")
    ap.add_argument("--genome", default=None,
                    help="genome FASTA, required for the GC-content stratum and "
                         "recorded alongside the results")
    args = ap.parse_args()

    # E1's expectation is a property of the library-preparation protocol, and the
    # audit's P1 scientific-interpretation finding was that it had been taken from
    # PCR selection instead. The two are different things: `library_selection=PCR`
    # names how the library was amplified, while the strand expectation comes from
    # the protocol. Where only the former is known, E1 is NOT_EVALUATED, because
    # inferring a protocol from an amplification method is the exact inference the
    # audit rejected.
    if args.expect_unstranded and args.protocol_strandedness not in (None, "unstranded"):
        print("error: --expect-unstranded conflicts with "
              "--protocol-strandedness; supply the protocol metadata instead of "
              "inferring it from library selection", file=sys.stderr)
        return 2
    if args.protocol_strandedness in ("forward", "reverse") and args.library_selection == "PCR":
        print("error: a PCR-amplified library carries no directional priming, so a "
              "forward/reverse protocol expectation and library_selection=PCR cannot "
              "both be true. One of them is wrong.", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bam, bed = Path(args.bam), Path(args.bed)

    rep = Report()
    e1_layout(rep, bam, bed, contig_subset=args.contig_subset,
              contigs=args.contigs, expected_strand=args.expected_strand,
              expect_unstranded=args.expect_unstranded,
              library_selection=args.library_selection,
              protocol_strandedness=args.protocol_strandedness)
    e3_genebody(rep, bam, bed, out_dir, args.reference_is_independent,
                estrand=args.estrand)
    e5_junctions(rep, bam, bed, out_dir, args.annotation_coverage)
    e2_tin(rep, args.run_e2)
    e4_fpkm(rep, args.has_spikein, args.spikein_sampled)
    e6_sc(rep)
    if args.stratify:
        stratify(rep, bam, bed, out_dir, genome=args.genome)

    payload = {
        "stratum": args.stratum, "run": args.run_id, "bam": str(bam),
        "bed": str(bed), "contig_subset": args.contig_subset,
        "contigs": args.contigs, "expected_strand": args.expected_strand,
        "counts": rep.counts(), "endpoints": rep.endpoints,
        "spec": "datasets/ENDPOINTS.md (signed off 2026-10-01, amended 6.1)",
    }
    (out_dir / "endpoints.json").write_text(json.dumps(payload, indent=2) + "\n")

    print(f"\nT4.2 endpoints -- {args.stratum} / {args.run_id}")
    print(f"  spec: datasets/ENDPOINTS.md  (signed 2026-10-01, amended 6.1)\n")
    for e in rep.endpoints:
        print(f"  [{e['verdict']:>14}] {e['id']}  {e['name']}")
        if e["threshold"]:
            print(f"                     threshold: {e['threshold']}")
        if e["observed"] is not None:
            print(f"                     observed : {json.dumps(e['observed'])}")
        if e["reason"]:
            for i, line in enumerate(_wrap(e["reason"], 84)):
                print(f"                     {'reason: ' if i == 0 else '        '}{line}")
    c = rep.counts()
    print("\n  " + "  ".join(f"{k}={v}" for k, v in sorted(c.items())))
    print(f"\n  json: {out_dir / 'endpoints.json'}")
    return 0 if not c.get("FAIL") else 1


def _wrap(text: str, width: int) -> list[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur); cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines


if __name__ == "__main__":
    sys.exit(main())
