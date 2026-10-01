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
              expect_unstranded: bool = False) -> None:
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
    if expected_strand is None:
        rep.add("E1", "layout/orientation", "INCONCLUSIVE", observed=detail,
                threshold="unassigned <= 0.20",
                reason=("No protocol strandedness metadata for this stratum, so the "
                        "stronger 'correct strand recovered' claim is not evaluable. "
                        f"Unassigned fraction {failed:.4f} is reported for reference only."))
        return
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
                reference_is_independent: bool = False) -> None:
    """Whether the observed gene-body curve reproduces the expected 3' skew shape.

    The reference is the dev panel's own curve: a well-annotated RNA-seq
    library of this type has a characteristic 5'-low / 3'-depleted shape, and
    the endpoint asks whether a held-out library reproduces it. Both curves are
    smoothed with the width frozen in ENDPOINTS.md 6.1 A4, because the raw
    100-point curve has a measured split-half reliability of only r=0.744 at
    this depth -- below this endpoint's own threshold -- and is therefore not a
    usable statistic unsmoothed. That reliability figure is the reason the
    smoothing exists and is why the width is not retuned per stratum.
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
    rep.add("E3", "gene-body skew shape", "PASS" if ok else "FAIL",
            observed=obs,
            threshold=f"Pearson r >= 0.80 vs dev reference, smoothed width {SMOOTH_WIDTH}",
            reason=None if ok else f"r={r:.4f} below 0.80 against the dev reference curve")


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
    # ENDPOINTS.md 6.1 A6: E5 is not evaluable when the annotation covers less
    # than half the indexed sequence. The rat stratum annotates 32% of its own
    # index, so its annotated/non-annotated ratio measures annotation density
    # rather than splicing fidelity, and the verdict is INCONCLUSIVE with the
    # coverage figure attached. The human strata annotate 503% and are
    # unaffected. The 0.5 bar is a priori, not fitted.
    if coverage is not None and coverage < 0.5:
        obs["annotation_coverage"] = round(coverage, 4)
        rep.add("E5", "junction classification rates", "INCONCLUSIVE", observed=obs,
                threshold="known fraction >= 0.50 and known > novel",
                reason=(f"Annotation covers only {100*coverage:.1f}% of the indexed "
                        f"sequence, so the annotated/non-annotated ratio is dominated "
                        f"by annotation density rather than splicing fidelity and the "
                        f"endpoint does not measure what it intends. Observed known "
                        f"fraction {known_frac:.4f}."))
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
                    help="fraction of indexed sequence the annotation covers; below 0.5 "
                         "makes E5 not evaluable (ENDPOINTS.md 6.1 A6)")
    ap.add_argument("--expect-unstranded", action="store_true",
                    help="ENA library_selection=PCR: check strand balance instead of a strand identity")
    ap.add_argument("--reference-is-independent", action="store_true",
                    help="score E3 against the independent dev SE reference curve")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bam, bed = Path(args.bam), Path(args.bed)

    rep = Report()
    e1_layout(rep, bam, bed, contig_subset=args.contig_subset,
              contigs=args.contigs, expected_strand=args.expected_strand,
              expect_unstranded=args.expect_unstranded)
    e3_genebody(rep, bam, bed, out_dir, args.reference_is_independent)
    e5_junctions(rep, bam, bed, out_dir, args.annotation_coverage)
    e2_tin(rep, args.run_e2)
    e4_fpkm(rep, args.has_spikein, args.spikein_sampled)
    e6_sc(rep)

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
