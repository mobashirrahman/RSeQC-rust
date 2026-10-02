#!/usr/bin/env python3
"""Verify that the rat reference preparation on disk is the one the manifest records.

The coordinate defect went undetected through two rounds of correction because the
*files* were trusted: a corrected annotation was written, a stale contig subset was
reused, and STAR reported success. Digests are the cheap check that catches the
recurrence, and this makes the manifest's recorded digests load-bearing rather than
decorative.

Also verifies the two claims that are easy to get wrong and hard to notice:

* the assembly is rn6 = Rnor_6.0, which is NOT mRatBN7.2 (UCSC serves mRatBN7.2 as rn7);
* every transcript's BED12 exon blocks reconstruct its own span exactly, which is the
  invariant an off-by-one frame error breaks.

Usage:
    python3 verification/check_rat_reference.py
    python3 verification/check_rat_reference.py --json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "datasets" / "manifest.yaml"
REF = REPO / "datasets" / "heldout" / "reference"
INDEX = REPO / "datasets" / "heldout" / "star_index"

# UCSC's assembly naming. rn6 is Rnor_6.0; mRatBN7.2 is rn7. Naming both at once
# would make every coordinate suspect, which is the defect class this checks for.
ASSEMBLY_FACTS = {
    "rn6": "Rnor_6.0",
    "rn7": "mRatBN7.2",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_manifest():
    """A minimal reader, so this check needs no YAML dependency.

    Only the fields this script needs are extracted, by locating them textually. A
    full parser would be better, but a check that cannot run without a dependency is
    a check that does not run in the environments that need it.
    """
    text = MANIFEST.read_text()
    out = {}
    for key in ("rn6_gtf_sha256", "rn6_bed12_sha256", "rn6_indexed_bed12_sha256"):
        marker = key + ":"
        i = text.find(marker)
        out[key] = text[i + len(marker):].split()[0] if i >= 0 else None
    i = text.find("index_annotation_sha256: >")
    out["index_annotation_sha256"] = None
    if i >= 0:
        j = text.find("--", i)
        if j > 0:
            out["index_annotation_sha256"] = text[i + len("index_annotation_sha256: >"):j].split()[0]
    # The rat stratum's assembly, not the first `assembly:` in the file -- the
    # human reference block comes first and says GRCh38, so a naive search reads the
    # wrong assembly and reports the rat naming as wrong.
    i = text.find("id: cross_organism")
    j = text.find("assembly: ", i) if i >= 0 else -1
    out["assembly"] = (text[j + len("assembly: "):].splitlines()[0].strip()
                       if j >= 0 else None)
    return out


def bed12_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        try:
            sizes = [int(x) for x in f[10].rstrip(",").split(",") if x != ""]
            starts = [int(x) for x in f[11].rstrip(",").split(",") if x != ""]
            rows.append((f[0], int(f[1]), int(f[2]), f[3], int(f[9]),
                         int(f[6]), int(f[7]), sizes, starts))
        except ValueError:
            continue
    return rows


def gtf_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 9:
            rows.append((f[0], f[2], int(f[3]), int(f[4]), f[6]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    results = []

    def check(name, ok, detail=""):
        results.append({"check": name, "pass": bool(ok), "detail": detail})
        print(f"  [{'ok  ' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
        return ok

    print("Digests recorded in datasets/manifest.yaml:")
    m = load_manifest()
    for key, fname in (("rn6_gtf_sha256", "rn6.gtf"),
                       ("rn6_bed12_sha256", "rn6.bed12"),
                       ("rn6_indexed_bed12_sha256", "rn6.indexed.bed12")):
        path = REF / fname
        if not path.exists():
            check(key, False, f"{path} missing")
            continue
        want = m.get(key)
        got = sha256(path)
        check(key, want == got,
              "matches the manifest" if want == got
              else f"manifest={want} on-disk={got}")

    # The junction database must be traceable to the SOURCE annotation, not only to
    # a copy that sits beside it. The index directory holds `annotation.gtf`, which
    # is the unpacked subset; a reader who digests that file proves only that the
    # copy is stable, not that it came from the corrected conversion. The subset is
    # therefore required to be exactly the exon records of the source GTF restricted
    # to the contigs the stamp names -- which is what the build script actually
    # does, so this asserts the property rather than a digest of a derived file.
    stamp_contigs = set()
    stamp_text = ""
    for candidate in sorted(INDEX.glob(".built-*")):
        stamp_text = candidate.read_text()
        for line in stamp_text.splitlines():
            if line.startswith("contigs="):
                stamp_contigs = set(line.split("=", 1)[1].split())
    source_gtf = REF / "rn6.gtf"
    index_gtf = INDEX / "annotation.subset.gtf"
    if source_gtf.exists() and index_gtf.exists() and stamp_contigs:
        def exon_set(path):
            out = set()
            for line in path.read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split("\t")
                if len(f) >= 9 and f[2] == "exon":
                    out.add((f[0], int(f[3]), int(f[4]), f[6], f[8]))
            return out

        source = exon_set(source_gtf)
        expected = {e for e in source if e[0] in stamp_contigs}
        actual = exon_set(index_gtf)
        check("index junction database is the corrected source annotation, "
              "restricted to the stamp's contigs",
              actual == expected and bool(expected),
              f"{len(actual)} index exon records, {len(expected)} expected from "
              f"{len(source)} source records over {sorted(stamp_contigs)}"
              if actual == expected else
              f"{len(actual - expected)} index records not in the source, "
              f"{len(expected - actual)} source records missing from the index")
    else:
        check("index junction database traceable to the source annotation", False,
              f"{source_gtf} or {index_gtf} missing, or the stamp names no contigs")

    idx_gtf = INDEX / "annotation.gtf"
    if idx_gtf.exists():
        want = m.get("index_annotation_sha256")
        got = sha256(idx_gtf)
        check("index_annotation_sha256", want == got,
              "the STAR junction database derives from the recorded annotation"
              if want == got else f"manifest={want} on-disk={got}")
        stamp = next((INDEX / n for n in ())
                     , None)
        for candidate in INDEX.glob(".built-*"):
            recorded = candidate.read_text()
            check("index stamp records the annotation digest",
                  got in recorded,
                  f"{candidate.name}")
    else:
        check("index annotation present", False, f"{idx_gtf} missing")

    print("\nAssembly naming:")
    assembly = m.get("assembly") or ""
    check("assembly is rn6 = Rnor_6.0",
          "rn6" in assembly and "Rnor_6.0" in assembly,
          assembly)
    check("assembly is NOT named as mRatBN7.2",
          "mRatBN" not in assembly,
          "mRatBN7.2 is rn7 in UCSC, a different assembly; naming it here would "
          "confuse two assemblies"
          if "mRatBN" in assembly else "")

    # The alignment must have been produced against THIS index. This is the check
    # that would have caught the state this file was written into: the STAR index
    # was rebuilt from the corrected annotation at 21:23, while the rat BAM on disk
    # was aligned at 11:30 against the previous index. Every digest above passed,
    # because the annotation on disk really was the corrected one -- what was stale
    # was the alignment. A re-run of the endpoint commands then reproduced the
    # earlier 0.585 figure while the recorded table said 0.617, and nothing in the
    # repository could say which was current.
    print("\nAlignment is bound to this index:")
    align_dir = REPO / "datasets" / "heldout" / "aligned" / "SRR1177982"
    align_json = align_dir / "SRR1177982.align.json"
    bam = align_dir / "SRR1177982.bam"
    if not align_json.exists() or not bam.exists():
        check("rat alignment present", False, f"{align_json} or {bam} missing")
    else:
        stamp_digest = None
        for line in stamp_text.splitlines():
            if line.startswith("subset_gtf_sha256="):
                stamp_digest = line.split("=", 1)[1].strip()
        try:
            record = json.loads(align_json.read_text())
        except json.JSONDecodeError as exc:
            record = {}
            check("align.json parses", False, str(exc))
        used = str(record.get("index_stamp", ""))
        check("the alignment records which annotation digest it used",
              bool(used) and "gtf_sha256=" in used,
              "the recorded index_stamp must name the annotation digest, otherwise "
              "'was this BAM aligned against this index?' is unanswerable"
              if used else "align.json has no index_stamp")
        if used and "gtf_sha256=" in used and stamp_digest:
            # Parse `key=value` pairs rather than splitting on whitespace: the stamp
            # contains a space-separated `contigs=` value, so `used.split()` yields
            # bare tokens like "chr2" alongside the real keys. A digest happens not to
            # contain a space so the comparison still worked, but the parse was
            # reading structure it did not understand.
            recorded = ""
            for token in used.split():
                if token.startswith("gtf_sha256="):
                    recorded = token.split("=", 1)[1]
            if recorded and not re.fullmatch(r"[0-9a-f]{64}", recorded):
                check("the alignment's recorded annotation digest is well formed",
                      False, f"gtf_sha256={recorded!r} is not a 64-character digest")
            check("the alignment used the index currently on disk",
                  recorded == stamp_digest,
                  f"align.json gtf_sha256={recorded[:16]}... "
                  f"index stamp gtf_sha256={stamp_digest[:16]}..."
                  if recorded != stamp_digest else
                  f"both {stamp_digest[:16]}...")
        if bam.exists() and align_json.exists() and stamp_digest:
            bam_mtime = bam.stat().st_mtime
            stamps = sorted(INDEX.glob(".built-*"))
            if stamps:
                stamp_mtime = max(s.stat().st_mtime for s in stamps)
                check("the BAM is newer than the index it claims to use",
                      bam_mtime >= stamp_mtime,
                      f"BAM {bam.stat().st_size} bytes, mtime {bam_mtime:.0f} vs "
                      f"index stamp {stamp_mtime:.0f}"
                      if bam_mtime < stamp_mtime else "")

    # The endpoint figures must agree across the invocations that produced them, and
    # the junction table must be the one both arms wrote. This is the check that makes
    # "reproduced rather than restated" an assertion rather than a claim: the
    # superseded endpoints.json next door carries 0.026 and stays on disk, so the only
    # way to know which JSON is current is to compare the artifact that both runs
    # must have produced.
    print("\nEndpoint results agree with the same-input comparison:")
    results_root = REPO / "datasets" / "heldout" / "endpoint_results"
    same_input = results_root / "SRR1177982_same_input" / "same-input-comparison.json"
    corrected = results_root / "SRR1177982_corrected"
    if not same_input.exists() or not corrected.exists():
        check("endpoint results present", False,
              f"{same_input} or {corrected} missing; the rat endpoints have not been "
              f"revalidated on the index-bound alignment")
    else:
        cmp_doc = json.loads(same_input.read_text())
        rates = cmp_doc.get("e5_rates", {})
        arms = cmp_doc.get("arms", {})
        up = arms.get("upstream", {}).get("junction_xls_sha256", "")
        rs = arms.get("rust", {}).get("junction_xls_sha256", "")
        check("both arms produced the same junction table", bool(up) and up == rs,
              f"upstream={up[:16]}... rust={rs[:16]}..." if up == rs
              else f"upstream={up[:16]}... rust={rs[:16]}... -- the arms disagree")
        for label, path in (("upstream", same_input.parent / "upstream.junction.xls"),
                            ("rust", same_input.parent / "rust.junction.xls"),
                            ("validator", corrected / "junc.junction.xls")):
            check(f"{label} junction table on disk matches the recorded digest",
                  path.exists() and sha256(path) == up,
                  f"{path.name} {'matches' if path.exists() and sha256(path) == up else 'differs'}")
        recorded_bam = cmp_doc.get("inputs", {}).get("bam_sha256", "")
        check("the comparison was run on the BAM currently on disk",
              bool(recorded_bam) and bam.exists() and sha256(bam) == recorded_bam,
              f"comparison={recorded_bam[:16]}... disk={sha256(bam)[:16]}..."
              if bam.exists() and recorded_bam else "")
        ej = corrected / "endpoints.json"
        if ej.exists():
            doc = json.loads(ej.read_text())
            e5 = next((e for e in doc.get("endpoints", []) if e["id"] == "E5"), {})
            obs = e5.get("observed", {})
            check("E5's total junctions match the same-input run",
                  obs.get("total") == rates.get("total_splicing_junctions"),
                  f"endpoints.json={obs.get('total')} "
                  f"same-input={rates.get('total_splicing_junctions')}")
            check("E5's annotated junctions match the same-input run",
                  obs.get("known") == rates.get("annotated_splicing_junctions"),
                  f"endpoints.json={obs.get('known')} "
                  f"same-input={rates.get('annotated_splicing_junctions')}")
        else:
            check("corrected endpoints.json present", False, f"{ej} missing")

    print("\nBED12 structure (the invariant an off-by-one frame error breaks):")
    for fname in ("rn6.bed12", "rn6.indexed.bed12"):
        path = REF / fname
        if not path.exists():
            check(f"{fname} present", False)
            continue
        rows = bed12_rows(path)
        if not rows:
            check(f"{fname} parses", False, "no rows")
            continue
        bad_count = bad_span = bad_order = bad_overlap = 0
        for _chrom, start, end, name, nblocks, tstart, tend, sizes, starts in rows:
            if len(sizes) != len(starts) or len(sizes) != nblocks:
                bad_count += 1
                continue
            cursor = 0
            ok = True
            for size, rel in zip(sizes, starts):
                if rel < cursor:
                    bad_order += 1
                    ok = False
                    break
                if rel + size > tend - tstart or start + rel < start:
                    bad_span += 1
                    ok = False
                    break
                cursor = rel + size
            # Blocks plus the trailing gap must consume the transcript exactly.
            if ok and tstart + cursor != tend:
                bad_span += 1
        check(f"{fname}: block counts match the column", bad_count == 0,
              f"{bad_count} rows" if bad_count else f"{len(rows)} transcripts")
        check(f"{fname}: blocks reconstruct each transcript's span", bad_span == 0,
              f"{bad_span} rows" if bad_span else "")
        check(f"{fname}: blocks ordered and non-overlapping", bad_order == 0,
              f"{bad_order} rows" if bad_order else "")
        check(f"{fname}: no zero-length exon block",
              all(s > 0 for r in rows for s in r[7]),
              "BED12 cannot represent an empty block")

    print("\nBED12 and GTF describe the same intervals:")
    bed = REF / "rn6.indexed.bed12"
    gtf = INDEX / "annotation.gtf"
    if bed.exists() and gtf.exists():
        bed_iv = {}
        for chrom, _s, _e, name, _n, tstart, _tend, sizes, starts in bed12_rows(bed):
            for size, rel in zip(sizes, starts):
                bed_iv.setdefault(chrom, set()).add((tstart + rel, tstart + rel + size))
        gtf_iv = {}
        for chrom, feat, s, e, _strand in gtf_rows(gtf):
            if feat == "exon":
                # GTF is 1-based inclusive; the half-open interval is [s-1, e).
                gtf_iv.setdefault(chrom, set()).add((s - 1, e))
        common = set(bed_iv) & set(gtf_iv)
        mismatched = [c for c in common if bed_iv[c] != gtf_iv[c]]
        check("BED12 and GTF exon intervals agree per contig",
              bool(common) and not mismatched,
              f"{len(common)} shared contigs"
              + (f", mismatched: {mismatched[:3]}" if mismatched else ""))
    else:
        check("BED12 and GTF both present", False,
              f"{bed} or {gtf} missing")

    # ------------------------------------------------------------------
    # Whole-genome stratum. The three-contig panel above turned out to retain only
    # 20.31% of uniquely mapped reads, because most rat reads fall outside
    # chr1/chr2/chr10 -- so a contig subset is a biased sample, not a smaller sample of
    # the same thing. The whole-genome alignment is now the primary substrate, and it
    # gets the same provenance checks plus one more: that the recorded byte-identity is
    # actually true of the files on disk. A stored hash of a claim is not the claim.
    # ------------------------------------------------------------------
    wg_index = Path("datasets/heldout/star_index/rn6-wholegenome-o100")
    wg_align = Path("datasets/heldout/aligned_wholegenome/SRR1177982")
    wg_evidence = Path("datasets/heldout/endpoint_results/SRR1177982_wholegenome")

    if wg_index.is_dir() and wg_align.is_dir():
        stamps = sorted(wg_index.glob(".built-*"))
        check("whole-genome index has a build record", bool(stamps),
              stamps[0].name if stamps else "no .built-* marker")
        stamp_text = stamps[0].read_text() if stamps else ""
        contig_field = next((l.split("=", 1)[1] for l in stamp_text.splitlines()
                             if l.startswith("contigs=")), "")
        check("whole-genome index records all 58 contigs",
              len(contig_field.split()) == 58,
              f"{len(contig_field.split())} contigs recorded")
        # The stamp's genome_sha256 is the digest of the UNPACKED genome.fa STAR
        # actually read, not of the rn6.fa.gz a user would re-fetch -- comparing it to
        # the .gz fails for a reason that has nothing to do with the index. So the
        # unpacked form is checked here, and the input form is checked against the
        # manifest above, which is where a re-fetchable digest belongs.
        for field, path in (("genome_sha256", wg_index / "genome.fa"),
                            ("gtf_sha256", REF / "rn6.gtf")):
            want = next((l.split("=", 1)[1] for l in stamp_text.splitlines()
                         if l.startswith(field + "=")), "")
            got = sha256(path) if path.is_file() else "(missing)"
            check(f"whole-genome index digest for {path.name} matches the stamp",
                  want == got, "matches" if want == got else f"stamp={want[:16]} disk={got[:16]}")

        ajson = wg_align / "SRR1177982.align.json"
        if ajson.is_file():
            aj = json.loads(ajson.read_text())
            check("whole-genome alignment records the index it was built against",
                  Path(aj.get("index_dir", "")).name == wg_index.name,
                  aj.get("index_dir", "(none)"))
            # align_run.sh flattens the marker into one space-separated line, so the
            # two are compared as PARSED key/value sets. A raw string comparison here
            # would fail on formatting alone, which is how a real binding difference
            # would get lost in a false positive.
            def parse_stamp(s):
                out = {}
                for part in s.split():
                    k, _, v = part.partition("=")
                    if v:
                        out[k] = v
                return out

            stamp_fields = parse_stamp(stamp_text)
            align_fields = parse_stamp(aj.get("index_stamp", ""))
            shared = set(stamp_fields) & set(align_fields)
            differing = sorted(k for k in shared if stamp_fields[k] != align_fields[k])
            check("whole-genome alignment is bound to this exact index build",
                  bool(shared) and not differing
                  and set(align_fields) == set(stamp_fields),
                  f"{len(align_fields)} fields compared"
                  + (f", DIFFERING: {differing}" if differing else "")
                  + (f", only in one side: "
                     f"{sorted(set(align_fields) ^ set(stamp_fields))}"
                     if set(align_fields) != set(stamp_fields) else ""))
            check("whole-genome alignment is indexed and coordinate-sorted",
                  (wg_align / "SRR1177982.bam.bai").is_file())

        if wg_evidence.is_dir() and (wg_evidence / "evidence.json").is_file():
            ev = json.loads((wg_evidence / "evidence.json").read_text())
            for name, rec in ev.get("byte_identical_artifacts", {}).items():
                up = wg_evidence / f"upstream.{name}"
                rs = wg_evidence / f"rust.{name}"
                if not (up.is_file() and rs.is_file()):
                    check(f"whole-genome {name} artifacts present", False,
                          "one side missing")
                    continue
                a, b = sha256(up), sha256(rs)
                check(f"whole-genome {name} is byte-identical as recorded", a == b,
                      f"{up.stat().st_size:,} bytes, sha {a[:16]}")
                check(f"whole-genome {name} digest matches the evidence record",
                      rec.get("upstream_sha256") == a,
                      "the recorded hash is the file's hash")
            # The scientific claim, checked rather than quoted: both arms must report
            # STAR's own spliced-junction count, so neither is inventing or dropping
            # junctions relative to the aligner.
            star_total = None
            final = wg_align / "SRR1177982.Log.final.out"
            if final.is_file():
                for line in final.read_text().splitlines():
                    if "Number of splices: Total" in line:
                        star_total = int(line.split("|")[-1].strip())
            claimed = ev.get("results", {}).get("upstream_total")
            check("both arms' junction total equals STAR's spliced-junction count",
                  star_total is not None and claimed == star_total
                  and ev.get("results", {}).get("rust_total") == star_total,
                  f"STAR={star_total} upstream={claimed} rust={ev.get('results',{}).get('rust_total')}")
    else:
        check("whole-genome rat stratum present", False,
              f"{wg_index} or {wg_align} missing")

    failed = [r for r in results if not r["pass"]]
    print()
    if args.json:
        print(json.dumps({"results": results, "failed": len(failed)}, indent=2))
    print(f"{len(results) - len(failed)}/{len(results)} rat reference checks passed")
    if failed:
        print("\nA mismatch here means the annotation on disk is not the one the "
              "manifest\nrecords. Either the files are stale or the manifest is "
              "wrong; both need\nresolving before any annotation-based result is "
              "quoted.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
