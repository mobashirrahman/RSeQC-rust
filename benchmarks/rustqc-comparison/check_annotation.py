#!/usr/bin/env python3
"""Card R2: compare BED12 and GTF transcript sets.

Compares the set of transcript IDs, and for each shared transcript its
chromosome, strand, start, end and exon blocks (BED is 0-based half-open,
GTF is 1-based inclusive).

Prints counts of: only in BED, only in GTF, shared and identical, shared
but different. Exits non-zero unless both "only" counts and the
"different" count are zero.

Usage:
    python3 check_annotation.py --bed <bed> --gtf <gtf> [--examples N]
    python3 check_annotation.py --workloads workloads.json [--examples N]
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_bed12(path: Path) -> dict:
    """Return {(chrom, transcript_id): (strand, start, end, ((s,e),...))}.

    Keyed by (chrom, transcript_id) like `load_transcripts` in
    benchmarks/generate_workload_real.py, not by bare transcript_id: the rat
    refGene GTF reuses one transcript_id on two contigs (e.g. NM_001000027
    on chr10 and chr10_KL568017v1_random), and BED12 itself repeats one name
    on two loci (e.g. NM_001000013 twice on chr10). All blocks sharing one
    (chrom, name) key are unioned, mirroring how the GTF loader merges exons
    of one (chrom, transcript_id) key.
    """
    blocks_by_key: dict[tuple, list] = {}
    strand_by_key: dict[tuple, str] = {}
    warned: set = set()
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                print(f"warning: {path}:{lineno}: only {len(cols)} columns, skipped",
                      file=sys.stderr)
                continue
            chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]
            strand = cols[5]
            nblock = int(cols[9])
            sizes = [int(x) for x in cols[10].rstrip(",").split(",") if x != ""]
            starts = [int(x) for x in cols[11].rstrip(",").split(",") if x != ""]
            if len(sizes) != nblock or len(starts) != nblock:
                print(f"warning: {path}:{lineno}: block count mismatch, skipped",
                      file=sys.stderr)
                continue
            key = (chrom, name)
            for s, z in zip(starts, sizes):
                blocks_by_key.setdefault(key, []).append((start + s, start + s + z))
            if key in strand_by_key and strand_by_key[key] != strand:
                if key not in warned:
                    print(f"warning: {path}:{lineno}: strand differs within {key}, kept first",
                          file=sys.stderr)
                    warned.add(key)
            else:
                strand_by_key.setdefault(key, strand)
    out = {}
    for key, blocks in blocks_by_key.items():
        blocks = sorted(set(blocks))
        # Merge overlapping blocks the same way the GTF loader does, so a
        # split BED representation of one locus still compares equal.
        merged = [blocks[0]]
        for s, e in blocks[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        start = min(s for s, _ in merged)
        end = max(e for _, e in merged)
        out[key] = (strand_by_key[key], start, end, tuple(merged))
    return out


def parse_gtf(path: Path) -> dict:
    """Return {(chrom, transcript_id): (strand, start, end, ((s,e),...))}.

    Keyed by (chrom, transcript_id) like `load_transcripts`: one
    transcript_id reused on two contigs yields two keys rather than one
    merged cross-contig transcript. Exon rows define the locus; transcripts
    with no exon rows fall back to their `transcript` row.
    """
    exons: dict[tuple, list] = {}
    strand_by_key: dict[tuple, str] = {}
    meta: dict[tuple, tuple] = {}
    warned: set = set()
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, feature, s, e, strand = cols[0], cols[2], cols[3], cols[4], cols[6]
            m = TX_RE.search(cols[8])
            if not m:
                continue
            tx = m.group(1)
            key = (chrom, tx)
            try:
                si, ei = int(s) - 1, int(e)  # to 0-based half-open
            except ValueError:
                continue
            if feature == "exon":
                exons.setdefault(key, []).append((si, ei))
                if key in strand_by_key and strand_by_key[key] != strand:
                    if key not in warned:
                        print(f"warning: {path}: strand differs within {key}, kept first",
                              file=sys.stderr)
                        warned.add(key)
                else:
                    strand_by_key.setdefault(key, strand)
            elif feature == "transcript" and key not in meta:
                meta[key] = (strand, si, ei)
    out = {}
    for key, blocks in exons.items():
        blocks = sorted(set(blocks))
        merged = [blocks[0]]
        for s, e in blocks[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        start = min(s for s, _ in merged)
        end = max(e for _, e in merged)
        out[key] = (strand_by_key[key], start, end, tuple(merged))
    for key, (strand, si, ei) in meta.items():
        if key not in out:
            out[key] = (strand, si, ei, ((si, ei),))
    return out


def compare(bed: dict, gtf: dict):
    bed_ids = set(bed)
    gtf_ids = set(gtf)
    only_bed = sorted(bed_ids - gtf_ids)
    only_gtf = sorted(gtf_ids - bed_ids)
    identical = []
    different = []  # (key, reason)
    for key in sorted(bed_ids & gtf_ids):
        if bed[key] == gtf[key]:
            identical.append(key)
        else:
            b, g = bed[key], gtf[key]
            reasons = []
            if b[0] != g[0]:
                reasons.append(f"strand {b[0]} vs {g[0]}")
            if b[1] != g[1]:
                reasons.append(f"start {b[1]} vs {g[1]}")
            if b[2] != g[2]:
                reasons.append(f"end {b[2]} vs {g[2]}")
            if b[3] != g[3]:
                reasons.append(f"exons {len(b[3])} blocks vs {len(g[3])} blocks")
            different.append((key, "; ".join(reasons) or "unknown"))
    return only_bed, only_gtf, identical, different


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default=None)
    ap.add_argument("--gtf", default=None)
    ap.add_argument("--workloads", default=None,
                    help="workloads.json; checks every workload's bed/gtf pair")
    ap.add_argument("--examples", type=int, default=5)
    args = ap.parse_args()

    pairs = []
    if args.workloads:
        import json
        data = json.loads(Path(args.workloads).read_text())
        for wid, w in data["workloads"].items():
            pairs.append((wid, Path(w["bed"]), Path(w["gtf"])))
    elif args.bed and args.gtf:
        pairs.append(("cli", Path(args.bed), Path(args.gtf)))
    else:
        ap.error("give --bed + --gtf, or --workloads")
        return 2

    overall_ok = True
    for label, bed_path, gtf_path in pairs:
        bed = parse_bed12(bed_path)
        gtf = parse_gtf(gtf_path)
        only_bed, only_gtf, identical, different = compare(bed, gtf)
        print(f"[{label}] bed={bed_path} ({len(bed)} (chrom, transcript) keys)")
        print(f"[{label}] gtf={gtf_path} ({len(gtf)} (chrom, transcript) keys)")
        print(f"[{label}] only in BED: {len(only_bed)}")
        print(f"[{label}] only in GTF: {len(only_gtf)}")
        print(f"[{label}] shared and identical: {len(identical)}")
        print(f"[{label}] shared but different: {len(different)}")
        for key in only_bed[: args.examples]:
            print(f"[{label}]   BED-only example: {key}")
        for key in only_gtf[: args.examples]:
            print(f"[{label}]   GTF-only example: {key}")
        for key, reason in different[: args.examples]:
            print(f"[{label}]   different example: {key} ({reason})")
        ok = (len(only_bed) == 0 and len(only_gtf) == 0 and len(different) == 0)
        print(f"[{label}] {'PASS' if ok else 'MISMATCH'}")
        overall_ok = overall_ok and ok
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
