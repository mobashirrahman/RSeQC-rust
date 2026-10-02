#!/usr/bin/env python3
"""Convert a UCSC `refGene.txt.gz` table into a GTF and a BED12.

Needed for the cross-organism held-out stratum. GENCODE's rat annotation was
not reachable from this network (both the release_47 and release_46 paths
404), and UCSC's rn6 `genes/` directory no longer serves a RefSeq GTF under any
name this project could resolve -- the same drift that forced GENCODE v47 for
human. `database/refGene.txt.gz` does resolve, so the table is converted
here instead.

UCSC refGene columns (per the UCSC schema):
    1 bin, 2 name, 3 chrom, 4 strand, 5 txStart, 6 txEnd, 7 cdsStart,
    8 cdsEnd, 9 exonCount, 10 exonStarts (0-based, comma-separated),
    11 exonEnds (1-based/exclusive, comma-separated), ...

Both outputs are emitted because two consumers need different things:
STAR's `--sjdbGTFfile` requires a GTF, and every RSeQC annotation-consuming
command requires BED12.

COORDINATE HANDLING -- this is the part worth being careful about, and this
module previously got it wrong twice.

UCSC refGene is *already* in the half-open 0-based genome-browser frame:
`txStart`/`exonStarts` are 0-based start offsets and `exonEnds`/`txEnd` are
exclusive end offsets. So a refGene exon is directly the half-open interval
[start, end) that BED12 stores verbatim:

  BED12 chromStart        = txStart                       (unchanged)
  BED12 exon size         = exonEnd - exonStart           (unchanged)
  BED12 exon start (rel)  = exonStart - txStart           (unchanged)
  GTF   exon start        = exonStart + 1                 (1-based inclusive)
  GTF   exon end          = exonEnd                       (already exclusive)

Sources: UCSC's documented counting systems and the refGene schema, which
defines exonStarts as 0-based and exonEnds as 1-based end offsets.

A worked probe (see datasets/test_refgene_to_gtf.py, which asserts it):

    refGene row  0  NM_TEST  chr1  +  100  400  100  400  2  100,300,  200,400,
    BED12         chr1 100 400 NM_TEST 0 + 100 400 255 2 100,100, 0,200,
    GTF exons     101-200 and 301-400

Two earlier revisions were wrong, and both errors were real defects rather
than harmless off-by-one noise:

  * Subtracting 1 from `txStart` while keeping the +1 exon-size rule treats the
    table as 1-based inclusive. That shifts the transcript one base left and
    makes every exon one base too long (start 99, sizes 101,101 above).
  * An earlier "+1 correction" that added 1 to the exon size to compensate for
    treating exonEnds as inclusive made every exon one base too long and, again,
    placed GTF starts one base early (100,300 above).

Either way the pair of errors did not cancel, because the BED and GTF consumers
disagree about inclusivity: a length computed with an inclusive end is wrong by
one base even when its start is right. The downstream symptom is junction
annotation reporting almost nothing as annotated (7 of 44,176 on the rat
stratum, 0.02%), which is what prompted the original investigation.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def parse(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            yield {
                "name": f[1],
                "chrom": f[2],
                "strand": f[3],
                "tx_start": int(f[4]),
                "tx_end": int(f[5]),
                "exon_starts": [int(x) for x in f[9].rstrip(",").split(",") if x],
                "exon_ends": [int(x) for x in f[10].rstrip(",").split(",") if x],
            }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refgene", required=True)
    ap.add_argument("--gtf-out", required=True)
    ap.add_argument("--bed-out", required=True)
    ap.add_argument("--contigs", nargs="+", default=None,
                    help="restrict to these contigs (the STAR index is built on a subset)")
    ap.add_argument("--min-exons", type=int, default=1)
    args = ap.parse_args()

    want = set(args.contigs) if args.contigs else None
    gtf = Path(args.gtf_out)
    bed = Path(args.bed_out)
    gtf.parent.mkdir(parents=True, exist_ok=True)

    n_gtf = n_bed = 0
    with gtf.open("w") as gf, bed.open("w") as bf:
        gf.write('##format: gtf\n')
        for t in parse(Path(args.refgene)):
            if want is not None and t["chrom"] not in want:
                continue
            if len(t["exon_starts"]) != len(t["exon_ends"]):
                continue
            if len(t["exon_starts"]) < args.min_exons:
                continue
            if any(e <= s for s, e in zip(t["exon_starts"], t["exon_ends"])):
                print(f"  skip {t['name']}: empty exon interval in refGene row",
                      file=sys.stderr)
                continue
            start0 = t["tx_start"]
            attr = f'gene_id "{t["name"]}"; transcript_id "{t["name"]}"; gene_name "{t["name"]}";'
            for s, e in zip(t["exon_starts"], t["exon_ends"]):
                # refGene exonStarts are 0-based, so GTF's 1-based inclusive start
                # is s + 1. exonEnds are already exclusive, which is exactly GTF's
                # inclusive end, so e passes through unchanged.
                gf.write(f'{t["chrom"]}\tsource\texon\t{s + 1}\t{e}\t.\t{t["strand"]}\t.\t{attr}\n')
                n_gtf += 1
            # refGene is already in the 0-based half-open frame BED12 uses, so the
            # exon size is the plain difference. See the module docstring.
            sizes = ",".join(str(e - s) for s, e in zip(t["exon_starts"], t["exon_ends"])) + ","
            starts = ",".join(str(s - t["tx_start"]) for s in t["exon_starts"]) + ","
            bf.write(
                f'{t["chrom"]}\t{start0}\t{t["tx_end"]}\t{t["name"]}\t0\t{t["strand"]}\t'
                f'{start0}\t{t["tx_end"]}\t255\t{len(t["exon_starts"])}\t{sizes}\t{starts}\n'
            )
            n_bed += 1

    print(f"gtf exon records : {n_gtf}  -> {gtf}")
    print(f"bed12 transcripts: {n_bed}  -> {bed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())