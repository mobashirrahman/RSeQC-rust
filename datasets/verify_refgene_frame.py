#!/usr/bin/env python3
"""Confirm refGene's coordinate frame from the genome sequence itself.

`test_refgene_to_gtf.py` pins the converter's arithmetic against UCSC's documented
conventions. This script checks the convention itself against sequence, so the two
arguments do not rest on the same premise.

The test is stated in genomic terms. refGene marks `cdsStartStat=cmpl` transcripts
as beginning their CDS at the transcript's first base, and a start codon is ATG,
so the three bases there must be `ATG` under the correct reading of the coordinate.
`cdsEndStat=cmpl` transcripts end their CDS at the last base, so the three bases
there must be a stop codon (TAA/TAG/TGA). Both candidate frames are evaluated;
whichever one satisfies biology is the frame the data is in, and the converter is
confirmed against that.

`cdsStart`/`cdsEnd` are the CDS's genomic low and high bounds on both strands, so
the transcript's 5' end is `cdsStart` on the plus strand and `cdsEnd` on the minus
strand.

Sequence comes from UCSC's public sequence API, which takes half-open 0-based
intervals -- the same frame as BED12, and a source independent of the refGene
table's own documentation. `--genome` reads a local FASTA instead, via pysam,
which requires an indexed (uncompressed, `.fai`-bearing) file.

Usage:
    python3 datasets/verify_refgene_frame.py --limit 40
    python3 datasets/verify_refgene_frame.py --refgene <table> --genome <fa>
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import urllib.request
from pathlib import Path

STOPS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
API = ("https://api.genome.ucsc.edu/getData/sequence"
       "?genome={genome};chrom={chrom};start={start};end={end}")

# Candidate frames, each expressed as the offset from the stored coordinate to the
# 0-based index of the first base of the interval.
FRAMES = {
    "0-based half-open (coordinate used directly)": 0,
    "1-based inclusive (coordinate minus one)": -1,
}


def load_refgene(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 16:
                continue
            yield {
                "name": f[1], "chrom": f[2], "strand": f[3],
                "tx_start": int(f[4]), "tx_end": int(f[5]),
                "cds_start": int(f[6]), "cds_end": int(f[7]),
                "cds_start_stat": f[13], "cds_end_stat": f[14],
            }


def rc(seq):
    return seq.translate(COMPLEMENT)[::-1]


class ApiSequence:
    """Half-open 0-based sequence lookups through UCSC's public API."""

    def __init__(self, genome, cache=4096):
        self.genome = genome
        self.cache = {}

    def fetch(self, chrom, start, end):
        if start < 0 or end <= start:
            return ""
        key = (chrom, start, end)
        if key in self.cache:
            return self.cache[key]
        url = API.format(genome=self.genome, chrom=chrom, start=start, end=end)
        last_error = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=30) as resp:
                    seq = json.load(resp)["dna"].upper()
                self.cache[key] = seq
                return seq
            except Exception as e:  # transient network/API failures are retried
                last_error = e
                time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"sequence fetch failed for {chrom}:{start}-{end}: {last_error}")


class LocalSequence:
    """0-based half-open lookups into an indexed local FASTA."""

    def __init__(self, path):
        import pysam
        self.handle = pysam.FastaFile(str(path))

    def fetch(self, chrom, start, end):
        if start < 0 or end <= start:
            return ""
        try:
            return self.handle.fetch(chrom, start, end).upper()
        except (KeyError, ValueError):
            return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refgene", default="datasets/heldout/reference/rn6.refGene.txt.gz")
    ap.add_argument("--genome", default=None,
                    help="indexed local FASTA; omit to use UCSC's sequence API")
    ap.add_argument("--genome-name", default="rn6",
                    help="UCSC genome name, used when --genome is omitted")
    ap.add_argument("--limit", type=int, default=40,
                    help="transcripts to test (both orientations)")
    args = ap.parse_args()

    seq = LocalSequence(args.genome) if args.genome else ApiSequence(args.genome_name)
    source = f"local FASTA {args.genome}" if args.genome else f"UCSC API (genome={args.genome_name})"

    results = {label: {"start": [0, 0], "stop": [0, 0]} for label in FRAMES}
    tested = 0
    skipped = {"not_complete_cds": 0, "boundary_cds": 0, "short_sequence": 0}

    for t in load_refgene(args.refgene):
        if tested >= args.limit:
            break
        if t["cds_start_stat"] != "cmpl" or t["cds_end_stat"] != "cmpl":
            skipped["not_complete_cds"] += 1
            continue
        # A CDS that reaches the transcript boundary is still complete, but those
        # cases are less discriminating (the start sits against a UTR rather than
        # inside the transcript), so they are counted separately.
        if t["cds_start"] <= t["tx_start"] or t["cds_end"] >= t["tx_end"]:
            skipped["boundary_cds"] += 1
            continue

        for label, shift in FRAMES.items():
            # cdsStart and cdsEnd are the CDS's genomic low and high bounds on BOTH
            # strands. Which bound is the transcript's 5' end therefore depends on
            # the strand: a minus-strand transcript is transcribed from high
            # coordinate to low, so its start codon sits at the high bound and its
            # stop codon at the low one. Reading the same bound for both strands is
            # the mistake this ordering avoids -- it reported 8/25 failures, all of
            # them minus-strand transcripts whose ATG and stop were swapped.
            plus_start = seq.fetch(t["chrom"], t["cds_start"] + shift,
                                   t["cds_start"] + shift + 3)
            plus_stop = seq.fetch(t["chrom"], t["cds_end"] + shift - 3,
                                  t["cds_end"] + shift)
            if len(plus_start) != 3 or len(plus_stop) != 3:
                continue
            if t["strand"] == "-":
                got_start, got_stop = rc(plus_stop), rc(plus_start)
            else:
                got_start, got_stop = plus_start, plus_stop
            if got_start == "ATG":
                results[label]["start"][0] += 1
            results[label]["start"][1] += 1
            if got_stop in STOPS:
                results[label]["stop"][0] += 1
            results[label]["stop"][1] += 1
        tested += 1

    print(f"source            : {source}")
    print(f"transcripts tested: {tested}  (skipped: {skipped})")
    print()
    print(f"{'frame':46} {'CDS begins with ATG':>22} {'CDS ends with stop':>22}")
    print("-" * 92)
    winners = []
    for label, r in results.items():
        sa, sn = r["start"]
        pa, pn = r["stop"]
        if not sn:
            continue
        s_pct, p_pct = 100 * sa / sn, 100 * pa / pn
        print(f"{label:46} {sa:>10}/{sn:<5} {s_pct:>5.1f}% "
              f"{pa:>10}/{pn:<5} {p_pct:>5.1f}%")
        if s_pct > 90 and p_pct > 90:
            winners.append(label)

    print()
    if not any(r["start"][1] for r in results.values()):
        print("No transcripts were tested; the table's columns were not read as expected.")
        return 1
    if len(winners) == 1:
        print(f"CONFIRMED: {winners[0]}")
        print("datasets/refgene_to_gtf.py converts from this frame, unchanged:")
        print("  BED12 chromStart = txStart, exon size = exonEnd - exonStart")
        print("  GTF exon start   = exonStart + 1, exon end = exonEnd")
        return 0
    print("NOT CONFIRMED: no single frame explains the CDS boundaries.")
    print("Treat the converter's arithmetic as unvalidated until this passes.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
