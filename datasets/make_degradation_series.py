#!/usr/bin/env python3
"""Build a controlled RNA-degradation series by truncating read 3' ends.

Endpoint E2 needs TIN measured against KNOWN sample conditions (testing.md
section 11.1, "experimentally degraded RNA or a measured integrity series").
Real degradation data would need a lab; this makes the perturbation explicit
and reversible instead.

Truncating the 3' end of every read is the standard proxy for a degraded
library: reads are shorter, and coverage along a transcript becomes more
uneven, which is exactly what TIN's entropy term responds to.

The truncation is from the 3' end only, keeping the 5' anchor, because TIN
depends on the spatial distribution of read starts and destroying the 5' end
would conflate the two ends.

Reads shorter than the requested length are dropped, and the drop is reported,
because silently keeping short reads would make a 'truncate to 50' series
really be a mixed-length series.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq1", required=True)
    ap.add_argument("--fastq2", default=None,
                    help="omit for single-end. TIN does not need pairs, and the "
                         "BAM-to-FASTQ path drops singletons asymmetrically, so a "
                         "paired series has to be built from a name-collated BAM "
                         "rather than two independently written files")
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-keep", type=int, default=0,
                    help="drop pairs where either mate is shorter than this")
    args = ap.parse_args()

    p1 = args.out_prefix + "_1.fastq.gz"
    p2 = args.out_prefix + "_2.fastq.gz" if args.fastq2 else None
    kept = dropped = 0
    dropped_len = 0
    f2 = gzip.open(args.fastq2, "rt") if args.fastq2 else None
    o2 = gzip.open(p2, "wt") if p2 else None
    with gzip.open(args.fastq1, "rt") as fin1, gzip.open(p1, "wt") as o1:
        while True:
            h1, s1, p1q, q1 = fin1.readline(), fin1.readline(), fin1.readline(), fin1.readline()
            if not h1:
                break
            if not (s1 and p1q and q1):
                sys.exit(f"error: {args.fastq1} ended mid-record")

            recs2 = None
            if f2 is not None:
                recs2 = (f2.readline(), f2.readline(), f2.readline(), f2.readline())
                if not all(recs2):
                    sys.exit(f"error: {args.fastq2} ended early; files are not paired in lockstep")

            shortest = len(s1) if recs2 is None else min(len(s1), len(recs2[1]))
            if shortest < args.length:
                dropped += 1
                dropped_len += shortest
                continue
            o1.write(f"{h1}{s1[:args.length]}\n+\n{q1[:args.length]}\n")
            if recs2 is not None:
                o2.write(f"{recs2[0]}{recs2[1][:args.length]}\n+\n{recs2[3][:args.length]}\n")
            kept += 1

    total = kept + dropped
    print(f"truncate to {args.length:>4}: kept {kept:>9,}  dropped {dropped:>9,} "
          f"({100*dropped/max(1,total):5.2f}%)  mean dropped length {dropped_len/max(1,dropped):.0f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
