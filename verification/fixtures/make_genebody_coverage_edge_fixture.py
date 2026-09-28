#!/usr/bin/env python3
"""Generate edge-case fixture for geneBody_coverage.py CIGAR D/N operations.

Covers: normal coverage (M), deletion-only D (should be visited/int 0),
skip-only N (should be visited/int 0), low-quality bases (should be int 0),
QC-fail flag (should be unvisited/float 0.0), duplicate flag (should be
unvisited/float 0.0). Tests Python duck-typing behavior where visited
positions (from pileup) render as int 0 even if no reads pass filters,
while unvisited positions stay float 0.0.

Run: oracle/venv/bin/python3 verification/fixtures/make_genebody_coverage_edge_fixture.py <outdir>
"""
import sys

import pysam


def build(out_dir: str) -> None:
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 5000}]}

    def rec(name, start, cigar, flag=0, seqlen=None, q="I"):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = 60
        a.cigarstring = cigar
        n = a.infer_query_length() if seqlen is None else seqlen
        a.query_sequence = "A" * n
        a.query_qualities = pysam.qualitystring_to_array(q * n)
        return a

    reads = [
        rec("norm", 1000, "50M"),            # normal coverage 1000-1049
        rec("del", 1100, "10M30D10M"),       # deletion-only 1110-1139
        rec("dup", 1200, "40M", flag=1024),   # duplicate-only 1200-1239
        rec("skip", 1300, "10M30N10M"),      # refskip-only 1310-1339
        rec("qcf", 1400, "40M", flag=512),    # qcfail-only
        rec("lowq", 1450, "40M", q="#"),      # base quality 2, below pysam default 13
        rec("norm2", 1000, "20M"),
    ]
    reads.sort(key=lambda r: r.reference_start)
    p = f"{out_dir}/edge.bam"
    with pysam.AlignmentFile(p, "wb", header=hdr) as o:
        for r in reads:
            o.write(r)
    pysam.index(p)
    # single-exon transcript chr1:990-1500 (>= mRNA length cutoff 100)
    open(f"{out_dir}/edge.bed12", "w").write("chr1\t990\t1500\tT1\t0\t+\t990\t1500\t0\t1\t510,\t0,\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_dir>", file=sys.stderr)
        sys.exit(1)
    build(sys.argv[1])
