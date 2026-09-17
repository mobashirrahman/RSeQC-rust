#!/usr/bin/env python3
"""Generate a small, tagged BAM fixture for differential verification of
sc_bamStat.py. Uses pysam so the file is guaranteed byte-compatible with
the real upstream tool; also readable by noodles-bam (standard BAM
format, no pysam-specific extensions used).

Covers: confident reads (xf bit0 set) spanning both region-type branches
(E/I), a confident read with an RE tag value that is neither E nor I
(exercises the upstream if/elif/else dead-code branch), sense vs.
antisense tags, a duplicate, a reverse-strand read, a read on the
mitochondrial contig, and a non-confident read (xf bit0 clear, excluded
from every downstream tally but still counted in total_reads_n).

Run: oracle/venv/bin/python3 verification/fixtures/make_sc_bamstat_fixture.py <out.bam>
Produces <out.bam> and <out.bam>.bai (sc_bamStat.py requires an index).
"""
import sys

import pysam


def build(out_path: str) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chrM", "LN": 200}],
    }

    with pysam.AlignmentFile(out_path, "wb", header=header) as out:
        def rec(name, ref_id, pos, flag, mapq, tags, cigar="20M", seq="A" * 20):
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = seq
            a.flag = flag
            a.reference_id = ref_id
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.set_tags(tags)
            return a

        # Confident, forward, non-dup, exonic, sense, with CB+UB.
        out.write(rec("r1", 0, 10, 0, 40, [("xf", 1, "i"), ("CB", "AAAA-1", "Z"), ("UB", "TTTT", "Z"), ("RE", "E", "A"), ("TX", 1, "i")]))

        # Confident, reverse, duplicate, intronic, antisense, no CB/UB.
        out.write(rec("r2", 0, 60, 0x10 | 0x400, 40, [("xf", 1, "i"), ("RE", "I", "A"), ("AN", 1, "i")]))

        # Confident, RE tag present but neither E nor I -- exercises the
        # upstream if/elif/else dead-code branch (silent no-op, not even
        # "other").
        out.write(rec("r3", 0, 110, 0, 40, [("xf", 1, "i"), ("RE", "N", "A"), ("TX", 1, "i")]))

        # Non-confident (xf bit0 clear): excluded from every confident-*
        # tally but still contributes to total_reads_n via its QNAME.
        out.write(rec("r5", 0, 150, 0, 40, [("xf", 0, "i")]))

        # Confident, on the mitochondrial contig (must sort after all
        # chr1 records: coordinate-sorted BAM requires ascending
        # reference_id, not just ascending position within a contig).
        out.write(rec("r4", 1, 5, 0, 40, [("xf", 1, "i"), ("RE", "E", "A"), ("TX", 1, "i")]))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "fixture.bam"
    build(out)
    pysam.index(out)
