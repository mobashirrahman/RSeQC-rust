#!/usr/bin/env python3
"""Build a minimal BAM fixture with a genuine deletion, for differential
verification of deletion_profile.py's actual counting logic. The
existing deletion_profile_no_deletions harness case only exercises the
zero-deletions early-exit path -- this fixture covers the command's
actual purpose.

One read, CIGAR 5M2D5M (M/S/I total = 10, matching read_align_length
10), with a real 2-base deletion starting right after the first 5
matched bases.

Run: oracle/venv/bin/python3 verification/fixtures/make_deletion_profile_fixture.py <out.bam>
"""
import sys

import pysam

HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000}],
}


def build(out_path: str) -> None:
    with pysam.AlignmentFile(out_path, "wb", header=HEADER) as fh:
        read = pysam.AlignedSegment()
        read.query_name = "del1"
        read.query_sequence = "AAAAAAAAAA"
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 10
        read.mapping_quality = 40
        read.cigarstring = "5M2D5M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 10)
        fh.write(read)
    pysam.index(out_path)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.bam")
