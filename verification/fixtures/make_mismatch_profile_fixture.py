#!/usr/bin/env python3
"""Build a minimal BAM fixture with a genuine mismatch (real MD/NM tags),
for differential verification of mismatch_profile.py's actual counting
logic. The existing mismatch_profile_no_mismatches harness case only
exercises the zero-mismatches early-exit path -- this fixture covers the
command's actual purpose.

One 10bp read, CIGAR 10M, with a single A->C mismatch at 0-based read
position 5 (MD tag "5A4", NM tag 1) -- an A2C genotype count in the
output table.

Run: oracle/venv/bin/python3 verification/fixtures/make_mismatch_profile_fixture.py <out.bam>
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
        read.query_name = "mm1"
        # Reference bases (per MD "5A4") are all 'A'; the query differs
        # only at 0-based position 5 (A -> C), matching MD/NM below.
        read.query_sequence = "AAAAACAAAA"
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 10
        read.mapping_quality = 40
        read.cigarstring = "10M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 10)
        read.set_tag("MD", "5A4")
        read.set_tag("NM", 1)
        fh.write(read)
    pysam.index(out_path)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.bam")
