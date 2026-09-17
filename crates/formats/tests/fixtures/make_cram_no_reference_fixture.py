#!/usr/bin/env python3
"""Regenerates cram_no_reference.cram: a real, pysam-written CRAM file
with no external reference configured. htslib itself falls back to
`embed_ref=2` (embedding the reference bases directly in the CRAM) in
this situation -- confirmed by the warning it prints when writing --
which is exactly the case rseqc_formats::open_alignments's CRAM
support (a default, empty noodles_cram reference-sequence repository)
can decode. One 4bp read, CIGAR 4M, no mismatches.

Run: oracle/venv/bin/python3 crates/formats/tests/fixtures/make_cram_no_reference_fixture.py
"""
import pysam

HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}


def build(out_path: str) -> None:
    out = pysam.AlignmentFile(out_path, "wc", header=HEADER)
    read = pysam.AlignedSegment()
    read.query_name = "r1"
    read.query_sequence = "ACGT"
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 10
    read.mapping_quality = 40
    read.cigarstring = "4M"
    read.query_qualities = pysam.qualitystring_to_array("IIII")
    out.write(read)
    out.close()


if __name__ == "__main__":
    import pathlib

    build(str(pathlib.Path(__file__).resolve().parent / "cram_no_reference.cram"))
