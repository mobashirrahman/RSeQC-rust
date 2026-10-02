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


def build_external_reference(out_path: str, ref_path: str) -> None:
    """The CONTRASTING case: a CRAM encoded against an external reference.

    This is the shape real CRAMs usually have -- `samtools view -C -T ref.fa` writes
    the reference's URI and M5 into the header and does NOT embed the bases. Decoding
    it needs a reference, and rseqc_formats resolves none by design, so the correct
    behaviour is a clean error naming the file and the reason. Before that was
    implemented, `noodles_cram` reached `.expect("invalid slice reference sequence
    name")` and the user saw a Rust panic with exit 101.

    A 1 kb reference of N's is enough: what is under test is the reference-resolution
    branch, not the sequence content. htslib requires the file to exist and be
    indexed at write time, which is why the bases are written first.
    """
    with open(ref_path, "w") as fh:
        fh.write(">chr1\n")
        fh.write("N" * 1000 + "\n")
    pysam.faidx(ref_path)
    out = pysam.AlignmentFile(out_path, "wc", header=HEADER,
                              reference_filename=ref_path)
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

    here = pathlib.Path(__file__).resolve().parent
    build(str(here / "cram_no_reference.cram"))
    build_external_reference(str(here / "cram_external_reference.cram"),
                             str(here / "cram_external_reference.fa"))
