#!/usr/bin/env python3
"""Generate a plain CRAM fixture for differential verification of
SAM-text/CRAM input support (DIV-0002/0004). Re-encodes the existing
bam_stat_basic.bam fixture as CRAM via pysam with NO external reference
configured, which makes htslib itself fall back to embedding the
reference in the CRAM file (embed_ref=2, confirmed by the warning it
prints) -- the case this port's open_alignments (a default, empty
noodles_cram reference-sequence repository) is designed to decode.

Deliberately reuses bam_stat_basic.bam's exact alignments (including
the unmapped1 read with an explicit MAPQ of 0) so any output difference
between running a command against the .bam vs the .cram is a genuine
format-handling bug, not a data difference -- this specific fixture is
what caught a real htslib/noodles-cram interop discrepancy in how
unmapped reads' MAPQ round-trips through CRAM (see
rseqc_formats::fix_unmapped_missing_mapping_quality's own doc comment).

Run: oracle/venv/bin/python3 verification/fixtures/make_bam_stat_cram_fixture.py <out.cram>
"""
import sys
from pathlib import Path

import pysam


def build(out_path: str) -> None:
    bam_path = Path(__file__).resolve().parent / "bam_stat_basic.bam"
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    out = pysam.AlignmentFile(out_path, "wc", header=bam.header)
    for read in bam:
        out.write(read)
    out.close()
    bam.close()


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.cram")
