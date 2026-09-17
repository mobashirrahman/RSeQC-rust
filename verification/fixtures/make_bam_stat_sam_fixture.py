#!/usr/bin/env python3
"""Generate a plain-text SAM fixture for differential verification of
SAM-text input support (DIV-0002/0004). Simply re-encodes the existing
bam_stat_basic.bam fixture as SAM text via pysam, so the two fixtures
describe the EXACT same alignments -- any output difference between
running a command against the .bam vs the .sam is a genuine format-
handling bug, not a data difference.

Run: oracle/venv/bin/python3 verification/fixtures/make_bam_stat_sam_fixture.py <out.sam>
"""
import sys
from pathlib import Path

import pysam


def build(out_path: str) -> None:
    bam_path = Path(__file__).resolve().parent / "bam_stat_basic.bam"
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    out = pysam.AlignmentFile(out_path, "wh", template=bam)
    for read in bam:
        out.write(read)
    out.close()
    bam.close()


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.sam")
