#!/usr/bin/env python3
"""Generate a plain-text SAM fixture for differential verification of
SAM-text input support (DIV-0002/0004) on RPKM_saturation.py. Re-encodes
the existing rpkm_saturation_basic.bam fixture as SAM text via pysam, so
both fixtures describe the EXACT same single alignment -- preserving the
"exactly one qualifying alignment" property that makes the whole
percentile-resampling table independent of random.shuffle's order (see
make_rpkm_saturation_fixture.py's own docstring for why that matters).

Run: oracle/venv/bin/python3 verification/fixtures/make_rpkm_saturation_sam_fixture.py <out.sam>
"""
import sys
from pathlib import Path

import pysam


def build(out_path: str) -> None:
    bam_path = Path(__file__).resolve().parent / "rpkm_saturation_basic.bam"
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    out = pysam.AlignmentFile(out_path, "wh", template=bam)
    for read in bam:
        out.write(read)
    out.close()
    bam.close()


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.sam")
