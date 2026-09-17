#!/usr/bin/env python3
"""Generate a minimal BAM fixture for differential verification of
RPKM_saturation.py. Deliberately contains exactly ONE qualifying
alignment (one exon block), so the whole percentile-resampling curve is
independent of `random.shuffle`'s order: with only one point in the
population, every percentile slice that includes it (all of them here)
produces the identical count regardless of shuffle order. This makes
the FULL saturation table byte-comparable, not just its 100% column
(the general case, since RPKM_saturation.py's percentile ranges are
cumulative across iterations -- see crates/commands/src/
rpkm_saturation.rs module docs -- only the final column is guaranteed
order-invariant with a larger population).

Run: oracle/venv/bin/python3 verification/fixtures/make_rpkm_saturation_fixture.py <out.bam>
Also writes <out_dir>/rpkm_saturation_model.bed12 alongside <out.bam>.
"""
import sys
from pathlib import Path

import pysam


def build(out_path: str) -> None:
    out = Path(out_path)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}

    with pysam.AlignmentFile(out_path, "wb", header=header) as fh:
        a = pysam.AlignedSegment()
        a.query_name = "r1"
        a.query_sequence = "A" * 20
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 40
        a.cigarstring = "20M"
        a.query_qualities = pysam.qualitystring_to_array("I" * 20)
        fh.write(a)
    pysam.index(out_path)

    bed_path = out.with_name("rpkm_saturation_model.bed12")
    bed_path.write_text("chr1\t50\t250\ttx1\t0\t+\t50\t250\t0\t1\t200,\t0,\n", encoding="utf-8")


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.bam")
