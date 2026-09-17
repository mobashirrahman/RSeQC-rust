#!/usr/bin/env python3
"""Generate tiny track fixtures for differential verification of the
BigWig/WIG command family: bam2wig.py, geneBody_coverage2.py,
normalize_bigwig.py, overlay_bigwig.py.

Covers: overlapping equal runs, negative values, explicit zeros versus
missing (gap) intervals, tail intervals, and a two-transcript BED12 model
(one single-exon 100bp boundary transcript, one two-exon transcript whose
first exon alone clears the legacy 100bp filter -- see DIV-0014).

Run: oracle/venv/bin/python3 verification/fixtures/make_track_fixtures.py <out_dir>
Writes: track_chrom.sizes, track_signal.bw, track_signal2.bw,
track_reads.bam(+.bai), track_model.bed12. All inputs are tiny and
human-auditable; hashes are validated by run_diff at case time.
"""
import sys
from pathlib import Path

import pysam


def build(out_dir: str) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    (out / "track_chrom.sizes").write_text("chr1\t1000\n", encoding="utf-8")

    import pyBigWig

    bw = pyBigWig.open(str(out / "track_signal.bw"), "w")
    bw.addHeader([("chr1", 1000)])
    # pyBigWig cannot store overlapping intervals, so the 150-160 sub-run
    # is written as its own piece: per-base values are 5.0 on 100-150 and
    # 2.5 on 150-160, with an explicit zero run and gaps (missing data).
    bw.addEntries(
        ["chr1"] * 5,
        [100, 150, 300, 500, 900],
        ends=[150, 160, 310, 505, 1000],
        values=[5.0, 2.5, -3.0, 0.0, 1.5],
    )
    bw.close()

    bw2 = pyBigWig.open(str(out / "track_signal2.bw"), "w")
    bw2.addHeader([("chr1", 1000)])
    bw2.addEntries(
        ["chr1"] * 4,
        [100, 150, 300, 700],
        ends=[150, 160, 310, 710],
        values=[1.0, 4.0, 2.0, 10.0],
    )
    bw2.close()

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    bam_path = str(out / "track_reads.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as fh:
        for name, pos in (("r1", 110), ("r2", 120), ("r3", 305), ("r4", 505), ("r5", 905)):
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "A" * 20
            a.flag = 0
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 40
            a.cigarstring = "20M"
            a.query_qualities = pysam.qualitystring_to_array("I" * 20)
            fh.write(a)
    pysam.index(bam_path)

    # BED12: tx1 single exon 100-200 (+); tx2 two exons 300-450,450-500 (+).
    (out / "track_model.bed12").write_text(
        "chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t1\t100,\t0,\n"
        "chr1\t300\t500\ttx2\t0\t+\t300\t500\t0\t2\t150,50,\t0,150,\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "track_fixtures")
