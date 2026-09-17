#!/usr/bin/env python3
"""Generate a small, diverse BAM fixture for differential verification
of bam_stat.py. Uses pysam so the file is guaranteed byte-compatible
with the real upstream tool; also readable by noodles-bam (standard
BAM format, no pysam-specific extensions used).

Covers: unique high-mapq forward reads, a low-mapq ("non-unique") read,
a duplicate, a QC-fail read, a secondary alignment, an unmapped read,
and a paired proper pair -- enough categories to exercise bam_stat.py's
real branch logic, not just a trivial single-read case.

Run: oracle/venv/bin/python3 verification/fixtures/make_bam_stat_fixture.py <out.bam>
"""
import sys

import pysam


def build(out_path: str) -> None:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 500}]}

    with pysam.AlignmentFile(out_path, "wb", header=header) as out:
        def rec(name, ref_id, pos, flag, mapq, cigar="20M", seq="A" * 20):
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = seq
            a.flag = flag
            a.reference_id = ref_id
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            return a

        # 3 unique, high-mapq, forward-strand reads on chr1.
        out.write(rec("u1", 0, 10, 0, 40))
        out.write(rec("u2", 0, 50, 0, 40))
        out.write(rec("u3", 0, 90, 0, 40))

        # 1 low-mapq ("non-unique" under the default -q 30 threshold).
        out.write(rec("lowmapq1", 0, 130, 0, 5))

        # 1 duplicate (flag 0x400) of a high-mapq read.
        out.write(rec("dup1", 0, 170, 0x400, 40))

        # 1 QC-fail read (flag 0x200).
        out.write(rec("qcfail1", 0, 210, 0x200, 40))

        # 1 secondary alignment (flag 0x100).
        out.write(rec("secondary1", 0, 250, 0x100, 40))

        # 1 unmapped read (flag 0x4); reference_id/pos are ignored by
        # convention for unmapped reads but pysam still needs valid
        # placeholder values when SO=coordinate.
        unmapped = rec("unmapped1", -1, -1, 0x4, 0)
        unmapped.reference_id = -1
        unmapped.reference_start = -1
        out.write(unmapped)

        # A proper pair on chr2: read1 forward + read2 reverse. Mate
        # fields must be set explicitly -- pysam defaults
        # next_reference_id to -1 (unset), which bam_stat.py's
        # "different chromosome" check would otherwise misread as a
        # cross-chromosome pair even though both mates are on chr2.
        r1 = rec("pair1", 1, 20, 0x1 | 0x2 | 0x40, 40)
        r1.next_reference_id = 1
        r1.next_reference_start = 60
        r1.template_length = 60
        out.write(r1)

        r2 = rec("pair1", 1, 60, 0x1 | 0x2 | 0x80 | 0x10, 40)
        r2.next_reference_id = 1
        r2.next_reference_start = 20
        r2.template_length = -60
        out.write(r2)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.bam")
