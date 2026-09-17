#!/usr/bin/env python3
"""Generate a small tagged BAM fixture for differential verification of
sc_editMatrix.py. Uses pysam so the file is guaranteed byte-compatible
with the real upstream tool; also readable by noodles-bam.

Deliberately includes both a fully-dense edit-count matrix cell shape
(every position observes the same substitution) and a sparse one, since
`sc_editMatrix.py`'s CSV float/int dtype is decided per matrix COLUMN
(no transpose in its pipeline, unlike sc_seqQual.py/sc_seqLogo.py) --
see the "Preserves upstream quirks" module docs in
crates/commands/src/sc_editmatrix.rs.

Run: oracle/venv/bin/python3 verification/fixtures/make_sc_editmatrix_fixture.py <out.bam>
"""
import sys

import pysam


def build(out_path: str) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}

    with pysam.AlignmentFile(out_path, "wb", header=header) as out:
        def rec(name, cr=None, cb=None, ur=None, ub=None):
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "ACGT"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 0
            a.mapping_quality = 40
            a.cigarstring = "4M"
            a.query_qualities = pysam.qualitystring_to_array("IIII")
            tags = []
            if cr is not None:
                tags.append(("CR", cr, "Z"))
            if cb is not None:
                tags.append(("CB", cb, "Z"))
            if ur is not None:
                tags.append(("UR", ur, "Z"))
            if ub is not None:
                tags.append(("UB", ub, "Z"))
            a.set_tags(tags)
            return a

        # CB: edited at position 1 (A->T), twice -- and again at
        # position 3 (T->A) once -- a sparse (non-uniform) matrix.
        out.write(rec("r1", cr="AAAA-1", cb="ATAA-1", ur="TTTT-1", ub="TTTT-1"))
        out.write(rec("r2", cr="CCCC", cb="CCCC", ur="GGGG", ub="GGGG"))
        out.write(rec("r3", cr="AAAA", cb="ATAA", ur=None, ub=None))
        out.write(rec("r4"))  # no tags at all -> both miss
        out.write(rec("r5", cr="GGGT", cb="GGGA", ur="AAAA", ub="AAAT"))


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "fixture.bam")
