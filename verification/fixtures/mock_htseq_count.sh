#!/bin/bash
# Mock htseq-count for differential verification of FPKM-UQ.py: the real
# htseq-count is not part of this sandbox's toolchain (a real per-alignment
# scan is not what's under test here -- calculate_fpkm's post-processing
# is). Ignores its real arguments and prints a small, fixed count table
# covering a zero-count gene, a non-protein-coding gene (excluded from
# the two normalization denominators but still rendered), and the
# htseq-count summary rows (`__no_feature`/`__ambiguous`) that both
# implementations must skip.
cat <<'COUNTS'
ENSG1	100
ENSG2	50
ENSG3	0
__no_feature	5
__ambiguous	2
COUNTS
