# DIV-0024: minimal reproducer for a geneBody_coverage divergence

`genebody_divergence_minimal.bed12` is a single BED12 row, committed so that the
gene-body coverage difference recorded in `compatibility/divergences.yaml` as
DIV-0024 cannot change silently, and so the eventual fix has a regression case.

```
chr1  256806475  256813678  NM_013162  0  -  ...  6 exons, 904 exon bases
```

## What it reproduces

Run both implementations against this model on the rat alignment
(`datasets/heldout/aligned/SRR1177982/SRR1177982.bam`):

```bash
oracle/venv/bin/python3 oracle/upstream-src/scripts/geneBody_coverage.py \
  -i <bam> -r verification/fixtures/genebody_divergence_minimal.bed12 \
  --out-prefix py --skip-plot
./target/release/geneBody_coverage \
  -i <bam> -r verification/fixtures/genebody_divergence_minimal.bed12 \
  --out-prefix rs --skip-plot
```

**76 of the 100 bins differ**, by up to 839 reads in a single bin, in both
directions. Bins 1-14 agree; the divergence starts at bin 15. Aggregated over the
full 5,359-transcript model, one bin differs by 32,697 reads.

## What has been ruled out, and what has

Ruled out by construction. A separate 100-base transcript was built in which each
10-base range exercises exactly one pileup knob -- plain reads, sub-threshold base
quality, duplicate flag, an overlapping proper pair with identical bases in the
overlap, an overlapping proper pair with mismatched bases, a paired-but-not-proper
orphan, a `20M20D20M` deletion, secondary and supplementary flags, and a plain tail.
BED12 base *i* maps to output bin *i*, so a disagreement localises to one mechanism.
**All nine agree exactly**, including the overlap rewriting the port implements as
htslib's `tweak_overlap_quality`, and including `min_base_quality=13`.

Also ruled out: the percentile arithmetic, the strand ordering, and the transcript
length filter.

**Identified: pysam/htslib `max_depth` semantics.** Uniform piles of 8,000 / 8,001 /
12,000 reads agree exactly at 8000, so the cap's value is not the defect. Three copies
of a single `20M20D20M` read agree byte-for-byte; 9,000 copies diverge at 50 of 100
bins. The port reimplements the cap as a per-position budget over the visited set, so
`is_del` pileups consume budget that upstream's buffer budget does not.

A synthetic version of that case is committed by `verification/synthetic_data.py` as
`depth_cap.bam` (8,100 `20M20D20M` reads) and asserted by the differential case
`genebody_coverage_depth_cap_divergence`.

## Why the differential suite does not see this

`verification/run_diff.py`'s `geneBody_coverage` fixture is small and shallow
enough that none of the candidate mechanisms produces a visible difference on it,
so the 90-case suite passes while this is wrong on real data. The benchmark
harness's structural gate caught it, and refuses to report a speedup for
`geneBody_coverage`. Adding a deep-coverage fixture to the differential suite is
one of the consequences recorded in DIV-0024.
