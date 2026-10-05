#!/usr/bin/env python3
"""Regenerate benchmarks/rustqc-comparison/gtf-subsets/ (orchestrator decision).

R2 found the source GTFs are whole-genome supersets of the 3-contig BEDs:
shared (chrom, transcript) keys are identical (5185/5185 rat, 50725/50725
human, 0 different), with 13465 rat / 334934 human GTF-only keys. RustQC
runs use a subset GTF containing exactly the BED's (chrom, transcript) keys
(header + exon/transcript/CDS/UTR lines for those keys + gene lines for
their genes), which passes check_annotation.py with 0/0/0.

Usage:
    python3 make_gtf_subsets.py   # writes gtf-subsets/rn6.3c.gtf + gencode.v47.3c.gtf
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from check_annotation import parse_bed12  # noqa: E402

TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')
GENE_RE = re.compile(r'gene_id\s+"([^"]+)"')


def subset_gtf(bed_path: Path, gtf_path: Path, out_path: Path):
    bed_keys = set(parse_bed12(bed_path))
    print(f"BED keys: {len(bed_keys)} from {bed_path}")
    opener_in = gzip.open if str(gtf_path).endswith(".gz") else open
    kept_genes: set[tuple] = set()
    with opener_in(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("exon", "transcript"):
                continue
            m = TX_RE.search(c[8])
            if m and (c[0], m.group(1)) in bed_keys:
                g = GENE_RE.search(c[8])
                if g:
                    kept_genes.add((c[0], g.group(1)))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_out = n_gene = 0
    with opener_in(gtf_path, "rt") as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            if c[2] in ("exon", "transcript"):
                m = TX_RE.search(c[8])
                if m and (c[0], m.group(1)) in bed_keys:
                    fout.write(line)
                    n_out += 1
            elif c[2] == "gene":
                g = GENE_RE.search(c[8])
                if g and (c[0], g.group(1)) in kept_genes:
                    fout.write(line)
                    n_gene += 1
            else:
                m = TX_RE.search(c[8])
                if m and (c[0], m.group(1)) in bed_keys:
                    fout.write(line)
                    n_out += 1
    print(f"wrote {n_out} feature lines + {n_gene} gene lines to {out_path}")


def main() -> int:
    repo = HERE.parent.parent
    subset_gtf(repo / "datasets/heldout/reference/rn6.indexed.bed12",
               repo / "datasets/heldout/reference/rn6.gtf",
               HERE / "gtf-subsets" / "rn6.3c.gtf")
    subset_gtf(repo / "datasets/heldout/reference/gencode.v47.indexed.bed12/model.bed12",
               repo / "datasets/reference/gencode.v47.annotation.gtf.gz",
               HERE / "gtf-subsets" / "gencode.v47.3c.gtf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
