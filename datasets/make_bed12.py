#!/usr/bin/env python3
"""Emit a BED12 gene model for the T4 panel, reusing the committed loader.

The annotation must not be produced by a second, divergent implementation:
if this script and `benchmarks/generate_workload_real.py` disagreed about
what a transcript is, a real-data failure could be blamed on the annotation
rather than the port. So the transcript loader is IMPORTED from the benchmark
generator rather than copied, and the BED12 line format is the generator's
verbatim.

Reads a GENCODE GTF and writes BED12 plus chrom.sizes for the given contigs.
Nothing is filtered by expression, subsampled, or reordered beyond the
generator's own stable sort by transcript id: the point is to test the port
against the annotation that was actually used to build the BAM's junctions.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
GEN_PATH = REPO / "benchmarks" / "generate_workload_real.py"


def load_generator():
    """Imports the benchmark's workload generator as a module.

    Imported by path rather than as a package: `benchmarks/` is a directory of
    standalone scripts, not a package, and adding an `__init__.py` to make it
    one would change how the committed benchmark is invoked.
    """
    spec = importlib.util.spec_from_file_location("generate_workload_real", GEN_PATH)
    if spec is None or spec.loader is None:
        sys.exit(f"error: cannot import {GEN_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--contigs", nargs="+", required=True)
    ap.add_argument("--out", required=True, help="output directory for model.bed12 + chrom.sizes")
    ap.add_argument("--min-exons", type=int, default=2,
                    help="transcripts with fewer exons are skipped (default 2, matching what "
                         "a gene-body/TIN analysis is meaningful for)")
    args = ap.parse_args()

    gen = load_generator()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    txs = gen.load_transcripts(Path(args.gtf), args.contigs, min_exons=args.min_exons)
    txs.sort(key=lambda t: t.tid)

    bed = out / "model.bed12"
    with bed.open("w") as fh:
        for t in txs:
            blocks = ",".join(str(e - s) for s, e in zip(t.starts, t.ends)) + ","
            starts = ",".join(str(s - t.starts[0]) for s in t.starts) + ","
            fh.write(
                f"{t.chrom}\t{t.starts[0]}\t{t.ends[-1]}\t{t.tid}\t0\t{t.strand}\t"
                f"{t.starts[0]}\t{t.ends[-1]}\t255\t{len(t.starts)}\t{blocks}\t{starts}\n"
            )

    # chrom.sizes from the GTF's own span, not from the BAM: the BAM carries
    # only the contigs STAR was given, and a gene model is about the
    # annotation, not about what happened to align.
    sizes: dict[str, int] = {}
    # The GTF may be gzipped (the committed reference is); the transcript
    # loader above already handles that via `gzip.open if str(gtf).endswith(".gz")`,
    # and this second pass must do the same or it dies on the gzip magic byte.
    opener = gzip.open if str(args.gtf).endswith(".gz") else open
    with opener(args.gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in set(args.contigs):
                continue
            sizes[f[0]] = max(sizes.get(f[0], 0), int(f[3]))
    with (out / "chrom.sizes").open("w") as fh:
        for c in sorted(sizes):
            fh.write(f"{c}\t{sizes[c]}\n")

    print(f"contigs   : {' '.join(args.contigs)}")
    print(f"transcripts: {len(txs)}  exons: {sum(len(t.starts) for t in txs)}")
    print(f"-> {bed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
