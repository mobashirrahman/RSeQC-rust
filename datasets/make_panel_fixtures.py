#!/usr/bin/env python3
"""Derive the real-data panel's fixture directory from the aligned runs.

The differential matrix's cases are written against a fixture directory with
specific names (`pe.bam`, `se.bam`, `model.bed12`, `chrom.sizes`) because
their expected outputs are derived from the BAM's stem -- a case comparing
`pe.tin.xls` needs the input to be called `pe.bam`. So pointing
RSEQC_REAL_DATA at the aligned runs directly does not work: the output names
would not match what the cases expect.

This builds that directory, and every derivation is recorded so it is
reproducible rather than a shell command someone ran once:

  pe.bam / se.bam      A windowed slice of each aligned development run, named
                       for the case that consumes it.
  model.bed12          Transcripts restricted to the same window.
  chrom.sizes          Assembly lengths taken from the BAM HEADER, which is
                       authoritative. Deliberately NOT taken from the GTF: the
                       GTF's last exon's end understates the contig length
                       (chr1: 248,936,715 from the annotation vs 248,956,422
                       in the assembly), and bam2wig would then emit tracks
                       against the wrong genome size.

`WINDOW` must be a gene-dense region. Choosing it badly is not a neutral
call: an early pick of chr17:1-30Mb (chr17's gene-poor short arm) produced
zero exonic fragments, so FPKM_count legitimately errored on both arms and
the case "failed" for a reason that had nothing to do with either port.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def samtools() -> str:
    """Path to samtools, preferring the pinned env over anything on PATH."""
    pinned = Path("/scratch/mdra00001/tmp/mamba/envs/t4star/bin/samtools")
    if pinned.is_file():
        return str(pinned)
    found = shutil.which("samtools")
    if not found:
        sys.exit("error: samtools not found; install the pinned t4star env or put samtools on PATH")
    return found


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True,
                    metavar="NAME=PATH", help="fixture name = aligned run BAM, e.g. pe=datasets/aligned/SRR/SRR.bam")
    ap.add_argument("--window", required=True, help="region, e.g. chr1:100000000-130000000")
    ap.add_argument("--bed", required=True, help="full BED12 to restrict")
    ap.add_argument("--out", required=True)
    ap.add_argument("--force", action="store_true", help="rebuild even if the output exists")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    st = samtools()

    for spec in args.run:
        name, _, run_bam = spec.partition("=")
        dest = out / f"{name}.bam"
        if dest.is_file() and not args.force:
            print(f"ok       {name}.bam (present)")
        else:
            # Regions are POSITIONAL in samtools >= 1.19. `view -r REGION` is
            # --read-group, and passing a region name there silently yields
            # every record rather than erroring -- which would produce a "slice"
            # that is a full copy. See datasets/README.md.
            subprocess.run([st, "view", "-b", "-o", str(dest), run_bam, args.window], check=True)
            subprocess.run([st, "index", str(dest)], check=True)
            n = subprocess.run([st, "view", "-c", str(dest)], check=True, capture_output=True, text=True)
            print(f"slice    {name}.bam  {n.stdout.strip()} reads  ({args.window})")
        # Fail loudly if the BAM is not sorted, rather than letting the
        # comparator discover it as an unexplained diff.
        check = subprocess.run([st, "quickcheck", "-v", str(dest)], capture_output=True, text=True)
        if check.returncode != 0:
            sys.exit(f"error: {dest} failed samtools quickcheck: {check.stderr}")

    contig, _, span = args.window.partition(":")
    lo, _, hi = span.partition("-")
    bed_out = out / "model.bed12"
    with open(args.bed) as src, bed_out.open("w") as dst:
        n = 0
        for line in src:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] != contig:
                continue
            # BED12 is half-open [start, end) in 0-based coordinates.
            if int(f[1]) < int(hi) and int(f[2]) > int(lo):
                dst.write(line)
                n += 1
    print(f"model    {bed_out}  {n} transcripts in {args.window}")

    # Assembly lengths from the BAM header, not the annotation: see module docstring.
    first = next(iter(args.run))[len(next(iter(args.run)).split("=", 1)[0]) + 1:]
    sizes = out / "chrom.sizes"
    header = subprocess.run([st, "view", "-H", first], check=True, capture_output=True, text=True).stdout
    with sizes.open("w") as fh:
        for line in header.splitlines():
            if not line.startswith("@SQ"):
                continue
            parts = dict(p.split(":", 1) for p in line.split("\t")[1:] if ":" in p)
            fh.write(f"{parts['SN']}\t{parts['LN']}\n")
    print(f"sizes    {sizes}  {len(sizes.read_text().splitlines())} contigs (from the BAM header)")

    print(f"\npanel ready: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
