#!/usr/bin/env python3
"""Independent truth cases for `refgene_to_gtf.py`.

Every expected value here is hand-derived from UCSC's documented counting
systems, not from any observed command output, so a converter that merely
reproduces the previous revision's numbers cannot pass. They are what the
readiness audit's P0 "scientific evidence" finding requires.

UCSC refGene is half-open 0-based:

    refGene txStart/exonStarts   0-based start offsets
    refGene txEnd/exonEnds       exclusive end offsets
    BED12                        0-based half-open    -> verbatim, unchanged
    GTF                          1-based inclusive    -> start + 1, end unchanged

The regression this pins is not a single off-by-one. The earlier revisions were
wrong in *both* directions at once (transcript start shifted -1 with exon
lengths +1; then lengths +1 with GTF starts unchanged), which is why the BED and
GTF outputs disagreed with each other about inclusivity instead of cancelling.

Runs on the standard library only: no oracle environment, no Rust build.
    python3 datasets/test_refgene_to_gtf.py
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "refgene_to_gtf.py"

COLUMNS = ("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart",
           "cdsEnd", "exonCount", "exonStarts", "exonEnds")


def refgene_row(name, chrom, strand, tx_start, tx_end, exon_starts, exon_ends,
                cds_start=None, cds_end=None, score=0):
    """Build one refGene line. Coordinates are half-open 0-based, as UCSC stores them."""
    return "\t".join([
        str(512), name, chrom, strand, str(tx_start), str(tx_end),
        str(tx_start if cds_start is None else cds_start),
        str(tx_end if cds_end is None else cds_end),
        str(len(exon_starts)),
        ",".join(str(x) for x in exon_starts) + ",",
        ",".join(str(x) for x in exon_ends) + ",",
        str(score), str(0), str(0), str(len(exon_starts)), "", "",
        f"{name},", "", "", "", "0", "",
    ])


def gtf_ranges(path):
    ranges = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        ranges.append((f[0], f[2], int(f[3]), int(f[4]), f[6]))
    return ranges


def bed_rows(path):
    return [line.split("\t") for line in path.read_text().splitlines() if line.strip()]


class RefGeneConversionTests(unittest.TestCase):
    def convert(self, rows, extra=()):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = Path(tmp.name)
        src = d / "refGene.txt"
        src.write_text("\n".join(rows) + "\n")
        cmd = [sys.executable, str(SCRIPT), "--refgene", str(src),
               "--gtf-out", str(d / "out.gtf"), "--bed-out", str(d / "out.bed12"),
               *extra]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        return d / "out.gtf", d / "out.bed12", proc.stderr

    def test_audit_probe_row(self):
        """The audit's probe input, with the audit's expected outputs.

        Input exons [100,200) and [300,400): correct BED12 start 100, sizes
        100,100; correct GTF starts 101 and 301. The pre-fix converter emitted
        start 99, sizes 101,101 and GTF starts 100,300.
        """
        row = refgene_row("NM_TEST", "chr1", "+", 100, 400, [100, 300], [200, 400])
        gtf, bed, _ = self.convert([row])

        self.assertEqual(bed_rows(bed), [[
            "chr1", "100", "400", "NM_TEST", "0", "+", "100", "400", "255", "2",
            "100,100,", "0,200,",
        ]])
        self.assertEqual(gtf_ranges(gtf), [
            ("chr1", "exon", 101, 200, "+"),
            ("chr1", "exon", 301, 400, "+"),
        ])

    def test_minus_strand_is_not_reversed(self):
        """Strand changes the frame label only; coordinates stay genomic and ascending.

        refGene stores exons in ascending genomic order on both strands, so a
        minus-strand row must not be reverse-complemented into descending order
        (a GTF reader sorts nothing for you).
        """
        row = refgene_row("NM_MINUS", "chrX", "-", 1000, 1500, [1000, 1400], [1100, 1500])
        gtf, bed, _ = self.convert([row])

        r = bed_rows(bed)[0]
        self.assertEqual(r[1], "1000")
        self.assertEqual(r[2], "1500")
        self.assertEqual(r[5], "-")
        self.assertEqual(r[9], "2")
        self.assertEqual(r[10], "100,100,")
        self.assertEqual(r[11], "0,400,")
        self.assertEqual(gtf_ranges(gtf), [
            ("chrX", "exon", 1001, 1100, "-"),
            ("chrX", "exon", 1401, 1500, "-"),
        ])

    def test_zero_coordinate_boundary(self):
        """A transcript at coordinate 0 is a valid case, not an off-by-one.

        Under the old `txStart - 1` form this produced a BED start of -1 and a
        relative first-exon start of +1; the correct answer keeps both at 0.
        """
        row = refgene_row("NM_ZERO", "chr1", "+", 0, 300, [0, 200], [100, 300])
        gtf, bed, _ = self.convert([row])

        r = bed_rows(bed)[0]
        self.assertEqual(r[1], "0", "chromStart must not go negative at coordinate 0")
        self.assertEqual(r[10], "100,100,")
        self.assertEqual(r[11], "0,200,")
        self.assertEqual(gtf_ranges(gtf), [
            ("chr1", "exon", 1, 100, "+"),
            ("chr1", "exon", 201, 300, "+"),
        ])

    def test_single_base_exon(self):
        """A one-base exon is size 1, and GTF places it at start..start, not start-1..start."""
        row = refgene_row("NM_ONE", "chr1", "+", 500, 500, [500], [501])
        gtf, bed, _ = self.convert([row])

        r = bed_rows(bed)[0]
        self.assertEqual(r[9], "1")
        self.assertEqual(r[10], "1,")
        self.assertEqual(r[11], "0,")
        self.assertEqual(gtf_ranges(gtf), [("chr1", "exon", 501, 501, "+")])

    def test_adjacent_exons_preserve_intron_length(self):
        """Exon lengths must sum with the introns to exactly txEnd - txStart.

        This is the invariant the old converter violated: its +1 per exon made
        the exons one base longer than the transcript they sit in, which is what
        destroyed junction matching downstream.
        """
        tx_start, tx_end = 1000, 5000
        starts, ends = [1000, 2000, 3000], [1500, 2500, 5000]
        row = refgene_row("NM_SUM", "chr1", "+", tx_start, tx_end, starts, ends)
        gtf, bed, _ = self.convert([row])

        r = bed_rows(bed)[0]
        sizes = [int(x) for x in r[10].rstrip(",").split(",")]
        rel = [int(x) for x in r[11].rstrip(",").split(",")]
        self.assertEqual(sum(sizes) + sum(
            (starts[i + 1] - ends[i]) for i in range(len(ends) - 1)
        ), tx_end - tx_start)
        # Relative starts must resolve to the refGene absolute starts.
        for size, offset, abs_start in zip(sizes, rel, starts):
            self.assertEqual(tx_start + offset, abs_start)
        # Walking the BED12 blocks must consume the transcript exactly: exons plus
        # the gaps the relative starts imply, with nothing left over and no overlap.
        self.assertEqual(rel[0], 0, "first relative start must be 0")
        for i in range(len(sizes) - 1):
            self.assertGreaterEqual(
                rel[i + 1], rel[i] + sizes[i],
                "BED12 blocks must not overlap; an off-by-one here is what the "
                "previous revision produced")
        last_end = rel[-1] + sizes[-1]
        self.assertEqual(tx_start + last_end, tx_end,
                         "blocks plus trailing intron must end exactly at txEnd")

        # GTF names the same bases one base up on the start side only: the
        # inclusive end is already the refGene exclusive end.
        self.assertEqual([(s - 1, e) for _, _, s, e, _ in gtf_ranges(gtf)],
                         list(zip(starts, ends)))

    def test_bed_and_gtf_describe_identical_intervals(self):
        """The strongest cross-consumer check: both outputs must name the same bases.

        BED12 and GTF disagree about inclusivity, so a converter can be
        self-consistent in each format and still describe different intervals.
        Reconstructing each exon's 0-based half-open interval and comparing them
        catches exactly that class of error.
        """
        rows = [
            refgene_row("NM_A", "chr2", "+", 100, 400, [100, 300], [200, 400]),
            refgene_row("NM_B", "chr2", "-", 900, 1200, [900, 1100], [1000, 1200]),
            refgene_row("NM_C", "chr2", "+", 0, 50, [0, 20], [20, 50]),
        ]
        gtf, bed, _ = self.convert(rows)

        by_name = {}
        for _chrom, _feat, start, end, strand in gtf_ranges(gtf):
            by_name.setdefault(strand, []).append((start - 1, end))
        bed_intervals = []
        for r in bed_rows(bed):
            strand, size_field, rel_field = r[5], r[10], r[11]
            sizes = [int(x) for x in size_field.rstrip(",").split(",")]
            rels = [int(x) for x in rel_field.rstrip(",").split(",")]
            base = int(r[1])
            for size, rel in zip(sizes, rels):
                bed_intervals.append((strand, base + rel, base + rel + size))
        self.assertEqual(sorted(by_name["+"] + by_name["-"]),
                         sorted((s, e) for _s, s, e in bed_intervals))

    def test_contig_filter_and_min_exons(self):
        rows = [
            refgene_row("NM_KEEP", "chr1", "+", 10, 60, [10, 40], [30, 60]),
            refgene_row("NM_SOLO", "chr1", "+", 100, 130, [100], [130]),
            refgene_row("NM_OTHER", "chr2", "+", 10, 60, [10, 40], [30, 60]),
        ]
        gtf, bed, _ = self.convert(rows, extra=["--contigs", "chr1", "--min-exons", "2"])
        self.assertEqual([r[3] for r in bed_rows(bed)], ["NM_KEEP"])
        self.assertEqual(len(gtf_ranges(gtf)), 2)

    def test_empty_exon_interval_is_skipped(self):
        """A zero-length exon is dropped with a note rather than emitted as size 0.

        BED12 has no representation for an empty block, so emitting one would put
        a malformed record into every downstream annotation consumer.
        """
        rows = [
            refgene_row("NM_BAD", "chr1", "+", 10, 60, [10, 40], [10, 60]),
            refgene_row("NM_OK", "chr1", "+", 100, 160, [100, 140], [130, 160]),
        ]
        gtf, bed, stderr = self.convert(rows)
        self.assertEqual([r[3] for r in bed_rows(bed)], ["NM_OK"])
        self.assertIn("NM_BAD", stderr)

    def test_mismatched_exon_lists_are_skipped(self):
        rows = [
            "\t".join(["512", "NM_MM", "chr1", "+", "10", "60", "10", "60", "2",
                       "10,40,", "30,", "0", "0", "0", "2", "", "", "NM_MM,", "", "", "",
                       "0", ""]),
            refgene_row("NM_OK", "chr1", "+", 100, 160, [100, 140], [130, 160]),
        ]
        gtf, bed, _ = self.convert(rows)
        self.assertEqual([r[3] for r in bed_rows(bed)], ["NM_OK"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
