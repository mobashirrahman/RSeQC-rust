#!/usr/bin/env python3
"""Tests for the data-preparation scripts' own decisions, not just their converters.

`datasets/refgene_to_gtf.py` had its own truth cases. `build_star_index.sh` had none,
and one of its decisions is load-bearing enough to break a build outright: the index
directory name embeds the contig set, so that two indexes differing only in contig
scope cannot collide. Spelled out in full, a whole-genome contig set exceeds the
255-byte filename limit -- which is exactly the scope the production pilot needs. The
failure surfaced as `File name too long` from `mkdir`, after the script had located
STAR and unpacked the genome.

So the rule is checked here: a short set stays readable, a long one is summarised
rather than dropped, the name stays unique per contig set, and the whole thing fits in a
filename. `build_star_index.sh --print-index-key` exposes the rule without building
anything, which is what makes it testable at all.
"""

from __future__ import annotations

import os
import subprocess
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BUILD_INDEX = REPO_ROOT / "datasets" / "build_star_index.sh"

# Linux NAME_MAX. A name longer than this fails at mkdir with ENAMETOOLONG.
NAME_MAX = 255


def index_key(contigs: str, genome: str = "rn6.fa.gz", gtf: str = "rn6.gtf") -> str:
    proc = subprocess.run(
        ["bash", str(BUILD_INDEX), "--print-index-key"],
        capture_output=True, text=True, timeout=120,
        env={**os.environ, "GENOME": genome, "GTF": gtf,
             "CONTIGS": contigs, "SJDB_OVERHANG": "100"},
    )
    if proc.returncode != 0:
        raise AssertionError(f"--print-index-key failed for {contigs!r}: {proc.stderr}")
    return proc.stdout.strip()


class IndexKeyNaming(unittest.TestCase):
    def test_short_contig_set_is_spelled_out(self):
        # A short set stays legible, because the directory name is how a reader tells
        # two indexes apart at a glance.
        key = index_key("chr1 chr2 chr10")
        self.assertIn("chr1_chr2_chr10", key)
        self.assertIn("-o100-", key, "the overhang must stay in the name")

    def test_a_three_contig_key_is_unchanged_by_the_bounding_rule(self):
        # The bounding rule must not perturb names that already fit, or every existing
        # index directory would become unreachable.
        key = index_key("chr1 chr2 chr10")
        self.assertEqual(key, "rn6.fa-rn6.gtf-o100-chr1_chr2_chr10")

    def test_whole_genome_set_fits_in_a_filename(self):
        contigs = " ".join(f"chr{i}" for i in range(1, 21))
        contigs += " " + " ".join(f"chrUn_KL5684{i:02d}v1" for i in range(100, 140))
        key = index_key(contigs)
        self.assertLessEqual(len(key.encode()), NAME_MAX,
                             f"index key is {len(key)} bytes, over NAME_MAX")
        # Spelled out in full this set is far past the limit, which is the bug.
        self.assertGreater(len(contigs), NAME_MAX)

    def test_a_long_set_records_its_size_and_stays_readable(self):
        contigs = " ".join(f"chrUn_KL5684{i:02d}v1" for i in range(100, 160))
        key = index_key(contigs)
        self.assertIn("plus60more-", key,
                      "a summarised set must say how much it summarises")

    def test_different_contig_sets_get_different_keys(self):
        # The whole reason the contig set is in the name. If summarising ever collided,
        # a whole-genome index could overwrite a subset one, or vice versa -- the same
        # silent-damage failure mode that destroyed the human rat index earlier.
        a = index_key(" ".join(f"chrUn_A{i:03d}v1" for i in range(100, 160)))
        b = index_key(" ".join(f"chrUn_B{i:03d}v1" for i in range(100, 160)))
        self.assertNotEqual(a, b)

    def test_the_same_contig_set_gets_the_same_key(self):
        contigs = " ".join(f"chrUn_KL5684{i:02d}v1" for i in range(100, 160))
        self.assertEqual(index_key(contigs), index_key(contigs))

    def test_assembly_and_overhang_still_separate_keys(self):
        contigs = "chr1 chr2 chr10"
        rn6 = index_key(contigs, genome="rn6.fa.gz", gtf="rn6.gtf")
        rn7 = index_key(contigs, genome="rn7.fa.gz", gtf="rn7.gtf")
        self.assertNotEqual(rn6, rn7,
                            "rn6 and rn7 must not share an index directory")
        self.assertIn("rn6", rn6)
        self.assertIn("rn7", rn7)


if __name__ == "__main__":
    unittest.main()