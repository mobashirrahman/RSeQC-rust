#!/usr/bin/env python3
"""Build small, deterministic BAM/BED fixtures for known edge cases.

The fixtures are intentionally minimal and are generated with pysam so that
the reference implementation and the Rust reader consume ordinary coordinate
sorted, indexed BAM files.  Regenerating them is idempotent:

    oracle/venv/bin/python3 verification/fixtures/make_regression_fixtures.py \
        verification/fixtures
"""
from __future__ import annotations

import pathlib
import sys
import tempfile

import pysam


HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 2000}],
}


def _record(
    name: str,
    start: int,
    flag: int = 0,
    cigar: str = "20M",
    mapq: int = 40,
    quality: int = 40,
    sequence_length: int | None = None,
    mate_start: int = -1,
    template_length: int = 0,
) -> pysam.AlignedSegment:
    read = pysam.AlignedSegment()
    read.query_name = name
    sequence_length = 20 if sequence_length is None else sequence_length
    read.query_sequence = "A" * sequence_length
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigarstring = cigar
    read.query_qualities = [quality] * sequence_length
    if flag & 0x1:
        read.next_reference_id = 0
        read.next_reference_start = mate_start
        read.template_length = template_length
    return read


def _write_indexed(path: pathlib.Path, records: list[pysam.AlignedSegment]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".bam", dir=path.parent, delete=False) as handle:
        unsorted = pathlib.Path(handle.name)
    try:
        with pysam.AlignmentFile(str(unsorted), "wb", header=HEADER) as bam:
            for record in records:
                bam.write(record)
        pysam.sort("-o", str(path), str(unsorted))
        pysam.index(str(path))
    finally:
        unsorted.unlink(missing_ok=True)


def _write_model(path: pathlib.Path) -> None:
    # One 100-bp, single-exon transcript.  Coordinates are BED half-open.
    path.write_text("chr1\t100\t200\ttx1\t0\t+\t100\t200\t0\t1\t100,\t0,\n")


def _write_splice_model(path: pathlib.Path) -> None:
    # The transcript starts inside the second aligned block of the spliced
    # read below.  Its first mate starts before tx_start, so fetch selection
    # must use CIGAR reference span rather than query length.
    path.write_text("chr1\t180\t220\ttx_splice\t0\t+\t180\t220\t0\t1\t40,\t0,\n")


def build(out_dir: pathlib.Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    model = out_dir / "regression_single_exon.bed12"
    _write_model(model)
    splice_model = out_dir / "regression_splice_fetch.bed12"
    _write_splice_model(splice_model)

    # FPKM denominator edge case: one pair has both mates overlapping the
    # transcript; the second has only read 1 overlapping it.  Upstream counts
    # the latter in total fragments but not as an exonic fragment.
    fpkm_records = [
        _record("inside", 110, 0x1 | 0x2 | 0x40, mate_start=140, template_length=50),
        _record("inside", 140, 0x1 | 0x2 | 0x80 | 0x10, mate_start=110, template_length=-50),
        _record("mate_outside", 120, 0x1 | 0x2 | 0x40, mate_start=500, template_length=400),
        _record("mate_outside", 500, 0x1 | 0x2 | 0x80 | 0x10, mate_start=120, template_length=-400),
    ]
    _write_indexed(out_dir / "regression_fpkm_mate_overlap.bam", fpkm_records)

    # The first mate starts before the transcript but its spliced reference
    # span reaches it.  htslib fetch returns this record; a query-length-only
    # end estimate incorrectly drops it before count_transcript can inspect
    # the mate start inside the exon.
    fetch_records = [
        _record("splice", 80, 0x1 | 0x2 | 0x40, cigar="20M100N20M", mate_start=190, template_length=150, sequence_length=40),
        _record("splice", 190, 0x1 | 0x2 | 0x80 | 0x10, mate_start=80, template_length=-150),
        _record("inside_fetch", 185, 0x1 | 0x2 | 0x40, mate_start=200, template_length=35),
        _record("inside_fetch", 200, 0x1 | 0x2 | 0x80 | 0x10, mate_start=185, template_length=-35),
    ]
    _write_indexed(out_dir / "regression_fpkm_fetch_span.bam", fetch_records)

    # CIGAR '=' consumes reference bases exactly like 'M'.  This catches
    # reference-span code that only handles M/D/N and therefore believes this
    # pair has zero transcript overlap.
    rna_records = [
        _record("equals", 100, 0x1 | 0x2 | 0x40, cigar="20=", mate_start=120, template_length=40),
        _record("equals", 120, 0x1 | 0x2 | 0x80 | 0x10, cigar="20=", mate_start=100, template_length=-40),
    ]
    _write_indexed(out_dir / "regression_rna_equals.bam", rna_records)

    # Overlapping paired ends: pysam pileup() defaults ignore_overlaps=True
    # and keeps one base (the higher-quality mate) in the shared interval.
    overlap_records = [
        _record("overlap", 100, 0x1 | 0x2 | 0x40, mate_start=110, template_length=30, quality=40),
        _record("overlap", 110, 0x1 | 0x2 | 0x80 | 0x10, mate_start=100, template_length=-30, quality=40),
    ]
    _write_indexed(out_dir / "regression_overlap_pair.bam", overlap_records)

    # pysam's default pileup max_depth is 8000.  8001 reads at one locus
    # expose an implementation that does not preserve that compatibility
    # limit.  The extra read prevents the fixture from being a one-coordinate
    # special case for transcript scaling.
    depth_records = [_record(f"depth{i}", 110) for i in range(8001)]
    depth_records.append(_record("tail", 150))
    _write_indexed(out_dir / "regression_genebody_depth.bam", depth_records)


if __name__ == "__main__":
    build(pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).parent)
