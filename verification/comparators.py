#!/usr/bin/env python3
"""Artifact comparators shared by the differential runner and the benchmark harness.

Both harnesses decide whether two independent implementations produced the same
result, so a weakness in a comparator is shared by both and becomes a false pass in
whichever one happens to be reporting. The audit found exactly that: the benchmark's
comparators accepted changed BAM quality scores, duplicate flags and NM tags, a
FASTQ with a trailing partial record, and finite text replaced by NaN -- all
reported as passing.

Three rules hold throughout this module:

  1. FAIL CLOSED. A comparator that cannot decide returns a failure with a reason.
     A false pass turns a silent regression into a reported speedup; a false
     failure costs one debugging session.
  2. COMPARE THE MEASUREMENT, NOT THE ENCODING. Compressed bytes, padding width and
     floating-point spelling are formatting, so they are normalised. Quality
     strings, flags, tags, coordinates, CIGARs and mate fields are measurements, so
     they are not.
  3. NONFINITE IS NOT A VALUE. NaN and infinity are rejected rather than compared,
     including when both arms produce the same one, because two arms both dividing
     by zero is not agreement.
"""
from __future__ import annotations

import re
from pathlib import Path

# A bare numeric literal, as opposed to a label or a filename.
NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
# Spellings float() accepts for a nonfinite value.
NONFINITE = re.compile(r"^[+-]?(inf(inity)?|nan)$", re.IGNORECASE)

SAMPLE_FIELDS = (
    "header", "sequence", "separator", "quality",
)


def is_nonfinite(token: str) -> bool:
    """True when `token` is a numeric cell that is not finite.

    A word that merely contains "nan", such as a label or a filename, is not a
    numeric cell and returns False.
    """
    t = token.strip()
    if not t:
        return False
    if NONFINITE.match(t):
        return True
    if not NUM.match(t):
        return False
    try:
        value = float(t)
    except ValueError:
        return False
    return value != value or value in (float("inf"), float("-inf"))


def _norm_line(line: str) -> str:
    """Strip and round bare numbers to a common precision.

    Two implementations may format the same value differently (`1.0` versus
    `1.000000`); that is encoding, not measurement.
    """
    parts = [p.strip() for p in line.rstrip().split("\t")]
    out = []
    for p in parts:
        if NUM.match(p):
            try:
                out.append(f"{float(p):.6g}")
                continue
            except ValueError:
                pass
        out.append(p)
    return "\t".join(out)


def compare_text(a: Path, b: Path, rtol: float = 1e-6, atol: float = 1e-9):
    """Line-wise text comparison with numeric tolerance.

    Returns (ok, reason). Any nonfinite cell fails, in either file, before any
    comparison: a metric that is not finite is a defect in the command that
    produced it, and symmetric nonfinite values are still not agreement.
    """
    la = [x for x in a.read_text(errors="replace").splitlines() if x.strip()]
    lb = [x for x in b.read_text(errors="replace").splitlines() if x.strip()]
    if len(la) != len(lb):
        return False, f"line count {len(la)} vs {len(lb)}"
    for i, (x, y) in enumerate(zip(la, lb)):
        for token in x.replace("\t", " ").split():
            if is_nonfinite(token):
                return False, f"line {i + 1}: nonfinite value {token!r} in first output"
        for token in y.replace("\t", " ").split():
            if is_nonfinite(token):
                return False, f"line {i + 1}: nonfinite value {token!r} in second output"
        if _norm_line(x) != _norm_line(y):
            # Fall back to a tolerance comparison for a line that differs only in
            # float spelling, and to whitespace-separated cells for the fixed-width
            # reports several commands print.
            if not _cells_within_tolerance(x, y, rtol, atol):
                return False, f"line {i + 1}: {x[:70]!r} vs {y[:70]!r}"
    return True, ""


def _cells_within_tolerance(x: str, y: str, rtol: float, atol: float) -> bool:
    fx, fy = x.split(), y.split()
    if len(fx) != len(fy):
        return False
    for u, v in zip(fx, fy):
        if u == v:
            continue
        if not (NUM.match(u) and NUM.match(v)):
            return False
        fu, fv = float(u), float(v)
        if not (fu == fu and fv == fv):
            return False
        if abs(fu - fv) > atol + rtol * abs(fv):
            return False
    return True


def read_fastq(path: Path):
    """Strictly parse FASTQ into (records, errors).

    A non-empty `errors` means the file is not valid FASTQ; the caller must treat
    that as a comparator failure. Stepping by four lines over `len - 3` silently
    discards a trailing partial record, which is how a truncated file compared
    equal to a complete one.
    """
    text = path.read_text(errors="replace")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        return [], [f"{path.name}: file is empty"]
    errors = []
    if len(lines) % 4 != 0:
        errors.append(f"{path.name}: {len(lines)} lines is not a multiple of 4 "
                      f"(trailing partial record)")
    records = []
    for i in range(0, len(lines) - 3, 4):
        header, seq, plus, qual = lines[i:i + 4]
        n = i // 4 + 1
        if not header.startswith("@"):
            errors.append(f"{path.name}: record {n}: header does not start with '@'")
            continue
        if not plus.startswith("+"):
            errors.append(f"{path.name}: record {n}: separator does not start with '+'")
            continue
        if len(seq) != len(qual):
            errors.append(f"{path.name}: record {n}: sequence length {len(seq)} != "
                          f"quality length {len(qual)}")
            continue
        records.append((header, seq, plus, qual))
    return records, errors


def compare_fastq(a: Path, b: Path):
    """Record-wise FASTQ comparison, refusing to discard a partial record."""
    ra, ea = read_fastq(a)
    rb, eb = read_fastq(b)
    if ea:
        return False, "; ".join(ea[:3])
    if eb:
        return False, "; ".join(eb[:3])
    if len(ra) != len(rb):
        return False, f"record count {len(ra)} vs {len(rb)}"
    for i, (x, y) in enumerate(zip(ra, rb)):
        if x != y:
            field = SAMPLE_FIELDS[next(k for k in range(4) if x[k] != y[k])]
            return False, f"record {i + 1} differs in {field}"
    return True, ""


def bam_signature(record) -> tuple:
    """Every field a scientific result can depend on, for one alignment.

    The earlier signature kept only `flag & 0xC0`, so a changed duplicate flag or a
    changed supplementary bit compared equal, and it omitted base qualities, mate
    coordinates and tags entirely -- which is why a BAM with corrupted qualities or
    NM tags passed the gate.

    Quality bytes are converted explicitly: pysam returns an `array.array`, which
    is unhashable and so cannot be part of a signature that gets counted.
    """
    return (
        record.query_name,
        record.flag,                       # every bit, not just 0xC0
        record.reference_name,
        record.reference_start,
        record.reference_end,
        record.cigarstring,
        record.mapping_quality,
        record.query_sequence,
        bytes(record.query_qualities) if record.query_qualities is not None else None,
        record.next_reference_id,
        record.next_reference_start,
        record.template_length,
        record.is_reverse,
        record.mate_is_reverse,
        record.mate_is_unmapped,
        tuple(sorted((tag, type(value).__name__, str(value))
                     for tag, value in record.get_tags())),
    )


def bam_header_signature(header: dict) -> dict:
    """Header facts that change what a consumer can do with the file.

    Only the reference dictionary and declared sort order. `@PG` lines legitimately
    name a different program in each implementation, so requiring them equal would
    fail a correct pair.
    """
    refs = []
    for entry in header.get("SQ", []) or []:
        if isinstance(entry, dict):
            refs.append((entry.get("SN"), entry.get("LN")))
        else:
            refs.append((getattr(entry, "name", None), getattr(entry, "length", None)))
    return {
        "references": tuple(refs),
        "sort_order": (header.get("HD") or {}).get("SO"),
    }


def _pysam():
    import pysam
    return pysam


def compare_bam(a: Path, b: Path, require_sorted: bool = False):
    """Compare BAMs as decoded records, not as compressed bytes.

    Compressed byte equality is meaningless across two writers: block boundaries,
    compression level and optional-field encoding all differ legitimately.

    `require_sorted=True` compares record order rather than a multiset. That is the
    right mode for a command whose contract is coordinate-sorted output; the default
    stays order-insensitive because several commands legitimately emit unsorted
    subsets, and for those per-file membership is not a meaningful claim.
    """
    pysam = _pysam()
    ka, kb, ha, hb = [], [], None, None
    for path, acc in ((a, ka), (b, kb)):
        with pysam.AlignmentFile(str(path)) as handle:
            if ha is None:
                ha = bam_header_signature(handle.header.to_dict())
            else:
                hb = bam_header_signature(handle.header.to_dict())
            for record in handle:
                acc.append(bam_signature(record))
    if ha is not None and hb is None:
        hb = ha
    if len(ka) != len(kb):
        return False, f"record count {len(ka)} vs {len(kb)}"

    if require_sorted:
        for i, (x, y) in enumerate(zip(ka, kb)):
            if x != y:
                return False, (f"record {i + 1} differs in sorted order: "
                               f"{x[:6]} vs {y[:6]}")
    else:
        from collections import Counter
        if Counter(ka) != Counter(kb):
            diff = Counter(ka) - Counter(kb)
            return False, (f"{sum(diff.values())} records differ; "
                           f"first: {list(diff)[0][:4]}")

    if ha != hb:
        if ha["references"] != hb["references"]:
            return False, (f"reference dictionary differs: "
                           f"{ha['references'][:2]} vs {hb['references'][:2]}")
        if ha["sort_order"] != hb["sort_order"]:
            return False, (f"declared sort order differs: "
                           f"{ha['sort_order']} vs {hb['sort_order']}")
    return True, ""


def numeric_cell_equal(left: str, right: str) -> bool:
    """Exact numeric equality for one table cell, refusing nonfinite on both sides.

    Decimal, not float: several commands print values that must match to every digit
    (`1.0` and `0.9999999999` are different results, not rounding noise), and float
    comparison would quietly accept them. The two implementations spell the same value
    identically in these tables, so any difference is worth reporting.

    A cell that does not parse as a number falls back to exact string equality, which
    keeps label columns working without a separate schema.
    """
    from decimal import Decimal, InvalidOperation

    try:
        lval, rval = Decimal(left), Decimal(right)
    except (InvalidOperation, ValueError):
        return left == right
    if not (lval.is_finite() and rval.is_finite()):
        return False
    return lval == rval


def compare_numeric_table(a: bytes, b: bytes):
    """Compare two whitespace-delimited numeric tables held in memory.

    Shared with the differential runner, which reports a diff line by line and names
    the file; this returns (ok, reason) so both harnesses decide equality the same way.
    Any nonfinite cell fails, in either table, including when both sides produce the
    same one.
    """
    # Accept str as well as bytes: the differential runner holds already-decoded text
    # after its own gzip/path normalisation, the benchmark holds raw bytes. Rejecting
    # one of them would push each caller back to its own private copy.
    text = (lambda v: v if isinstance(v, str) else v.decode("utf-8", errors="replace"))
    a_rows = text(a).splitlines()
    b_rows = text(b).splitlines()
    if len(a_rows) != len(b_rows):
        return False, f"row count {len(a_rows)} vs {len(b_rows)}"
    for i, (ra, rb) in enumerate(zip(a_rows, b_rows)):
        ca, cb = ra.split(), rb.split()
        if len(ca) != len(cb):
            return False, f"row {i + 1}: cell count {len(ca)} vs {len(cb)}"
        for token in ca + cb:
            if is_nonfinite(token):
                return False, f"row {i + 1}: nonfinite value {token!r}"
        for x, y in zip(ca, cb):
            if not numeric_cell_equal(x, y):
                return False, f"row {i + 1}: {x!r} vs {y!r}"
    return True, ""


def write_bam(path: Path, records, header=None) -> None:
    """Write a BAM from plain dicts. Used by the tests to build known-good inputs."""
    pysam = _pysam()
    header = header or {"HD": {"VN": "1.6", "SO": "coordinate"},
                        "SQ": [{"SN": "chr1", "LN": 100000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in records:
            a = pysam.AlignedSegment()
            a.query_name = r["name"]
            a.query_sequence = r.get("seq", "ACGTACGTAC")
            a.flag = r.get("flag", 0)
            a.reference_id = r.get("ref_id", 0)
            a.reference_start = r.get("start", 100)
            a.mapping_quality = r.get("mapq", 60)
            a.cigar = r.get("cigar", [(0, 10)])
            a.query_qualities = pysam.qualitystring_to_array(r.get("qual", "IIIIIIIIII"))
            a.next_reference_id = r.get("next_ref", -1)
            a.next_reference_start = r.get("next_start", -1)
            a.template_length = r.get("tlen", 0)
            for tag, value in r.get("tags", {}).items():
                a.set_tag(tag, value, value_type="i")
            out.write(a)
