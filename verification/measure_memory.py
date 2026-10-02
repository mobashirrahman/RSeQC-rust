#!/usr/bin/env python3
"""Measure memory against input size for each cost driver the audit named.

The audit's P1 production-envelope finding has three parts, and this script
addresses the first two by measurement and the third from the source:

1. *The CRAM branch collects the entire file before processing.* -- a source fact,
   reported by `source_facts()` rather than inferred from timing, because a 10-record
   CRAM fixture cannot demonstrate anything.
2. *Duplicate counting retains unique sequences and positions.* -- measured.
3. *WIG generation retains covered positions.* -- measured.

What is measured is peak RSS against alignment size, for the command whose algorithm
retains that structure. The point is to find which commands grow with their retained
set rather than with the input stream, and to record the measured slope instead of
assuming one. A `streams` driver is included as a control: if the streaming command's
RSS is not flat against input size, the streaming claim is wrong.

Two sizes come from one fixture by restricting to a contig, using the project's own
`bam2wig` output to build a chrom.sizes file, and reading back the BAM's index with
pysam to fetch a range. Sizes are counted by records so the denominator is exact.

Each measurement runs `/usr/bin/time -v`, whose "Maximum resident set size" is the
largest SINGLE child's RSS (Linux documents `RUSAGE_CHILDREN.ru_maxrss` that way).
These commands are single-process, so it is the quantity of interest; a helper-heavy
workflow needs the separate aggregate measure noted in benchmarks/protocol-v2.md §11.

Usage:
    oracle/venv/bin/python3 verification/measure_memory.py
    oracle/venv/bin/python3 verification/measure_memory.py --driver duplication --quick
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RELEASE = REPO / "target" / "release"
TINY = REPO / "datasets" / "aligned" / "tiny"

# Drivers whose cost follows TRANSCRIPT count rather than alignment-record count.
#
# These are measured separately because restricting an alignment by coordinate does not
# change how many transcripts a command builds state for: `geneBody_coverage` computes
# a percentile list per transcript, `tin` holds an exon index per transcript, and
# `junction_annotation` holds the reference's intron sets. Running them through the
# per-record sweep would report a slope of nearly zero and a reader would conclude they
# stream, which is a statement about the sweep rather than about the command.
#
# Swept by varying the BED12, holding the alignment fixed.
TRANSCRIPT_DRIVERS = {
    "genebody_coverage": {
        "command": "geneBody_coverage",
        "retains": "a percentile list per transcript",
        "interpretation": (
            "Per-transcript state, so the cost should follow the transcript count and "
            "be nearly flat in alignment size. A rising slope in records here would "
            "mean the command reads the alignment more than once, or retains per-read "
            "state the sliding window was meant to avoid."),
    },
    "tin": {
        "command": "tin",
        "retains": "an exon block index per transcript",
        "interpretation": (
            "Per-transcript state. The read-side window was already fixed to be "
            "streaming (see CHANGELOG), so this isolates the annotation cost that "
            "remains regardless of how many reads arrive."),
    },
    "junction_annotation": {
        "command": "junction_annotation",
        "retains": "the reference model's intron start/end sets per contig",
        "interpretation": (
            "Cost follows the reference annotation's intron count, not the alignment. "
            "The record-side loop streams, so this measures the model side alone."),
    },
}

DRIVERS = {
    "streams": {
        "command": "bam_stat",
        "retains": "nothing beyond the current record",
        "interpretation": (
            "Control case. A streaming command's peak RSS should be flat against "
            "alignment size. A rising slope here would mean the streaming claim is "
            "false, whatever the code comments say."),
    },
    "duplication": {
        "command": "read_duplication",
        "retains": "unique sequences and their positions",
        "interpretation": (
            "Duplicate detection keeps one entry per unique sequence/position, so RSS "
            "should track the number of DISTINCT reads rather than the alignment "
            "size. A slope near one would mean every read is retained."),
    },
    "wig": {
        "command": "bam2wig",
        "retains": "covered positions",
        "interpretation": (
            "WIG generation keeps a per-position value for every covered base, so RSS "
            "is proportional to covered positions. This is why bam2wig cannot process "
            "a whole-genome alignment at ordinary coverage within a modest memory "
            "budget, and it is the driver that bounds the declared envelope."),
    },
}


# How many times each sweep point is measured, so the instrument's own repeatability
# can be established per command rather than assumed. `/proc/meminfo` has a ~90 MB
# noise floor on this class of machine (measured by `measure_concurrency.py`), but
# that is a SYSTEM-WIDE figure; a single process's own RSS via `/usr/bin/time -v` is
# far more repeatable, and a fixed threshold derived from the wrong instrument's noise
# would discard a real monotone signal. Measuring the repeatability directly is the
# only honest way to set the bar.
SWEEP_REPEATS = 3

# A fitted slope is trusted only when the quantity it was fitted from moves by more
# than this multiple of the sweep's own observed repeatability. A monotone, repeatable
# rise is signal; a rise inside the noise is the denominator's range expressed as a
# number.
MIN_SIGNAL_TO_NOISE = 3.0


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_time_v(stderr: str) -> dict:
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", stderr)
    if m:
        secs = 0.0
        for p in m.group(1).split(":"):
            secs = secs * 60 + float(p)
        out["wall_s"] = secs
    m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
    if m:
        out["peak_rss_mb"] = int(m.group(1)) / 1024.0
    m = re.search(r"User time \(seconds\):\s*([\d.]+)", stderr)
    if m:
        out["user_s"] = float(m.group(1))
    return out


def measure(argv, cwd, timeout=5400):
    proc = subprocess.run(["/usr/bin/time", "-v"] + [str(a) for a in argv],
                          cwd=str(cwd), capture_output=True, text=True,
                          timeout=timeout)
    m = parse_time_v(proc.stderr)
    m["exit_code"] = proc.returncode
    if proc.returncode != 0:
        m["stderr_tail"] = proc.stderr[-300:]
    return m


def record_count(path):
    """Exact alignment-record count, via htslib."""
    import pysam
    with pysam.AlignmentFile(str(path), "rb") as f:
        return sum(1 for _ in f.fetch(until_eof=True))


def covered_bases(path):
    """Distinct reference positions covered by any alignment."""
    import pysam
    with pysam.AlignmentFile(str(path), "rb") as f:
        covered = set()
        for r in f.fetch(until_eof=True):
            if r.is_unmapped or r.cigartuples is None:
                continue
            ref = r.reference_start
            for op, n in r.cigartuples:
                if op in (0, 2, 7, 8):  # M, D, =, X consume reference
                    covered.update(range(ref, ref + n))
                    ref += n
    return len(covered)


def restrict_to_range(src, dst, contig, start, end):
    """A real BAM covering `[start, end)` of one contig, written with pysam.

    Written here rather than shelled out to a project command, because the commands
    under test must not also be the tool that builds their own input. The index makes
    this cheap, and restricting by coordinate rather than by contig matters: the tiny
    fixture carries every record on one contig, so a per-contig split would yield one
    empty subset and measure nothing.
    """
    import pysam
    with pysam.AlignmentFile(str(src), "rb") as s, \
            pysam.AlignmentFile(str(dst), "wb", template=s) as d:
        n = 0
        for r in s.fetch(contig, start, end):
            d.write(r)
            n += 1
    return n


def chrom_sizes_for(bam, dst):
    import pysam
    with pysam.AlignmentFile(str(bam), "rb") as f:
        lines = [f"{e.get('SN')}\t{e.get('LN')}"
                 for e in f.header.to_dict().get("SQ", [])]
    dst.write_text("\n".join(lines) + "\n")
    return dst


def build_argv(command, exe, bam, bed, sizes, work):
    argv = [exe]
    if command == "bam2wig":
        # bam2wig needs all three of input, chrom sizes and an output prefix. Omitting
        # the prefix makes clap exit 2 immediately, which measures argument parsing
        # rather than memory.
        return argv + ["-i", str(bam), "-s", str(sizes),
                       "-o", str(work / "w")]
    argv += ["-i", str(bam)]
    if command in ("geneBody_coverage", "junction_annotation", "tin"):
        argv += ["-r", str(bed)]
    if command == "tin":
        argv += ["-n", "50", "-o", str(work)]
    elif command in ("read_duplication", "geneBody_coverage", "junction_annotation"):
        argv += ["--out-prefix", str(work / "out")]
    return argv


def slope(points, key="peak_rss_mb", denom="records"):
    """Least-squares slope of `key` against `denom`, and the observed range.

    Reported with both endpoints rather than as a bare slope, because a slope over a
    narrow range of one size point says nothing.
    """
    xs = [p[denom] for p in points if p.get(key) and p.get(denom)]
    ys = [p[key] for p in points if p.get(key) and p.get(denom)]
    if len(xs) < 2:
        return None
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    denom_v = sum((x - mx) ** 2 for x in xs)
    if denom_v == 0:
        return None
    m = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom_v
    return {
        "mb_per_unit": round(m, 9),
        "bytes_per_record": round(m * 1024 * 1024, 3),
        "bytes_per_unit": round(m * 1024 * 1024, 3),
        "denominator": denom,
        "records_range": [min(xs), max(xs)],
        "rss_range_mb": [round(min(ys), 1), round(max(ys), 1)],
    }


def transcript_sweep(name: str, spec: dict, bam: Path, bed: Path, work: Path,
                     min_span: float) -> dict:
    """Measure a per-transcript driver while varying the TRANSCRIPT count.

    The BED is subsampled by taking a prefix of its rows, and the transcript count of
    each sample is COUNTED from the resulting file rather than assumed from the
    fraction requested. That distinction is the whole point: the audit's protocol
    section 3 warns that an earlier scaling file labelled a point "20,000 transcripts"
    when it contained 9,179, and a sweep that repeats that would produce the same class
    of mislabelled figure. Each row here carries the measured count.
    """
    exe = RELEASE / spec["command"]
    if not exe.exists():
        return {"status": "SKIPPED", "reason": f"{exe} missing"}
    rows = [line for line in bed.read_text().splitlines()
            if line.strip() and not line.startswith(("#", "track", "browser"))]
    if not rows:
        return {"status": "SKIPPED", "reason": f"{bed} has no rows"}
    points = []
    for frac in (0.0625, 0.125, 0.25, 0.5, 1.0):
        keep = max(1, int(len(rows) * frac))
        d = work / f"{name}-t{keep}"
        d.mkdir(parents=True, exist_ok=True)
        subset = d / "model.bed12"
        subset.write_text("\n".join(rows[:keep]) + "\n")
        measured = len([l for l in subset.read_text().splitlines() if l.strip()])
        argv = build_argv(spec["command"], exe, bam, subset, d / "sizes", d)
        runs = []
        for rep in range(SWEEP_REPEATS):
            out = d / f"rep{rep}"
            out.mkdir(parents=True, exist_ok=True)
            runs.append(measure(build_argv(spec["command"], exe, bam, subset,
                                           d / "sizes", out), out))
        peaks = [r.get("peak_rss_mb", 0.0) for r in runs if r["exit_code"] == 0]
        # The PEAK across repeats, not the mean: peak RSS is what the OS is asked to
        # guarantee, and an allocator that occasionally spikes is a real property of
        # the command. The spread across repeats is the noise floor.
        peak = max(peaks) if peaks else 0.0
        spread = (max(peaks) - min(peaks)) if len(peaks) > 1 else 0.0
        best = runs[0] if runs else {"exit_code": None}
        points.append({"transcripts_requested_fraction": frac,
                       "transcripts": measured,
                       "repeats": len(runs),
                       "peak_rss_mb_across_repeats": [round(p, 1) for p in peaks],
                       "repeat_spread_mb": round(spread, 1),
                       "peak_rss_mb": peak,
                       "exit_code": best["exit_code"],
                       "wall_s": min((r.get("wall_s", 0.0) for r in runs
                                      if r["exit_code"] == 0), default=0.0)})
        print(f"  {name:22} transcripts={measured:>7,} "
              f"rss={peak:8.1f} MB (spread {spread:5.1f}) "
              f"exit={best['exit_code']}")
    entry = {"command": spec["command"], "retains": spec["retains"],
             "interpretation": spec["interpretation"],
             "alignment": str(bam), "alignment_records": record_count(bam),
             "points": points}
    counts = [pt["transcripts"] for pt in points]
    span = (max(counts) / max(1, min(counts))) if counts else 0
    entry["transcript_span_factor"] = round(span, 3)
    entry["all_runs_succeeded"] = all(pt["exit_code"] == 0 for pt in points)
    if not entry["all_runs_succeeded"]:
        # A cost per transcript fitted over runs that never produced a result is
        # exactly the class of figure this script was written to stop publishing. It
        # happened: every point failed instantly and the sweep still reported
        # 1234.5 bytes/transcript.
        entry["slope_status"] = (
            "NOT REPORTED: at least one run failed, so no per-transcript cost is "
            f"computed from these points ({sum(1 for pt in points if pt['exit_code'])} "
            f"of {len(points)} failed). The failures are: "
            + "; ".join(pt.get("stderr_tail", "")[-90:].replace("\n", " ")
                        for pt in points if pt["exit_code"] != 0)[:400]
        )
        return entry
    if len(points) >= 2 and span >= min_span:
        s = slope(points, "peak_rss_mb", "transcripts")
        # A per-unit cost is only interpretable if the quantity it was fitted from
        # actually MOVED. A 16x change in transcript count that moves peak RSS by
        # 0.6 MB produces a slope that is arithmetically fine and evidentially
        # worthless: the line is being drawn through measurement noise, and the
        # number's precision reflects the denominator's range rather than the
        # numerator's. The threshold is a multiple of the instrument's own noise floor
        # (measured, not assumed, in the concurrency tool) so "moved" means "moved by
        # more than the measurement can resolve".
        movement = s["rss_range_mb"][1] - s["rss_range_mb"][0]
        s["rss_movement_mb"] = round(movement, 1)
        noise = max((pt.get("repeat_spread_mb", 0.0) for pt in points), default=0.0)
        s["repeat_noise_mb"] = round(noise, 1)
        s["signal_to_noise"] = round(abs(movement) / noise, 1) if noise else None
        entry["repeat_noise_mb"] = round(noise, 1)
        if abs(movement) < MIN_SIGNAL_TO_NOISE * noise:
            entry["slope_status"] = (
                f"WITHHELD: the transcript counts span {span:.2f}x but peak RSS moves "
                f"only {abs(movement):.1f} MB, against a repeat-to-repeat spread of "
                f"{noise:.1f} MB measured by running each point "
                f"{SWEEP_REPEATS} times. The signal is "
                f"{abs(movement) / noise:.1f}x the noise where "
                f"{MIN_SIGNAL_TO_NOISE}x is required, so a per-transcript cost fitted "
                f"to it would express the denominator's range rather than a measured "
                f"cost. This command's per-transcript state is smaller than this "
                f"measurement can resolve on this alignment."
            )
            entry["slope_vs_transcripts_raw"] = s
        elif s["mb_per_unit"] < 0:
            # A negative per-unit cost is not a cost. It means two effects are moving
            # in opposite directions and the larger one is measured: with a SMALL
            # gene model, the commands that hold a read-side window (geneBody_coverage
            # and tin) cannot retire a read until every transcript it might belong to
            # has been scored, so a smaller model retains MORE reads. That read
            # retention dominates the annotation cost the sweep is trying to measure.
            #
            # Reporting the fitted number would publish "memory falls as the
            # annotation grows" as though it were a property of the annotation. The
            # interaction is real and worth recording; the slope is not interpretable
            # on its own, so it is withheld and the confound named.
            # Whether the read-side window explains the sign is a claim about the
            # command, so it is only made for the commands that actually hold one.
            # `junction_annotation` does not: it scores junctions in a single streaming
            # pass, so its annotation cost should dominate and a small movement with
            # model size is the annotation cost being small in absolute terms. The
            # generic explanation would have been wrong for it, which is why the text
            # is per-command rather than shared.
            WINDOW_HOLDERS = {"geneBody_coverage", "tin"}
            if spec["command"] in WINDOW_HOLDERS:
                why = (
                    "The read-side window these commands hold cannot retire a read "
                    "until every transcript it might belong to is scored, so a "
                    "SMALLER model retains MORE reads. That effect runs opposite to "
                    "the annotation cost and dominates it here, so the fitted value "
                    "measures the confound, not the annotation. Separating the two "
                    "needs the read stream held fixed as well as the model."
                )
            else:
                why = (
                    "This command scores its records in a single streaming pass and "
                    "holds no read-side window, so the annotation cost should "
                    "dominate. Peak RSS is essentially flat across a 16x change in "
                    "transcript count, which says the per-transcript state is small "
                    "next to the fixed cost of the input -- not that there is no "
                    "per-transcript cost. Isolating it needs a range of model sizes "
                    "large enough for the annotation cost to exceed the constant."
                )
            entry["slope_status"] = (
                f"WITHHELD: the fitted slope is negative "
                f"({s['bytes_per_unit']:.1f} bytes/transcript), which is not a cost. "
                f"Peak RSS moves from {s['rss_range_mb'][1]:.0f} MB to "
                f"{s['rss_range_mb'][0]:.0f} MB as the model grows from "
                f"{s['records_range'][0]:,} to {s['records_range'][1]:,} transcripts. "
                f"{why}"
            )
            entry["slope_vs_transcripts_raw"] = s
        else:
            entry["slope_vs_transcripts"] = s
    elif len(points) >= 2:
        entry["slope_status"] = (
            f"NOT REPORTED: the transcript counts span only {span:.2f}x "
            f"({min(counts):,}..{max(counts):,}), below the {min_span}x minimum")
    return entry


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", action="append", default=None, choices=sorted(DRIVERS))
    ap.add_argument("--quick", action="store_true",
                    help="one input size per driver instead of two")
    ap.add_argument("--json", default=None)
    ap.add_argument("--bam", default=None,
                    help="measure against this alignment instead of the tiny "
                         "fixture. Required for the production pilot: the fixture's "
                         "621,783 records cannot show what a command does at 10M or "
                         "50M pairs, and a slope fitted over a 22%% range is an "
                         "extrapolation dressed as a measurement.")
    ap.add_argument("--bed", default=None,
                    help="BED12 model matching --bam (default: the fixture's)")
    ap.add_argument("--fractions", type=float, nargs="+", default=None,
                    help="coordinate-window fractions of the busiest contig to "
                         "measure (default: 0.25 0.5 1.0). On a real alignment the "
                         "reads are spread along the contig, so many small windows "
                         "give genuinely distinct sizes; on the tiny fixture they "
                         "do not, which is why the default sweep cannot produce a "
                         "slope there.")
    ap.add_argument("--transcript-drivers", nargs="*", default=None,
                    choices=sorted(TRANSCRIPT_DRIVERS),
                    help="measure these per-transcript drivers by varying the BED12 "
                         "while holding the alignment fixed, instead of (or in "
                         "addition to) the per-record sweep")
    ap.add_argument("--min-span", type=float, default=2.0,
                    help="refuse to report a per-unit cost unless the measured inputs "
                         "span at least this factor (default 2.0). A slope fitted over "
                         "a 1.2x range cannot distinguish a line from a step, and "
                         "reporting one anyway is how the first version of this script "
                         "produced a confident 135.6 bytes/record from two points 22%% "
                         "apart.")
    args = ap.parse_args()

    if not RELEASE.exists():
        print("target/release not found; run: cargo build --workspace --release --locked")
        return 2
    # Resolved, because every command runs with its cwd set to a per-run scratch
    # directory. A relative `--bam` therefore resolved against that scratch directory
    # and not the invocation directory: the per-record sweep was unaffected because it
    # writes its windows to absolute paths first, but the per-transcript sweep passes
    # the alignment through unchanged, and every one of its runs failed instantly with
    # "No such file or directory" -- which the slope function then fitted, producing a
    # confident 1234.5 bytes/transcript from fifteen failed runs.
    bam = (Path(args.bam) if args.bam else TINY / "tiny.bam").resolve()
    bed = (Path(args.bed) if args.bed else TINY / "model.bed12").resolve()
    if not bam.exists():
        print(f"alignment missing: {bam}")
        return 2
    if not bed.exists():
        print(f"annotation missing: {bed}")
        return 2

    results = {
        "source_facts": source_facts(),
        "transcript_drivers": {},
        "input": {
            "bam": str(bam),
            "bed": str(bed),
            "bam_sha256": file_sha256(bam),
            "bed_sha256": file_sha256(bed),
        },
        "drivers": {},
    }

    import pysam
    with pysam.AlignmentFile(str(bam)) as f:
        contigs = [(e.get("SN"), e.get("LN")) for e in f.header.to_dict().get("SQ", [])]
    results["input"]["total_records"] = record_count(bam)
    results["input"]["total_covered_bases"] = covered_bases(bam)
    print(f"input: {bam}")
    print(f"  sha256 {results['input']['bam_sha256']}")
    print(f"  {results['input']['total_records']:,} records, "
          f"{results['input']['total_covered_bases']:,} covered bases")
    print()

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)

        # The per-transcript sweep first, and outside the per-record machinery, because
        # it varies a different thing and must not inherit the contig-windowing.
        if args.transcript_drivers:
            specs = {k: v for k, v in TRANSCRIPT_DRIVERS.items()
                     if k in args.transcript_drivers}
            print(f"per-transcript sweep: holding {bam.name} "
                  f"({results['input']['total_records']:,} records) fixed, varying "
                  f"the model")
            for name, spec in specs.items():
                entry = transcript_sweep(name, spec, bam, bed, work, args.min_span)
                results["transcript_drivers"][name] = entry
                s = entry.get("slope_vs_transcripts")
                if s:
                    print(f"  {name:22} {s['bytes_per_unit']:8.1f} bytes/transcript "
                          f"over {s['records_range'][0]:,}..{s['records_range'][1]:,} "
                          f"transcripts ({s['rss_range_mb']} MB)")
                elif entry.get("slope_status"):
                    print(f"  {name:22} slope {entry['slope_status']}")
                elif entry.get("status") == "SKIPPED":
                    print(f"  {name:22} SKIPPED ({entry['reason']})")
            print()

        drivers = args.driver or sorted(DRIVERS)
        with pysam.AlignmentFile(str(bam), "rb") as f:
            pass
        # Input sizes from the one fixture, by restricting a coordinate range on the
        # contig that actually carries reads. The tiny panel puts every record on a
        # single contig, so splitting by contig would give one empty subset; the
        # resulting empty-BAM error would then be measured as if it were memory.
        # The contig that actually carries records, chosen by COUNTING them rather
        # than by length: the tiny panel's longest contig is chr1 (249 Mb) while every
        # read sits on chr17, so choosing by length would restrict to an empty range
        # and measure the empty-input path three times.
        busy, best, length = None, -1, 0
        for name, ln in contigs:
            with pysam.AlignmentFile(str(bam), "rb") as probe:
                n = sum(1 for _ in probe.fetch(name))
            if n > best:
                busy, best, length = (name, ln), n, ln
        contig, length = busy[0], busy[1]
        if args.fractions:
            points_frac = sorted(set(args.fractions))
        elif args.quick:
            points_frac = [0.25]
        else:
            points_frac = [0.25, 0.5, 1.0]
        inputs = []
        for frac in points_frac:
            end_coord = int(length * frac)
            dst = work / f"{contig}-{frac}.bam"
            restrict_to_range(bam, dst, contig, 0, end_coord)
            inputs.append({
                "label": f"{contig}:0-{end_coord:,} ({frac:.0%})",
                "path": dst,
                "records": record_count(dst),
                "covered_bases": covered_bases(dst),
            })
        inputs.sort(key=lambda x: x["records"])

        # Windows that turn out to contain the same number of records are the same
        # measurement, and a slope fitted through a repeated point is not a fit at
        # all. The tiny fixture produces exactly this: reads are concentrated, so
        # the 50% and 100% windows both hold all 621,783 records and the earlier
        # version of this script reported a confident bytes/record figure over the
        # resulting 1.22x "range". Duplicates are dropped here rather than in the
        # slope function so the printed table cannot show the same point twice.
        deduped = []
        for i in inputs:
            if deduped and i["records"] == deduped[-1]["records"]:
                print(f"  note: dropping {i['label']} -- it holds the same "
                      f"{i['records']:,} records as {deduped[-1]['label']}, so it is "
                      f"a repeat of that measurement, not a new size")
                continue
            deduped.append(i)
        inputs = deduped
        sizes = chrom_sizes_for(bam, work / "chrom.sizes")

        widest = max((i["records"] for i in inputs), default=0)
        for i in inputs:
            i["nested"] = i["records"] < widest
        print(f"{'input':22} {'records':>10} {'covered bases':>15}")
        for i in inputs:
            note = " (prefix of a larger window)" if i.get("nested") else ""
            print(f"{i['label']:22} {i['records']:>10,} {i['covered_bases']:>15,}{note}")
        if inputs and all(i.get("nested") or len(inputs) == 1 for i in inputs):
            print("\n  note: every window is a prefix of a larger one, so these are "
                  "nested\n        measurements of the same alignment rather than "
                  "independent sizes. A slope across nested\n        windows is a "
                  "valid cost-per-record estimate, but it says nothing about "
                  "behaviour\n        that changes with alignment size, and it is "
                  "labelled as nested above.")
        print()

        for name in drivers:
            d = DRIVERS[name]
            exe = RELEASE / d["command"]
            if not exe.exists():
                results["drivers"][name] = {"status": "SKIPPED",
                                            "reason": f"{exe} missing"}
                print(f"{name}: SKIPPED ({exe} missing)")
                continue
            points = []
            for i in inputs:
                out_dir = work / f"{name}-{i['records']}"
                out_dir.mkdir(parents=True, exist_ok=True)
                runs = []
                for rep in range(SWEEP_REPEATS):
                    out = out_dir / f"rep{rep}"
                    out.mkdir(parents=True, exist_ok=True)
                    runs.append(measure(
                        build_argv(d["command"], exe, i["path"], bed, sizes, out),
                        out))
                peaks = [r.get("peak_rss_mb", 0.0) for r in runs
                         if r["exit_code"] == 0]
                peak = max(peaks) if peaks else 0.0
                spread = (max(peaks) - min(peaks)) if len(peaks) > 1 else 0.0
                first = runs[0] if runs else {"exit_code": None}
                points.append({"label": i["label"], "records": i["records"],
                               "covered_bases": i["covered_bases"],
                               "repeats": len(runs),
                               "peak_rss_mb_across_repeats":
                                   [round(x, 1) for x in peaks],
                               "repeat_spread_mb": round(spread, 1),
                               "peak_rss_mb": peak,
                               "exit_code": first["exit_code"],
                               "wall_s": min((r.get("wall_s", 0.0) for r in runs
                                              if r["exit_code"] == 0),
                                             default=0.0)})
                print(f"  {name:12} {i['label']:18} records={i['records']:>9,} "
                      f"rss={peak:8.1f} MB (spread {spread:5.1f}) "
                      f"exit={first['exit_code']}")
            entry = {
                "command": d["command"],
                "retains": d["retains"],
                "interpretation": d["interpretation"],
                "points": points,
            }
            recs = [p["records"] for p in points]
            span = (max(recs) / max(1, min(recs))) if recs else 0
            entry["input_span_factor"] = round(span, 3)
            entry["all_runs_succeeded"] = all(p["exit_code"] == 0 for p in points)
            if not entry["all_runs_succeeded"]:
                entry["slope_status"] = (
                    "NOT REPORTED: at least one run failed, so no per-record cost is "
                    f"computed from these points ({sum(1 for p in points if p['exit_code'])} "
                    f"of {len(points)} failed)"
                )
                results["drivers"][name] = entry
                continue
            if len(points) >= 2 and span >= args.min_span:
                s = slope(points, "peak_rss_mb", "records")
                movement = s["rss_range_mb"][1] - s["rss_range_mb"][0]
                noise = max((p.get("repeat_spread_mb", 0.0) for p in points),
                            default=0.0)
                s["rss_movement_mb"] = round(movement, 1)
                s["repeat_noise_mb"] = round(noise, 1)
                s["signal_to_noise"] = (round(abs(movement) / noise, 1)
                                        if noise else None)
                entry["repeat_noise_mb"] = round(noise, 1)
                if abs(movement) < MIN_SIGNAL_TO_NOISE * noise:
                    entry["slope_status"] = (
                        f"WITHHELD: peak RSS moves {abs(movement):.1f} MB across the "
                        f"range against a repeat-to-repeat spread of {noise:.1f} MB, "
                        f"below the {MIN_SIGNAL_TO_NOISE}x signal-to-noise this "
                        f"requires. The command is flat to within this measurement's "
                        f"own repeatability over these inputs; a per-record cost "
                        f"fitted to it would be the denominator's range."
                    )
                    entry["slope_vs_records_raw"] = s
                elif s["mb_per_unit"] < 0:
                    # Same rule as the per-transcript sweep: a negative per-record
                    # cost is a measurement that does not support a cost claim.
                    entry["slope_status"] = (
                        f"WITHHELD: fitted slope is negative "
                        f"({s['bytes_per_record']:.1f} bytes/record). A streaming "
                        f"command can drift slightly below zero across a narrow range; "
                        f"that is flatness within noise, not a cost, and it is "
                        f"reported as such rather than as a number."
                    )
                    entry["slope_vs_records_raw"] = s
                else:
                    entry["slope_vs_records"] = s
                    entry["slope_vs_covered_bases"] = slope(
                        points, "peak_rss_mb", "covered_bases")
                    print(f"  {name:12} slope fitted over {span:.2f}x "
                          f"({min(recs):,}..{max(recs):,} records)")
            elif len(points) >= 2:
                entry["slope_status"] = (
                    f"NOT REPORTED: inputs span only {span:.2f}x "
                    f"({min(recs):,}..{max(recs):,} records), below the "
                    f"{args.min_span}x minimum. A slope over a range this narrow "
                    f"cannot distinguish a linear cost from a step or a saturating "
                    f"one, so no bytes/record figure is published for it."
                )
                print(f"  {name:12} slope NOT reported ({span:.2f}x < "
                      f"{args.min_span}x minimum)")
            else:
                entry["slope_status"] = ("NOT REPORTED: fewer than two completed "
                                         "sizes")
            results["drivers"][name] = entry

    print()
    print("Source facts (not measurements):")
    for k, v in results["source_facts"].items():
        print(f"  {k}\n    {v}")
    print("\nMeasured slopes:")
    for name, e in results["drivers"].items():
        s = e.get("slope_vs_records")
        if s:
            print(f"  {name:12} {s['bytes_per_record']:.1f} bytes/record over "
                  f"{s['records_range'][0]:,}..{s['records_range'][1]:,} records "
                  f"({s['rss_range_mb']} MB)")
        elif "status" in e:
            print(f"  {name:12} {e['status']}")
        elif e.get("slope_status"):
            print(f"  {name:12} {e['slope_status']}")

    if results["transcript_drivers"]:
        print("\nPer-transcript costs (alignment held fixed, model varied):")
        for name, e in results["transcript_drivers"].items():
            s = e.get("slope_vs_transcripts")
            if s:
                print(f"  {name:22} {s['bytes_per_unit']:.1f} bytes/transcript over "
                      f"{s['records_range'][0]:,}..{s['records_range'][1]:,} "
                      f"transcripts ({s['rss_range_mb']} MB)")
            elif e.get("slope_status"):
                print(f"  {name:22} {e['slope_status']}")
            elif e.get("status") == "SKIPPED":
                print(f"  {name:22} SKIPPED ({e['reason']})")

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2) + "\n")
        print(f"\nwritten to {args.json}")

    print("\nThis is a starting measurement, not a declared envelope. A published "
          "memory limit\nrequires the measured slope at a named scale plus the "
          "observed failure behaviour beyond it;\nneither is established by a single "
          "sweep over one fixture.")
    return 0


def source_facts():
    """Facts read from the source, for what a timing run cannot demonstrate."""
    facts = {}
    lib = (REPO / "crates/formats/src/lib.rs").read_text()
    if "CramBuffered" in lib:
        facts["cram_decode_is_whole_file"] = (
            "`AlignmentRecords::CramBuffered` collects the entire file before "
            "processing: `noodles-cram` 0.99 exposes record iteration only as "
            "`records(&header)`, which is single-use, so re-entering a drained reader "
            "yields a spurious decode error instead of EOF. This bounds CRAM input by "
            "available memory rather than by alignment size. It is a code fact read "
            "from crates/formats/src/lib.rs, not an inference from timing: the CRAM "
            "fixture is 10 records and could not demonstrate it.")
    dup = (REPO / "crates/commands/src/read_duplication.rs").read_text()
    if "HashSet" in dup or "HashMap" in dup:
        facts["duplication_retains_unique_entries"] = (
            "read_duplication keeps one entry per unique sequence/position in a hash "
            "container, so its memory is proportional to DISTINCT reads. The measured "
            "slope below is the evidence for how that behaves on a real alignment.")
    wig = (REPO / "crates/commands/src/bam2wig.rs").read_text()
    if "Vec<" in wig or "HashMap" in wig:
        facts["wig_retains_covered_positions"] = (
            "bam2wig materialises a per-position value map over covered bases, so its "
            "memory is proportional to covered positions and cannot be streamed down. "
            "The measured slope below quantifies it on this fixture.")
    return facts


if __name__ == "__main__":
    sys.exit(main())
