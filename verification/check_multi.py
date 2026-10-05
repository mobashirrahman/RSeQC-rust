#!/usr/bin/env python3
"""Single-pass driver check (card C1).

Runs each registered command alone and through `rseqc_multi` on the same
input and compares every file and stream byte for byte. Fails when a file
is missing on either side, when any byte differs, or when exit codes
differ. Each side runs in its own directory with a relative `-o out` so
prefix-embedding artifacts (e.g. the `.GC_plot.r` pdf path) compare
exactly.

Usage:
    python3 verification/check_multi.py --bin-dir target/release --work DIR

Exit 0: all comparisons identical. Exit 1: any mismatch (first bytes
shown), with the failing comparison named.
"""

import argparse
import difflib
import random
import shutil
import subprocess
import sys
from pathlib import Path

COMMANDS = ("bam_stat", "read_GC", "read_NVC", "read_quality")


def build_fixture(path: Path, n_reads: int = 300, seed: int = 7) -> None:
    """Deterministic small BAM: 101 bp reads stepping along chr1."""
    import pysam

    rng = random.Random(seed)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        pos = 100
        for i in range(n_reads):
            read = pysam.AlignedSegment()
            read.query_name = f"r{i}"
            read.query_sequence = "".join(rng.choice("ACGT") for _ in range(101))
            read.flag = 0
            read.reference_id = 0
            read.reference_start = pos
            read.mapping_quality = 40
            read.cigarstring = "101M"
            read.query_qualities = pysam.qualitystring_to_array("I" * 101)
            out.write(read)
            pos += 150


def run(cmd, cwd: Path):
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True)
    return proc.returncode, proc.stdout, proc.stderr


def check_equal(label: str, solo: bytes, multi: bytes, failures: list) -> None:
    if solo != multi:
        failures.append(label)
        print(f"MISMATCH: {label} ({len(solo)} vs {len(multi)} bytes)")
        a = solo.decode("utf-8", "replace").splitlines()[:20]
        b = multi.decode("utf-8", "replace").splitlines()[:20]
        for line in difflib.unified_diff(a, b, "solo", "multi", lineterm=""):
            print(f"  {line}")
    else:
        print(f"identical: {label} ({len(solo)} bytes)")


def check_files(label: str, solo_dir: Path, multi_dir: Path, names: list, failures: list) -> None:
    for name in names:
        solo_file, multi_file = solo_dir / name, multi_dir / name
        tag = f"{label}/{name}"
        if not solo_file.exists():
            failures.append(tag)
            print(f"MISSING on solo side: {tag}")
            continue
        if not multi_file.exists():
            failures.append(tag)
            print(f"MISSING on multi side: {tag}")
            continue
        check_equal(tag, solo_file.read_bytes(), multi_file.read_bytes(), failures)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--mapq", type=int, default=20)
    args = parser.parse_args()

    bin_dir = Path(args.bin_dir).resolve()
    work = Path(args.work).resolve()
    if work.exists():
        shutil.rmtree(work)
    solo_dir = work / "solo"
    multi_dir = work / "multi"
    solo_dir.mkdir(parents=True)
    multi_dir.mkdir(parents=True)

    bam = work / "input.bam"
    build_fixture(bam)
    q = str(args.mapq)
    failures: list = []

    # Standalone bam_stat (stdout-only command: streams are the output).
    rc, solo_bs_out, solo_bs_err = run(
        [str(bin_dir / "bam_stat"), "-i", str(bam), "-q", q], solo_dir
    )
    print(f"solo bam_stat exit={rc}")
    if rc != 0:
        failures.append("solo bam_stat exit code")

    # Standalone read_GC (files + streams).
    rc, solo_gc_out, solo_gc_err = run(
        [str(bin_dir / "read_GC"), "-i", str(bam), "-o", "out", "-q", q, "--skip-plot"],
        solo_dir,
    )
    print(f"solo read_GC exit={rc}")
    if rc != 0:
        failures.append("solo read_GC exit code")

    # Standalone read_NVC (files + streams).
    rc, solo_nvc_out, solo_nvc_err = run(
        [str(bin_dir / "read_NVC"), "-i", str(bam), "-o", "out", "-q", q, "--skip-plot"],
        solo_dir,
    )
    print(f"solo read_NVC exit={rc}")
    if rc != 0:
        failures.append("solo read_NVC exit code")

    # Standalone read_quality (R script only + streams; no data table).
    rc, solo_rq_out, solo_rq_err = run(
        [str(bin_dir / "read_quality"), "-i", str(bam), "-o", "out", "-q", q, "--skip-plot"],
        solo_dir,
    )
    print(f"solo read_quality exit={rc}")
    if rc != 0:
        failures.append("solo read_quality exit code")

    # Through the driver, same flags, same relative prefix.
    rc, multi_own_out, multi_own_err = run(
        [
            str(bin_dir / "rseqc_multi"),
            "-i", str(bam),
            "-o", "out",
            "--run", ",".join(COMMANDS),
            "-q", q,
            "--skip-plot",
        ],
        multi_dir,
    )
    print(f"multi exit={rc}")
    if rc != 0:
        failures.append("multi exit code")
    if multi_own_out != b"" or multi_own_err != b"":
        failures.append("multi own streams (expected silence on success)")
        print(f"multi own stdout={multi_own_out!r} stderr={multi_own_err!r}")

    # Streams, byte for byte: (label, solo bytes, multi stream file).
    for label, solo_bytes, stream_name in (
        ("bam_stat.stdout", solo_bs_out, "out.bam_stat.stdout"),
        ("bam_stat.stderr", solo_bs_err, "out.bam_stat.stderr"),
        ("read_GC.stdout", solo_gc_out, "out.read_GC.stdout"),
        ("read_GC.stderr", solo_gc_err, "out.read_GC.stderr"),
        ("read_NVC.stdout", solo_nvc_out, "out.read_NVC.stdout"),
        ("read_NVC.stderr", solo_nvc_err, "out.read_NVC.stderr"),
        ("read_quality.stdout", solo_rq_out, "out.read_quality.stdout"),
        ("read_quality.stderr", solo_rq_err, "out.read_quality.stderr"),
    ):
        stream_file = multi_dir / stream_name
        if not stream_file.exists():
            failures.append(f"multi {stream_name} missing")
            print(f"MISSING on multi side: {stream_name}")
            continue
        check_equal(label, solo_bytes, stream_file.read_bytes(), failures)

    # Data files, byte for byte, missing on either side is a failure.
    check_files(
        "files",
        solo_dir,
        multi_dir,
        ["out.GC.xls", "out.GC_plot.r", "out.NVC.xls", "out.NVC_plot.r", "out.qual.r"],
        failures,
    )

    # --nx forwarding: solo -x vs multi --nx, same relative prefix.
    nx_solo_dir = work / "nx_solo"
    nx_multi_dir = work / "nx_multi"
    nx_solo_dir.mkdir()
    nx_multi_dir.mkdir()
    rc, nx_solo_out, nx_solo_err = run(
        [str(bin_dir / "read_NVC"), "-i", str(bam), "-o", "outnx", "-q", q, "-x", "--skip-plot"],
        nx_solo_dir,
    )
    if rc != 0:
        failures.append("nx solo read_NVC -x exit code")
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outnx",
         "--run", "read_NVC", "-q", q, "--nx", "--skip-plot"],
        nx_multi_dir,
    )
    if rc != 0:
        failures.append("nx multi --nx exit code")
    check_equal(
        "nx read_NVC.stdout",
        nx_solo_out,
        (nx_multi_dir / "outnx.read_NVC.stdout").read_bytes()
        if (nx_multi_dir / "outnx.read_NVC.stdout").exists()
        else b"<MISSING>",
        failures,
    )
    if not (nx_multi_dir / "outnx.read_NVC.stdout").exists():
        failures.append("nx multi outnx.read_NVC.stdout missing")
    check_equal(
        "nx read_NVC.stderr",
        nx_solo_err,
        (nx_multi_dir / "outnx.read_NVC.stderr").read_bytes()
        if (nx_multi_dir / "outnx.read_NVC.stderr").exists()
        else b"<MISSING>",
        failures,
    )
    if not (nx_multi_dir / "outnx.read_NVC.stderr").exists():
        failures.append("nx multi outnx.read_NVC.stderr missing")
    check_files(
        "nx files",
        nx_solo_dir,
        nx_multi_dir,
        ["outnx.NVC.xls", "outnx.NVC_plot.r"],
        failures,
    )

    # --reduce forwarding: solo -r vs multi --reduce (long-only on the
    # driver: -r is already the reference BED there), same relative prefix.
    red_solo_dir = work / "reduce_solo"
    red_multi_dir = work / "reduce_multi"
    red_solo_dir.mkdir()
    red_multi_dir.mkdir()
    rc, red_solo_out, red_solo_err = run(
        [str(bin_dir / "read_quality"), "-i", str(bam), "-o", "outred", "-q", q, "-r", "1000", "--skip-plot"],
        red_solo_dir,
    )
    if rc != 0:
        failures.append("reduce solo read_quality -r exit code")
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outred",
         "--run", "read_quality", "-q", q, "--reduce", "1000", "--skip-plot"],
        red_multi_dir,
    )
    if rc != 0:
        failures.append("reduce multi --reduce exit code")
    check_equal(
        "reduce read_quality.stdout",
        red_solo_out,
        (red_multi_dir / "outred.read_quality.stdout").read_bytes()
        if (red_multi_dir / "outred.read_quality.stdout").exists()
        else b"<MISSING>",
        failures,
    )
    if not (red_multi_dir / "outred.read_quality.stdout").exists():
        failures.append("reduce multi outred.read_quality.stdout missing")
    check_equal(
        "reduce read_quality.stderr",
        red_solo_err,
        (red_multi_dir / "outred.read_quality.stderr").read_bytes()
        if (red_multi_dir / "outred.read_quality.stderr").exists()
        else b"<MISSING>",
        failures,
    )
    if not (red_multi_dir / "outred.read_quality.stderr").exists():
        failures.append("reduce multi outred.read_quality.stderr missing")
    check_files(
        "reduce files",
        red_solo_dir,
        red_multi_dir,
        ["outred.qual.r"],
        failures,
    )

    # Subset selection: --run bam_stat alone produces only bam_stat streams.
    sub_dir = work / "subset"
    sub_dir.mkdir()
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "out", "--run", "bam_stat"],
        sub_dir,
    )
    if rc != 0:
        failures.append("subset --run bam_stat exit code")
        print("subset run failed")
    else:
        for name in ("out.bam_stat.stdout", "out.bam_stat.stderr"):
            if not (sub_dir / name).exists():
                failures.append(f"subset {name} missing")
                print(f"MISSING in subset run: {name}")
        for name in ("out.GC.xls", "out.GC_plot.r", "out.read_GC.stdout"):
            if (sub_dir / name).exists():
                failures.append(f"subset unexpected {name}")
                print(f"UNEXPECTED in subset run: {name}")
        print("subset --run bam_stat: only bam_stat streams present")

    if failures:
        print(f"\nFAILED: {len(failures)} comparison(s): {failures}")
        return 1
    print("\nAll multi-driver comparisons identical.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
