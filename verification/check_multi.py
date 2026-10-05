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

COMMANDS = ("bam_stat", "read_GC", "read_NVC", "read_quality", "clipping_profile", "insertion_profile")


def build_bed12(path: Path) -> None:
    """Minimal BED12 covering the fixture's reads, for the gene-model commands.

    The fixture's reads all start at or after position 100 on chr1, so two
    wide transcripts are enough to give these commands something real to
    match against. Not a real annotation: check_multi compares the port
    against ITSELF through two invocation paths, so the model's contents are
    irrelevant as long as both sides see the same file.
    """
    path.write_text(
        "chr1\t100\t40000\tgene1\t0\t+\t100\t40000\t0,1,2,3\t0,10000,20000,30000\t0,0,0,0\n"
        "chr1\t100\t40000\tgene2\t0\t-\t100\t40000\t0,1,2,3\t0,10000,20000,30000\t0,0,0,0\n"
    )


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
    bed = work / "model.bed12"
    build_bed12(bed)
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

    # Standalone clipping_profile SE (files + streams; driver defaults to SE).
    rc, solo_cp_out, solo_cp_err = run(
        [str(bin_dir / "clipping_profile"), "-i", str(bam), "-o", "out", "-s", "SE", "-q", q, "--skip-plot"],
        solo_dir,
    )
    print(f"solo clipping_profile exit={rc}")
    if rc != 0:
        failures.append("solo clipping_profile exit code")

    # `infer_experiment` is checked separately below: it needs `--reference-bed`,
    # which the main matrix run does not pass, so it cannot be part of COMMANDS.

    # Standalone insertion_profile SE (files + streams; driver defaults to SE).
    rc, solo_ip_out, solo_ip_err = run(
        [str(bin_dir / "insertion_profile"), "-i", str(bam), "-o", "out", "-s", "SE", "-q", q, "--skip-plot"],
        solo_dir,
    )
    print(f"solo insertion_profile exit={rc}")
    if rc != 0:
        failures.append("solo insertion_profile exit code")

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
        ("clipping_profile.stdout", solo_cp_out, "out.clipping_profile.stdout"),
        ("clipping_profile.stderr", solo_cp_err, "out.clipping_profile.stderr"),
        ("insertion_profile.stdout", solo_ip_out, "out.insertion_profile.stdout"),
        ("insertion_profile.stderr", solo_ip_err, "out.insertion_profile.stderr"),
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
        ["out.GC.xls", "out.GC_plot.r", "out.NVC.xls", "out.NVC_plot.r", "out.qual.r",
         "out.clipping_profile.xls", "out.clipping_profile.r",
         "out.insertion_profile.xls", "out.insertion_profile.r"],
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

        # PE sequencing forwarding: solo -s PE vs multi -s PE, same prefix.
    pe_solo_dir = work / "pe_solo"
    pe_multi_dir = work / "pe_multi"
    pe_solo_dir.mkdir()
    pe_multi_dir.mkdir()
    rc, pe_solo_out, pe_solo_err = run(
        [str(bin_dir / "clipping_profile"), "-i", str(bam), "-o", "outpe", "-s", "PE", "-q", q, "--skip-plot"],
        pe_solo_dir,
    )
    if rc != 0:
        failures.append("pe solo clipping_profile -s PE exit code")
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outpe",
         "--run", "clipping_profile", "-q", q, "-s", "PE", "--skip-plot"],
        pe_multi_dir,
    )
    if rc != 0:
        failures.append("pe multi -s PE exit code")
    check_equal(
        "pe clipping_profile.stdout",
        pe_solo_out,
        (pe_multi_dir / "outpe.clipping_profile.stdout").read_bytes()
        if (pe_multi_dir / "outpe.clipping_profile.stdout").exists()
        else b"<MISSING>",
        failures,
    )
    if not (pe_multi_dir / "outpe.clipping_profile.stdout").exists():
        failures.append("pe multi outpe.clipping_profile.stdout missing")
    check_equal(
        "pe clipping_profile.stderr",
        pe_solo_err,
        (pe_multi_dir / "outpe.clipping_profile.stderr").read_bytes()
        if (pe_multi_dir / "outpe.clipping_profile.stderr").exists()
        else b"<MISSING>",
        failures,
    )
    if not (pe_multi_dir / "outpe.clipping_profile.stderr").exists():
        failures.append("pe multi outpe.clipping_profile.stderr missing")
    check_files(
        "pe files",
        pe_solo_dir,
        pe_multi_dir,
        ["outpe.clipping_profile.xls", "outpe.clipping_profile.r"],
        failures,
    )

    # infer_experiment: needs -r (gene model) and the shared header, and
    # writes no output file (its report IS stdout). Checked apart from the
    # main matrix because it requires --reference-bed.
    ie_solo_dir = work / "ie_solo"
    ie_multi_dir = work / "ie_multi"
    ie_solo_dir.mkdir()
    ie_multi_dir.mkdir()
    rc, ie_solo_out, ie_solo_err = run(
        [str(bin_dir / "infer_experiment"), "-i", str(bam), "-r", str(bed), "-q", q],
        ie_solo_dir,
    )
    if rc != 0:
        failures.append("ie solo infer_experiment exit code")
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outie", "-r", str(bed),
         "--run", "infer_experiment", "-q", q],
        ie_multi_dir,
    )
    if rc != 0:
        failures.append("ie multi --reference-bed exit code")
    for label, solo_bytes, stream in (
        ("ie stdout", ie_solo_out, "outie.infer_experiment.stdout"),
        ("ie stderr", ie_solo_err, "outie.infer_experiment.stderr"),
    ):
        path = ie_multi_dir / stream
        if not path.exists():
            failures.append(f"multi {stream} missing")
            print(f"MISSING on multi side: {stream}")
            continue
        check_equal(label, solo_bytes, path.read_bytes(), failures)

    # --sample-size forwarding: long-only on the driver (-s is sequencing).
    ss_solo_dir = work / "ss_solo"
    ss_multi_dir = work / "ss_multi"
    ss_solo_dir.mkdir()
    ss_multi_dir.mkdir()
    rc, ss_solo_out, ss_solo_err = run(
        [str(bin_dir / "infer_experiment"), "-i", str(bam), "-r", str(bed), "-s", "5", "-q", q],
        ss_solo_dir,
    )
    if rc != 0:
        failures.append("ss solo infer_experiment -s 5 exit code")
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outss", "-r", str(bed),
         "--run", "infer_experiment", "-q", q, "--sample-size", "5"],
        ss_multi_dir,
    )
    if rc != 0:
        failures.append("ss multi --sample-size 5 exit code")
    for label, solo_bytes, stream in (
        ("ss stdout", ss_solo_out, "outss.infer_experiment.stdout"),
        ("ss stderr", ss_solo_err, "outss.infer_experiment.stderr"),
    ):
        path = ss_multi_dir / stream
        if not path.exists():
            failures.append(f"multi {stream} missing")
            print(f"MISSING on multi side: {stream}")
            continue
        check_equal(label, solo_bytes, path.read_bytes(), failures)
    # The sample-size warning must survive the driver, not be dropped.
    if b"below 1,000" not in (ss_multi_dir / "outss.infer_experiment.stderr").read_bytes():
        failures.append("multi dropped the sub-1000 sample-size warning")
        print("UNEXPECTED: sample-size warning missing from multi stderr")

    # Refusal check: selecting infer_experiment without --reference-bed
    # must exit 2 (usage error), not run against an empty gene model.
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outie2",
         "--run", "infer_experiment", "-q", q],
        ie_multi_dir,
    )
    if rc != 2:
        failures.append(f"multi infer_experiment without -r exit code (expected 2, got {rc})")
        print(f"UNEXPECTED: missing --reference-bed gave exit {rc}, not 2")

    # deletion_profile: its own -l/-n, plus the required-flag refusal.
    del_solo_dir = work / "del_solo"
    del_multi_dir = work / "del_multi"
    del_solo_dir.mkdir()
    del_multi_dir.mkdir()
    rc, del_solo_out, del_solo_err = run(
        [str(bin_dir / "deletion_profile"), "-i", str(bam), "-o", "outdel",
         "-l", "101", "-n", "100000", "-q", q, "--skip-plot"],
        del_solo_dir,
    )
    if rc != 0:
        failures.append("del solo deletion_profile exit code")
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outdel",
         "--run", "deletion_profile", "--read-align-length", "101",
         "--read-num", "100000", "-q", q, "--skip-plot"],
        del_multi_dir,
    )
    if rc != 0:
        failures.append("del multi --read-align-length exit code")
    check_equal(
        "del deletion_profile.stdout",
        del_solo_out,
        (del_multi_dir / "outdel.deletion_profile.stdout").read_bytes()
        if (del_multi_dir / "outdel.deletion_profile.stdout").exists()
        else b"<MISSING>",
        failures,
    )
    if not (del_multi_dir / "outdel.deletion_profile.stdout").exists():
        failures.append("del multi outdel.deletion_profile.stdout missing")
    check_equal(
        "del deletion_profile.stderr",
        del_solo_err,
        (del_multi_dir / "outdel.deletion_profile.stderr").read_bytes()
        if (del_multi_dir / "outdel.deletion_profile.stderr").exists()
        else b"<MISSING>",
        failures,
    )
    if not (del_multi_dir / "outdel.deletion_profile.stderr").exists():
        failures.append("del multi outdel.deletion_profile.stderr missing")
    check_files(
        "del files",
        del_solo_dir,
        del_multi_dir,
        ["outdel.deletion_profile.txt", "outdel.deletion_profile.r"],
        failures,
    )

    # Refusal check: selecting deletion_profile without --read-align-length
    # must exit 2 (usage error), not run with an invented length.
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outdel2",
         "--run", "deletion_profile", "-q", q, "--skip-plot"],
        del_multi_dir,
    )
    if rc != 2:
        failures.append(f"multi deletion_profile without -l exit code (expected 2, got {rc})")
        print(f"UNEXPECTED: missing --read-align-length gave exit {rc}, not 2")

    # mismatch_profile: its own -l/-n, plus the shared required-flag refusal.
    mm_solo_dir = work / "mm_solo"
    mm_multi_dir = work / "mm_multi"
    mm_solo_dir.mkdir()
    mm_multi_dir.mkdir()
    rc, mm_solo_out, mm_solo_err = run(
        [str(bin_dir / "mismatch_profile"), "-i", str(bam), "-o", "outmm",
         "-l", "101", "-n", "100000", "-q", q, "--skip-plot"],
        mm_solo_dir,
    )
    if rc != 0:
        failures.append("mm solo mismatch_profile exit code")
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outmm",
         "--run", "mismatch_profile", "--read-align-length", "101",
         "--read-num", "100000", "-q", q, "--skip-plot"],
        mm_multi_dir,
    )
    if rc != 0:
        failures.append("mm multi --read-align-length exit code")
    check_equal(
        "mm mismatch_profile.stdout",
        mm_solo_out,
        (mm_multi_dir / "outmm.mismatch_profile.stdout").read_bytes()
        if (mm_multi_dir / "outmm.mismatch_profile.stdout").exists()
        else b"<MISSING>",
        failures,
    )
    if not (mm_multi_dir / "outmm.mismatch_profile.stdout").exists():
        failures.append("mm multi outmm.mismatch_profile.stdout missing")
    check_equal(
        "mm mismatch_profile.stderr",
        mm_solo_err,
        (mm_multi_dir / "outmm.mismatch_profile.stderr").read_bytes()
        if (mm_multi_dir / "outmm.mismatch_profile.stderr").exists()
        else b"<MISSING>",
        failures,
    )
    if not (mm_multi_dir / "outmm.mismatch_profile.stderr").exists():
        failures.append("mm multi outmm.mismatch_profile.stderr missing")
    check_files(
        "mm files",
        mm_solo_dir,
        mm_multi_dir,
        ["outmm.mismatch_profile.xls", "outmm.mismatch_profile.r"],
        failures,
    )

    # Refusal check: selecting mismatch_profile without --read-align-length
    # must exit 2 (usage error), same as deletion_profile.
    rc, _, _ = run(
        [str(bin_dir / "rseqc_multi"), "-i", str(bam), "-o", "outmm2",
         "--run", "mismatch_profile", "-q", q, "--skip-plot"],
        mm_multi_dir,
    )
    if rc != 2:
        failures.append(f"multi mismatch_profile without -l exit code (expected 2, got {rc})")
        print(f"UNEXPECTED: missing --read-align-length gave exit {rc}, not 2")

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
