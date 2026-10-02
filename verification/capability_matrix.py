#!/usr/bin/env python3
"""Per-command capability matrix, derived by executing the binaries.

The audit's Stage A asks for a "capability manifest" that records, per command, the
inputs it accepts, which outputs depend on a helper program, whether it is
stochastic, and which upstream quirks it deliberately reproduces. `compatibility/`
records those things in prose across three files, which is fine for a reader and
unusable as a machine-readable claim: nothing verifies that a command which says it
invokes `Rscript` actually does, or that a command which says it takes a `.cram`
actually parses one.

So this derives the matrix by running each command, not by reading documentation:

* **input formats** -- for each declared format, a real file of that type is fed to
  the command and the result is recorded. A command that exits 0 on a `.cram` is
  recorded as accepting CRAM; a command that exits non-zero is recorded as not
  accepting it, whatever its `--help` text says.
* **helper dependencies** -- the command is re-run with a `PATH` that contains no
  `Rscript` at all. The difference between that run and a normal one is what
  establishes whether the plot is native or delegated, and it also catches the
  failure mode where a command silently skips plotting and still exits 0.
* **stochastic behaviour** -- run twice with identical inputs and compared
  byte-for-byte. Inherently random commands (upstream seeds `random.Random`) will
  differ; that is recorded as a property to disclose, not as a defect to fix.
* **skip-plot contract** -- whether `--skip-plot` exists and, when it does, whether
  it actually suppresses the R invocation rather than being accepted and ignored.

The output feeds the archive's build/capability manifest. Claims that are asserted
in documentation but contradicted here are the point of running it.

Usage:
    cargo build --workspace --release --locked
    oracle/venv/bin/python3 verification/capability_matrix.py
    oracle/venv/bin/python3 verification/capability_matrix.py --json capabilities.json
    oracle/venv/bin/python3 verification/capability_matrix.py --emit-yaml
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RELEASE = REPO / "target" / "release"
ORACLE_SCRIPTS = REPO / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHON = REPO / "oracle" / "venv" / "bin" / "python3"

TINY = REPO / "datasets" / "aligned" / "tiny"
INTEROP = REPO / "verification" / "fixtures"

# Minimal PATH contents: everything a shell needs to run a command, and nothing that
# provides Rscript, wigToBigWig or htseq-count. Built by symlinking the resolved
# paths rather than by editing a PATH string, because PATH search order and the
# presence of /usr/bin in a default PATH are exactly what made the first version of
# this probe meaningless -- it "removed" Rscript and still found it.
NEEDED_TOOLS = [
    "bash", "sh", "cat", "ls", "rm", "mv", "cp", "mkdir", "dirname", "basename",
    "uname", "grep", "sed", "awk", "sort", "head", "tail", "wc", "tr", "cut",
    "python3", "mktemp", "date", "env", "printf", "test", "true", "false", "id",
    "gzip", "gunzip", "chmod", "find", "touch", "expr", "readlink", "which",
]

# Commands and the arguments that make a well-formed minimal invocation. `plot` marks
# the upstream `--skip-plot`/`--rscript` family; those need an R invocation unless
# skipped, which is precisely what the helper probe is for.
# The probe order. `COMMAND_SPECS` below is the single source of truth for what each
# command needs; this list only fixes the order, and it is validated against the specs
# and against the release archive's allowlist at run time so a new command cannot be
# added to one and forgotten in the other.
PROBE_ORDER = [
    "bam_stat", "bam2fq", "divide_bam", "split_bam", "split_paired_bam",
    "read_GC", "read_NVC", "read_quality", "read_duplication",
    "clipping_profile", "insertion_profile", "deletion_profile", "mismatch_profile",
    "read_hexamer", "infer_experiment", "read_distribution", "inner_distance",
    "RNA_fragment_size", "junction_annotation", "junction_saturation", "tin",
    "FPKM_count", "FPKM_UQ", "RPKM_saturation", "bam2wig", "geneBody_coverage",
    "geneBody_coverage2", "normalize_bigwig", "overlay_bigwig", "sc_bamStat",
    "sc_editMatrix", "sc_seqQual", "sc_seqLogo",
]

# The input containers an alignment-consuming command is tested against. The matrix
# derives each command's real support from running it; this is the list to try.
ALIGNMENT_INPUTS = ["bam", "sam", "cram"]

# SAM and CRAM fixtures are derived from the tiny BAM at probe time rather than
# committed, because a CRAM of the whole tiny panel is a large binary file that changes
# whenever the panel does. Deriving them here also means the format probe tests the
# same records in all three containers, so a difference is a difference in the format
# handling and not in the data.
def derive_alignment_fixtures(bam: Path, work: Path) -> dict:
    """Write tiny.sam and tiny.cram beside the probe's inputs; return what exists."""
    import pysam
    out = {"bam": bam}
    sam = work / (bam.stem + ".sam")
    cram = work / (bam.stem + ".cram")
    if not sam.exists():
        # `template=` needs an AlignmentFile, not a header object; passing
        # `src.header` raises a TypeError inside pysam. Writing the header explicitly
        # works with both, and the records are copied verbatim either way.
        with pysam.AlignmentFile(str(bam)) as src:
            header = src.header
            records = list(src.fetch(until_eof=True))
        with pysam.AlignmentFile(str(sam), "w", header=header) as dst:
            for r in records:
                dst.write(r)
    out["sam"] = sam
    # A CRAM needs a reference file to exist at write time -- htslib opens it to build
    # its reference cache, so a missing or unindexed FASTA fails with "failure when
    # setting reference filename" before a single record is written. A FASTA of N's at
    # the header's contig lengths is enough: what is under test is the port's CRAM
    # container handling, not the reference sequence, and the records are copied
    # verbatim from the BAM.
    if not cram.exists():
        ref = work / "probe_reference.fa"
        with pysam.AlignmentFile(str(bam)) as src:
            for entry in src.header.to_dict()["SQ"]:
                with open(ref, "a") as fh:
                    fh.write(f">{entry['SN']}\n")
                    fh.write("N" * entry["LN"] + "\n")
        pysam.faidx(str(ref))
        with pysam.AlignmentFile(str(bam)) as src:
            header = src.header
            records = list(src.fetch(until_eof=True))
        with pysam.AlignmentFile(str(cram), "wc", header=header,
                                 reference_filename=str(ref)) as dst:
            for r in records:
                dst.write(r)
    out["cram"] = cram
    return out


NON_ALIGNMENT_INPUTS = ["fastq", "bigwig", "fa"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def make_minimal_path(dirpath: Path) -> Path:
    """A directory of symlinks providing shell basics and no Rscript."""
    dirpath.mkdir(parents=True, exist_ok=True)
    for tool in NEEDED_TOOLS:
        found = shutil.which(tool)
        if found:
            try:
                (dirpath / tool).symlink_to(found)
            except FileExistsError:
                pass
    return dirpath


# Per-command required options and the kind of primary input each command reads.
#
# This table is what makes the input-format probe meaningful. The first version built
# every argv as `-i <bam> ...`, which is correct for maybe two thirds of the commands
# and silently wrong for the rest, so those failed for a reason that had nothing to do
# with the format under test: `bam2wig` needs `-s chrom.sizes`, `overlay_bigwig` and
# `normalize_bigwig` read BigWig rather than alignments at all, `geneBody_coverage2`
# reads BigWig, and `sc_bamStat` needs a single-cell BAM with CB/UB tags. Every one of
# them was recorded as "REJECTED bam", which is a statement about the harness.
#
# `input_kind` is what the command actually consumes:
#   alignment  -- BAM/SAM/CRAM; the format probe applies
#   bigwig     -- a coverage track; the format probe does not apply
#   fastq      -- a barcode/quality FASTQ
# `extra` are the options the command genuinely requires, beyond input and output.
# `stdout_only` commands have no output-file option at all; adding one makes the probe
# fail on an unknown flag and report the command as not accepting its own input format.
COMMAND_SPECS = {
    "bam_stat": dict(input_kind="alignment", extra=[], stdout_only=True),
    "infer_experiment": dict(input_kind="alignment", extra=["-r", "@BED@"],
                             stdout_only=True),
    "bam2fq": dict(input_kind="alignment", extra=[]),
    # `-n` is required (the help text says so); without it the command exits 2 on a
    # missing argument and the probe recorded that as "did not succeed".
    "divide_bam": dict(input_kind="alignment", extra=["-n", "2"]),
    "split_bam": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    "split_paired_bam": dict(input_kind="alignment", extra=[]),
    "read_GC": dict(input_kind="alignment", extra=[]),
    "read_NVC": dict(input_kind="alignment", extra=[]),
    "read_quality": dict(input_kind="alignment", extra=[]),
    "read_duplication": dict(input_kind="alignment", extra=[]),
    # `--sequencing` is a single value, SE or PE -- the layout, not a chemistry and a
    # read length. The probe passed "pa" then "100" as two values, which clap rejected
    # against the enum, so both commands reported "did not succeed" for a flag they do
    # accept.
    "clipping_profile": dict(input_kind="alignment",
                             extra=["--sequencing", "PE"]),
    "insertion_profile": dict(input_kind="alignment",
                              extra=["--sequencing", "PE"]),
    "deletion_profile": dict(input_kind="alignment",
                             extra=["--read-align-length", "100"]),
    "mismatch_profile": dict(input_kind="alignment",
                             extra=["--read-align-length", "100"]),
    # read_hexamer takes comma-separated FASTA/FASTQ READ files, not an alignment.
    # `@READS@` is filled in with the single-cell FASTQ fixture.
    "read_hexamer": dict(input_kind="fastq", extra=["-i", "@READS@"],
                         input_flag=None),
    "read_distribution": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    "inner_distance": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    # `-r` is required here, which the help text states but the previous entry did
    # not reflect -- the probe then failed on a missing required argument and recorded
    # the command as not accepting any format.
    "RNA_fragment_size": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    "junction_annotation": dict(input_kind="alignment",
                                extra=["-r", "@BED@", "--skip-bed", "--skip-interact"]),
    "junction_saturation": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    "tin": dict(input_kind="alignment", extra=["-r", "@BED@", "-n", "50"]),
    "FPKM_count": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    # FPKM-UQ's flags are entirely different: --bam/--gtf/--info/--output, no -i and
    # no -r. A probe that assumed the common shape could not invoke it at all.
    "FPKM_UQ": dict(input_kind="alignment",
                    extra=["--bam", "@BAM@", "--gtf", "@GTF@", "--info", "@INFO@"],
                    input_flag=None),
    # Upstream re-parses the gene model on every percentile iteration; the port parses
    # it once (DIV-0012). The R-script data columns still come from unseeded
    # resampling, so this command is stochastic by design like junction_saturation.
    "RPKM_saturation": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    "bam2wig": dict(input_kind="alignment",
                    extra=["-s", "@CHROMSIZES@", "-o", "@PREFIX@w"]),
    "geneBody_coverage": dict(input_kind="alignment", extra=["-r", "@BED@"]),
    # These three read a coverage track, not an alignment. Probing them with a BAM
    # measures nothing about their input support, so they are excluded from the
    # format probe and recorded as taking BigWig.
    # @BED@ resolves to a BED12 whose chromosomes match the BigWig actually being
    # read, derived in _fixture_context. Pairing pyBigWig's chr1-based test track with
    # the tiny panel's chr17 model fails with "no valid BED12 records matched
    # chromosomes in the BigWig", which is a statement about the fixture pairing.
    "geneBody_coverage2": dict(input_kind="bigwig", extra=["-r", "@TRACKBED@"]),
    "normalize_bigwig": dict(input_kind="bigwig", extra=[]),
    "overlay_bigwig": dict(input_kind="bigwig", extra=[]),
    # Single-cell: needs a BAM carrying CB/UB/RE tags, which the tiny fixture is not.
    # The single-cell commands need fixtures carrying CB/UB/RE tags and a barcode
    # FASTQ; the tiny panel has none of those, so they get the dedicated fixtures
    # `verification/fixtures/make_sc_bamstat_fixture.py` and
    # `make_sc_editmatrix_fixture.py` already build for the differential suite.
    "sc_bamStat": dict(input_kind="single_cell_alignment", extra=[]),
    "sc_editMatrix": dict(input_kind="single_cell_alignment", extra=[]),
    "sc_seqQual": dict(input_kind="fastq", extra=[]),
    "sc_seqLogo": dict(input_kind="fastq", extra=[]),
}


# Commands whose output argument is a DIRECTORY rather than a file name or prefix.
# Upstream `tin.py` requires its `-o` directory to exist and errors out if it does not,
# which is a documented contract of that command and must be honoured by anything
# invoking it. Every other command writes files, so pre-creating a directory named
# `out` made them fail with "Is a directory (os error 21)" -- a statement about the
# probe, recorded as if it were a statement about the command.
OUTDIR_COMMANDS = {"tin"}


def _make_output_dir(name: str, run_dir: Path) -> None:
    """Create whatever the command's output argument must already exist.

    Two cases, and both were found by the probe failing for a reason unrelated to what
    it was testing:

    * `tin`'s `-o` IS a directory and must exist (upstream errors out if it does not),
      so it is created with `exist_ok=True` -- the determinism probe runs each command
      twice in the same parent, so a plain `mkdir` would raise FileExistsError inside
      the harness.
    * Every other command's output argument is a file PREFIX whose PARENT must exist.
      A handful (RPKM_saturation among them) surface a missing parent as a bare
      `No such file or directory (os error 2)` after doing all the work, so the parent
      is created for them too.
    """
    if name in OUTDIR_COMMANDS:
        (run_dir / "out").mkdir(parents=True, exist_ok=True)
    else:
        (run_dir / "out").parent.mkdir(parents=True, exist_ok=True)


def resolve_extra(value: str, ctx: dict) -> str:
    """Substitute the placeholders in a spec's `extra` list."""
    return {
        "@BED@": ctx["bed"],
        "@TRACKBED@": ctx.get("track_bed") or ctx["bed"],
        "@CHROMSIZES@": ctx["chrom_sizes"],
        "@BAM@": ctx["bam"],
        "@GTF@": ctx["gtf"],
        "@INFO@": ctx["info"],
        "@READS@": ctx.get("sc_reads", ""),
        "@PREFIX@": "out",
    }.get(value, value)


# Which context key holds each command's PRIMARY input. A command that reads BigWig or
# a barcode FASTQ was previously handed the tiny BAM, so it failed for a reason
# unrelated to the property being probed and its row said nothing either way.
PRIMARY_INPUT_KEY = {
    "bigwig": "bigwig",
    "fastq": "sc_reads",
    "single_cell_alignment": "sc_bam",
    "single_cell": "sc_edit_bam",
}


def primary_input(name: str, bam: Path, ctx: dict) -> Path:
    """The file this command actually reads, given the command's declared input kind."""
    kind = COMMAND_SPECS[name]["input_kind"]
    key = PRIMARY_INPUT_KEY.get(kind, "bam")
    value = ctx.get(key)
    if not value:
        return bam
    return Path(value)


def build_argv(name: str, exe: Path, bam: Path, bed: Path, work: Path,
               ext: str, primary: Path | None = None, ctx: dict | None = None) -> list[str]:
    """A minimal, well-formed invocation for `name`, per COMMAND_SPECS."""
    spec = COMMAND_SPECS.get(name)
    if spec is None:
        raise KeyError(f"{name} has no COMMAND_SPECS entry; the capability matrix "
                       f"would be measuring a guessed invocation")
    ctx = ctx or {"bed": str(bed), "chrom_sizes": str(work / "chrom.sizes"),
                  "bam": str(bam), "gtf": "", "info": "", "reads": ""}
    if primary is not None:
        src = Path(primary)
    elif spec["input_kind"] in PRIMARY_INPUT_KEY:
        src = primary_input(name, bam, ctx)
    else:
        src = Path(ctx["bam"]).with_suffix("." + ext)
    help_text = subprocess.run([str(exe), "--help"], capture_output=True,
                               text=True, timeout=60).stdout
    # Detect the output flag from the usage line, not from a bare count of "-o":
    # "-o" appears inside "--out-prefix" and inside help prose, so a count matched
    # commands that have no such option and made their probes fail on an unknown flag.
    # `-o` is advertised in two spellings across the 33 commands: `-o, --output <X>`
    # (clap's short-then-long form) and a bare `-o <X>`. Both are matched. A count of
    # "-o" occurrences was the earlier approach and it matched `--out-prefix` too, so
    # commands with no such option were invoked with one and failed on it.
    out_flag = "--out-prefix" if "--out-prefix" in help_text else (
        "-o" if re.search(r"^\s*-o(?:,|\s)", help_text, re.M) else None)
    # The output argument is a RELATIVE name, because every probe runs with the
    # command's cwd set to its own scratch directory. Deriving it from `work` instead
    # -- the obvious first version -- produced "out/out" whenever the caller passed a
    # relative `work`, because `work / "out"` is relative too. That string is a valid
    # path whose parent does not exist, so RPKM_saturation did all its work and then
    # failed with a bare ENOENT, which the determinism probe recorded as "command did
    # not succeed".
    out_value = "out"

    if name == "overlay_bigwig":
        # Two BigWig inputs plus an action; `@PRIMARY@` stands in for the first.
        # "Average", not "mean": upstream's own vocabulary for this flag, which the
        # binary lists in its error message. A plausible synonym fails the same way an
        # unsupported format does, which is a misleading reason to record.
        return [str(exe), "-i", str(src), "-j", str(src), "-a", "Average",
                "-o", out_value, "--overwrite"]
    if name == "normalize_bigwig":
        return [str(exe), "-i", str(src), "-o", out_value, "--overwrite"]
    if name in ("bam2wig",):
        # Absolute, because every probe runs with its cwd set to the run's scratch
        # directory. A relative path here would write chrom.sizes into whichever
        # directory the probe happened to be launched from -- which, during a manual
        # run of this script, is the repository root.
        sizes = Path(ctx["chrom_sizes"]).resolve()
        if not sizes.exists():
            import pysam
            with pysam.AlignmentFile(str(bam)) as f:
                sizes.parent.mkdir(parents=True, exist_ok=True)
                sizes.write_text("\n".join(
                    f"{e['SN']}\t{e['LN']}"
                    for e in f.header.to_dict()["SQ"]) + "\n")
        return [str(exe), "-i", str(src), "-s", str(sizes), "-o", out_value + "w"]

    # `input_flag: None` means the command names its input some other way entirely
    # (FPKM-UQ's `--bam`), and everything needed is already in `extra`.
    argv = [str(exe)] if spec.get("input_flag", "-i") is None else [
        str(exe), spec.get("input_flag", "-i"), str(src)]
    for x in spec["extra"]:
        argv.append(resolve_extra(x, ctx))
    if out_flag and not spec.get("stdout_only"):
        argv += [out_flag, out_value]
    return argv


def run(argv: list[str], cwd: Path, env_path: str | None = None,
        timeout: int = 300) -> dict:
    env = dict(os.environ)
    if env_path:
        env["PATH"] = env_path
    else:
        env["PATH"] = "/usr/local/bin:/usr/bin:/bin"
    try:
        p = subprocess.run([str(a) for a in argv], cwd=str(cwd), env=env,
                           capture_output=True, text=True, timeout=timeout)
        return {"exit": p.returncode, "stderr": p.stderr[-300:],
                "stdout": p.stdout[-300:]}
    except subprocess.TimeoutExpired:
        return {"exit": None, "stderr": "TIMEOUT", "stdout": ""}


def probe_formats(name: str, exe: Path, bam: Path, bed: Path, work: Path,
                  fixtures: dict, ctx: dict) -> dict:
    """Feed a real file of each declared format and record what happens.

    Only meaningful for commands whose primary input IS an alignment. A command that
    reads BigWig or a single-cell BAM is recorded with `applies: false` and the kind of
    input it actually takes, rather than being handed a BAM and reported as not
    supporting BAM.
    """
    kind = COMMAND_SPECS[name]["input_kind"]
    if kind != "alignment":
        return {
            "applies": False,
            "primary_input_kind": kind,
            "note": "not an alignment-consuming command; its input support is not "
                    "covered by feeding it a BAM",
        }
    # Per-command format restrictions that are REAL and are inherited from upstream.
    # `tin` resolves its inputs through upstream's own `getBamFiles.get_bam_files`,
    # which globs `*.bam`, so a `.sam` or `.cram` path is accepted on the command line
    # and then reported as "No BAM files found". Reproducing that is correct for a port:
    # the alternative is a command that works on inputs upstream cannot. It is recorded
    # here with its cause so the restriction is a disclosed property rather than an
    # unexplained rejection, and so a future change that alters it is visible.
    # Each entry is a real, inherited restriction with the mechanism that causes it, so
    # the matrix reports a disclosed property rather than an unexplained failure. The
    # eight BAM-only commands are not a port limitation: upstream's own `validate_args`
    # rejects a non-`.bam` extension with `parser.error` before reading anything, and
    # the port reproduces that. They are recorded individually because each names a
    # different upstream script.
    BAM_ONLY_UPSTREAM = {
        "tin", "divide_bam", "split_bam", "split_paired_bam", "deletion_profile",
        "mismatch_profile", "RNA_fragment_size", "FPKM_count", "bam2wig",
    }
    FORMAT_RESTRICTIONS = {
        name: ("BAM only, inherited from upstream: the command's own `validate_args` "
               "rejects a non-`.bam` extension with parser.error before reading the "
               "file. Reproduced deliberately; a port that accepted more than upstream "
               "here would disagree on an input the user was told was invalid."
               ) for name in sorted(BAM_ONLY_UPSTREAM)
    }
    results = {"applies": True, "primary_input_kind": kind}
    if name in FORMAT_RESTRICTIONS:
        results["known_restriction"] = FORMAT_RESTRICTIONS[name]
    for ext in ALIGNMENT_INPUTS:
        src = fixtures.get(ext)
        if src is None or not Path(src).exists():
            results[ext] = {"status": "NO_FIXTURE",
                            "detail": f"no {ext} fixture available"}
            continue
        if name in ("read_hexamer", "sc_seqLogo", "sc_editMatrix",
                    "sc_seqQual"):
            results[ext] = {
                "status": "NOT_APPLICABLE",
                "detail": "command reads FASTQ/single-cell input, not an alignment",
            }
            continue
        # The per-run directory must be keyed by COMMAND as well as format. Sharing
        # one `fmt-bam/` across commands meant an earlier command that writes a FILE
        # named `out` (RNA_fragment_size's `-o out`) left a file where `tin` needs a
        # directory, and `tin` then failed with FileExistsError inside the harness --
        # a probe defect that would have read as "tin does not accept BAM".
        d = work / f"fmt-{name}-{ext}"
        d.mkdir(parents=True, exist_ok=True)
        _make_output_dir(name, d)
        argv = build_argv(name, exe, bam, bed, d, ext, primary=Path(src), ctx=ctx)
        r = run(argv, d)
        if r["exit"] == 0:
            results[ext] = {"status": "ACCEPTED", "exit": 0, "detail": ""}
        elif _mentions_missing_helper(r["stderr"]):
            # Not a format verdict: the command never reached the point of reading the
            # input, so nothing about format support was tested.
            results[ext] = {
                "status": "NOT_MEASURED",
                "exit": r["exit"],
                "detail": "requires an external helper that is absent, so input "
                          "support was not tested",
            }
        else:
            results[ext] = {"status": "REJECTED", "exit": r["exit"],
                            "detail": r["stderr"].strip()[-120:]}
    return results


def probe_helper(name: str, exe: Path, bam: Path, bed: Path, work: Path,
                 minpath: str, ctx: dict) -> dict:
    """Does the command need a helper, and does it say so when the helper is absent?"""
    d = work / f"helper-{name}"
    d.mkdir(parents=True, exist_ok=True)
    _make_output_dir(name, d)
    argv = build_argv(name, exe, bam, bed, d, "bam", ctx=ctx)
    r = run(argv, d, env_path=minpath)
    pdfs = list(d.glob("*.pdf"))
    said = r["stderr"].strip()
    mentions_helper = any(k in said for k in
                          ("Rscript", "Rscript executable not found",
                           "wigToBigWig", "htseq-count", "not found"))
    if r["exit"] == 0 and not pdfs and not said:
        return {"needs_helper": "unknown", "requires_rscript": False,
                "detail": "exit 0, no output, no message"}
    if r["exit"] == 0 and not pdfs:
        return {"needs_helper": "no", "requires_rscript": False,
                "detail": "succeeds with no helper, no PDF",
                "plot_artifact": "absent"}
    if r["exit"] != 0 and mentions_helper:
        # WHICH helper is absent, not merely that one is. Deriving "requires Rscript"
        # from a single generic needs_helper flag labels FPKM-UQ as Rscript-dependent
        # when its helper is htseq-count, which is a measured fact turned into a false
        # one for exactly the command that does not need R. The detail names the
        # helper, so record it and let the two questions stay separate.
        named = None
        for candidate in ("htseq-count", "Rscript", "wigToBigWig", "wigToBigWigHelper"):
            if candidate in said:
                named = candidate
                break
        return {"needs_helper": "yes", "reports_missing_helper": True,
                "helper": named,
                "requires_rscript": named in ("Rscript",),
                "detail": said[-120:]}
    if r["exit"] != 0:
        return {"needs_helper": "unknown", "requires_rscript": False,
                "reports_missing_helper": False, "detail": said[-120:]}
    return {"needs_helper": "no", "requires_rscript": False,
            "plot_artifact": "present even without a helper"}


def probe_skip_plot(name: str, exe: Path, bam: Path, bed: Path, work: Path,
                    ctx: dict) -> dict:
    """--skip-plot must exist AND suppress the helper, not be accepted and ignored."""
    help_text = subprocess.run([str(exe), "--help"], capture_output=True,
                               text=True, timeout=60).stdout
    has_flag = "--skip-plot" in help_text
    entry = {"has_skip_plot": has_flag}
    if not has_flag:
        entry["suppresses_helper"] = None
        return entry
    d = work / f"skipplot-{name}"
    d.mkdir(parents=True, exist_ok=True)
    _make_output_dir(name, d)
    argv = build_argv(name, exe, bam, bed, d, "bam", ctx=ctx) + ["--skip-plot"]
    r = run(argv, d)
    entry["exit"] = r["exit"]
    entry["pdfs_written"] = len(list(d.glob("*.pdf")))
    entry["rscripts_written"] = len(list(d.glob("*.r")))
    # A --skip-plot that is accepted, produces no PDF, and mentions no helper means
    # the plot was genuinely skipped rather than attempted-and-failed.
    entry["suppresses_helper"] = (
        r["exit"] == 0 and not list(d.glob("*.pdf"))
        and "Rscript" not in r["stderr"]
    )
    return entry


# R writes /CreationDate and /ModDate into every PDF it renders, so a plot is never
# byte-identical between two runs even when its data is. Comparing PDFs byte-for-byte
# therefore reports every R-backed command as non-deterministic, which is noise rather
# than a finding. The data artifacts (the .xls/.txt tables a consumer actually reads)
# are compared exactly, and the PDFs are compared with their timestamps normalised, so
# a genuine difference in a rendered plot is still caught.
PDF_DATE = re.compile(rb"/(?:Creation|Mod)Date \(D:[0-9]+\)")

# Phrases a command uses when it cannot find a helper it needs. Matched so a
# helper-absent run is recorded as a capability fact rather than as a failure, and so
# no hand-maintained list of helper-dependent commands is needed.
HELPER_MISSING_PHRASES = (
    "executable not found",
    "cannot find htseq-count",
    "Install it with",            # an R package is missing, e.g. pheatmap
    "could not find wigToBigWig",
)


def _mentions_missing_helper(stderr: str) -> bool:
    return any(p in stderr for p in HELPER_MISSING_PHRASES)


def digest_outputs(d: Path) -> dict:
    out = {}
    for p in sorted(d.rglob("*")):
        if p.is_file():
            data = p.read_bytes()
            if p.suffix == ".pdf":
                data = PDF_DATE.sub(b"/Date(D:0)", data)
            out[str(p.relative_to(d))] = (p.stat().st_size, hashlib.sha256(data).hexdigest())
    return out


def probe_stochastic(name: str, exe: Path, bam: Path, bed: Path,
                     work: Path, ctx: dict) -> dict:
    """Two identical runs; do they agree?

    Data artifacts are compared byte-for-byte. PDFs are compared with R's embedded
    timestamps normalised, and the two verdicts are reported separately, because "the
    numbers are reproducible but the plot carries a render time" is a different claim
    from "the command is stochastic" and a release note needs the first.
    """
    # A command that needs an external helper cannot produce a determinism verdict
    # here, because the helper is absent. That is a real, reportable capability fact
    # -- FPKM-UQ needs htseq-count, sc_editMatrix needs R's pheatmap -- and it is
    # recorded as such rather than as "did not succeed", which reads like a defect in
    # the command. Detected from the diagnostic the command itself produced, so no
    # per-command list of helper-dependent commands has to be maintained here.
    pre = work / f"preflight-{name}"
    pre.mkdir(parents=True, exist_ok=True)
    _make_output_dir(name, pre)
    probe = run(build_argv(name, exe, bam, bed, Path("out"), "bam", ctx=ctx), pre)
    if probe["exit"] != 0 and _mentions_missing_helper(probe["stderr"]):
        last = probe["stderr"].strip().splitlines()
        return {
            "data_deterministic": None,
            "plot_deterministic": None,
            "stdout_deterministic": None,
            "detail": "not measured: requires an external helper that is absent "
                      f"({last[-1][:140] if last else 'no diagnostic'})",
        }

    data_digests, plot_digests, stdout_digests = [], [], []
    for i in (1, 2):
        d = work / f"stoch-{name}-{i}"
        d.mkdir(parents=True, exist_ok=True)
        # The SAME relative output path is used in both runs, differing only in the
        # parent directory. A generated .r script embeds its own output path
        # (`pdf("o.GC_plot.pdf")`), so giving each run a different --out-prefix made
        # the scripts differ by construction and every command with a plot contract
        # looked stochastic. That was a defect in the probe, not in the commands.
        # The `out` parent must exist: a few commands require their -o directory to
        # be present (upstream tin.py errors if it is not), and its absence would make
        # both runs fail identically and be reported as "did not succeed".
        _make_output_dir(name, d)
        argv = build_argv(name, exe, bam, bed, Path("out"), "bam", ctx=ctx)
        r = run(argv, d)
        if r["exit"] != 0:
            return {"data_deterministic": None, "plot_deterministic": None,
                    "stdout_deterministic": None,
                    "detail": "command did not succeed",
                    "stderr": r["stderr"][-160:]}
        digests = digest_outputs(d)
        data_digests.append({k: v for k, v in digests.items()
                             if not k.endswith(".pdf")})
        plot_digests.append({k: v for k, v in digests.items() if k.endswith(".pdf")})
        # A stdout-only command (bam_stat, infer_experiment) writes no data file, so
        # comparing artifacts alone would report it as "nothing to compare" and skip
        # the determinism question entirely -- which is the question that matters most
        # for a command whose entire result is its stdout. The report is part of the
        # result.
        stdout_digests.append((r["exit"], hashlib.sha256(
            r["stdout"].encode("utf-8", "replace")).hexdigest()))
    entry = {
        "stdout_deterministic": (stdout_digests[0] == stdout_digests[1]
                                 if stdout_digests[0][0] == 0 else None),
        "data_deterministic": (data_digests[0] == data_digests[1]
                               if data_digests[0] else None),
        "plot_deterministic_after_timestamp_normalisation":
            plot_digests[0] == plot_digests[1] if plot_digests[0] else None,
        "data_artifacts": len(data_digests[0]),
        "plot_artifacts": len(plot_digests[0]),
    }
    if not data_digests[0] and entry["stdout_deterministic"] is None:
        entry["detail"] = "no data artifacts and no successful stdout to compare"
    return entry


def _check_coverage(names: list[str]) -> None:
    """Fail loudly if the probe does not cover exactly the shipped command set.

    The previous probe iterated a hand-written list that had drifted from both the
    per-command spec table and the release archive's allowlist: `RPKM_saturation` was
    absent from it and `read_hexamer` appeared twice, so 32 of the 33 shipped commands
    were probed and nobody could tell from the output that one was missing. Cross-
    checking both sources is what makes "every advertised capability" checkable.
    """
    missing_specs = sorted(set(PROBE_ORDER) - set(COMMAND_SPECS))
    extra_specs = sorted(set(COMMAND_SPECS) - set(PROBE_ORDER))
    if missing_specs or extra_specs:
        raise SystemExit(
            f"PROBE_ORDER and COMMAND_SPECS disagree: "
            f"in the order but not the specs {missing_specs}; "
            f"in the specs but not the order {extra_specs}"
        )
    script = REPO / "scripts" / "build-release-archive.sh"
    if script.is_file():
        text = script.read_text()
        m = re.search(r'ALLOWLIST="(.*?)"', text, re.S)
        if m:
            allowlist = set(m.group(1).split())
            gap = sorted(allowlist - set(PROBE_ORDER))
            if gap:
                raise SystemExit(
                    f"the release archive ships {gap} but the capability matrix does "
                    f"not probe them; every shipped command needs a capability record"
                )


# Context keys that name a fixture, as opposed to a derived value. Used only for the
# "extra fixtures available" log line.
_CONTEXT_FIXTURE_KEYS = {
    "sc_bam", "sc_edit_bam", "sc_reads", "sc_logo_reads", "gtf", "info", "bigwig",
}


TRACK_MODEL = REPO / "verification" / "fixtures" / "track" / "track_model.bed12"


def _bed_for_track(bigwig: Path, work: Path) -> Path:
    """A minimal BED12 covering a BigWig's own intervals, with matching contig names."""
    import pyBigWig
    out = work / "track_derived.bed12"
    rows = []
    with pyBigWig.open(str(bigwig)) as bw:
        for chrom, length in bw.chroms().items():
            intervals = bw.intervals(chrom) or []
            if not intervals:
                # No intervals means no coverage to intersect against; a gene body
                # spanning the whole contig is the only thing that can match.
                rows.append((chrom, 0, int(length)))
                continue
            for start, _end, _value in intervals:
                rows.append((chrom, int(start), int(start) + 10))
    with open(out, "w") as fh:
        for i, (chrom, start, end) in enumerate(rows):
            span = end - start
            fh.write(f"{chrom}\t{start}\t{end}\ttx{i}\t0\t+\t{start}\t{end}\t0"
                     f"\t1\t{span},\t0,\n")
    return out


def _fixture_context(bam: Path, bed: Path, work: Path) -> dict:
    """Resolve every input path the probe may need, once.

    Each is optional: a command whose fixture cannot be produced is still probed, and
    the missing input is reported in its row rather than aborting the whole matrix.
    That distinction matters -- aborting would leave a partially-written JSON that
    reads like a complete capability record.
    """
    from verification_fixture_paths import (
        fpkm_uq_gtf, fpkm_uq_info, sc_bamstat, sc_editmatrix, sc_seqlogo,
        sc_seqqual, track_bigwig,
    )

    def attempt(fn, *args):
        try:
            return str(fn(*args))
        except Exception as exc:
            print(f"note: {fn.__name__} unavailable: {exc}")
            return ""

    # A BED12 whose chromosomes match the BigWig the probe reads.
    # `geneBody_coverage2` intersects its gene model against the track's chromosomes,
    # so pairing pyBigWig's chr1 test track with the tiny panel's chr17 model fails
    # with "no valid BED12 records matched chromosomes in the BigWig" -- a statement
    # about the fixture pairing, not about the command. The track fixtures under
    # `verification/fixtures/track/` already carry a matching model; when only the
    # pyBigWig copy is present its contig names are read and a minimal model covering
    # the track's own intervals is written.
    track_bed = TRACK_MODEL
    try:
        bigwig = track_bigwig()
    except Exception:
        bigwig = None
    if bigwig is not None and not (track_bed and track_bed.is_file()):
        track_bed = _bed_for_track(bigwig, work)

    return {
        "bed": str(bed),
        "track_bed": str(track_bed or ""),
        "chrom_sizes": str(work / "chrom.sizes"),
        "bam": str(bam),
        "gtf": attempt(fpkm_uq_gtf),
        "info": attempt(fpkm_uq_info),
        "reads": attempt(sc_seqqual),
        "sc_reads": attempt(sc_seqqual),
        "sc_logo_reads": attempt(sc_seqlogo),
        "sc_bam": attempt(sc_bamstat),
        "sc_edit_bam": attempt(sc_editmatrix),
        "bigwig": attempt(track_bigwig),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--commands", nargs="*", default=None)
    ap.add_argument("--skip-formats", action="store_true",
                    help="skip the per-format input probe (the slowest part)")
    args = ap.parse_args()

    if not RELEASE.exists():
        print("target/release not found; run: cargo build --workspace --release --locked")
        return 2
    bam = TINY / "tiny.bam"
    bed = TINY / "model.bed12"
    if not bam.exists():
        print(f"fixture panel missing: {bam}")
        return 2

    names = args.commands or list(PROBE_ORDER)
    _check_coverage(names)
    results = {
        "schema": 1,
        "generated_from": "executed binaries; every field is observed, not declared",
        "fixture": {"bam": str(bam.relative_to(REPO)),
                    "bam_sha256": sha256(bam),
                    "bed": str(bed.relative_to(REPO))},
        "commands": {},
    }

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        minpath = str(make_minimal_path(work / "minpath"))
        # A chrom.sizes file for the commands that need one.
        import pysam
        with pysam.AlignmentFile(str(bam)) as f:
            (work / "chrom.sizes").write_text(
                "\n".join(f"{e['SN']}\t{e['LN']}"
                          for e in f.header.to_dict()["SQ"]) + "\n")
        fixtures = (derive_alignment_fixtures(bam, work)
                    if not args.skip_formats else {"bam": bam})
        print("derived SAM and CRAM fixtures from the tiny panel: "
              f"{[k for k in fixtures if k in ('sam', 'cram')]}\n")

        # Placeholders for the commands whose input is not the tiny panel. The
        # single-cell fixtures are built by the generators the differential suite
        # already uses, so there is one definition of each rather than two.
        ctx = _fixture_context(bam, bed, work)
        print("extra fixtures available: "
              + ", ".join(sorted(k for k, v in ctx.items()
                                 if v and k in _CONTEXT_FIXTURE_KEYS)) + "\n")
        print(f"minimal PATH has no Rscript: "
              f"{shutil.which('Rscript', path=minpath) is None}\n")

        for name in names:
            exe = RELEASE / name
            if not exe.exists():
                results["commands"][name] = {"status": "MISSING_BINARY"}
                print(f"{name}: MISSING")
                continue
            entry: dict = {"binary": name, "binary_sha256": sha256(exe)}
            if not args.skip_formats:
                entry["input_formats"] = probe_formats(name, exe, bam, bed, work,
                                                       fixtures, ctx)
            entry["helper_dependency"] = probe_helper(name, exe, bam, bed, work,
                                                      minpath, ctx)
            entry["skip_plot_contract"] = probe_skip_plot(name, exe, bam, bed, work, ctx)
            entry["determinism"] = probe_stochastic(name, exe, bam, bed, work, ctx)
            results["commands"][name] = entry

            raw_fmts = entry.get("input_formats") or {}
            # `applies` and `primary_input_kind` are metadata about the probe, not
            # results; only the per-format entries have a "status".
            fmts = {k: v for k, v in raw_fmts.items() if isinstance(v, dict)}
            if raw_fmts.get("applies") is False:
                kind = raw_fmts.get("primary_input_kind")
                accepted, rejected = [f"n/a ({kind})"], []
            else:
                accepted = [k for k, v in fmts.items()
                            if v.get("status") == "ACCEPTED"]
                rejected = [k for k, v in fmts.items()
                            if v.get("status") in ("REJECTED", "NO_FIXTURE",
                                                   "NOT_APPLICABLE")]
            det = entry["determinism"]
            print(f"{name:22} formats accepted={accepted or '-'} "
                  f"rejected={rejected or '-'} "
                  f"helper={entry['helper_dependency'].get('needs_helper')} "
                  f"skip-plot={entry['skip_plot_contract'].get('suppresses_helper')} "
                  f"stdout-det={det.get('stdout_deterministic')} "
                  f"data-det={det.get('data_deterministic')} "
                  f"plot-det={det.get('plot_deterministic_after_timestamp_normalisation')}")

    # A summary block, because the release notes need the aggregate facts and
    # reading 33 rows to count them is how a claim goes stale.
    helper_yes = sorted(n for n, e in results["commands"].items()
                        if e.get("helper_dependency", {}).get("needs_helper") == "yes")
    nondet = sorted(n for n, e in results["commands"].items()
                    if e.get("determinism", {}).get("data_deterministic") is False
                    or e.get("determinism", {}).get("stdout_deterministic") is False)
    nondet_plots = sorted(
        n for n, e in results["commands"].items()
        if e.get("determinism", {}).get(
            "plot_deterministic_after_timestamp_normalisation") is False)
    format_support: dict[str, list[str]] = {}
    non_alignment: dict[str, list[str]] = {}
    for n, e in results["commands"].items():
        fmts = e.get("input_formats") or {}
        if fmts.get("applies") is False:
            non_alignment.setdefault(
                fmts.get("primary_input_kind", "unknown"), []).append(n)
            continue
        for fmt, v in fmts.items():
            if isinstance(v, dict) and v.get("status") == "ACCEPTED":
                format_support.setdefault(fmt, []).append(n)
    # Commands whose output legitimately varies between runs are a property to
    # DISCLOSE, not a defect to fix -- but only after confirming upstream varies too.
    # `junction_saturation` and `RPKM_saturation` resample reads with an unseeded
    # `random.shuffle`, so their per-percentile intermediate counts differ on every run
    # while the underlying data does not. The port reproduces that deliberately (see
    # `crates/commands/src/junction_saturation.rs`'s own module comment: "random.shuffle
    # has no seed anywhere upstream, so the ORIGINAL command's output already varies
    # run-to-run; this port shuffles with ..."). Asserting byte-equality for these
    # commands would fail for a correct implementation, and asserting byte-inequality
    # for a deterministic one would hide a regression, so the property is named
    # explicitly and the two cases are separated.
    # `divide_bam` is the third: it assigns reads to output files with an RNG that is
    # UNSEEDED unless `--seed` is given, which is upstream's own behaviour
    # (DIV-0017 -- upstream seeds `random.Random(seed)` and calls `rng.randrange(n)`
    # per distinct query name, and the port does the same with a different generator,
    # so the two do not agree on WHICH reads land in which file, only on the properties
    # that are contractual: paired mates stay together and subset sizes are roughly
    # equal). Deterministic when `--seed` is supplied. The probe runs without it, which
    # is the default a user gets.
    STOCHASTIC_BY_DESIGN = {"junction_saturation", "RPKM_saturation", "divide_bam"}
    stochastic_disclosed = sorted(set(nondet) & STOCHASTIC_BY_DESIGN)
    stochastic_unexplained = sorted(set(nondet) - STOCHASTIC_BY_DESIGN)

    # Commands that can only be judged with an external helper installed, recorded as
    # their own list because "helper absent" and "capability absent" are different
    # answers and a release note needs to give the first.
    helper_absent = sorted(
        n for n, e in results["commands"].items()
        if "requires an external helper" in str(
            e.get("determinism", {}).get("detail", "")))

    results["summary"] = {
        "commands_probed": len(results["commands"]),
        "not_judged_because_a_helper_is_absent": helper_absent,
        "requires_helper_at_runtime": helper_yes,
        "non_deterministic_by_design": stochastic_disclosed,
        "non_deterministic_unexplained": stochastic_unexplained,
        "nondeterministic_data": nondet,
        "nondeterministic_plots": nondet_plots,
        "accepts_format": {k: sorted(v) for k, v in sorted(format_support.items())},
        "primary_input_kind": {k: sorted(v) for k, v in sorted(non_alignment.items())},
    }
    print(f"\n{len(helper_yes)} commands require a helper program at runtime: {helper_yes}")
    print(f"non-deterministic across two identical runs, by design "
          f"(unseeded resampling, as upstream): {stochastic_disclosed}")
    if stochastic_unexplained:
        print(f"  UNEXPLAINED non-determinism, which is a finding: "
              f"{stochastic_unexplained}")
    print(f"non-deterministic PLOTS (timestamps normalised): {nondet_plots}")
    if helper_absent:
        print(f"no determinism verdict, because a helper they need is absent: "
              f"{helper_absent}")

    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2) + "\n")
        print(f"written to {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
