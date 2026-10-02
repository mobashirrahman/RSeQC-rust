#!/usr/bin/env python3
"""Where each specialised verification fixture lives, built on demand.

`verification/capability_matrix.py` needs an input per command, and the single-cell
commands, FPKM-UQ and the BigWig readers need inputs the tiny panel cannot provide:
a BAM carrying CB/UB/RE/xf tags, a barcode FASTQ, a GTF plus the gene-information file
FPKM-UQ parses, and a real BigWig track.

Those fixtures already exist for the differential suite. Re-declaring their paths here
would create a second place to update when one moves, so the existing generator scripts
are invoked and their documented outputs reused. The single-cell BAMs are rebuilt here
into a temporary location rather than written into the repository, because
`verification/fixtures/` is committed and the capability matrix is a probe, not a
differential input.

Each accessor returns a `Path`, or raises if the input cannot be produced -- callers
report that rather than silently substituting the wrong input, which is how a
"no determinism verdict" row should be read.
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FIXTURES = REPO / "verification" / "fixtures"
TRACK = FIXTURES / "track"

# The oracle interpreter is the one that has pysam. Falling back to sys.executable
# keeps this usable under a system python that also has pysam, and the error otherwise
# names the interpreter that could not build the fixture.
def oracle_python() -> str:
    candidate = REPO / "oracle" / "venv" / "bin" / "python3"
    return str(candidate if candidate.exists() else Path(sys.executable))


def _build(generator: str, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    script = FIXTURES / generator
    if not script.is_file():
        raise FileNotFoundError(f"{script} not found")
    proc = subprocess.run([oracle_python(), str(script), str(out)],
                          capture_output=True, text=True, timeout=600)
    if proc.returncode != 0 or not out.exists():
        raise RuntimeError(
            f"{generator} failed (exit {proc.returncode}): "
            f"{(proc.stderr or proc.stdout).strip()[-200:]}"
        )
    return out


# Scratch for the single-cell BAMs. Module-level so the path is stable for the process
# lifetime; the caller decides when to clean it up.
_SCRATCH = Path(tempfile.mkdtemp(prefix="rseqc-capability-fixtures-"))


def sc_bamstat() -> Path:
    """A BAM with CB/UB/RE/xf tags, as sc_bamStat requires."""
    return _build("make_sc_bamstat_fixture.py", _SCRATCH / "sc_bamstat_basic.bam")


def sc_editmatrix() -> Path:
    """A BAM with CR/CB/UR/UB tags, as sc_editMatrix requires."""
    return _build("make_sc_editmatrix_fixture.py", _SCRATCH / "sc_editmatrix_basic.bam")


def sc_seqqual() -> Path:
    """A barcode FASTQ, as sc_seqQual requires."""
    path = FIXTURES / "regression_sc_seqqual.fq"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; it is committed, so this "
                                f"indicates an incomplete checkout")
    return path


def sc_seqlogo() -> Path:
    """A read FASTA/FASTQ, as sc_seqLogo requires."""
    path = FIXTURES / "regression_sc_seqlogo.fa"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; it is committed, so this "
                                f"indicates an incomplete checkout")
    return path


def fpkm_uq_gtf() -> Path:
    """The GTF FPKM-UQ parses for its exon model."""
    path = FIXTURES / "regression_fpkm_uq_dummy.gtf"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found")
    return path


def fpkm_uq_info() -> Path:
    """The gene-information file FPKM-UQ requires alongside the GTF."""
    path = FIXTURES / "regression_fpkm_uq_genes.info.txt"
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found")
    return path


def track_bigwig() -> Path:
    """A real BigWig track, for the three BigWig-reading commands.

    pyBigWig's own test file, copied unmodified (see
    `crates/formats/tests/fixtures/README.md`), so the track is independently
    verifiable rather than one this project generated. The track fixtures under
    `verification/fixtures/track/` are preferred when present, since they are the ones
    the differential suite uses.
    """
    for candidate in (TRACK / "track_signal.bw",
                      REPO / "crates" / "formats" / "tests" / "fixtures" /
                      "pybigwig_test.bw"):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"no BigWig fixture found; expected {TRACK / 'track_signal.bw'} or "
        f"crates/formats/tests/fixtures/pybigwig_test.bw"
    )


def chrom_sizes_for(bam: Path, out: Path) -> Path:
    """Write a two-column chromosome-size file for `bam2wig`."""
    import pysam
    with pysam.AlignmentFile(str(bam)) as f:
        out.write_text("\n".join(f"{e['SN']}\t{e['LN']}"
                                 for e in f.header.to_dict()["SQ"]) + "\n")
    return out
