#!/usr/bin/env bash
# Builds a standalone, self-contained distribution archive: an ALLOWLISTED set
# of release binaries, the original upstream RSeQC script-name aliases, the
# documentation and licence metadata a redistribution needs, and a
# machine-readable build/capability manifest. Packaged as a checksummed
# .tar.gz and then smoke-tested on a real workload OUTSIDE the source tree.
#
# Usage: scripts/build-release-archive.sh [output-dir]
#   output-dir   Where to write the archive and checksum file.
#                Default: dist/ relative to the repo root.
#
# WHAT THIS DOES NOT DO, disclosed rather than silently skipped. Each of these
# was an audit finding; the first three are now closed, the rest remain recorded.
#   - No multi-platform build: HOST target triple only. macOS/Windows/ARM
#     artifacts ship only after installed-artifact validation on each platform
#     (release Stage C).
#   - The smoke test runs on this host with the repository's oracle environment
#     on PATH. It exercises the EXTRACTED archive against a real BAM workload
#     rather than only `--help`, which is the audit's requirement, but it is not
#     a clean room: wigToBigWig and htseq-count remain absent here, so rows that
#     need them are labelled rather than claimed.
#   - This project's own licence variant is still blocked on DIV-0003 (upstream's
#     own licence metadata is internally inconsistent: GPL-3.0-or-later in the
#     README vs GPLv3 text vs a GPLv2 classifier). The repository's existing
#     LICENSE file is bundled as-is and the pending decision is recorded in the
#     manifest, so a redistributor sees the ambiguity rather than inheriting it
#     silently.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-$REPO_ROOT/dist}"

VERSION="$(grep -m1 '^version' "$REPO_ROOT/Cargo.toml" | sed -E 's/.*"([^"]+)".*/\1/')"
TARGET_TRIPLE="$(rustc -vV | sed -n 's/^host: //p')"
ARCHIVE_NAME="rseqc-rust-${VERSION}-${TARGET_TRIPLE}"
STAGE_DIR="$(mktemp -d)"
SMOKE_DIR=""
cleanup() { rm -rf "$STAGE_DIR" ${SMOKE_DIR:+"$SMOKE_DIR"}; }
trap cleanup EXIT

# ---------------------------------------------------------------------------------
# The allowlist
# ---------------------------------------------------------------------------------
# An earlier version staged every executable file it found at the top of
# target/release, which is how build scripts, example binaries and unrelated
# targets would silently end up in a published archive. The list below is the
# complete set of names this archive claims to ship, and the staging step fails
# if any of them is missing rather than publishing a partial archive.
#
# Derived from scripts/install-aliases.sh's PAIRS, which is the single source of
# truth for binary-to-alias naming. Kept as an explicit list rather than parsed
# out of that script so that adding a command to the alias table does not
# silently widen a published archive's contents.
ALLOWLIST="bam_stat split_paired_bam bam2fq divide_bam read_GC read_NVC
read_quality read_duplication clipping_profile insertion_profile
deletion_profile mismatch_profile infer_experiment RNA_fragment_size
read_distribution inner_distance junction_annotation split_bam
junction_saturation tin FPKM_count RPKM_saturation FPKM_UQ read_hexamer
bam2wig geneBody_coverage geneBody_coverage2 normalize_bigwig
overlay_bigwig sc_bamStat sc_editMatrix sc_seqQual sc_seqLogo"

# Documentation and metadata a redistribution of this project's code requires.
# Cargo.toml is bundled because it is the crate's own canonical metadata: it carries
# the repository/homepage fields a `cargo`-aware reader looks up, and it is one of the
# two files check_release_metadata.py refuses to publish with a placeholder URL.
DOC_FILES="LICENSE README.md CHANGELOG.md CITATION.cff Cargo.toml"
OPTIONAL_DOC_FILES="testing.md CONTRIBUTING.md"
# Files that must be present for the archive to be compliant.
REQUIRED_DOC_FILES="LICENSE"

echo "Building release binaries ..." >&2
# The build invocation is part of the artifact, not an implementation detail.
# `cargo build -p rseqc-cli --release` and `cargo build --workspace --release`
# produce DIFFERENT bytes for the same binary on the same source (verified: bam_stat
# hashes d9841652... under --workspace and c497e29f... under -p). So a capability
# matrix probed after a per-package build describes different bytes than this archive
# ships, even though the source revision is identical.
#
# That is why --workspace, --release and --locked are all spelled out here, and why
# scripts/check_release_metadata.py compares each command's recorded binary hash
# against the hash of the binary actually staged here. Provenance that binds to a
# revision but not to the bytes is not provenance.
(cd "$REPO_ROOT" && cargo build --workspace --release --locked)

PKG_DIR="$STAGE_DIR/$ARCHIVE_NAME"
mkdir -p "$PKG_DIR/bin" "$PKG_DIR/doc" "$PKG_DIR/examples"

missing=0
staged=0
for bin in $ALLOWLIST; do
    src="$REPO_ROOT/target/release/$bin"
    if [ ! -x "$src" ]; then
        echo "error: allowlisted binary missing: $bin (refusing to publish a partial archive)" >&2
        missing=$((missing + 1))
        continue
    fi
    cp "$src" "$PKG_DIR/bin/"
    staged=$((staged + 1))
done
if [ "$missing" -gt 0 ]; then
    echo "error: $missing allowlisted binary/binaries absent; archive not published" >&2
    exit 1
fi
echo "Staged $staged allowlisted binaries ..." >&2

# Nothing unexpected may appear in bin/: the staging step copied by name, but a
# future edit that switched to a glob would reintroduce the original defect, so
# the invariant is asserted rather than assumed.
staged_count="$(find "$PKG_DIR/bin" -maxdepth 1 -type f | wc -l | tr -d ' ')"
expected_count="$(echo $ALLOWLIST | wc -w | tr -d ' ')"
if [ "$staged_count" != "$expected_count" ]; then
    echo "error: bin/ has $staged_count files, allowlist expects $expected_count" >&2
    exit 1
fi

echo "Installing original-name aliases ..." >&2
"$REPO_ROOT/scripts/install-aliases.sh" "$PKG_DIR/bin" --copy

for f in $REQUIRED_DOC_FILES; do
    if [ ! -f "$REPO_ROOT/$f" ]; then
        echo "error: required file missing from the repository: $f" >&2
        exit 1
    fi
done
for f in $DOC_FILES $OPTIONAL_DOC_FILES; do
    if [ -f "$REPO_ROOT/$f" ]; then
        cp "$REPO_ROOT/$f" "$PKG_DIR/doc/$f"
    elif echo "$REQUIRED_DOC_FILES" | grep -qw "$f"; then
        echo "error: required doc file missing: $f" >&2
        exit 1
    else
        echo "note: optional doc file absent, not bundled: $f" >&2
    fi
done

# Third-party notices: this project's own dependency licences, so a
# redistributor is not left to assemble them.
echo "Recording dependency licences ..." >&2
{
    echo "Third-party notices for $ARCHIVE_NAME"
    echo "Generated $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo
    echo "This binary distribution embeds the Rust crates listed below as static"
    echo "code. Their licence terms apply to those crates, not to this project's"
    echo "own licence (see doc/LICENSE)."
    echo
    if [ -f "$REPO_ROOT/target/release/../release/deps" ]; then :; fi
    (cd "$REPO_ROOT" && cargo metadata --format-version 1 --locked 2>/dev/null) \
        | python3 -c '
import json, sys
try:
    meta = json.load(sys.stdin)
except Exception:
    print("  (cargo metadata unavailable; run cargo fetch then rebuild)")
    sys.exit(0)
seen = {}
for node in meta.get("packages", []):
    if node.get("source") is None:
        continue  # workspace member
seen = {}
for node in meta.get("packages", []):
    name = node["name"]
    if name in seen:
        continue
    seen[name] = node.get("license") or "(not declared)"
for name in sorted(seen):
    print(f"  {name}: {seen[name]}")
' || echo "  (could not enumerate dependency licences)"
} > "$PKG_DIR/doc/THIRD_PARTY_NOTICES.txt"

# Small example inputs/outputs, so an extracted archive can be validated without
# a real dataset in hand.
if [ -d "$REPO_ROOT/datasets/aligned/tiny" ]; then
    cp "$REPO_ROOT/datasets/aligned/tiny"/* "$PKG_DIR/examples/" 2>/dev/null || true
fi

# Refuse to publish metadata with unresolved placeholders. The audit found
# `TBD` in both Cargo.toml and CITATION.cff, which would ship a URL that resolves to
# nobody. The permanent URL is a maintainer decision (release Stage A); this check
# makes it impossible to reach the publish step without making it.
#
# ALLOW_UNRESOLVED_URL=1 downgrades this to a warning. It exists so the archive can
# still be BUILT and smoke-tested on a developer machine before the URL decision is
# made, and it is deliberately refused in CI: a tag push has no excuse, and the
# release workflow never sets it.
# The check runs against the *staged* copy of the metadata, not the source tree, so
# what is verified is the file that will actually be unpacked by a user. Checking the
# source tree instead would pass on a build that had already dropped or rewritten the
# file after staging.
if [ "${ALLOW_UNRESOLVED_URL:-0}" = "1" ]; then
  echo "Checking release metadata (placeholders downgraded to warnings) ..." >&2
  python3 "$REPO_ROOT/scripts/check_release_metadata.py" --root "$PKG_DIR/doc" || true
else
  echo "Checking release metadata ..." >&2
  python3 "$REPO_ROOT/scripts/check_release_metadata.py" --strict --root "$PKG_DIR/doc" \
    || { echo "error: unresolved placeholders in release metadata (see above)" >&2
         echo "       set the permanent repository URL, or set ALLOW_UNRESOLVED_URL=1" >&2
         echo "       to build a local test archive" >&2
         exit 1; }
fi

# ---------------------------------------------------------------------------------
# Build / capability manifest
# ---------------------------------------------------------------------------------
echo "Writing build/capability manifest ..." >&2
ARCHIVE_NAME="$ARCHIVE_NAME" TARGET_TRIPLE="$TARGET_TRIPLE" VERSION="$VERSION" \
PKG_DIR="$PKG_DIR" REPO_ROOT="$REPO_ROOT" python3 - <<'PY' > "$PKG_DIR/manifest.json"
import json, os, subprocess, sys, hashlib
from pathlib import Path

pkg = Path(os.environ["PKG_DIR"])
repo = Path(os.environ["REPO_ROOT"])


def per_command_capabilities(repo):
    """Per-command support, from the measured matrix rather than from prose.

    The audit's Stage A asks for a capability manifest that records, per command, the
    inputs it accepts, which outputs depend on a helper, whether it is stochastic, and
    which upstream quirks it reproduces. This function embeds
    benchmarks/capability-matrix.json, which `verification/capability_matrix.py`
    produces by EXECUTING each command -- feeding it a real file of each input format,
    re-running it with a PATH containing no Rscript, and running it twice to compare
    artifacts.

    The matrix is embedded rather than regenerated here because a build must not depend
    on probe fixtures being present, and because regenerating would mean the archive's
    capability claims could differ from the ones that were reviewed. A missing or stale
    matrix is therefore recorded as such, and `check_release_metadata.py` requires the
    per-command block to be present and to cover every shipped binary.
    """
    import json as _json
    matrix_path = repo / "benchmarks" / "capability-matrix.json"
    if not matrix_path.is_file():
        return {
            "source": "unavailable",
            "detail": f"{matrix_path} is absent, so no per-command capability "
                      f"record is claimed. Run verification/capability_matrix.py.",
            "commands": {},
        }
    matrix = _json.loads(matrix_path.read_text())
    cmds = matrix.get("commands", {})
    per = {"source": "benchmarks/capability-matrix.json",
           "source_note": matrix.get("generated_from", ""),
           "fixture": matrix.get("fixture", {}),
           "summary": matrix.get("summary", {}),
           "commands": {}}
    for name, entry in sorted(cmds.items()):
        formats = entry.get("input_formats") or {}
        accepted = sorted(k for k, v in formats.items()
                          if isinstance(v, dict) and v.get("status") == "ACCEPTED")
        det = entry.get("determinism") or {}
        per["commands"][name] = {
            "input_formats_accepted": accepted,
            "input_formats_measured": sorted(
                k for k, v in formats.items() if isinstance(v, dict)),
            "primary_input_kind": formats.get("primary_input_kind"),
            "known_input_restriction": formats.get("known_restriction"),
            # From the MEASURED helper name, not from needs_helper: 16 commands need
            # some external helper, but FPKM-UQ's is htseq-count, so the generic flag
            # would claim it requires Rscript. The release metadata checker now fails
            # that specific mislabelling.
            "requires_rscript": (entry.get("helper_dependency") or {}).get(
                "requires_rscript") is True,
            "missing_helper": (entry.get("helper_dependency") or {}).get("helper"),
            "reports_a_missing_helper": (entry.get("helper_dependency") or {}).get(
                "reports_missing_helper"),
            "has_skip_plot": (entry.get("skip_plot_contract") or {}).get(
                "has_skip_plot"),
            "skip_plot_suppresses_rscript": (entry.get("skip_plot_contract") or {}).get(
                "suppresses_helper"),
            "reproducible": det.get("stdout_deterministic") if det.get(
                "stdout_deterministic") is not None else det.get("data_deterministic"),
            "reproducibility_note": det.get("detail"),
            "binary_sha256": entry.get("binary_sha256"),
        }
    return per

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def sh(*cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception:
        return "unknown"

bins = {}
for p in sorted((pkg / "bin").iterdir()):
    # Aliases are byte-identical copies of their binaries, so listing both would
    # double the count and suggest 66 distinct commands where there are 33.
    if p.is_file() and not p.name.endswith(".py"):
        bins[p.name] = {"sha256": sha256(p), "bytes": p.stat().st_size}
aliases = {p.name for p in (pkg / "bin").iterdir() if p.name.endswith(".py")}

docs = {p.name: sha256(p) for p in sorted((pkg / "doc").iterdir()) if p.is_file()}

def lock_digest():
    lock = repo / "compatibility" / "upstream.lock"
    return sha256(lock) if lock.exists() else None

manifest = {
    "name": os.environ["ARCHIVE_NAME"],
    "version": os.environ["VERSION"],
    "schema": 1,
    "target": {
        "triple": os.environ["TARGET_TRIPLE"],
        "family": "unix",
        "os": "linux",
        "arch": "x86_64",
        "abi_floor": {
            "glibc": "2.28",
            "note": "Rust's stable x86_64-unknown-linux-gnu baseline. Newer "
                    "binaries linked against a newer glibc will not run on "
                    "older distributions; this floor is the minimum the archive "
                    "is claimed to support.",
        },
        "musl": False,
    },
    "built_from": {
        "git_commit": sh("git", "-C", str(repo), "rev-parse", "HEAD"),
        "git_dirty": bool(sh("git", "-C", str(repo), "status", "--porcelain")),
        "rustc": sh("rustc", "--version"),
        "cargo": sh("cargo", "--version"),
    },
    "binaries": bins,
    "aliases": sorted(aliases),
    "docs": docs,
    "capabilities": {
        "commands": len(bins),
        "upstream_aliases": "each binary is also present under its original "
                            "upstream script name (bam_stat.py, ...); aliases are "
                            "byte-identical copies, so the archive holds "
                            f"{len(bins) + len(aliases)} files for {len(bins)} commands",
        "native_plots": "most plotting is implemented natively; every command that "
                        "generates an R script still invokes Rscript unless "
                        "--skip-plot is given, and the per-command capability matrix "
                        "below names which",
        "python_api": "NOT IMPLEMENTED in this release; the rseqc-python crate is a "
                      "stub. Do not announce Python API compatibility.",
        "helper_programs_required": {
            "Rscript": "optional; required by the commands listed in "
                       "per_command.requires_rscript",
            "wigToBigWig": "optional; bam2wig/normalize_bigwig write bedGraph "
                           "without it and cannot produce BigWig",
            "htseq-count": "optional; required by FPKM_UQ only",
        },
        "input_formats": ["BAM", "CRAM", "SAM", "BED12", "BigWig", "FASTQ", "FASTA"],
    },
    "per_command": per_command_capabilities(repo),
    "known_limitations": [
        "Python API is a stub and is excluded from this release's claims.",
        "Commands documented as experimental are included but unqualified; their "
        "numerical results have not been through the endpoint suite.",
        "Upstream's own licence metadata is internally inconsistent (GPL-3.0-or-later "
        "README vs GPLv3 text vs GPLv2 classifier); doc/LICENSE is the repository's "
        "existing file and the variant decision (DIV-0003) is still open.",
        "Build provenance is bound to the git commit above. A dirty tree is flagged, "
        "not rejected, so a locally built archive from uncommitted changes is "
        "distinguishable from a clean one.",
        "geneBody_coverage.py does NOT currently reproduce upstream's gene-body "
        "coverage curve on real data (compatibility/divergences.yaml DIV-0024, open): "
        "on an 8.2M-record rat alignment 76 of 100 bins differ, by up to 839 reads, "
        "in both directions, bisected to a single 904-base transcript. Identified "
        "cause: pysam/htslib max_depth semantics -- every other pileup filter has been "
        "verified equivalent, one mechanism at a time, on a fixture built so that each "
        "base range exercises exactly one. The fix is not yet written. No speedup or "
        "scientific claim is made for this command.",
        "The permanent repository URL is a maintainer decision and is not yet set "
        "(release Stage A). An archive built before that decision is a local test "
        "artefact, not a publication; the release workflow refuses to publish one.",
    ],
    "reference_environment": {
        "oracle_lock_sha256": lock_digest(),
        "note": "Digest of compatibility/upstream.lock, for attributing differential "
                "results to a pinned oracle.",
    },
}
print(json.dumps(manifest, indent=2))
PY

# ---------------------------------------------------------------------------------
# Archive + checksum
# ---------------------------------------------------------------------------------
mkdir -p "$OUT_DIR"
ARCHIVE_PATH="$OUT_DIR/${ARCHIVE_NAME}.tar.gz"

echo "Creating archive: $ARCHIVE_PATH" >&2
tar -C "$STAGE_DIR" -czf "$ARCHIVE_PATH" "$ARCHIVE_NAME"

echo "Computing checksum ..." >&2
(cd "$OUT_DIR" && sha256sum "$(basename "$ARCHIVE_PATH")" > "${ARCHIVE_NAME}.tar.gz.sha256")

# ---------------------------------------------------------------------------------
# Smoke test: the EXTRACTED archive, outside the source tree, on a real workload
# ---------------------------------------------------------------------------------
# The audit's requirement, and the reason this is no longer just `--help`:
# validating the archive that will actually be installed, with real inputs, in a
# directory that is not the repository, so a missing runtime dependency or a path
# assumption baked in at build time shows up here.
echo "Smoke-testing the extracted archive outside the source tree ..." >&2
# Extracted into a directory whose name contains a space, because that is where an
# unquoted path in a wrapper script or a Makefile breaks, and because the audit
# requires paths with spaces to be exercised rather than assumed.
SMOKE_DIR="$(mktemp -d)"
SMOKE_PARENT="$SMOKE_DIR/a directory with spaces"
mkdir -p "$SMOKE_PARENT"
tar -C "$SMOKE_PARENT" -xzf "$ARCHIVE_PATH"
PKG="$SMOKE_PARENT/$ARCHIVE_NAME"

# A fresh PATH containing only the archive's bin/ plus the system basics, so
# nothing in the repository's target/ or oracle/ can satisfy a missing file.
export PATH="$PKG/bin:/usr/local/bin:/usr/bin:/bin"

for required in LICENSE manifest.json; do
    if [ ! -e "$PKG/doc/$required" ] && [ ! -e "$PKG/$required" ]; then
        echo "error: extracted archive is missing $required" >&2
        exit 1
    fi
done

"$PKG/bin/bam_stat.py" --help >/dev/null
"$PKG/bin/infer_experiment.py" --help >/dev/null

EX_BAM="$PKG/examples/tiny.bam"
if [ -f "$EX_BAM" ]; then
    "$PKG/bin/bam_stat.py" -i "$EX_BAM" > "$SMOKE_DIR/bam_stat.out" 2>&1 \
        || { echo "error: bam_stat failed on the bundled example" >&2
             cat "$SMOKE_DIR/bam_stat.out" >&2; exit 1; }
    grep -qi 'total' "$SMOKE_DIR/bam_stat.out" \
        || { echo "error: bam_stat produced no recognisable report" >&2
             cat "$SMOKE_DIR/bam_stat.out" >&2; exit 1; }

    # An annotation-consuming command, which is the path most likely to depend on
    # something the archive failed to bundle.
    "$PKG/bin/infer_experiment.py" -i "$EX_BAM" -r "$PKG/examples/model.bed12" \
        > "$SMOKE_DIR/infer.out" 2>&1 \
        || { echo "error: infer_experiment failed on the bundled example" >&2
             cat "$SMOKE_DIR/infer.out" >&2; exit 1; }

    # Writing files, not just printing: proves the process can create output in
    # a directory the build never touched. bam2fq names paired output <p>.R1.fastq
    # and <p>.R2.fastq, matching upstream, so accept either spelling rather than
    # asserting one and failing a correct build.
    ( cd "$SMOKE_DIR" && "$PKG/bin/bam2fq.py" -i "$EX_BAM" -o "$SMOKE_DIR/r1" \
        >/dev/null 2>&1 ) \
        || { echo "error: bam2fq failed to write output" >&2; exit 1; }
    if [ ! -s "$SMOKE_DIR/r1.R1.fastq" ] && [ ! -s "$SMOKE_DIR/r1.fastq" ]; then
        echo "error: bam2fq produced no FASTQ (looked for r1.R1.fastq and r1.fastq)" >&2
        exit 1
    fi
    echo "  real-workload smoke test passed (bam_stat, infer_experiment, bam2fq)" >&2
else
    echo "  no example BAM bundled; only --help was exercised" >&2
    echo "  (this is weaker than the audit requires; rebuild with" >&2
    echo "   datasets/aligned/tiny present to close it)" >&2
fi

# A command invoked by its upstream alias, from a working directory that is not the
# source tree and whose name contains a space, writing its output there.
( cd "$SMOKE_PARENT" && "$PKG/bin/infer_experiment.py" -i "$PKG/examples/tiny.bam" \
    -r "$PKG/examples/model.bed12" > "$SMOKE_PARENT/spaces.out" 2>&1 ) \
  || { echo "error: upstream alias failed from a path containing spaces" >&2
       cat "$SMOKE_PARENT/spaces.out" >&2; exit 1; }
grep -q 'Fraction of reads' "$SMOKE_PARENT/spaces.out" \
  || { echo "error: no report from the spaces run" >&2; exit 1; }
echo "  spaces-in-path invocation passed (infer_experiment.py)" >&2

# An output path containing spaces. read_hexamer takes FASTQ/FASTA rather than a
# BAM, so the fixture is derived from the bundled example rather than passing it a
# file it cannot read -- otherwise this check would fail for the wrong reason and
# report a spaces defect that is not one.
( cd "$SMOKE_PARENT" && "$PKG/bin/bam2fq.py" -i "$PKG/examples/tiny.bam" \
    -o "$SMOKE_PARENT/reads" >/dev/null 2>&1 ) \
  || { echo "error: could not write FASTQ to a spaced path" >&2; exit 1; }
FQ="$SMOKE_PARENT/reads.R1.fastq"
[ -s "$FQ" ] || { echo "error: no FASTQ produced at a spaced path" >&2; exit 1; }
( cd "$SMOKE_PARENT" && "$PKG/bin/read_hexamer.py" -i "$FQ" \
    -o "$SMOKE_PARENT/hex out.txt" >/dev/null 2>&1 ) \
  && [ -s "$SMOKE_PARENT/hex out.txt" ] \
  || { echo "error: could not write an output file whose path contains spaces" >&2; exit 1; }
echo "  spaces-in-output-path write passed (read_hexamer.py)" >&2

echo "Extraction + alias + real-workload + spaces smoke test passed" >&2

echo "$ARCHIVE_PATH"
