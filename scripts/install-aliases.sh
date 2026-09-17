#!/usr/bin/env bash
# Creates the original upstream RSeQC script-name aliases (e.g.
# `bam_stat.py`, `FPKM-UQ.py`) pointing at this port's compiled
# binaries, so scripts/pipelines invoking commands by their original
# name work unmodified. Cargo bin names can't contain '.' or match
# upstream's exact casing/hyphenation in every case, so this is a
# separate step rather than something Cargo itself can produce.
#
# Usage: scripts/install-aliases.sh [bin-dir] [--copy]
#   bin-dir   Directory containing the compiled rseqc-cli binaries.
#             Default: target/release relative to the repo root.
#   --copy    Copy instead of symlink. Use this when distributing a
#             standalone archive, since symlinks may not survive
#             extraction/transport to another machine.
#
# Run after building, e.g.:
#   cargo build --workspace --release
#   scripts/install-aliases.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="$REPO_ROOT/target/release"
MODE="link"

for arg in "$@"; do
    case "$arg" in
        --copy) MODE="copy" ;;
        *) BIN_DIR="$arg" ;;
    esac
done

if [ ! -d "$BIN_DIR" ]; then
    echo "error: bin directory not found: $BIN_DIR" >&2
    echo "Build first, e.g.: cargo build --workspace --release" >&2
    exit 1
fi

# bin_name:py_name pairs, one per RSeQC command (33 total). Keep in
# sync with the "packaging aliases this to ..." comments in
# crates/cli/Cargo.toml -- that file is the source of truth.
PAIRS="
bam_stat:bam_stat.py
split_paired_bam:split_paired_bam.py
bam2fq:bam2fq.py
divide_bam:divide_bam.py
read_GC:read_GC.py
read_NVC:read_NVC.py
read_quality:read_quality.py
read_duplication:read_duplication.py
clipping_profile:clipping_profile.py
insertion_profile:insertion_profile.py
deletion_profile:deletion_profile.py
mismatch_profile:mismatch_profile.py
infer_experiment:infer_experiment.py
RNA_fragment_size:RNA_fragment_size.py
read_distribution:read_distribution.py
inner_distance:inner_distance.py
junction_annotation:junction_annotation.py
split_bam:split_bam.py
junction_saturation:junction_saturation.py
tin:tin.py
FPKM_count:FPKM_count.py
RPKM_saturation:RPKM_saturation.py
FPKM_UQ:FPKM-UQ.py
read_hexamer:read_hexamer.py
bam2wig:bam2wig.py
geneBody_coverage:geneBody_coverage.py
geneBody_coverage2:geneBody_coverage2.py
normalize_bigwig:normalize_bigwig.py
overlay_bigwig:overlay_bigwig.py
sc_bamStat:sc_bamStat.py
sc_editMatrix:sc_editMatrix.py
sc_seqQual:sc_seqQual.py
sc_seqLogo:sc_seqLogo.py
"

created=0
missing=0

for pair in $PAIRS; do
    bin_name="${pair%%:*}"
    py_name="${pair##*:}"
    src="$BIN_DIR/$bin_name"
    dst="$BIN_DIR/$py_name"

    if [ ! -x "$src" ]; then
        echo "warning: binary not found, skipping: $src" >&2
        missing=$((missing + 1))
        continue
    fi

    rm -f "$dst"
    if [ "$MODE" = "copy" ]; then
        cp "$src" "$dst"
    else
        ln -s "$bin_name" "$dst"
    fi
    created=$((created + 1))
done

echo "Created $created alias(es) in $BIN_DIR" >&2
if [ "$missing" -gt 0 ]; then
    echo "warning: $missing binary(ies) were missing -- build may be incomplete" >&2
    exit 1
fi
