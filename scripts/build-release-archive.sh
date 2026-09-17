#!/usr/bin/env bash
# Builds a standalone, self-contained distribution archive: release
# binaries plus the original upstream RSeQC script-name aliases (via
# install-aliases.sh --copy, since symlinks may not survive extraction/
# transport), packaged as a checksummed .tar.gz. Part of PORTING_PLAN.md
# Step 10 ("Produce versioned source archives, supported binary
# distributions...").
#
# Usage: scripts/build-release-archive.sh [output-dir]
#   output-dir   Where to write the archive and checksum file.
#                Default: dist/ relative to the repo root.
#
# What this does NOT do yet (disclosed, not silently skipped):
#   - No LICENSE file is bundled: this project's own release license
#     metadata is blocked on DIV-0003 (upstream's own license metadata
#     is internally inconsistent -- GPL-3.0-or-later vs GPLv3 text vs a
#     GPLv2 classifier -- resolving which variant THIS project releases
#     under is a user policy decision, not something a build script
#     should decide unilaterally). Do not add a LICENSE file to this
#     archive until that decision is made and recorded.
#   - No README is bundled (none exists in the repo yet -- PORTING_PLAN
#     Step 11 item 1).
#   - Single-platform only: builds for the HOST target triple. Add
#     cross-compilation (PORTING_PLAN's "supported binary distributions"
#     plural, multiple OS/ABI targets) as a separate, later step.
#   - Does not test the archive in a clean-room environment without
#     Python/R/wigToBigWig present (PORTING_PLAN Step 10's explicit
#     requirement) -- only confirms extraction + alias invocation here
#     on the SAME machine the archive was built on.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-$REPO_ROOT/dist}"

VERSION="$(grep -m1 '^version' "$REPO_ROOT/Cargo.toml" | sed -E 's/.*"([^"]+)".*/\1/')"
TARGET_TRIPLE="$(rustc -vV | sed -n 's/^host: //p')"
ARCHIVE_NAME="rseqc-rust-${VERSION}-${TARGET_TRIPLE}"
STAGE_DIR="$(mktemp -d)"
trap 'rm -rf "$STAGE_DIR"' EXIT

echo "Building release binaries ..." >&2
(cd "$REPO_ROOT" && cargo build --workspace --release --locked)

PKG_DIR="$STAGE_DIR/$ARCHIVE_NAME"
mkdir -p "$PKG_DIR/bin"

echo "Staging binaries ..." >&2
find "$REPO_ROOT/target/release" -maxdepth 1 -type f -executable -exec cp {} "$PKG_DIR/bin/" \;

echo "Installing original-name aliases ..." >&2
"$REPO_ROOT/scripts/install-aliases.sh" "$PKG_DIR/bin" --copy

mkdir -p "$OUT_DIR"
ARCHIVE_PATH="$OUT_DIR/${ARCHIVE_NAME}.tar.gz"

echo "Creating archive: $ARCHIVE_PATH" >&2
tar -C "$STAGE_DIR" -czf "$ARCHIVE_PATH" "$ARCHIVE_NAME"

echo "Computing checksum ..." >&2
(cd "$OUT_DIR" && sha256sum "$(basename "$ARCHIVE_PATH")" > "${ARCHIVE_NAME}.tar.gz.sha256")

echo "Smoke-testing extraction on this machine ..." >&2
SMOKE_DIR="$(mktemp -d)"
trap 'rm -rf "$STAGE_DIR" "$SMOKE_DIR"' EXIT
tar -C "$SMOKE_DIR" -xzf "$ARCHIVE_PATH"
"$SMOKE_DIR/$ARCHIVE_NAME/bin/bam_stat.py" --help >/dev/null
echo "Extraction + alias smoke test passed" >&2

echo "$ARCHIVE_PATH"
