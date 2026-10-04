#!/usr/bin/env python3
"""Verify the frozen upstream RSeQC source tree against its recorded digest.

`compatibility/upstream.lock` carries a `source_tree_sha256_scripts_and_src`
field so a recorded differential result can be attributed to a specific
upstream tree. Nothing checked it, which is why it drifted: the digest it
recorded was taken over `tar -cf - -C oracle/upstream-src .`, which sweeps up
`__pycache__`, `*.pyc` and the `qcmodule.egg-info` that `pip install -e`
generates. Those are build artifacts of *this host*, so the digest changed the
first time anything was imported and could never be reproduced by anyone --
including the person who wrote it.

This recomputes the digest over tracked source only, deterministically:

  * paths relative to `oracle/upstream-src`, `scripts/` and `src/` only
  * build artifacts excluded: `__pycache__/`, `*.pyc`, `*.egg-info/`
  * `LC_ALL=C sort`, so ordering does not depend on locale
  * per-file sha256, then sha256 of the `<digest>  <path>` listing

Prints the digest it computed, so a refresh can copy the value straight into
the lock. Exits 1 on mismatch or on a missing tree.

    python3 verification/check_oracle_source_digest.py
    python3 verification/check_oracle_source_digest.py --update-lock
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TREE = REPO / "oracle/upstream-src"
LOCK = REPO / "compatibility/upstream.lock"
FIELD = "source_tree_sha256_scripts_and_src"

EXCLUDED_DIRS = {"__pycache__", "qcmodule.egg-info"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}


def tracked_files() -> list[Path]:
    files: list[Path] = []
    for top in ("scripts", "src"):
        base = TREE / top
        if not base.is_dir():
            sys.exit(f"error: {base} is missing; the upstream checkout is incomplete")
        for path in base.rglob("*"):
            if not path.is_file():
                continue
            if EXCLUDED_DIRS & set(path.relative_to(TREE).parts):
                continue
            if path.suffix in EXCLUDED_SUFFIXES:
                continue
            files.append(path)
    return sorted(files, key=lambda p: str(p.relative_to(TREE)))


def tree_digest() -> tuple[str, int] | None:
    """Digest of the tracked upstream source, or None if there is no checkout.

    `oracle/upstream-src` is a gitignored working copy, so a fresh CI checkout
    has no oracle tree. That is a skip, not a pass and not a failure: there is
    nothing to attest. Any host that DOES have the tree gets the real check,
    which is what makes the field trustworthy rather than decorative.
    """
    if not TREE.is_dir():
        return None
    files = tracked_files()
    outer = hashlib.sha256()
    for path in files:
        rel = path.relative_to(TREE).as_posix()
        inner = hashlib.sha256(path.read_bytes()).hexdigest()
        outer.update(f"{inner}  {rel}\n".encode())
    return outer.hexdigest(), len(files)


def recorded_digest() -> str | None:
    m = re.search(rf"^\s*{FIELD}:\s*([0-9a-f]{{64}})\s*$", LOCK.read_text(), re.M)
    return m.group(1) if m else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--update-lock", action="store_true",
                    help="write the computed digest into the lock file in place")
    args = ap.parse_args()

    computed = tree_digest()
    if computed is None:
        print(f"SKIP: no upstream checkout at {TREE} (gitignored working copy).\n"
              "Nothing to attest. Populate oracle/upstream-src at the pinned commit\n"
              "to enable this check.")
        return 0
    digest, n = computed
    recorded = recorded_digest()
    print(f"tracked files:  {n}")
    print(f"computed:       {digest}")
    print(f"recorded:       {recorded}")

    if args.update_lock:
        if recorded is None:
            sys.exit(f"error: no {FIELD} field in {LOCK}")
        LOCK.write_text(LOCK.read_text().replace(recorded, digest))
        print(f"updated {LOCK}")
        return 0

    if recorded is None:
        sys.exit(f"error: no {FIELD} field in {LOCK}")
    if digest != recorded:
        print("MISMATCH: the upstream source tree differs from the frozen lock.\n"
              "If this is an intended refresh, re-pin with --update-lock and say why\n"
              "in the lock's scope-change comment; every recorded differential and\n"
              "benchmark row is attributed to this digest.", file=sys.stderr)
        return 1
    print("OK: upstream source tree matches the frozen lock.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())