#!/usr/bin/env python3
"""Check release metadata for unresolved placeholders and unsupported claims.

A published archive carrying `https://github.com/TBD/rseqc-rust` is worse than no
archive: the URL resolves to somebody else's page, or to nothing, and a reader
cannot tell which. The audit found these placeholders in `Cargo.toml` and
`CITATION.cff`; this check makes them impossible to ship unnoticed.

The second job here is the one that matters more. The same archive's manifest
declares a Python API that does not exist and plot support that still needs
Rscript. Those are legitimate later milestones, but announcing them while the
Python crate is a stub is a false claim about software, and the audit called it
out specifically. So the check verifies that the manifest *says so*, rather than
trusting that someone remembered.

Exit codes: 0 clean, 1 problems found. `scripts/build-release-archive.sh` calls
this in `--strict` mode before publishing, where a placeholder is a hard failure;
without `--strict` the placeholders are reported as warnings, because a local
build is not a publication.

Usage:
    python3 scripts/check_release_metadata.py
    python3 scripts/check_release_metadata.py --strict
    python3 scripts/check_release_metadata.py --manifest dist/.../manifest.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Files whose public-facing metadata must not ship a placeholder.
# The Bioconda recipe is here too, and for the same reason it is here for the other
# two: its `about.home`, `about.doc_url` and `source.url` ship verbatim into a package
# index. A `TBD` in Cargo.toml and a `TBD` in a submitted recipe are the same defect --
# a published URL that resolves to somebody else's page -- so they are reported
# together instead of the recipe's going unnoticed until a maintainer reads it.
METADATA_FILES = ["Cargo.toml", "CITATION.cff", "recipes/rseqc-rust/meta.yaml"]

# What the archive must disclose, and the evidence that it does.
REQUIRED_MANIFEST_CLAIMS = {
    "python_api_unimplemented": (
        lambda m: "NOT IMPLEMENTED" in m["capabilities"].get("python_api", ""),
        "capabilities.python_api must state the Python API is not implemented",
    ),
    "known_limitations_present": (
        lambda m: len(m.get("known_limitations", [])) >= 3,
        "known_limitations must list what the release does not cover",
    ),
    "abi_floor_declared": (
        lambda m: bool(m["target"]["abi_floor"].get("glibc")),
        "target.abi_floor must name a minimum glibc version",
    ),
    "license_decision_recorded": (
        lambda m: any("DIV-0003" in x for x in m.get("known_limitations", [])),
        "known_limitations must record that the licence variant decision is open",
    ),
    "binaries_digested": (
        lambda m: all(len(v.get("sha256", "")) == 64
                      for v in m["binaries"].values()),
        "every binary must carry a full SHA256, not a prefix",
    ),
    "docs_digested": (
        lambda m: all(len(v) == 64 for v in m.get("docs", {}).values()),
        "every bundled document must carry a full SHA256",
    ),
    "provenance_bound": (
        lambda m: bool(m.get("built_from", {}).get("git_commit")),
        "built_from.git_commit must record the revision the binaries came from",
    ),
    "oracle_lock_bound": (
        lambda m: bool(m["reference_environment"].get("oracle_lock_sha256")),
        "reference_environment must bind the differential oracle's lock digest",
    ),
    # The per-command capability record. Checked for coverage, not just presence: a
    # manifest that carries a capability block naming 30 of 33 shipped commands would
    # otherwise read as though the other three are unqualified-but-present, which is
    # a different and much weaker statement.
    "per_command_capabilities_present": (
        lambda m: isinstance(m.get("per_command"), dict)
        and isinstance(m["per_command"].get("commands"), dict)
        and m["per_command"].get("commands") != {},
        "manifest must carry a per-command capability record, not just a summary "
        "capability block",
    ),
    "per_command_capabilities_cover_every_binary": (
        lambda m: set(m["per_command"]["commands"]) == set(m["binaries"]),
        "the per-command capability record must cover exactly the binaries the "
        "archive ships, so no shipped command is unrecorded",
    ),
    "per_command_capabilities_name_helper_use": (
        lambda m: all(isinstance(v.get("requires_rscript"), bool)
                      for v in m["per_command"]["commands"].values()),
        "each command's record must state whether it requires Rscript, rather than "
        "leaving it unstated",
    ),
    "per_command_capabilities_state_reproducibility": (
        lambda m: all("reproducible" in v
                      for v in m["per_command"]["commands"].values()),
        "each command's record must state whether its output is reproducible",
    ),
    # A recorded hash that no longer matches the shipped binary is worse than no
    # hash: it is a claim that the evidence belongs to these bytes. Coverage alone
    # cannot catch this -- all 33 records existed while all 33 hashes were stale
    # after a rebuild, and the checks above all passed.
    "per_command_capability_hashes_match_binaries": (
        lambda m: all(
            v.get("binary_sha256") == (m.get("binaries", {}).get(name) or {}).get("sha256")
            for name, v in m["per_command"]["commands"].items()
        ),
        "each command's recorded binary hash must equal the hash the manifest "
        "records for that binary, so a capability record cannot be carried over "
        "from different bytes",
    ),
    "per_command_capabilities_distinguish_rscript_from_other_helpers": (
        lambda m: all(
            # FPKM_UQ's helper is htseq-count, not Rscript. The generic probe
            # answers "needs a helper absent from this environment", and labelling
            # every such command as requiring Rscript turns a measured fact into a
            # false one for exactly the command that does not need R.
            not (v.get("requires_rscript") and v.get("missing_helper") == "htseq-count")
            for v in m["per_command"]["commands"].values()
        ),
        "a command whose absent helper is htseq-count must not be recorded as "
        "requiring Rscript",
    ),
}


def find_placeholders(root: Path | None = None) -> list[tuple[str, int, str]]:
    """Placeholder mentions in the metadata files, including explanatory comments."""
    root = root or REPO
    hits = []
    for name in METADATA_FILES:
        path = root / name
        if not path.exists():
            hits.append((name, 0, "file is missing"))
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if PLACEHOLDER_TOKEN.search(line):
                hits.append((name, i, line.strip()[:100]))
    return hits


# A placeholder in an *assigned value* is what actually ships, so that is what counts
# as a hard failure. A comment or prose block explaining that the placeholder is
# deliberate must not itself be reported, or every future edit that documents the
# placeholder re-fails the check and the check gets disabled.
#
# Both syntaxes are matched because the two metadata files disagree on them:
# `Cargo.toml` uses TOML (`repository = "..."`) and `CITATION.cff` uses YAML
# (`repository-code: "..."`). The earlier version of this pattern accepted only `=`,
# which is the same failure in miniature as the defect it exists to catch: it read
# as if it covered the published metadata while silently exempting CITATION.cff, the
# file a reader's citation tooling is most likely to read. `\s*[:=]\s*` covers both.
PLACEHOLDER_TOKEN = re.compile(
    r"\bTBD\b|example\.com|<your[- ]|\bYOUR[-_ ]?(?:USERNAME|ORG|ACCOUNT|REPO)\b",
    re.IGNORECASE,
)

VALUE_ASSIGNMENT = re.compile(
    r"^\s*-?\s*"
    r"(?:repository|repository-code|url|homepage|documentation|doi|repository_url)"
    r"\s*[:=]\s*"
    r"[\"']?(?P<value>[^\"'#]*?)[\"']?\s*$"
)


def find_placeholder_values(root: Path | None = None) -> list[tuple[str, int, str]]:
    """Placeholders in *assigned values* -- the ones that would actually ship."""
    root = root or REPO
    hits = []
    for name in METADATA_FILES:
        path = root / name
        if not path.exists():
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            m = VALUE_ASSIGNMENT.match(line)
            if m and PLACEHOLDER_TOKEN.search(m.group("value")):
                hits.append((name, i, m.group("value").strip()))
    return hits


def check_manifest(path: Path) -> list[str]:
    problems = []
    try:
        m = json.loads(path.read_text())
    except Exception as e:
        return [f"cannot read manifest {path}: {e}"]
    for key, (check, message) in REQUIRED_MANIFEST_CLAIMS.items():
        try:
            ok = check(m)
        except (KeyError, TypeError) as e:
            problems.append(f"{message} (manifest error: {e})")
            continue
        if not ok:
            problems.append(message)
    return problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict", action="store_true",
                    help="treat unresolved placeholders as failures (used when publishing)")
    ap.add_argument("--manifest", default=None,
                    help="path to a built manifest.json to check as well")
    ap.add_argument("--root", default=None,
                    help="directory holding Cargo.toml and CITATION.cff "
                         "(default: the source checkout). Pointing this at an "
                         "extracted release archive checks the metadata that "
                         "actually ships, not the tree it was built from.")
    args = ap.parse_args()
    root = Path(args.root).resolve() if args.root else REPO

    problems = []
    # A placeholder in an assigned VALUE is what ships, so it always counts.
    # A placeholder mentioned in a comment is only informational: this file and
    # Cargo.toml both explain that the value is deliberately undecided, and that
    # explanation must not be indistinguishable from the defect.
    values = find_placeholder_values(root)
    for name, line, value in values:
        msg = f"unresolved placeholder in {name}:{line}: {value}"
        if args.strict:
            problems.append(msg)
        else:
            print(f"  warning: {msg}")
    for name, line, text in find_placeholders(root):
        if not any(n == name and l == line for n, l, _ in values):
            print(f"  note: {name}:{line} mentions a placeholder in prose: {text}")

    if args.manifest:
        for problem in check_manifest(Path(args.manifest)):
            problems.append(problem)
    else:
        print("  note: no --manifest given; only the metadata files were checked")

    if problems:
        print("\nrelease metadata problems:")
        for p in problems:
            print(f"  {p}")
        print("\nThe permanent repository URL is a maintainer decision (release Stage A). "
              "Set it in\nCargo.toml and CITATION.cff before publishing. It is not guessed "
              "here on purpose:\na wrong URL in published metadata resolves to somebody "
              "else's page.")
        return 1

    print("release metadata: clean"
          + (f" (root {root})" if args.root else "")
          + (f" (manifest {args.manifest} checked)" if args.manifest else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
