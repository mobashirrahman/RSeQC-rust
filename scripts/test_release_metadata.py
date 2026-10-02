#!/usr/bin/env python3
"""Executable tests for `scripts/check_release_metadata.py`.

Why this file exists: the check's whole job is to fail when published metadata still
carries a placeholder, so a version of it that does not fail is worse than no check
at all. The first version of `VALUE_ASSIGNMENT` matched only TOML's `key = "value"`
and not YAML's `key: "value"`, so it read as though it covered `CITATION.cff` -- the
file citation tooling actually reads -- while exempting that file completely. It
still reported a problem, because `Cargo.toml` uses `=`, so the bug was invisible in
the one test a person would write by hand (run the check, see it fail, stop).

Every test below therefore asserts on a *defective* input, and the reference-POSITIVE
test is the one that caught it: with `CITATION.cff` deliberately holding the only
placeholder, the check must still report it.

Run: `python3 scripts/test_release_metadata.py`
"""
from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import check_release_metadata as crm  # noqa: E402

SCRIPT = Path(__file__).resolve().parent / "check_release_metadata.py"

CLEAN_CARGO = """\
[package]
name = "rseqc-rust"
version = "0.1.0"
repository = "https://github.com/example-org/rseqc-rust"
homepage = "https://example-org.github.io/rseqc-rust"
"""

CLEAN_CFF = """\
cff-version: 1.2.0
message: cite this
title: RSeQC-rust
type: software
version: 0.1.0
repository-code: "https://github.com/example-org/rseqc-rust"
url: "https://github.com/example-org/rseqc-rust"
license: GPL-3.0-or-later
"""

CLEAN_RECIPE = """\
{% set name = "rseqc-rust" %}
{% set version = "0.1.0" %}
package:
  name: {{ name|lower }}
  version: {{ version }}
source:
  url: https://github.com/example-org/rseqc-rust/archive/v{{ version }}.tar.gz
  sha256: 0000000000000000000000000000000000000000000000000000000000000000
about:
  home: https://github.com/example-org/rseqc-rust
"""

# A recipe whose placeholders are UNRESOLVED. The same defect as in Cargo.toml and
# CITATION.cff: `about.home` and `source.url` ship verbatim into a package index, and
# `sha256: TBD` is not a digest. Including it here proves the recipe is checked rather
# than present-but-ignored, which is what it was before this file was added.
REAL_RECIPE = """\
{% set name = "rseqc-rust" %}
{% set version = "0.1.0" }
package:
  name: {{ name|lower }}
  version: {{ version }}
source:
  url: https://github.com/TBD/rseqc-rust/archive/v{{ version }}.tar.gz
  sha256: TBD
about:
  home: https://github.com/TBD/rseqc-rust
"""

# The exact shape currently in the tree, including the explanatory comments that must
# NOT be reported as shipped values.
REAL_CARGO = """\
# UNRESOLVED PLACEHOLDER. The permanent repository URL is a maintainer decision
# (release Stage A) and is deliberately left as TBD rather than guessed: a wrong
# URL in published metadata resolves to somebody else's page. The check at
# scripts/check_release_metadata.py --strict
# fails a release build while this is still TBD, so the placeholder cannot ship
# unnoticed.
repository = "https://github.com/TBD/rseqc-rust"
"""

REAL_CFF = """\
cff-version: 1.2.0
title: RSeQC-rust
# NOTE: this repository has not yet been published under a public URL -- see
# compatibility/divergences.yaml's DIV-0003 for why this project's own release license (and, by
# extension, a permanent citable URL) is still pending a decision. Update `repository-code` and
# `url` below once that happens; Cargo.toml's own `repository` field uses the same placeholder.
repository-code: "https://github.com/TBD/rseqc-rust"
url: "https://github.com/TBD/rseqc-rust"
license: GPL-3.0-or-later
"""


def _line_of(text: str, needle: str) -> int:
    """1-based line number of `needle` in `text` (asserted present)."""
    for i, line in enumerate(text.splitlines(), 1):
        if needle in line:
            return i
    raise AssertionError(f"{needle!r} not found in fixture")


class _Repo:
    """Point crm.REPO at a temporary directory for the duration of a test."""

    def __init__(self, cargo: str, cff: str, recipe: str = CLEAN_RECIPE) -> None:
        self.cargo = cargo
        self.cff = cff
        self.recipe = recipe
        self._orig = crm.REPO

    def __enter__(self) -> "_Repo":
        self._tmp = tempfile.TemporaryDirectory()
        root = Path(self._tmp.name)
        (root / "Cargo.toml").write_text(self.cargo)
        (root / "CITATION.cff").write_text(self.cff)
        # The recipe is written too, and this matters: a metadata file listed in
        # METADATA_FILES but absent from the temporary repository is reported as
        # "file is missing", so every clean-metadata test would fail for a reason
        # that has nothing to do with what it checks.
        recipe_path = root / "recipes" / "rseqc-rust" / "meta.yaml"
        recipe_path.parent.mkdir(parents=True, exist_ok=True)
        recipe_path.write_text(self.recipe)
        crm.REPO = root
        return self

    def __exit__(self, *exc: object) -> None:
        crm.REPO = self._orig
        self._tmp.cleanup()


class PlaceholderValues(unittest.TestCase):
    """A placeholder in an assigned value is what ships, so it must be found."""

    def test_yaml_cff_placeholder_is_found(self) -> None:
        """The regression: CFF's `key: "value"` must be matched, not just `key = "value"`."""
        with _Repo(CLEAN_CARGO, REAL_CFF):
            hits = crm.find_placeholder_values()
        self.assertEqual(
            sorted((n, l) for n, l, _ in hits),
            [("CITATION.cff", _line_of(REAL_CFF, "repository-code:")),
             ("CITATION.cff", _line_of(REAL_CFF, 'url: "'))],
            "both CFF placeholder values must be reported",
        )

    def test_toml_placeholder_is_found(self) -> None:
        with _Repo(REAL_CARGO, CLEAN_CFF):
            hits = crm.find_placeholder_values()
        self.assertEqual([(n, l) for n, l, _ in hits],
                         [("Cargo.toml", _line_of(REAL_CARGO, "repository ="))])

    def test_real_tree_placeholders_found_in_both_files(self) -> None:
        """The actual checked-in metadata must be detected in both files."""
        with _Repo(REAL_CARGO, REAL_CFF):
            hits = crm.find_placeholder_values()
        files = {n for n, _, _ in hits}
        self.assertEqual(files, {"Cargo.toml", "CITATION.cff"})
        self.assertEqual(len(hits), 3)

    def test_only_cff_placeholder_still_fails(self) -> None:
        """The case the hand-written test misses: CFF alone must be enough to fail."""
        with _Repo(CLEAN_CARGO, REAL_CFF):
            values = crm.find_placeholder_values()
            prose = crm.find_placeholders()
        self.assertTrue(values, "a CFF-only placeholder must still be reported")
        # Every prose line here is also an assigned value, so strict mode fails.
        self.assertEqual(
            sorted(l for _, l, _ in values), sorted(l for _, l, _ in prose)
        )

    def test_only_recipe_placeholder_still_fails(self) -> None:
        """The Bioconda recipe alone must be enough to fail.

        `about.home` and `source.url` are copied verbatim into a package index, so a
        `TBD` there is the same published-wrong-URL defect as in Cargo.toml. Before
        this file was in METADATA_FILES the recipe's placeholders were invisible to
        the release gate, which meant the very place most likely to publish a URL
        nobody had decided on was the one place nobody checked.
        """
        with _Repo(CLEAN_CARGO, CLEAN_CFF, REAL_RECIPE):
            values = crm.find_placeholder_values()
            prose = crm.find_placeholders()
        recipe_values = {n for n, _, _ in values if n.endswith("meta.yaml")}
        self.assertTrue(recipe_values,
                        "a recipe-only placeholder must still be reported")
        self.assertEqual(len(recipe_values), 1)
        # And the prose report must name the recipe too, so it is visible in the
        # non-strict summary rather than only in strict mode.
        self.assertTrue({n for n, _, _ in prose if n.endswith("meta.yaml")})

    def test_clean_metadata_is_clean(self) -> None:
        with _Repo(CLEAN_CARGO, CLEAN_CFF):
            self.assertEqual(crm.find_placeholder_values(), [])
            self.assertEqual(crm.find_placeholders(), [])

    def test_comment_prose_is_not_a_shipped_value(self) -> None:
        """Documenting the placeholder must not itself re-fail the check."""
        with _Repo(REAL_CARGO, REAL_CFF):
            prose = {(n, l) for n, l, _ in crm.find_placeholders()}
            values = {(n, l) for n, l, _ in crm.find_placeholder_values()}
        explained = ("Cargo.toml", _line_of(REAL_CARGO, "deliberately left as TBD"))
        self.assertIn(explained, prose, "the comment must still be reported as prose")
        self.assertNotIn(explained, values, "a comment is not a shipped value")
        self.assertLess(len(values), len(prose))
        # CFF's own explanatory comment mentions the decision but not the token, so it
        # is not a placeholder hit at all; the CFF *values* are, which is the point.
        self.assertEqual({n for n, _ in values}, {"Cargo.toml", "CITATION.cff"})

    def test_upstream_reference_url_is_not_flagged(self) -> None:
        """RSeQC's own real URL must pass: over-matching is also a broken check."""
        cff = CLEAN_CFF + (
            "references:\n"
            '  - type: software\n'
            "    title: RSeQC\n"
            '    repository-code: "https://github.com/liguowang/RSeQC"\n'
        )
        with _Repo(CLEAN_CARGO, cff):
            self.assertEqual(crm.find_placeholder_values(), [])

    def test_other_placeholder_tokens(self) -> None:
        for bad in ("https://example.com/o/r", "https://github.com/<your-org>/r",
                    "https://github.com/YOUR_ORG/r"):
            with self.subTest(bad=bad):
                cff = CLEAN_CFF.replace("example-org/rseqc-rust", bad)
                with _Repo(CLEAN_CARGO, cff):
                    self.assertTrue(crm.find_placeholder_values())


def _good_manifest() -> dict:
    """A manifest that passes every claim, including the per-command record.

    Kept in one place and reused by the per-command tests as a base, so a new claim
    does not need to be added to several fixtures -- which is how a check ends up
    asserted in one test and satisfied by no fixture in another.
    """
    m = {
        "capabilities": {"python_api": "NOT IMPLEMENTED"},
        "known_limitations": [
            "DIV-0003: release licence variant decision is open",
            "Linux x86_64 only",
            "CRAM is buffered whole-file",
        ],
        "target": {"abi_floor": {"glibc": "2.28"}},
        "binaries": {"bam_stat": {"sha256": "a" * 64}},
        "docs": {"LICENSE": "b" * 64, "README.md": "c" * 64},
        "built_from": {"git_commit": "deadbeef"},
        "reference_environment": {"oracle_lock_sha256": "d" * 64},
        "per_command": {
            "source": "benchmarks/capability-matrix.json",
            "commands": {
                # binary_sha256 is part of a real record and must equal the hash the
                # manifest records for the same binary. All 33 records existed while all
                # 33 hashes were stale after a rebuild, and every coverage check passed,
                # so a fixture without the field would have let that through.
                "bam_stat": {"requires_rscript": False, "reproducible": True,
                             "binary_sha256": "a" * 64},
            },
        },
    }
    return m


class ManifestClaims(unittest.TestCase):
    """The manifest must disclose the things the release does not do."""

    def _check(self, m: dict) -> list[str]:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "manifest.json"
            p.write_text(json.dumps(m))
            return crm.check_manifest(p)

    def test_good_manifest_passes(self) -> None:
        self.assertEqual(self._check(_good_manifest()), [])

    def test_python_api_claim_must_admit_stub(self) -> None:
        m = _good_manifest()
        m["capabilities"]["python_api"] = "supported"
        self.assertTrue(self._check(m))

    def test_short_sha_is_rejected(self) -> None:
        m = _good_manifest()
        m["binaries"]["bam_stat"]["sha256"] = "a" * 12
        self.assertTrue(self._check(m))

    def test_short_doc_sha_is_rejected(self) -> None:
        m = _good_manifest()
        m["docs"]["LICENSE"] = "b" * 8
        self.assertTrue(self._check(m))

    def test_abi_floor_must_be_declared(self) -> None:
        m = _good_manifest()
        m["target"]["abi_floor"] = {}
        self.assertTrue(self._check(m))

    def test_licence_decision_must_be_recorded(self) -> None:
        m = _good_manifest()
        m["known_limitations"] = ["a", "b", "c"]
        self.assertTrue(self._check(m))

    def test_provenance_must_be_bound(self) -> None:
        m = _good_manifest()
        m["built_from"] = {"git_commit": ""}
        self.assertTrue(self._check(m))

    def test_oracle_lock_must_be_bound(self) -> None:
        m = _good_manifest()
        m["reference_environment"] = {"oracle_lock_sha256": ""}
        self.assertTrue(self._check(m))

    def test_known_limitations_must_exist(self) -> None:
        m = _good_manifest()
        m["known_limitations"] = []
        self.assertTrue(self._check(m))

    def test_unreadable_manifest_reports_error(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "manifest.json"
            p.write_text("{not json")
            self.assertTrue(crm.check_manifest(p))


class PerCommandCapabilities(unittest.TestCase):
    """The manifest's per-command capability record must be present and complete.

    The audit's Stage A asks for a capability manifest recording, per command, the
    inputs it accepts, which outputs need a helper, whether it is stochastic, and which
    upstream quirks it reproduces. The manifest previously carried a *summary*
    capability block only ("16 commands require Rscript"), and nothing checked that the
    summary matched reality or covered the commands actually shipped -- so a manifest
    naming 30 of 33 commands would have read as a complete record.
    """

    def _good(self):
        """The shared valid manifest, with a second binary so coverage is meaningful.

        Derived from `_good_manifest()` rather than written out again: a hand-built
        copy is a second place to forget a claim, and the failure mode is a test that
        passes for reasons unrelated to what it names.
        """
        m = copy.deepcopy(_good_manifest())
        m["binaries"]["tin"] = {"sha256": "b" * 64}
        m["per_command"]["commands"]["tin"] = {
            "requires_rscript": False, "reproducible": True,
            "binary_sha256": "b" * 64}
        return m

    def test_good_manifest_passes(self):
        self.assertEqual(self._check(self._good()), [])

    def test_stale_capability_hash_is_rejected(self):
        """A record carried over from different bytes is a false claim.

        Coverage checks cannot see this: the record is present, names the right
        command, and states the right capabilities. Only comparing the hash to the
        binary the archive actually ships notices.
        """
        m = self._good()
        m["per_command"]["commands"]["bam_stat"]["binary_sha256"] = "0" * 64
        problems = self._check(m)
        self.assertTrue(any("recorded binary hash must equal" in p for p in problems),
                        problems)

    def test_absent_capability_hash_is_rejected(self):
        m = self._good()
        m["per_command"]["commands"]["bam_stat"].pop("binary_sha256")
        problems = self._check(m)
        self.assertTrue(any("recorded binary hash must equal" in p for p in problems),
                        problems)

    def test_htseq_helper_must_not_be_labelled_rscript(self):
        """FPKM-UQ's helper is htseq-count, not Rscript.

        16 commands need some external helper, and deriving `requires_rscript` from
        that single generic flag labels this one as Rscript-dependent -- turning a
        measured fact into a false one for exactly the command that does not need R.
        """
        m = self._good()
        m["binaries"]["FPKM_UQ"] = {"sha256": "e" * 64}
        m["per_command"]["commands"]["FPKM_UQ"] = {
            "requires_rscript": True, "reproducible": False,
            "binary_sha256": "e" * 64, "missing_helper": "htseq-count"}
        problems = self._check(m)
        self.assertTrue(any("htseq-count" in p for p in problems), problems)

        m["per_command"]["commands"]["FPKM_UQ"]["requires_rscript"] = False
        self.assertEqual(self._check(m), [])

    def test_per_command_block_required(self):
        m = self._good()
        m.pop("per_command")
        self.assertTrue(self._check(m))

    def test_every_shipped_binary_must_have_a_record(self):
        m = self._good()
        m["per_command"]["commands"].pop("tin")
        problems = self._check(m)
        self.assertTrue(any("cover exactly the binaries" in p for p in problems),
                        problems)

    def test_extra_records_are_also_a_mismatch(self):
        m = self._good()
        m["per_command"]["commands"]["ghost"] = {
            "requires_rscript": False, "reproducible": True}
        self.assertTrue(any("cover exactly the binaries" in p for p in self._check(m)))

    def test_helper_use_must_be_stated(self):
        m = self._good()
        for v in m["per_command"]["commands"].values():
            v.pop("requires_rscript")
        self.assertTrue(any("requires Rscript" in p for p in self._check(m)))

    def test_reproducibility_must_be_stated(self):
        m = self._good()
        for v in m["per_command"]["commands"].values():
            v.pop("reproducible")
        self.assertTrue(any("reproducible" in p for p in self._check(m)))

    def test_empty_per_command_block_is_not_a_record(self):
        m = self._good()
        m["per_command"]["commands"] = {}
        self.assertTrue(self._check(m))

    def _check(self, m):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "manifest.json"
            p.write_text(json.dumps(m))
            return crm.check_manifest(p)


class EndToEnd(unittest.TestCase):
    """The published CLI contract, exercised as the release script calls it."""

    def _run(self, cargo: str, cff: str, *args: str) -> subprocess.CompletedProcess:
        with _Repo(cargo, cff) as repo:
            return subprocess.run(
                [sys.executable, str(SCRIPT), "--root", str(repo._tmp.name), *args],
                capture_output=True, text=True,
            )

    def test_strict_exits_1_on_real_tree_placeholders(self) -> None:
        r = self._run(REAL_CARGO, REAL_CFF, "--strict")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertIn("Cargo.toml", r.stdout)
        self.assertIn("CITATION.cff", r.stdout)

    def test_non_strict_exits_0_but_warns(self) -> None:
        r = self._run(REAL_CARGO, REAL_CFF)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("warning: unresolved placeholder", r.stdout)

    def test_strict_exits_0_on_clean_tree(self) -> None:
        r = self._run(CLEAN_CARGO, CLEAN_CFF, "--strict")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("clean", r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
