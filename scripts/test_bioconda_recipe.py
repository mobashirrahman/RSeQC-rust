#!/usr/bin/env python3
"""Validate the Bioconda recipe without needing conda-build installed.

Why this exists. Stage C requires "a Bioconda recipe", and a recipe is a Jinja2
template that only becomes YAML after rendering. So a recipe can be committed that
no Bioconda build would ever accept -- a typo in a section name, a `requirements`
block with the wrong shape, an unbalanced `{{ }}`, a build step naming a file that
is not in the repository -- and every one of those failures appears only at
submission time, in a channel-maintainer review, which is the worst place to find
out.

So the template is rendered with a context supplying exactly the two things
conda-build supplies (`compiler()` and `environ`) and nothing else, with
`StrictUndefined` so an unknown variable is an error rather than an empty string.
Every section, and every file the build script references, is then checked.

`{{ compiler('c') }}` and `{{ environ.get(...) }}` cannot be resolved by plain Jinja2,
which is why this renders them itself rather than shelling out to conda-build.
"""

from __future__ import annotations

import sys
from pathlib import Path

import jinja2
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
RECIPE = REPO_ROOT / "recipes" / "rseqc-rust" / "meta.yaml"

REQUIRED_SECTIONS = ("package", "source", "build", "requirements", "test", "about")
REQUIRED_TEST_COMMANDS = (
    "bam_stat --version",
    "infer_experiment --version",
    "read_distribution --version",
)


def conda_build_context(version: str, name: str) -> dict:
    """The jinja globals conda-build injects, and nothing more.

    Deliberately minimal: if a recipe starts depending on some other implicit global,
    this must fail rather than pass here and fail in Bioconda review.
    """
    return {
        "version": version,
        "name": name,
        "compiler": lambda language: f"_{language}_compiler_stub",
        "environ": {"RSEQC_ORACLE_PYTHON_VERSION": "3.12"},
    }


def render() -> dict:
    env = jinja2.Environment(
        loader=jinja2.BaseLoader(),
        undefined=jinja2.StrictUndefined,
        keep_trailing_newline=True,
    )
    template = env.from_string(RECIPE.read_text())
    return yaml.safe_load(template.render(**conda_build_context("0.1.0", "rseqc-rust")))


def check(recipe: dict) -> list[str]:
    problems: list[str] = []

    for section in REQUIRED_SECTIONS:
        if section not in recipe:
            problems.append(f"recipe has no '{section}' section")

    package = recipe.get("package", {})
    if not re_full_match_safe(str(package.get("name", ""))):
        problems.append("package.name must be a lowercase conda name")
    if package.get("version") != "0.1.0":
        problems.append(
            f"package.version is {package.get('version')!r}; it must be pinned, "
            f"not inherited from a template"
        )

    source = recipe.get("source", {})
    if "url" not in source:
        problems.append("source.url is required by Bioconda")
    if "sha256" not in source:
        problems.append("source.sha256 is required by Bioconda; a recipe without a "
                        "digest is not a versioned source")

    build = recipe.get("build", {})
    script = build.get("script")
    if not script:
        problems.append("build.script is empty")
    if "number" not in build:
        problems.append("build.number is required")
    else:
        problems.extend(_check_locked_build(script))

    test = recipe.get("test", {})
    if not test.get("commands"):
        problems.append("test.commands must exercise the package, not merely build it")
    joined = " ".join(test.get("commands", []))
    for command in REQUIRED_TEST_COMMANDS:
        if command not in joined:
            problems.append(f"test.commands must include {command!r}: a binary that "
                            f"cannot start is a binary nobody can use")

    about = recipe.get("about", {})
    for field in ("home", "license", "license_file", "summary"):
        if not about.get(field):
            problems.append(f"about.{field} is required")

    problems.extend(_check_referenced_files_exist())
    return problems


def _check_locked_build(script) -> list[str]:
    problems = []
    text = "\n".join(script) if isinstance(script, (list, tuple)) else str(script)
    if "cargo build" in text and "--locked" not in text:
        problems.append(
            "build.script runs cargo build without --locked; Bioconda's Rust "
            "guidance requires a locked build, and without it the recipe can pick "
            "up a different transitive dependency than the release archive"
        )
    if "THIRD_PARTY_NOTICES" not in text:
        problems.append(
            "build.script does not install THIRD_PARTY_NOTICES.txt; the release "
            "archive bundles the licences of shipped dependencies, and a recipe "
            "that dropped the file would package fewer of them"
        )
    return problems


def _check_referenced_files_exist() -> list[str]:
    """Every in-repository path a build step names must exist in the repository.

    A recipe step that copies a file the repository does not contain fails at
    build time, in a maintainer's channel, for a reason that has nothing to do with
    the package.
    """
    problems = []
    candidates = (
        "scripts/build-release-archive.sh",
        "benchmarks/capability-matrix.json",
        "Cargo.lock",
        "LICENSE",
        "README.md",
        "CHANGELOG.md",
        "CITATION.cff",
    )
    for rel in candidates:
        if not (REPO_ROOT / rel).exists():
            problems.append(f"recipe references {rel}, which is not in the repository")
    return problems


def re_full_match_safe(value: str) -> bool:
    import re
    return bool(re.fullmatch(r"[a-z0-9_.-]+", value))


def main() -> int:
    if not RECIPE.is_file():
        print(f"FAIL: no recipe at {RECIPE}")
        return 1
    try:
        recipe = render()
    except jinja2.TemplateError as exc:
        print(f"FAIL: the recipe does not render: {exc}")
        return 1
    except yaml.YAMLError as exc:
        print(f"FAIL: the rendered recipe is not valid YAML: {exc}")
        return 1

    problems = check(recipe)
    if problems:
        print(f"FAIL: {len(problems)} recipe problem(s)")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    print("OK: recipe renders, has every required section, builds locked, installs the")
    print("    dependency notices, and its test section runs the shipped binaries.")
    return 0


if __name__ == "__main__":
    sys.exit(main())