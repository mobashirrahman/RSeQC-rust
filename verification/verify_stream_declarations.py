#!/usr/bin/env python3
"""Check every `EXPECTED_STREAMS` declaration against what the command actually does.

Why this exists. The benchmark gate rejects a command whose declared streams
disagree with reality, and that check is only as good as the declarations. Two real
commands were gated on a rule that did not describe them, and in both cases the two
arms agreed byte-for-byte while the gate failed them:

  * `bam2wig` prints a one-line `wigToBigWig` command to stdout; the declaration said
    no stdout, so a byte-identical pair was rejected with "unexpected stdout content
    for a command that declares none".
  * `RNA_fragment_size`'s entire per-transcript table IS stdout; the declaration said
    no stdout, so 369,382 identical bytes were rejected.

The `RNA_fragment_size` declaration had been written from an observation, and the
observation was wrong in an instructive way: a probe looked for `label value` pairs,
found none in a table whose *header row* carries the column names, and concluded
there was no stdout. "No labelled stdout" was read as "no stdout".

So the declarations are not trusted here; they are re-derived from a live run. For
each command the script runs both arms once against a small workload and compares:

  * whether stdout was produced, against the `stdout` declaration;
  * whether stdout carries the declared labels, when the declaration uses labels;
  * whether stdout carries the declared substrings, when it uses `stdout_mode: text`;
  * which files the command produced, against the declared artifacts.

A mismatch is a failure. This is not a substitute for the equivalence gate -- it
checks the *declaration*, not the two implementations' agreement -- but a gate built
on a wrong declaration cannot report anything trustworthy about the commands it
gates.

Usage:
    verification/verify_stream_declarations.py [--workload DIR] [--output-dir DIR]
                                               [--commands NAME ...]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "bench_decl", REPO_ROOT / "benchmarks" / "bench.py")
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)

ORACLE_PY = bench.VENV_PY
UPSTREAM = bench.UPSTREAM_SRC
RELEASE = bench.RELEASE


def run(argv: list[str], cwd: Path, timeout: int = 900) -> tuple[int, str, str]:
    proc = subprocess.run(argv, cwd=cwd, capture_output=True, text=True,
                          timeout=timeout, env=dict(os.environ))
    return proc.returncode, proc.stdout, proc.stderr


def observed_artifacts(d: Path) -> set[str]:
    """Files the command produced, mapped back to the declared-artifact vocabulary.

    A declared artifact like `nvc.NVC.xls` is a suffix of a produced file name, since
    the run directory name is chosen by the harness. Matching on suffix keeps the
    check meaningful without hard-coding the harness's directory naming into it.
    """
    produced = {p.name for p in d.iterdir() if p.is_file()}
    log_dir = d / "logs"
    if log_dir.is_dir():
        produced |= {p.name for p in log_dir.iterdir() if p.is_file()}
    return produced


def check_command(name: str, workload: Path, work: Path) -> tuple[list[str], list[str]]:
    script, binary, argfn, _eclass = bench.COMMANDS[name]
    spec = bench.EXPECTED_STREAMS[name]
    problems: list[str] = []
    unverified: list[str] = []

    results = {}
    argvs = {}
    for arm in ("py", "rs"):
        d = work / f"{name}_{arm}"
        d.mkdir(parents=True, exist_ok=True)
        # The output directory handed to the command IS this arm's own directory, so
        # the artifacts can be found by looking where they were written. Sharing one
        # output directory across arms would also make each arm's files look like the
        # other's, which is precisely the mistake that makes an artifact check
        # meaningless.
        ctx = bench.Ctx(workload, d)
        if arm == "py":
            argvs[arm] = [str(ORACLE_PY), str(UPSTREAM / script)] + [
                str(a) if isinstance(a, Path) else a for a in argfn(ctx)]
        else:
            argvs[arm] = [str(RELEASE / binary)] + [
                str(a) if isinstance(a, Path) else a for a in argfn(ctx)]
        code, out, err = run(argvs[arm], d)
        results[arm] = {"exit": code, "stdout": out, "stderr": err, "dir": d}

    for arm in ("py", "rs"):
        r = results[arm]
        if r["exit"] != 0:
            # UNVERIFIED, not wrong. A command whose input shape this workload does
            # not provide cannot be observed, and saying "the declaration is wrong"
            # would be a different and unsupported claim. It is still worth printing,
            # because a declaration nobody has ever checked is exactly the failure
            # this script exists to prevent -- three commands here were gated on
            # declarations that were wrong, and all three had been "observed" by a
            # probe that could not actually see their stdout.
            unverified.append(
                f"{name}/{arm}: exited {r['exit']} on this workload, so its "
                f"declaration is UNVERIFIED here: "
                f"{(r['stderr'] or r['stdout']).strip()[:150]}")
            continue

        declared_stdout = spec.get("stdout")
        actual_stdout = r["stdout"].strip()
        if declared_stdout is False and actual_stdout:
            problems.append(
                f"{name}/{arm}: declared no stdout but wrote "
                f"{len(r['stdout'])} characters -- first line "
                f"{r['stdout'].splitlines()[0][:90]!r}")
        if declared_stdout is True:
            if not actual_stdout:
                problems.append(
                    f"{name}/{arm}: declares stdout but wrote none")
            elif spec.get("stdout_mode") == "text":
                for needle in spec.get("stdout_contains", ()):
                    if needle not in r["stdout"]:
                        problems.append(
                            f"{name}/{arm}: declared stdout_contains {needle!r} "
                            f"but it is absent from the stream")
            else:
                labels = bench.extract_stdout_metrics(bench.normalise_stdout(r["stdout"], None))
                for label in sorted(spec.get("labels") or ()):
                    if label not in labels:
                        problems.append(
                            f"{name}/{arm}: declares label {label!r} but the stream "
                            f"carries {sorted(labels)[:6]}")

        for artifact in spec.get("artifacts") or ():
            produced = observed_artifacts(r["dir"])
            if not any(name_of == artifact or name_of.endswith("." + artifact)
                       for name_of in produced):
                problems.append(
                    f"{name}/{arm}: declares artifact {artifact!r} but the run "
                    f"directory holds {sorted(produced)[:6]}")

    return problems, unverified


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workload", type=Path,
                    default=Path("/tmp/opencode/wl-rat"),
                    help="a workload directory the harness can drive")
    ap.add_argument("--commands", nargs="*", default=None)
    # `action="append"` is load-bearing. With the default store action and nargs="*",
    # each occurrence REPLACES the previous list, so
    #   --workload-for a=X --workload-for b=Y --workload-for c=Z
    # silently keeps only c. Every override but the last is dropped, and the run then
    # reports those commands as "UNVERIFIED because the workload does not supply their
    # input shape" -- which reads as a property of the commands rather than a bug in the
    # invocation. It was exactly that, and it is why the sweep sat at 24/29 while every
    # individual override worked when passed alone.
    ap.add_argument("--workload-for", action="append", nargs="*", default=[],
                    metavar="CMD=DIR",
                    help="per-command workload directory, as the harness takes; "
                         "repeatable, and several CMD=DIR values may follow one flag")
    args = ap.parse_args()

    if not args.workload.is_dir():
        print(f"skip: no workload at {args.workload}")
        return 0

    names = args.commands or [
        n for n in bench.COMMANDS if n not in bench.EXCLUDED]
    overrides = {}
    for group in args.workload_for:
        for spec_str in group:
            cmd, _, path = spec_str.partition("=")
            if not cmd or not path:
                print(f"skip: malformed --workload-for {spec_str!r}; expected CMD=DIR")
                continue
            overrides[cmd] = Path(path).resolve()
    if overrides:
        print(f"per-command workloads: {', '.join(sorted(overrides))}")

    all_problems: list[str] = []
    all_unverified: list[str] = []
    checked = 0
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for name in names:
            problems, unverified = check_command(
                name, overrides.get(name, args.workload), work)
            all_problems.extend(problems)
            all_unverified.extend(unverified)
            checked += 0 if unverified else 1
            if problems:
                status = "BAD"
            elif unverified:
                status = "unv"
            else:
                status = "ok "
            print(f"  [{status}] {name}")
            for problem in problems:
                print(f"          {problem}")
            for note in unverified:
                print(f"          {note}")

    print()
    if all_problems:
        print(f"FAIL: {len(all_problems)} declaration/reality mismatch(es)")
        print("A gate built on a wrong declaration cannot report anything trustworthy")
        print("about the commands it gates. Fix EXPECTED_STREAMS, or the command.")
        return 1
    print(f"OK: {checked} of {len(names)} declarations verified against a live run")
    if all_unverified:
        print(f"     {len(all_unverified)} arm(s) unverified because the workload "
              f"does not supply their input shape.")
        print("     Pass --workload-for CMD=DIR to cover those too.")
    return 0


if __name__ == "__main__":
    sys.exit(main())