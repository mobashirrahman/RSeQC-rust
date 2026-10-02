#!/usr/bin/env python3
"""Derive per-command expected streams and artifacts by running both arms once.

The audit's P0 benchmark-gate finding was that the gate compared only two file
trees: stdout was never compared, and nothing declared what a command was
supposed to produce, so two empty directories passed and a stdout-only metric
received a gate pass with no metric checking. Writing those expectations by hand
from reading the sources is exactly the kind of guess the audit is about.

This script instead runs each command once in each arm and records what actually
appears: which files, which stdout labels, which exit statuses. Its output is a
draft for `EXPECTED_STREAMS` that a human then reviews -- it observes behaviour,
it does not certify correctness, so nothing it emits is trusted as an oracle.

Usage:
    oracle/venv/bin/python3 benchmarks/derive_expected_streams.py \
        --workload <workload dir> --output benchmarks/expected-streams.observed.json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location(
    "bench", REPO / "benchmarks" / "bench.py")
assert SPEC and SPEC.loader
bench = importlib.util.module_from_spec(SPEC)
sys.modules["bench"] = bench
SPEC.loader.exec_module(bench)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workload", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--commands", nargs="*", default=sorted(bench.COMMANDS))
    ap.add_argument("--timeout", type=int, default=300)
    args = ap.parse_args()

    workload = Path(args.workload).resolve()
    env = dict(os.environ)
    env.update(bench.PINNED_ENV)

    observed = {}
    for name in args.commands:
        if name not in bench.COMMANDS:
            print(f"!! unknown command {name}")
            continue
        script, binary, argfn, eclass = bench.COMMANDS[name]
        entry = {"upstream_script": script, "rust_binary": binary,
                 "experiment_class": eclass}
        for arm in ("py", "rs"):
            tmp = Path(tempfile.mkdtemp(prefix=f"derive-{name}-{arm}-"))
            try:
                argv = ([str(bench.VENV_PY), str(bench.UPSTREAM_SRC / script)]
                        if arm == "py" else [str(bench.RELEASE / binary)])
                argv += [str(a) if isinstance(a, Path) else a
                         for a in argfn(bench.Ctx(workload, tmp))]
                m, so, se = bench.run_arm(argv, tmp, env, args.timeout)
                files = sorted(
                    str(p.relative_to(tmp)) for p in tmp.rglob("*")
                    if p.is_file() and p.name not in ("py.code", "rs.code", "time.txt"))
                metrics = bench.extract_stdout_metrics(so)
                entry[arm] = {
                    "exit_code": m["exit_code"],
                    "timed_out": m["timed_out"],
                    "files": [{"path": f, "bytes": (tmp / f).stat().st_size}
                              for f in files],
                    "stdout_labels": sorted(metrics),
                    "stdout_line_count": len(so.splitlines()),
                    "stderr_tail": se[-300:],
                }
            finally:
                shutil.rmtree(tmp, ignore_errors=True)
        both = set(f["path"] for f in entry.get("py", {}).get("files", []))
        only_rs = set(f["path"] for f in entry.get("rs", {}).get("files", [])) - both
        only_py = both - set(f["path"] for f in entry.get("rs", {}).get("files", []))
        entry["asymmetry"] = {
            "only_upstream": sorted(only_py),
            "only_rust": sorted(only_rs),
            "stdout_labels_only_upstream":
                sorted(set(entry.get("py", {}).get("stdout_labels", []))
                       - set(entry.get("rs", {}).get("stdout_labels", []))),
            "stdout_labels_only_rust":
                sorted(set(entry.get("rs", {}).get("stdout_labels", []))
                       - set(entry.get("py", {}).get("stdout_labels", []))),
        }
        observed[name] = entry
        print(f"{name}: py exit {entry['py']['exit_code']} "
              f"({len(entry['py']['files'])} files, "
              f"{len(entry['py']['stdout_labels'])} stdout labels) | "
              f"rs exit {entry['rs']['exit_code']} "
              f"({len(entry['rs']['files'])} files, "
              f"{len(entry['rs']['stdout_labels'])} stdout labels)")

    Path(args.output).write_text(json.dumps(observed, indent=2))
    print(f"\nobserved -> {args.output}")
    print("Review before using as EXPECTED_STREAMS: this records behaviour, not "
          "correctness.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
