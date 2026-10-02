#!/usr/bin/env python3
"""Check the live oracle environment against `compatibility/upstream.lock`.

The lock is only useful if something checks it. Before this, the pinned file
recorded upstream's source commit but not the environment that executed it, which
is the audit's P1 finding: a recorded differential result could not be attributed
to a specific reference environment, because the environment was never pinned and
so never compared.

Three outcomes, all recorded rather than silently tolerated:

  * MATCH   every checkable item agrees.
  * DRIFT   the environment differs from the lock. A drifted oracle can still be
            scientifically fine, but its numbers are not the ones the lock claims,
            so this exits non-zero and says what differs.
  * ABSENT  an external helper the lock marks absent has appeared, or one marked
            present has gone. Recorded because an absent wigToBigWig is a
            documented cause of two benchmark rows measuring a failure path.

Usage:
    oracle/venv/bin/python3 verification/check_oracle_env.py
    oracle/venv/bin/python3 verification/check_oracle_env.py --json
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
LOCK = REPO / "compatibility" / "upstream.lock"
VENV_PY = REPO / "oracle" / "venv" / "bin" / "python3"


def parse_lock(path: Path):
    """Minimal reader for this one file's shape.

    A YAML parser is not assumed to be available: this check must run in a bare
    environment, because an environment too broken to parse its own lock is
    exactly the situation worth reporting. The subset needed is
    `key:` / `  key: value` / `  - name: x` blocks.
    """
    section = {}
    current_map = None
    current_list = None
    for raw in path.read_text().splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip())
        line = raw.strip()
        if indent == 0:
            current_map, current_list = None, None
            key, _, value = line.partition(":")
            section[key] = value.strip() or {}
            continue
        if line.startswith("- "):
            item = {}
            key, _, value = line[2:].partition(":")
            item[key.strip()] = value.strip()
            current_list = current_list if current_list is not None else []
            current_list.append(item)
            section.setdefault("__list__", []).append(current_list)
            continue
        key, _, value = line.partition(":")
        value = value.strip()
        key = key.strip()
        if current_list is not None and indent >= 4:
            current_list[-1][key] = value
        else:
            current_map = current_map if current_map is not None else {}
            parent = next(reversed(section.get("__order__", []))) if False else None
            current_map[key] = value
    return section


def pip_freeze(python):
    try:
        out = subprocess.run([str(python), "-m", "pip", "freeze"],
                             capture_output=True, text=True, timeout=120).stdout
    except Exception as e:
        return None, f"pip freeze failed: {e}"
    pkgs = {}
    for line in out.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("-e "):
            # An editable install of the upstream tree. Its name and version come
            # from pyproject metadata, not from the pinned commit, so it is
            # recorded separately in the lock rather than treated as a drift.
            spec = line[3:].strip()
            name, _, version = spec.partition("@")
            if version and "==" in version:
                name, _, version = version.partition("==")
                pkgs[name.strip().lower()] = version.strip()
            continue
        name, _, version = line.partition("==")
        if version:
            pkgs[name.strip().lower()] = version.strip()
    return pkgs, None


def sh(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True,
                              timeout=60).stdout.strip()
    except Exception:
        return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--lock", default=str(LOCK))
    args = ap.parse_args()

    text = Path(args.lock).read_text()
    findings = []
    status = "MATCH"

    def drift(msg):
        nonlocal status
        findings.append(("DRIFT", msg))
        status = "DRIFT"

    def absent(msg):
        nonlocal status
        findings.append(("ABSENT", msg))
        status = "ABSENT"

    # --- Python interpreter -----------------------------------------------------
    want_python = re.search(r'version:\s*"([\d.]+)"', text)
    got_python = sh([sys.executable, "--version"])
    if want_python:
        pinned = want_python.group(1)
        if pinned not in got_python:
            drift(f"python: lock pins {pinned}, interpreter reports {got_python or 'unknown'}")
    else:
        drift("lock records no python version")

    # --- pip packages -----------------------------------------------------------
    freeze_block = text.split("pip_freeze:", 1)
    if len(freeze_block) == 2:
        block = freeze_block[1].split("importable_version_checks:")[0]
        pinned = {}
        for line in block.splitlines():
            m = re.match(r"\s*([A-Za-z0-9_.-]+):\s*\"?([^\"\s]+)\"?", line)
            if m and not line.strip().startswith("#"):
                pinned[m.group(1).lower()] = m.group(2)
        live, err = pip_freeze(VENV_PY)
        if err:
            drift(f"could not read the oracle environment: {err}")
        else:
            for name, version in sorted(pinned.items()):
                actual = live.get(name)
                if actual is None:
                    drift(f"{name}: lock pins {version}, package is absent")
                elif actual != version:
                    drift(f"{name}: lock pins {version}, installed {actual}")
            extra = sorted(set(live) - set(pinned))
            if extra:
                # An extra package is only a problem if the oracle can import it,
                # since an unused wheel cannot change behaviour. Recorded either way.
                findings.append(("DRIFT",
                                 f"present but not in the lock (check whether the oracle "
                                 f"imports them): {', '.join(extra)}"))
                status = "DRIFT" if status == "MATCH" else status
    else:
        drift("lock has no pip_freeze block")

    # --- external helper programs ----------------------------------------------
    for block in re.findall(r"- name: (\S+)(.*?)(?=\n  - name:|\n# ---|\Z)",
                            text, re.S):
        name, body = block
        m = re.search(r"status:\s*(\S+)", body)
        declared = m.group(1) if m else "unknown"
        path_m = re.search(r"path:\s*(\S+)", body)
        declared_path = path_m.group(1) if path_m and path_m.group(1) != "null" else None
        if declared == "present" and declared_path:
            if not Path(declared_path).exists():
                absent(f"{name}: lock records it present at {declared_path}, "
                       f"which does not exist")
        elif declared == "absent":
            found = sh(["which", name])
            if found:
                absent(f"{name}: lock records it absent, but {found} is on PATH")

    report = {"status": status, "findings": findings,
              "lock": str(LOCK), "python": got_python}
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"oracle environment check: {status}")
        print(f"  lock    : {LOCK}")
        print(f"  python  : {got_python}")
        for kind, msg in findings:
            print(f"  {kind}: {msg}")
        if not findings:
            print("  every checkable item matches the lock")
        if status != "MATCH":
            print("\nThe oracle environment differs from its lock. Differential and "
                  "benchmark\nnumbers produced here are not the pinned baseline; "
                  "re-pin before\nciting them, or restore the locked environment.")
    return 0 if status == "MATCH" else 1


if __name__ == "__main__":
    sys.exit(main())
