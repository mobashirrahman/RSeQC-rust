#!/usr/bin/env python3
"""Synthetic sweep: comprehensive byte-for-byte comparison of RSeQC commands.

Runs upstream Python and Rust binary on various synthetic workloads
across multiple commands and flag combinations, reporting mismatches.

Usage:
    oracle/venv/bin/python3 verification/synthetic_sweep.py

Output: summary of mismatches found, if any.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
ORACLE_PYTHON = REPO_ROOT / "oracle" / "venv" / "bin" / "python3"
RUST_BIN_DIR = REPO_ROOT / "target" / "release"
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"
ORACLE_PYTHONPATH = str(REPO_ROOT / "oracle" / "upstream-src" / "src")
ORACLE_SCRIPTS = REPO_ROOT / "oracle" / "upstream-src" / "scripts"

# Workload configurations to test: (size, seed)
WORKLOADS = [
    (200, 42),
    (200, 99),
]

# Commands to test with flag variants
COMMAND_TESTS = [
    ("bam_stat", [
        (["-i", "reads.bam", "-q", "30"],
         ["-i", "reads.bam", "-q", "30"]),
        (["-i", "reads.bam", "-q", "0"],
         ["-i", "reads.bam", "-q", "0"]),
    ]),
    ("read_GC", [
        (["-i", "reads.bam", "-o", "out", "-q", "30", "--skip-plot"],
         ["-i", "reads.bam", "-o", "out", "-q", "30", "--skip-plot"]),
    ]),
    ("read_NVC", [
        (["-i", "reads.bam", "-o", "out", "-q", "30", "--skip-plot"],
         ["-i", "reads.bam", "-o", "out", "-q", "30", "--skip-plot"]),
    ]),
    ("read_quality", [
        (["-i", "reads.bam", "-o", "out", "-q", "30", "-r", "1", "--skip-plot"],
         ["-i", "reads.bam", "-o", "out", "-q", "30", "-r", "1", "--skip-plot"]),
    ]),
    ("read_duplication", [
        (["-i", "reads.bam", "-o", "out", "-q", "30", "--skip-plot"],
         ["-i", "reads.bam", "-o", "out", "-q", "30", "--skip-plot"]),
    ]),
    ("bam2fq", [
        (["-i", "reads.bam", "-o", "out"],
         ["-i", "reads.bam", "-o", "out"]),
    ]),
    ("split_paired_bam", [
        (["-i", "reads.bam", "-o", "out"],
         ["-i", "reads.bam", "-o", "out"]),
    ]),
    ("read_distribution", [
        (["-i", "reads.bam", "-r", "model.bed12"],
         ["-i", "reads.bam", "-r", "model.bed12"]),
    ]),
    ("infer_experiment", [
        (["-i", "reads.bam", "-r", "model.bed12", "-q", "30"],
         ["-i", "reads.bam", "-r", "model.bed12", "-q", "30"]),
    ]),
    ("read_hexamer", [
        (["-i", "reads_1.fastq"],
         ["-i", "reads_1.fastq"]),
    ]),
]


def run_command(argv: list[str], *, cwd: Path, timeout_s: float = 60.0,
               env: dict | None = None) -> tuple[int | None, str, str]:
    """Run command and return (exit_code, stdout, stderr)."""
    try:
        run_env = dict(os.environ)
        if env:
            run_env.update(env)
        proc = subprocess.Popen(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=str(cwd),
            env=run_env,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
            return proc.returncode, stdout, stderr
        except subprocess.TimeoutExpired:
            proc.kill()
            return None, "", f"[timeout after {timeout_s}s]"
    except Exception as e:
        return None, "", f"[error running command: {e}]"


def generate_workload(size: int, seed: int, output_dir: Path) -> bool:
    """Generate a synthetic workload. Returns True on success."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(ORACLE_PYTHON),
        str(BENCHMARKS_DIR / "generate_workload.py"),
        "--size", str(size),
        "--seed", str(seed),
        "--output-dir", str(output_dir),
    ]
    code, stdout, stderr = run_command(cmd, cwd=REPO_ROOT)
    if code != 0:
        print(f"WARNING: Failed to generate workload size={size} seed={seed}", file=sys.stderr)
        return False
    return True


def test_command(command: str, py_flags: list[str], rust_flags: list[str],
                workload_dir: Path) -> dict:
    """Test a single command. Returns result dict."""
    result = {
        "command": command,
        "success": True,
        "mismatches": [],
    }

    # Separate output dirs per side: sharing one lets the first run's files
    # trip the second run's overwrite guard and makes file comparison moot.
    py_outdir = workload_dir / f"out_py_{command}"
    rust_outdir = workload_dir / f"out_rust_{command}"
    for d in (py_outdir, rust_outdir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir()

    # Replace relative paths with absolute paths
    def make_absolute(flags: list[str], outdir: Path) -> list[str]:
        result_flags = []
        i = 0
        while i < len(flags):
            flag = flags[i]
            result_flags.append(flag)
            # Check if this flag takes an argument
            if flag in ["-i", "--input", "--input-file", "-r", "--refgene", "--refgenome", "-g"]:
                if i + 1 < len(flags):
                    i += 1
                    filename = flags[i]
                    result_flags.append(str(workload_dir / filename))
            elif flag in ["-o", "--output", "--out-prefix", "--outfile"]:
                if i + 1 < len(flags):
                    i += 1
                    result_flags.append(str(outdir / flags[i]))
            i += 1
        return result_flags

    py_flags_abs = make_absolute(py_flags, py_outdir)
    rust_flags_abs = make_absolute(rust_flags, rust_outdir)

    # Run Python version
    py_cmd = [str(ORACLE_PYTHON), str(ORACLE_SCRIPTS / f"{command}.py")] + py_flags_abs
    env = {"PYTHONPATH": ORACLE_PYTHONPATH}
    py_code, py_out, py_err = run_command(py_cmd, cwd=py_outdir, env=env, timeout_s=60.0)
    result["py_exit"] = py_code

    # Run Rust version
    rust_bin = RUST_BIN_DIR / command
    if not rust_bin.exists():
        result["mismatches"].append(f"Rust binary not found: {rust_bin}")
        result["success"] = False
        return result

    rust_cmd = [str(rust_bin)] + rust_flags_abs
    rust_code, rust_out, rust_err = run_command(rust_cmd, cwd=rust_outdir, timeout_s=60.0)
    result["rust_exit"] = rust_code

    # Compare exit codes
    if py_code != rust_code:
        result["mismatches"].append(f"Exit mismatch: py={py_code}, rust={rust_code}")
        if rust_err:
            result["mismatches"].append(f"rust stderr: {rust_err.strip()[-200:]}")
        result["success"] = False
        return result

    def norm(data: bytes, outdir: Path) -> bytes:
        return data.replace(str(outdir).encode(), b"<OUT>")

    if norm(py_out.encode(), py_outdir) != norm(rust_out.encode(), rust_outdir):
        result["mismatches"].append("stdout mismatch")
        result["success"] = False

    # log.txt is upstream's unscoped CWD side effect (DIV-0022), not an output.
    py_files = {p.name for p in py_outdir.iterdir() if p.is_file()} - {"log.txt"}
    rust_files = {p.name for p in rust_outdir.iterdir() if p.is_file()}
    for name in sorted(py_files ^ rust_files):
        side = "python" if name in py_files else "rust"
        result["mismatches"].append(f"file only on {side} side: {name}")
        result["success"] = False
    for name in sorted(py_files & rust_files):
        if name.endswith((".bam", ".bai", ".gz")):
            continue  # container bytes differ by design (BGZF/gzip headers)
        py_bytes = norm((py_outdir / name).read_bytes(), py_outdir)
        rust_bytes = norm((rust_outdir / name).read_bytes(), rust_outdir)
        if py_bytes != rust_bytes:
            result["mismatches"].append(f"file differs: {name}")
            result["success"] = False

    return result


def main():
    print("RSeQC Synthetic Sweep", file=sys.stderr)
    print(f"Testing {len(WORKLOADS)} workloads × {len(COMMAND_TESTS)} commands", file=sys.stderr)

    all_results = []

    for size, seed in WORKLOADS:
        print(f"\nWorkload: size={size}, seed={seed}", file=sys.stderr)
        with tempfile.TemporaryDirectory() as tmpdir:
            workload_dir = Path(tmpdir) / "workload"
            if not generate_workload(size, seed, workload_dir):
                continue

            print(f"  Workload created in {workload_dir}", file=sys.stderr)

            # Test each command
            for command, flag_variants in COMMAND_TESTS:
                for py_flags, rust_flags in flag_variants:
                    result = test_command(command, py_flags, rust_flags, workload_dir)
                    all_results.append(result)

                    status = "OK" if result["success"] else "FAIL"
                    print(f"  [{status}] {command}", file=sys.stderr)
                    if result["mismatches"]:
                        for m in result["mismatches"]:
                            print(f"      {m}", file=sys.stderr)

    # Summary
    print("\n=== SUMMARY ===", file=sys.stderr)
    total = len(all_results)
    passed = sum(1 for r in all_results if r["success"])
    failed = total - passed

    print(f"Total: {total}, Passed: {passed}, Failed: {failed}", file=sys.stderr)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
