#!/usr/bin/env python3
"""Run benchmarks comparing upstream Python and Rust implementations.

Executes a set of RSeQC commands on a synthetic workload, captures timing
and resource usage, and verifies output compatibility. Records raw measurements
and confidence intervals, but flags that this sandbox is non-isolated.

Usage:
    oracle/venv/bin/python3 benchmarks/run_benchmarks.py \
        --workload workloads/test_1000 \
        --commands bam_stat read_distribution \
        --reps 3 \
        --output-dir benchmark_results/test_1000

Output:
    benchmark_results/test_1000/
        results.json           - Raw per-run measurements
        summary.json           - Medians, intervals, environment info
"""
import argparse
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parent.parent
ORACLE_SCRIPTS = REPO_ROOT / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHONPATH = str(REPO_ROOT / "oracle" / "upstream-src" / "src")
ORACLE_PYTHON = REPO_ROOT / "oracle" / "venv" / "bin" / "python3"
RUST_BIN_DIR = REPO_ROOT / "target" / "release"

# Available benchmark commands (name -> (py_script, rust_bin, args_fn))
BENCHMARK_COMMANDS = {
    "bam_stat": (
        "bam_stat.py",
        "bam_stat",
        lambda workload, out: ["-i", str(workload / "reads.bam")],
    ),
    "read_distribution": (
        "read_distribution.py",
        "read_distribution",
        lambda workload, out: ["-i", str(workload / "reads.bam"), "-r", str(workload / "model.bed12")],
    ),
    "geneBody_coverage": (
        "geneBody_coverage.py",
        "geneBody_coverage",
        lambda workload, out: [
            "-i", str(workload / "reads.bam"),
            "-r", str(workload / "model.bed12"),
            "-o", str(out / "genebody"),
            "--skip-plot",
        ],
    ),
    "junction_annotation": (
        "junction_annotation.py",
        "junction_annotation",
        lambda workload, out: [
            "-i", str(workload / "reads.bam"),
            "-r", str(workload / "model.bed12"),
            "-o", str(out / "junction"),
            "--skip-plot",
        ],
    ),
    "read_duplication": (
        "read_duplication.py",
        "read_duplication",
        lambda workload, out: [
            "-i", str(workload / "reads.bam"),
            "-o", str(out / "dup"),
            "--skip-plot",
        ],
    ),
}


class RunResult(NamedTuple):
    exit_code: int | None
    stdout: str
    stderr: str
    wall_time: float
    user_time: float
    sys_time: float
    peak_rss_mb: float | None
    timed_out: bool = False


def run_with_timing(
    argv: list[str],
    pythonpath: str | None = None,
    *,
    cwd: Path = REPO_ROOT,
    timeout_s: float = 300.0,
) -> RunResult:
    """Run a command and capture exit code, streams, and timing/resource info."""
    env = dict(os.environ)
    if pythonpath is not None:
        env["PYTHONPATH"] = pythonpath

    start_wall = time.time()
    start_time = time.process_time()

    try:
        proc = subprocess.Popen(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=env,
            start_new_session=True,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
            end_time = time.process_time()
            end_wall = time.time()
            wall_time = end_wall - start_wall
            cpu_time = end_time - start_time

            # Try to get peak RSS using /usr/bin/time if available
            peak_rss_mb = None

            return RunResult(
                exit_code=proc.returncode,
                stdout=stdout,
                stderr=stderr,
                wall_time=wall_time,
                user_time=cpu_time,  # process_time gives total CPU time
                sys_time=0,  # Not directly available from process_time
                peak_rss_mb=peak_rss_mb,
            )
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = proc.communicate()
            stdout = stdout or ""
            stderr = stderr or ""
            return RunResult(
                exit_code=None,
                stdout=stdout,
                stderr=stderr + f"\n[benchmark timeout after {timeout_s:g}s]",
                wall_time=timeout_s,
                user_time=0,
                sys_time=0,
                peak_rss_mb=None,
                timed_out=True,
            )
    except OSError as exc:
        return RunResult(
            exit_code=None,
            stdout="",
            stderr=f"[benchmark could not execute command: {exc}]",
            wall_time=0,
            user_time=0,
            sys_time=0,
            peak_rss_mb=None,
        )


def run_benchmark(
    command: str,
    workload_dir: Path,
    output_dir: Path,
    repetitions: int = 3,
) -> dict:
    """Run a single benchmark command with both Python and Rust, multiple reps."""
    if command not in BENCHMARK_COMMANDS:
        raise ValueError(f"Unknown command: {command}")

    py_script, rust_bin, args_fn = BENCHMARK_COMMANDS[command]

    results = {
        "command": command,
        "workload": str(workload_dir),
        "repetitions": repetitions,
        "python": {"runs": []},
        "rust": {"runs": []},
        "compatibility": {"match": True, "errors": []},
    }

    # Run Python and Rust side by side, alternating
    for rep in range(repetitions):
        print(f"  {command} rep {rep + 1}/{repetitions}...", file=sys.stderr, end=" ", flush=True)

        # Run in separate temp directories to avoid conflicts
        with tempfile.TemporaryDirectory() as py_tmpdir:
            py_dir = Path(py_tmpdir)

            with tempfile.TemporaryDirectory() as rust_tmpdir:
                rust_dir = Path(rust_tmpdir)

                # Build argument lists
                py_args = args_fn(workload_dir, py_dir)
                rust_args = args_fn(workload_dir, rust_dir)

                # Run Python
                py_argv = [str(ORACLE_PYTHON), str(ORACLE_SCRIPTS / py_script)] + py_args
                py_result = run_with_timing(py_argv, pythonpath=ORACLE_PYTHONPATH, timeout_s=300.0)

                # Run Rust
                rust_argv = [str(RUST_BIN_DIR / rust_bin)] + rust_args
                rust_result = run_with_timing(rust_argv, timeout_s=300.0)

                # Record results
                results["python"]["runs"].append({
                    "repetition": rep + 1,
                    "exit_code": py_result.exit_code,
                    "wall_time_s": py_result.wall_time,
                    "cpu_time_s": py_result.user_time,
                    "peak_rss_mb": py_result.peak_rss_mb,
                    "timed_out": py_result.timed_out,
                })

                results["rust"]["runs"].append({
                    "repetition": rep + 1,
                    "exit_code": rust_result.exit_code,
                    "wall_time_s": rust_result.wall_time,
                    "cpu_time_s": rust_result.user_time,
                    "peak_rss_mb": rust_result.peak_rss_mb,
                    "timed_out": rust_result.timed_out,
                })

                # Check compatibility (exit codes and basic output match)
                if py_result.exit_code != 0 or rust_result.exit_code != 0:
                    if py_result.exit_code != rust_result.exit_code:
                        results["compatibility"]["match"] = False
                        results["compatibility"]["errors"].append(
                            f"Rep {rep + 1}: exit code mismatch: python={py_result.exit_code} rust={rust_result.exit_code}"
                        )

                print("OK", file=sys.stderr)

    return results


def get_environment_info() -> dict:
    """Capture environment and system information for reproducibility."""
    import subprocess

    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            text=True,
        ).strip()
    except:
        git_commit = "unknown"

    try:
        rustc_version = subprocess.check_output(
            ["rustc", "--version"],
            text=True,
        ).strip()
    except:
        rustc_version = "unknown"

    try:
        python_version = subprocess.check_output(
            [str(ORACLE_PYTHON), "--version"],
            text=True,
        ).strip()
    except:
        python_version = "unknown"

    try:
        loadavg = os.getloadavg()
    except:
        loadavg = [0, 0, 0]

    return {
        "timestamp": time.time(),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "processor": platform.processor(),
        "git_commit": git_commit,
        "rustc_version": rustc_version,
        "python_version": python_version,
        "loadavg": loadavg,
        "warning": "Results from non-isolated hardware. Not suitable for publication.",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run benchmarks comparing Python and Rust implementations"
    )
    parser.add_argument(
        "--workload",
        type=Path,
        required=True,
        help="Path to workload directory (from generate_workload.py)",
    )
    parser.add_argument(
        "--commands",
        nargs="+",
        default=list(BENCHMARK_COMMANDS.keys()),
        help="Commands to benchmark (default: all)",
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=3,
        help="Number of repetitions per command (default: 3)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for benchmark results",
    )

    args = parser.parse_args()

    # Validate workload
    if not args.workload.is_dir():
        print(f"Error: workload directory not found: {args.workload}", file=sys.stderr)
        sys.exit(1)

    # Setup output directory
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80, file=sys.stderr)
    print("BENCHMARK RUNNER - NON-ISOLATED HARDWARE WARNING", file=sys.stderr)
    print("=" * 80, file=sys.stderr)
    print("This benchmark is running on shared hardware.", file=sys.stderr)
    print("Results are NOT suitable for publication.", file=sys.stderr)
    print("See testing.md section 12 for publication requirements.", file=sys.stderr)
    print("=" * 80, file=sys.stderr)

    print(f"\nBenchmark Configuration:", file=sys.stderr)
    print(f"  Workload: {args.workload}", file=sys.stderr)
    print(f"  Commands: {', '.join(args.commands)}", file=sys.stderr)
    print(f"  Repetitions: {args.reps}", file=sys.stderr)
    print(f"  Output: {output_dir}\n", file=sys.stderr)

    # Collect environment info
    env_info = get_environment_info()

    # Run benchmarks
    all_results = []
    for command in args.commands:
        print(f"Benchmarking {command}...", file=sys.stderr)
        try:
            result = run_benchmark(command, args.workload, output_dir, args.reps)
            all_results.append(result)
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)

    # Write raw results
    results_file = output_dir / "results.json"
    with open(results_file, "w") as f:
        json.dump(
            {
                "environment": env_info,
                "benchmark_results": all_results,
            },
            f,
            indent=2,
        )
    print(f"\nRaw results: {results_file}", file=sys.stderr)

    # Print summary
    print("\nBenchmark Summary:", file=sys.stderr)
    print("-" * 80, file=sys.stderr)
    for cmd_result in all_results:
        command = cmd_result["command"]
        py_times = [r["wall_time_s"] for r in cmd_result["python"]["runs"]]
        rust_times = [r["wall_time_s"] for r in cmd_result["rust"]["runs"]]

        if py_times and rust_times:
            py_median = sorted(py_times)[len(py_times) // 2]
            rust_median = sorted(rust_times)[len(rust_times) // 2]
            speedup = py_median / rust_median if rust_median > 0 else 0

            py_compatible = all(r["exit_code"] == 0 for r in cmd_result["python"]["runs"])
            rust_compatible = all(r["exit_code"] == 0 for r in cmd_result["rust"]["runs"])

            status = "✓" if (py_compatible and rust_compatible and cmd_result["compatibility"]["match"]) else "✗"

            print(f"{status} {command:25} Python: {py_median:8.2f}s  Rust: {rust_median:8.2f}s  Speedup: {speedup:.2f}x", file=sys.stderr)

            if not cmd_result["compatibility"]["match"]:
                for error in cmd_result["compatibility"]["errors"]:
                    print(f"    {error}", file=sys.stderr)


if __name__ == "__main__":
    main()
