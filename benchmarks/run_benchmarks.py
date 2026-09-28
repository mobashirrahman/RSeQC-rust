#!/usr/bin/env python3
"""Run benchmarks comparing upstream Python and Rust implementations.

Executes a set of RSeQC commands on a synthetic workload with:
- Deterministic seeded interleaving and randomization of repetitions
- Warmup runs (discarded from analysis)
- Peak RSS measurement via /usr/bin/time
- Paired wall-time ratios and bootstrap confidence intervals
- Output file comparison where applicable

Usage:
    oracle/venv/bin/python3 benchmarks/run_benchmarks.py \
        --workload workloads/test_1000 \
        --commands bam_stat read_distribution \
        --reps 3 \
        --warmup 1 \
        --output-dir benchmark_results/test_1000
"""
import argparse
import json
import os
import platform
import random
import re
import signal
import subprocess
import sys
import tempfile
import time
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parent.parent
ORACLE_SCRIPTS = REPO_ROOT / "oracle" / "upstream-src" / "scripts"
ORACLE_PYTHONPATH = str(REPO_ROOT / "oracle" / "upstream-src" / "src")
ORACLE_PYTHON = REPO_ROOT / "oracle" / "venv" / "bin" / "python3"
RUST_BIN_DIR = REPO_ROOT / "target" / "release"

# Benchmark commands: name -> (py_script, rust_bin, args_fn, compare_files)
BENCHMARK_COMMANDS = {
    "bam_stat": (
        "bam_stat.py",
        "bam_stat",
        lambda workload, out: ["-i", str(workload / "reads.bam")],
        (),
    ),
    "read_distribution": (
        "read_distribution.py",
        "read_distribution",
        lambda workload, out: ["-i", str(workload / "reads.bam"), "-r", str(workload / "model.bed12")],
        (),
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
        ("genebody.geneBodyCoverage.txt",),
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
        ("junction.junction.xls",),
    ),
    "read_duplication": (
        "read_duplication.py",
        "read_duplication",
        lambda workload, out: [
            "-i", str(workload / "reads.bam"),
            "-o", str(out / "dup"),
            "--skip-plot",
        ],
        ("dup.dup.xls",),
    ),
}


class RunResult(NamedTuple):
    exit_code: int | None
    stdout: str
    stderr: str
    wall_time: float
    peak_rss_mb: float | None
    timed_out: bool = False


def run_with_timing(
    argv: list[str],
    pythonpath: str | None = None,
    *,
    cwd: Path = REPO_ROOT,
    timeout_s: float = 300.0,
) -> RunResult:
    """Run command and capture exit code, streams, timing, and peak RSS."""
    env = dict(os.environ)
    if pythonpath is not None:
        env["PYTHONPATH"] = pythonpath

    start_wall = time.time()
    
    # Try using /usr/bin/time if available
    use_time_cmd = Path("/usr/bin/time").exists()
    
    if use_time_cmd:
        # Use /usr/bin/time -v to capture peak RSS
        time_cmd = ["/usr/bin/time", "-v"]
        full_argv = time_cmd + argv
    else:
        full_argv = argv

    try:
        proc = subprocess.Popen(
            full_argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            env=env,
            start_new_session=True,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout_s)
            end_wall = time.time()
            wall_time = end_wall - start_wall
            
            # Parse peak RSS from /usr/bin/time output if available
            peak_rss_mb = None
            if use_time_cmd:
                # Look for "Maximum resident set size" line in stderr
                match = re.search(r"Maximum resident set size[^:]*:\s+(\d+)", stderr)
                if match:
                    # Value is in KB, convert to MB
                    peak_rss_mb = int(match.group(1)) / 1024.0
            
            return RunResult(
                exit_code=proc.returncode,
                stdout=stdout,
                stderr=stderr,
                wall_time=wall_time,
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
                peak_rss_mb=None,
                timed_out=True,
            )
    except OSError as exc:
        return RunResult(
            exit_code=None,
            stdout="",
            stderr=f"[benchmark could not execute command: {exc}]",
            wall_time=0,
            peak_rss_mb=None,
        )


def numeric_equal(left: str, right: str) -> bool:
    """Compare values with numeric semantics."""
    try:
        lval = Decimal(left)
        rval = Decimal(right)
    except InvalidOperation:
        return left == right
    if lval.is_nan() or rval.is_nan():
        return lval.is_nan() and rval.is_nan()
    return lval == rval


def numeric_table_equal(left_bytes: bytes, right_bytes: bytes) -> bool:
    """Compare whitespace-delimited tables with numeric cell semantics."""
    left_rows = left_bytes.decode("utf-8", errors="replace").splitlines()
    right_rows = right_bytes.decode("utf-8", errors="replace").splitlines()
    if len(left_rows) != len(right_rows):
        return False
    
    for left_row, right_row in zip(left_rows, right_rows):
        left_cells = left_row.split()
        right_cells = right_row.split()
        if len(left_cells) != len(right_cells):
            return False
        for left_cell, right_cell in zip(left_cells, right_cells):
            try:
                Decimal(left_cell)
                Decimal(right_cell)
            except InvalidOperation:
                if left_cell != right_cell:
                    return False
            else:
                if not numeric_equal(left_cell, right_cell):
                    return False
    return True


def compare_files(py_dir: Path, rust_dir: Path, file_list: tuple[str, ...]) -> tuple[bool, list[str]]:
    """Compare output files between Python and Rust runs."""
    errors = []
    for rel_path in file_list:
        py_file = py_dir / rel_path
        rust_file = rust_dir / rel_path
        if not py_file.is_file() or not rust_file.is_file():
            errors.append(
                f"file '{rel_path}': python_exists={py_file.is_file()} rust_exists={rust_file.is_file()}"
            )
            continue
        
        py_bytes = py_file.read_bytes()
        rust_bytes = rust_file.read_bytes()
        
        # Use numeric table comparison for text files
        if not numeric_table_equal(py_bytes, rust_bytes):
            # Report the first difference
            py_lines = py_bytes.decode("utf-8", errors="replace").splitlines()
            rust_lines = rust_bytes.decode("utf-8", errors="replace").splitlines()
            for i, (py_line, rust_line) in enumerate(zip(py_lines, rust_lines)):
                if py_line != rust_line:
                    errors.append(f"file '{rel_path}' line {i+1}: python={py_line[:50]!r} rust={rust_line[:50]!r}")
                    break
            else:
                if len(py_lines) != len(rust_lines):
                    errors.append(f"file '{rel_path}': python_lines={len(py_lines)} rust_lines={len(rust_lines)}")
                else:
                    errors.append(f"file '{rel_path}': content differs")
    
    return len(errors) == 0, errors


def bootstrap_ci(values: list[float], seed: int = 42, resamples: int = 10000) -> tuple[float, float]:
    """Compute 95% confidence interval via bootstrap."""
    if not values or len(values) < 2:
        return 0, 0
    
    rng = random.Random(seed)
    bootstrap_medians = []
    
    for _ in range(resamples):
        sample = [rng.choice(values) for _ in range(len(values))]
        bootstrap_medians.append(sorted(sample)[len(sample) // 2])
    
    sorted_medians = sorted(bootstrap_medians)
    lower_idx = int(0.025 * len(sorted_medians))
    upper_idx = int(0.975 * len(sorted_medians))
    return sorted_medians[lower_idx], sorted_medians[upper_idx]


def run_benchmark(
    command: str,
    workload_dir: Path,
    output_dir: Path,
    num_reps: int = 3,
    num_warmup: int = 1,
    rng_seed: int = 42,
) -> dict:
    """Run a single benchmark command with both Python and Rust."""
    if command not in BENCHMARK_COMMANDS:
        raise ValueError(f"Unknown command: {command}")
    
    py_script, rust_bin, args_fn, compare_files_list = BENCHMARK_COMMANDS[command]
    
    results = {
        "command": command,
        "workload": str(workload_dir),
        "num_reps": num_reps,
        "num_warmup": num_warmup,
        "python": {"runs": []},
        "rust": {"runs": []},
        "compatibility": {"match": True, "errors": []},
        "statistics": {},
    }
    
    # Generate randomized interleaved order of (rep_index, side) tuples
    rng = random.Random(rng_seed)
    run_order = []
    
    # Warmup runs (discarded)
    for _ in range(num_warmup):
        run_order.append(("warmup", rng.choice(["python", "rust"])))
    
    # Measured runs - interleaved
    for rep in range(num_reps):
        run_order.append((rep, "python"))
        run_order.append((rep, "rust"))
    
    # Randomize the order
    rng.shuffle(run_order)
    
    # Track measurements per side and per rep
    py_times = {}
    rust_times = {}
    py_dirs = {}
    rust_dirs = {}
    
    # Execute in randomized order
    for order_item, side in run_order:
        is_warmup = order_item == "warmup"
        rep = 0 if is_warmup else order_item
        
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            
            # Build argument list
            args = args_fn(workload_dir, run_dir)
            
            if side == "python":
                argv = [str(ORACLE_PYTHON), str(ORACLE_SCRIPTS / py_script)] + args
                result = run_with_timing(argv, pythonpath=ORACLE_PYTHONPATH, cwd=run_dir, timeout_s=300.0)
            else:
                argv = [str(RUST_BIN_DIR / rust_bin)] + args
                result = run_with_timing(argv, cwd=run_dir, timeout_s=300.0)
            
            if not is_warmup:
                # Record measured run
                run_rec = {
                    "repetition": rep + 1,
                    "exit_code": result.exit_code,
                    "wall_time_s": result.wall_time,
                    "peak_rss_mb": result.peak_rss_mb,
                    "timed_out": result.timed_out,
                }
                
                if side == "python":
                    results["python"]["runs"].append(run_rec)
                    if result.wall_time > 0:
                        py_times[rep] = result.wall_time
                    py_dirs[rep] = run_dir
                else:
                    results["rust"]["runs"].append(run_rec)
                    if result.wall_time > 0:
                        rust_times[rep] = result.wall_time
                    rust_dirs[rep] = run_dir
                
                # Check exit codes
                if result.exit_code != 0:
                    results["compatibility"]["match"] = False
                    results["compatibility"]["errors"].append(
                        f"Rep {rep + 1} ({side}): exit code {result.exit_code}"
                    )
    
    # Compare output files across reps (when both sides completed)
    for rep in range(num_reps):
        if rep in py_dirs and rep in rust_dirs:
            files_ok, file_errors = compare_files(py_dirs[rep], rust_dirs[rep], compare_files_list)
            if not files_ok:
                results["compatibility"]["match"] = False
                for err in file_errors:
                    results["compatibility"]["errors"].append(f"Rep {rep + 1}: {err}")
    
    # Compute statistics (only if we have matching repetitions)
    all_reps = set(py_times.keys()) & set(rust_times.keys())
    if len(all_reps) > 0:
        py_times_list = [py_times[rep] for rep in sorted(all_reps)]
        rust_times_list = [rust_times[rep] for rep in sorted(all_reps)]
        
        # Paired time ratios (Python / Rust)
        time_ratios = [py_times_list[i] / rust_times_list[i] if rust_times_list[i] > 0 else 0 
                      for i in range(len(py_times_list))]
        
        py_median = sorted(py_times_list)[len(py_times_list) // 2]
        rust_median = sorted(rust_times_list)[len(rust_times_list) // 2]
        
        py_ci = bootstrap_ci(py_times_list, seed=rng_seed)
        rust_ci = bootstrap_ci(rust_times_list, seed=rng_seed)
        ratio_ci = bootstrap_ci(time_ratios, seed=rng_seed + 1)
        
        results["statistics"] = {
            "python_median_s": py_median,
            "python_ci_95": {"lower": py_ci[0], "upper": py_ci[1]},
            "rust_median_s": rust_median,
            "rust_ci_95": {"lower": rust_ci[0], "upper": rust_ci[1]},
            "time_ratio_median": py_median / rust_median if rust_median > 0 else 0,
            "time_ratio_ci_95": {"lower": ratio_ci[0], "upper": ratio_ci[1]},
        }
    
    return results


def get_environment_info() -> dict:
    """Capture environment for reproducibility."""
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
        help="Number of measured repetitions per command (default: 3)",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=1,
        help="Number of warmup runs (discarded, default: 1)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for benchmark results",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible interleaving (default: 42)",
    )
    
    args = parser.parse_args()
    
    # Validate workload
    if not args.workload.is_dir():
        print(f"Error: workload directory not found: {args.workload}", file=sys.stderr)
        sys.exit(1)
    
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80, file=sys.stderr)
    print("BENCHMARK RUNNER - NON-ISOLATED HARDWARE WARNING", file=sys.stderr)
    print("=" * 80, file=sys.stderr)
    print("Results are NOT suitable for publication.", file=sys.stderr)
    print("See testing.md section 12 for publication requirements.", file=sys.stderr)
    print("=" * 80, file=sys.stderr)
    
    print(f"\nBenchmark Configuration:", file=sys.stderr)
    print(f"  Workload: {args.workload}", file=sys.stderr)
    print(f"  Commands: {', '.join(args.commands)}", file=sys.stderr)
    print(f"  Measured reps: {args.reps}, Warmup: {args.warmup}", file=sys.stderr)
    print(f"  Output: {output_dir}\n", file=sys.stderr)
    
    env_info = get_environment_info()
    
    all_results = []
    for command in args.commands:
        print(f"Benchmarking {command}...", file=sys.stderr)
        try:
            result = run_benchmark(command, args.workload, output_dir, args.reps, args.warmup, args.seed)
            all_results.append(result)
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
    
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
    print(f"Results: {results_file}", file=sys.stderr)
    
    # Print summary
    print("\nBenchmark Summary:", file=sys.stderr)
    print("-" * 80, file=sys.stderr)
    for cmd_result in all_results:
        command = cmd_result["command"]
        match_status = "✓" if cmd_result["compatibility"]["match"] else "✗"
        
        if cmd_result["statistics"]:
            stats = cmd_result["statistics"]
            ratio = stats["time_ratio_median"]
            ratio_ci = stats["time_ratio_ci_95"]
            print(
                f"{match_status} {command:25} Speedup: {ratio:6.2f}x  95% CI: [{ratio_ci['lower']:6.2f}, {ratio_ci['upper']:6.2f}]",
                file=sys.stderr
            )
        else:
            print(f"{match_status} {command:25} (no statistics)", file=sys.stderr)
        
        if cmd_result["compatibility"]["errors"]:
            for err in cmd_result["compatibility"]["errors"]:
                print(f"    {err}", file=sys.stderr)


if __name__ == "__main__":
    main()
