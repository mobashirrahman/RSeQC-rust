#!/usr/bin/env python3
"""Card R5 tests: timeout leaves harness alive + writes row; failure
recorded with exit status; block order reproducible from seed."""
import os
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from timing import build_schedule, parse_time_v, run_timed


def test_timeout_writes_row_and_harness_alive():
    with tempfile.TemporaryDirectory() as td:
        m, out, err, code, to = run_timed(
            ["sleep", "30"], cwd=td, env=dict(os.environ), timeout_s=0.5)
        assert to is True, "expected timed_out=True"
        # Harness alive: can still run another command.
        m2, _, _, code2, to2 = run_timed(
            ["true"], cwd=td, env=dict(os.environ), timeout_s=5)
        assert to2 is False and code2 == 0, "harness did not survive timeout"
        assert "wall_s" in m, "row lacks wall measurement"
        print(f"PASS timeout row + harness alive (wall={m['wall_s']:.2f}s)")


def test_failed_process_recorded():
    with tempfile.TemporaryDirectory() as td:
        m, out, err, code, to = run_timed(
            ["bash", "-c", "exit 3"], cwd=td, env=dict(os.environ), timeout_s=5)
        assert to is False, "false timeout"
        assert code == 3, f"expected exit 3, got {code}"
        print("PASS failed process exit=3 recorded")


def test_block_order_reproducible():
    s1 = build_schedule(3, 2, seed=42)
    s2 = build_schedule(3, 2, seed=42)
    s3 = build_schedule(3, 2, seed=7)
    assert s1 == s2, "same seed gave different schedules"
    assert s1 != s3, "different seeds gave identical schedules"
    # U-seq appears only in blocks < reps_upstream
    for item in s1:
        if item["config"] == "U-seq":
            assert item["block"] < 2, "U-seq beyond reps_upstream"
    # Every block's orders are a permutation starting at 0
    from collections import defaultdict
    by_block = defaultdict(list)
    for item in s1:
        by_block[item["block"]].append(item["order"])
    for block, orders in by_block.items():
        assert sorted(orders) == list(range(len(orders))), f"block {block} orders broken"
    print(f"PASS reproducible schedule ({len(s1)} runs, seed 42)")


def test_parse_time_v():
    sample = ("Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.23\n"
              "User time (seconds): 0.90\n"
              "System time (seconds): 0.10\n"
              "Maximum resident set size (kbytes): 102400\n")
    m = parse_time_v(sample)
    assert abs(m["peak_rss_mb"] - 100.0) < 1e-9, m
    assert abs(m["user_s"] - 0.90) < 1e-9, m
    print("PASS parse_time_v")


if __name__ == "__main__":
    test_timeout_writes_row_and_harness_alive()
    test_failed_process_recorded()
    test_block_order_reproducible()
    test_parse_time_v()
    print("all timing tests passed")
