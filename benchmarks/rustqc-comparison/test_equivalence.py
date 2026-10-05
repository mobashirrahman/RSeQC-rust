#!/usr/bin/env python3
"""Card R4 tests: failure demonstrations for equivalence.verdict_file/text."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from equivalence import verdict_file, verdict_text


def _write(td: Path, name: str, content: str) -> Path:
    p = td / name
    p.write_text(content)
    return p


def test_changed_digit_gives_differs():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        a = _write(td, "a.xls", "tx\t10.0\ntx2\t20.0\n")
        b = _write(td, "b.xls", "tx\t10.0\ntx2\t20.1\n")
        v, d = verdict_file(a, b)
        assert v == "differs", f"expected differs, got {v}: {d}"
        print("PASS changed digit -> differs")


def test_deleted_file_gives_missing():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        a = _write(td, "a.xls", "x\n")
        b = td / "absent.xls"
        v, d = verdict_file(a, b)
        assert v == "missing", f"expected missing, got {v}: {d}"
        print("PASS deleted file -> missing")


def test_two_empty_dirs_no_pass():
    # Two empty output trees must not yield a passing verdict: model as two
    # empty files plus two empty streams.
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        a = _write(td, "a.xls", "")
        b = _write(td, "b.xls", "")
        v, _ = verdict_file(a, b)
        assert v not in ("byte-identical", "numerically-identical", "within-tolerance"), \
            f"empty files must not pass, got {v}"
        v2, _ = verdict_text("", "")
        assert v2 not in ("byte-identical", "numerically-identical", "within-tolerance"), \
            f"empty streams must not pass, got {v2}"
        print(f"PASS two empty -> no pass ({v} / {v2})")


def test_nan_gives_differs():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        a = _write(td, "a.xls", "tx\t1.5\n")
        b = _write(td, "b.xls", "tx\tNaN\n")
        v, d = verdict_file(a, b)
        assert v == "differs", f"expected differs for NaN, got {v}: {d}"
        # symmetric NaN is still not agreement
        c = _write(td, "c.xls", "tx\tNaN\n")
        v2, _ = verdict_file(b, c)
        assert v2 == "differs", f"symmetric NaN must differ, got {v2}"
        print("PASS NaN -> differs (including symmetric NaN)")


def test_byte_identical_and_numeric():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        a = _write(td, "a.xls", "tx\t1.0\tfoo\n")
        b = _write(td, "b.xls", "tx\t1.0\tfoo\n")
        v, _ = verdict_file(a, b)
        assert v == "byte-identical", f"expected byte-identical, got {v}"
        c = _write(td, "c.xls", "tx\t1.000000\tfoo\n")
        v2, _ = verdict_file(a, c)
        assert v2 == "numerically-identical", f"expected numerically-identical, got {v2}"
        print("PASS byte-identical + numerically-identical")


if __name__ == "__main__":
    test_changed_digit_gives_differs()
    test_deleted_file_gives_missing()
    test_two_empty_dirs_no_pass()
    test_nan_gives_differs()
    test_byte_identical_and_numeric()
    print("all equivalence tests passed")
