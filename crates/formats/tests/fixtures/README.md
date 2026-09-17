`pybigwig_test.bw` is copied unmodified from pyBigWig's own test suite
(`pyBigWigTest/test.bw` in the `pyBigWig` PyPI package, MIT licensed).
Used as a real, independently-verifiable BigWig fixture for
`crates/formats/src/bigwig.rs` tests: ground truth for its contents was
captured via `pyBigWig.open(...)` directly (see
`crates/formats/src/bigwig.rs` test module for the exact expected
values). Two chromosomes: `1` (195471971bp) and `10` (130694993bp), a
handful of small intervals near the start of each.
