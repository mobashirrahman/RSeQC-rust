# Compatibility levels and baseline (Step 1 / G0)

Baseline: see `compatibility/upstream.lock` (commit `59a24c5a0c4b0c86b7ac3f5babf92e4eb69c2f24`,
declared package version `5.0.5`, GitHub chosen as primary baseline pending an
independent PyPI cross-check — open item).

## Compatibility levels

This project reports compatibility per command/API item at one or more of these
levels. A command may pass some levels and not others; that is recorded
explicitly rather than collapsed into a single "compatible" claim.

1. **Invocation compatibility** — original executable name (including `.py`
   suffix on PATH), flags, defaults, `--help` text intent, exit codes, and
   where output is written match upstream.
2. **Deterministic textual compatibility** — byte-identical stdout/output files
   for a given input, modulo explicitly named normalization (paths, timestamps,
   or other documented unstable fields — normalization must be named per field,
   never applied broadly).
3. **Numeric equivalence** — computed metrics match within a metric-specific,
   justified absolute/relative tolerance, with explicit NaN/Inf/zero/rounding
   rules.
4. **Binary-file semantic equivalence** — decoded BAM/BigWig/etc. content
   (headers, records, values) matches; raw compressed bytes may differ.
5. **Visual equivalence** — plot data and, where the interface promises it,
   generated scripts/rendered output match within documented limits.
6. **Python API equivalence** — imports, signatures, types, values, exceptions,
   mutability, and iterator behavior match for the `qcmodule` surface.

## Status

No command has passed a differential verification run yet — the verification
harness (`verification/`, PORTING_PLAN.md Step 4) does not exist yet. The
`bam_stat` vertical slice (Step 5) is in progress: a from-scratch Rust
implementation exists and has its own independent unit tests (hand-computed
expected values from synthetic fixtures), but it has NOT been differentially
compared against the oracle. Do not read "has a Rust implementation" as "is
compatibility-verified" anywhere in this repo until Step 4/5 gates close.

See `compatibility/commands.yaml` and `compatibility/api.yaml` for the full
inventory, and `compatibility/divergences.yaml` for the known-divergence
ledger.
