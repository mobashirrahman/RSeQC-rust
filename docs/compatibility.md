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

As reviewed on 2026-09-17 at implementation commit
`b52114315e98e84b201dc769ddf750125570b1a0`, all 33 commands have implementation
attempts, the locked release build succeeds, and 224 Rust unit tests pass.
The [differential runner](../verification/run_diff.py) exists and its five basic
cases pass against the local oracle: bam_stat, read_NVC, read_GC,
read_duplication, and read_quality. These cases reuse one small BAM and check
selected outputs; they do not establish full compatibility for those commands.

The deeper audit reproduced false-pass paths in the runner and scientific
mismatches in FPKM_count, RNA_fragment_size, and geneBody_coverage. Native
rendering, standalone execution of all commands, and Python API compatibility
remain incomplete. No publication, full-compatibility, or speedup claim is
established by the current test results.

The [scientific testing plan](../testing.md) records the audit evidence,
reproducer recipes, all-command test matrix, and separate acceptance gates for
scientific validity, compatibility, performance, and released artifacts.

Since that historical audit, the working tree has an executable first regression
gate: 11 differential cases pass against the pinned local oracle, including
focused FPKM, CIGAR-span, gene-body, and TIN overlap/depth fixtures; the
runner's eight unit tests pass; and strict workspace Clippy is clean. These are
restricted regression results, not evidence that the full 33-command/API,
standalone, real-data, or publication gates are complete. See the working-tree
gate note in [testing.md](../testing.md) and retain the historical revision
numbers above when comparing reports.

See `compatibility/commands.yaml` and `compatibility/api.yaml` for the full
inventory, and `compatibility/divergences.yaml` for the known-divergence
ledger. Newly audited issues in `testing.md` must also be reviewed and migrated
into the versioned issue/divergence records as the plan is implemented; the
existing ledger is not yet an exhaustive account of known risks.
