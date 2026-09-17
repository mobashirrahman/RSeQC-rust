# Contributing

Thanks for your interest in RSeQC-rust. This project is a from-scratch Rust reimplementation of
[RSeQC](https://github.com/liguowang/RSeQC), so its overriding priority is **behavioral
compatibility with a frozen upstream baseline**, not just "does it look right." Please read
[`docs/PORTING_PLAN.md`](docs/PORTING_PLAN.md) and [`testing.md`](testing.md) before making
substantial changes — they're the authoritative plan and verification standard this project holds
itself to.

## Before you start

- Check [`compatibility/divergences.yaml`](compatibility/divergences.yaml) for known, already-
  documented differences from upstream. If your change touches one of these, update the entry
  rather than leaving it stale.
- Check [`compatibility/commands.yaml`](compatibility/commands.yaml) for the command inventory and
  [`compatibility/upstream.lock`](compatibility/upstream.lock) for the exact upstream commit this
  port targets.

## Development setup

```sh
git clone <this-repo>
cd RSeQC-rust
cargo build --workspace
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

The differential verification harness (`verification/run_diff.py`) additionally needs a real
upstream Python environment, which is gitignored and not part of a normal clone — see
`docs/PORTING_PLAN.md` for how to set up `oracle/venv` and `oracle/upstream-src`. Without it, you
can still build, run unit tests, and lint; you just can't run the harness that diffs this port's
output against the real upstream CLI.

## Compatibility-preservation philosophy

This is the single most important thing to internalize before contributing code:

- **Never silently "fix" an upstream bug, typo, or quirk.** If upstream's console output has a
  literal typo (e.g. "Totoal reads used" instead of "Total"), this port reproduces it exactly.
  Silently correcting it would make behavior diverge from what real users of upstream actually
  see.
- **Any accepted difference from upstream must be disclosed**, not just handled in code. Add or
  update an entry in `compatibility/divergences.yaml` explaining what differs, why, and what
  `status` applies (`open`, `accepted`, or `fixed-working-tree`).
- **Verify against the real upstream CLI, not just your own expectations.** A `python3 -c` probe
  of upstream's literal behavior, or better, a `verification/run_diff.py` case that runs the real
  upstream script and this port's compiled binary against the same fixture and diffs their actual
  output, is worth far more than a unit test asserting what you *believe* the correct output to
  be. This project's own history has repeatedly found real bugs this way that looked "obviously
  correct" in isolation.

## Making a change

1. If you're fixing a compatibility bug, first reproduce it: build a minimal fixture, run both the
   real upstream script and this port's binary against it, and confirm they actually differ.
2. Make the fix. Prefer editing existing command/format modules over introducing new
   abstractions — see `crates/commands/src/*.rs` for the existing per-command pattern (a pure
   `compute_*`/`render_*` function, no I/O, tested independently of the CLI layer).
3. Re-run the differential harness case for the affected command (or add one if none exists yet —
   see any `ensure_*_fixture` function in `verification/run_diff.py` for the established pattern).
4. Run `cargo test --workspace` and `cargo clippy --workspace --all-targets -- -D warnings`; both
   must stay clean.
5. If your change touches a documented divergence, update `compatibility/divergences.yaml`.

## Code style

- No premature abstraction: three similar lines are better than a speculative helper function
  built for hypothetical future reuse.
- Comments should explain *why*, not *what* — especially for anything preserving a non-obvious
  upstream quirk. A future reader needs to know *why* this port does something that looks wrong at
  a glance.
- `cargo fmt` drift currently exists across parts of the codebase predating any formatting policy
  (see the CI workflow's own note on this) — don't mix a repo-wide reformat into an unrelated
  change; that belongs in its own dedicated PR.

## Reporting a compatibility bug

Open an issue with: the exact command and flags, the input file's shape (or a minimal
reproduction), the real upstream output, and this port's actual output. If you can also link the
exact line(s) in `oracle/upstream-src` (or your own local upstream checkout) responsible for the
behavior, that speeds up triage significantly.
