## What does this change

## Why

<!-- Especially for a compatibility fix: what did you observe differing from real upstream, and
how did you confirm it (a verification/run_diff.py case is the strongest evidence)? -->

## Verification

- [ ] `cargo test --workspace` passes
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` passes
- [ ] If this changes a command's output: a `verification/run_diff.py` case exists covering it,
      and passes against real upstream
- [ ] If this closes or changes a known divergence: `compatibility/divergences.yaml` is updated

## Checklist

- [ ] I have not silently "fixed" an upstream bug/quirk without disclosing it (see CONTRIBUTING.md)
- [ ] I have not introduced new abstractions beyond what this change needs
