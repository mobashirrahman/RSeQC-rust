# ADR 0001: Rust with standalone CLI artifacts

Status: accepted target direction; dependency and platform selections pending investigation.

## Decision

Implement the RSeQC port in Rust. Prioritize performance and standalone command binaries, as selected by the user on 2026-09-16.

All scientific algorithms live in Rust. The primary CLI distribution must execute all 33 commands, including plotting and format conversion, without separately installed Python, R, or helper executables. Preserve original executable names and CLI behavior. Linking/bundling native libraries is allowed subject to dependency review and verified artifact portability; standalone does not imply pure Rust or fully static linkage.

Provide optional Python bindings and launchers for the existing `qcmodule` namespace and Python script invocation. Their runtime requirements do not apply to the standalone CLI. Full-package compatibility remains incomplete until this API surface passes its own verification.

## Architecture and consequences

- Cargo workspace: core, formats, commands, rendering, CLI, and Python bindings.
- One shared Rust implementation powers both interfaces; no production fallback to upstream Python.
- Pin the Rust toolchain and dependencies; record target triples and minimum platform/ABI requirements.
- Native plotting is in scope. Compare underlying data, legacy script artifacts where required, and rendered output; document any unavoidable visual differences.
- Use upstream Python/R only in isolated reference environments and development verification.
- Start with a verified single-thread implementation. Add bounded parallelism only after profiling and preserve ordering, numerical, and sampling semantics.
- Validate format-library and rendering feasibility early; select exact dependencies after source/documentation review and interoperability probes.
- Audit standalone artifacts in environments without reference runtimes. Verify Python packages separately.

## Release criteria

The standalone CLI milestone does not by itself establish full Python-package compatibility. The complete release requires all compatibility, scientific verification, benchmark, packaging, and publication gates in [the execution plan](../PORTING_PLAN.md).

No speedup factor is promised before measurement. Publish per-command performance results, including regressions and remaining limitations.
