# 0040 — Versioned C ABI boundary

- **Status:** Accepted
- **Date:** 2026-09-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amended by:** [0105](./0105-drawings-join-the-discrete-and-application-bundles.md): the C ABI's scope includes drawings; C ABI symbols are never feature-gated.

## Context

Native downstream applications need Axiolid without adopting Rust types or Rust's unstable ABI. The boundary must survive hostile scalar values, make exact and approximate results distinguishable, prevent unwinding, and assign every allocation to exactly one owner. The existing architecture forbids unsafe Rust in facade and contract crates; dereferencing caller-provided C buffers necessarily requires one narrow, audited exception.

## Decision

We will expose ABI version 0.4 exclusively from `axiolid-capi`, over the supported `axiolid::application` facade.

- Every exported symbol is prefixed `axiolid_v0_4_`; the version query reports ABI and crate versions independently.
- Contexts own all mesh and operation-result handles. Handles are non-zero, globally unique integer tokens and never expose Rust addresses.
- Callers own import/export and error-message buffers. No Rust allocation crosses the ABI.
- Every operation result reports `TriangleMesh` or `ExactBrep` before ownership can be consumed. Unsupported exact operations return `UnsupportedExact`; they never tessellate silently.
- All exports contain panics with `catch_unwind`. The crate denies unsafe operations in unsafe functions, while every other facade crate continues to forbid unsafe code.
- C integer inputs are validated as integers rather than accepted as Rust enum discriminants. Nulls, lengths, finite values, limits, handles, and result kinds are checked before use.
- A generated C11 header, Rust boundary tests, and a compiled C smoke test are mandatory gates.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Export Rust structs or trait objects | Rust layout and vtables are not a stable C ABI and would leak implementation ownership. |
| Return Rust-owned buffers with a free function | More allocation provenance and double-free states than caller-sized copies; unnecessary for the v0.4 surface. |
| Use raw pointer handles | Enables accidental dereference/lifetime coupling and weakens stale-handle detection. |
| Generate bindings from handwritten C declarations | Allows the declarations to drift from exported Rust symbols. |
| Catch panics only in callers | Rust unwinding across a foreign frame is forbidden and cannot be delegated safely. |

## Consequences

**Positive**

- C and C++ consumers get a fixed-width, allocation-neutral, provider-hidden surface.
- Ownership, refusal, panic, concurrency, and buffer behavior are executable tests.
- New ABI versions can coexist under different symbol prefixes.

**Negative / costs**

- The global registry serializes calls in v0.4; independent contexts are thread-safe but not yet parallel.
- Exact B-rep results can currently be classified and destroyed but not serialized through this ABI.
- Every new pointer-bearing export requires an explicit unsafe audit and generated-header update.

**Follow-ups / risks to watch**

- Native packaging and CMake import targets are tracked separately; this ADR does not claim cross-platform binary distribution.
- Benchmark registry contention before replacing the simple lock; do not claim a throughput improvement without evidence.
- Preserve ABI 0.4 symbols when adding future versions.

## Relation to existing code

- `crates/facade/axiolid-capi/`
- `crates/facade/axiolid/src/application.rs`
- `tools/xtask/src/ffi.rs`
- `scripts/check-capi.sh`
