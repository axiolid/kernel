# 0105 — Drawings join the `discrete` and `application` bundles

- **Status:** Accepted
- **Date:** 2026-11-07
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0090](./0090-mesh-hidden-line-drawings.md), [0040](./0040-versioned-c-abi-boundary.md)

## Context

Issue #371. ADR 0090 kept `mesh-hidden-lines`, `dispatch-mesh-hidden-lines`
and `hidden-lines-provider` out of `discrete`, `application` and `full`, and
left open whether the application boundary and the C ABI carry hidden lines
"for when a provider exists". Two providers exist now: `ScalarHiddenLines`
(ADR 0090) and `ScalarVisibleFaces` (ADR 0094). A native consumer, which
reaches Axiolid only through `axiolid::application` and the C ABI over it
(ADR 0040), cannot draw while they are opt-in, and the C ABI cannot expose
them: its header must equal its library for every build
(`scripts/check-capi.sh`, ADR 0041), so it cannot change with Cargo features.

Measured on origin/main `e1d941f4`: the two crates and the contract add
exactly `axiolid-hidden-lines` and `axiolid-mesh-hidden-lines-contract` to
the closures of `rust-facade-application` and `c-abi-profile`. No third-party
crate is added (`num-bigint`, `exact-arrangement` and `spatial` are already
there). The face suite takes about 0.45 s, the line suite about 12 ms;
`Application::portable()` takes about 1 ms and the C ABI calls it on every
`context_create`.

## Decision

**Hidden lines and visible faces join the bundles together.**

- `discrete` gains `mesh-hidden-lines` (the contract). `application` gains
  `dispatch-mesh-hidden-lines` and `hidden-lines-provider`. `full` therefore
  gains the contract only. The leaf features stay, so a Rust caller that only
  draws still avoids `application`.
- `Application` gains `hidden_lines`, `audit_hidden_lines` and
  `visible_faces`, and `ApplicationBuilder` gains
  `with_portable_hidden_lines` and `with_portable_visible_faces`;
  `Application::portable()` selects both. Each failure carries a `CallContext`
  with the operation, the selected provider (`axiolid-application` when none
  is selected) and the tolerance. A refusal that names a mesh keeps it as data
  in `ApplicationErrorSource::HiddenLines`.
- The descriptor advertises `MESH_HIDDEN_LINES` (`scalar-hidden-lines`) and
  `MESH_VISIBLE_FACES` (`scalar-visible-faces`) with output
  `Representation::Profile2d` and `Exactness::ToleranceBounded`: the output is
  binary64 and collapsed pieces are counted, not exact. A new
  `Representation::Drawing2d` would break `axiolid-contracts`; it waits for
  that crate's next breaking release.
- **Conformance is verified lazily, once per process.** The two builders do
  not run the suites (unlike `register_conformant` for booleans and sections);
  the first call that reaches a provider does, and keeps the verdict in a
  `OnceLock`. The providers are stateless and the suites deterministic, so the
  verdict is the same for every application. The one-time suite runs outside
  the caller's memory budget and cancellation: the call checks cancellation
  first, then pays the suite once. A failed verdict refuses every
  call with `ApplicationErrorSource::ProviderConformance` and no output.
  Verifying eagerly inside `portable()` would make context creation about 450
  times slower.
- `ApplicationErrorSource` and `ProviderConformanceReport` become
  `#[non_exhaustive]` and gain the new variants. This is breaking for a
  `match` on either enum; it rides the pending, unpublished `axiolid` 0.5.0,
  where it is free.

**The joining rule for later operations.** Later variants of these operations
(perspective, oblique clips) are request options, not features: they reach
`application` as soon as the bundled provider claims them. A new operation,
such as exact B-rep hidden lines in its own provider crate, starts as an
opt-in leaf feature and joins `application` by an amendment of this ADR once
it has a conformant portable provider and a reviewed closure diff. C ABI
symbols are never feature-gated.

## Alternatives considered

| Option | Why not |
| --- | --- |
| An `application-drawing` feature and a Cargo feature on `axiolid-capi` | Creates a matrix of `Application` shapes and C libraries; the header no longer matches every build. |
| Stay opt-in | Native consumers could never draw, and the C ABI for faces (#386) stays blocked. |
| Verify eagerly in `with_portable_*` | The face suite costs about 0.45 s on every `context_create`. |
| Skip the suites | A provider is advertised only with conformance evidence. |
| Add `Representation::Drawing2d` | Breaking for `axiolid-contracts` for no information the exactness label does not carry. |

## Consequences

**Positive**

- A Rust application and, in a later step, a C consumer draw through one
  boundary with typed, mesh-naming refusals.
- Leaf features keep the pay-for-what-you-use route.

**Negative / costs**

- The closures of `rust-facade-application` and `c-abi-profile` grow by two
  crates; `full` by one.
- The first drawing call in a process pays the suite once (about 0.5 s for
  faces).
- A `match` on the two enums needs a wildcard arm.

**Follow-ups / risks to watch**

- The C ABI for hidden lines and visible faces (ADR 0040 scope), a separate
  step.
- Exact B-rep and smooth-silhouette drawings join by amendment of this ADR.

## Relation to existing code

- `crates/facade/axiolid/Cargo.toml`: features `discrete`, `application`.
- `crates/facade/axiolid/src/application/drawing.rs`: builder calls, methods,
  descriptor rows and the `verified` helper.
- `crates/facade/axiolid/src/application.rs`: the two non-exhaustive enums.
- `docs/architecture/closure-profiles.toml`: `rust-facade-application`,
  `c-abi-profile`, `full`.
