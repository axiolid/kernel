# axiolid/kernel

A standalone, pure-Rust, format-agnostic geometry kernel: neutral geometry
values, explicit operation contracts, and replaceable execution providers.
The `axiolid-model` DAG is the input seam, and applications select providers.
This file is the only `AGENTS.md`. Before editing a crate, read its
`README.md` and the `//!` docs of the module you touch (ADR 0078).

## Layout

- `crates/`: publishable packages, nested by ownership: `foundation/`
  (`axiolid-core`, the dependency root), `representations/`, `contracts/`,
  `algorithms/`, `providers/`, `execution/`, `facade/` (`axiolid`, and
  `axiolid-capi`, the only unsafe C ABI boundary).
- `tools/`: `xtask` (architecture, closure, context, FFI and ledger
  checks), `benchmark`, `oracle`. Local-only, never published.
- `tests/`: black-box consumer fixtures and downstream/native probes.
- `native/`: CMake source-build and installed-package integration.
- `docs/`: the VitePress site. ADRs live in `docs/adr/`.
  `docs/architecture/` holds the machine-checked declarations:
  `closure-profiles.toml` lists minimal downstream closures,
  `capability-ledger.toml` grades every OCCT/CGAL capability against
  Axiolid, and `semver-exceptions.toml` lists accepted breaking changes.
  Its generated maps come from `cargo xtask architecture docs` and are
  never edited by hand.

## Dependency rule

Production direction is
`foundation <- representations <- contracts <- algorithms/providers <- execution <- facade`.
This is a role DAG, not a licence to depend on every earlier layer: each
crate's exact internal edges are allowlisted in its
`[package.metadata.axiolid]`, and `cargo xtask architecture check` enforces
them together with naming (ADR 0064), placement, unsafe policy and
format neutrality. Contracts never depend on providers or dispatch.
Algorithms do not select execution policy. Upward edges are dev-only and
limited to conformance tests. A closure change is an API change: update
`expected_internal` deliberately and record why in an ADR, never to
silence the gate. ADR 0035 owns the package topology.

## Behaviour rules

- No IFC, file-format, vendor, renderer or GPU-API types in `crates/`.
  No C++ dependency path.
- `unsafe` is forbidden everywhere except `axiolid-capi`, which denies
  unsafe operations inside unsafe functions.
- Refuse with a typed error or diagnostic. Never return substitute
  geometry, a guessed default or a silently degraded result. Broad-phase
  candidates are never labelled exact.
- A provider advertises only what it implements. It lands after its typed
  contract, refusal behaviour and conformance suite, and it needs a portable
  scalar correctness oracle before claiming an operation trait. Concrete
  providers stay optional. Never tessellate exact intent silently.
- A capability claim needs an implementation, typed refusal and
  conformance evidence; package metadata is not one. A performance claim
  needs a committed benchmark confirmed in wall clock, not only in
  instruction counts.
- Output order is deterministic. CPU dispatch is chosen at runtime, never
  by `target-cpu=native`.
- Public values implement `Debug` and `Clone`; add other standard traits
  only when they are semantically valid. Split a module before unrelated
  data, validation and algorithms grow together; add no placeholder files.

## Open work

Open work lives in GitHub issues (ADR 0078). A marker in code names its
issue as `TODO(#N)`. Do not add plans, checklists, progress logs or
nested `AGENTS.md` files: `cargo xtask context check` rejects them.
To choose capability work, start from the ledger: `cargo xtask gaps` lists
what is ready by priority, and `cargo xtask gaps show <row|issue>` gives the
evidence and the OCCT/CGAL packages to read. A landing that changes a
row's level edits that row in the same commit. `scoped` rows and
`needs_decision` issues are maintainer decisions: ask, don't start coding.

## Gate

Iterate with focused crate checks, then run the full gate before landing
workspace-wide changes. Judge it by exit code:

```bash
scripts/gate.sh
```

It needs `cargo-semver-checks`, pinned because newer releases need a
newer rustc: `cargo install cargo-semver-checks --version 0.44.0 --locked`.
The main steps can run alone:

```bash
cargo xtask architecture check      # after metadata changes: cargo xtask architecture docs
cargo xtask architecture closure check
cargo xtask context check
cargo xtask gaps check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
scripts/geometry-feature-matrix.sh
scripts/check-capi.sh
```

Mutation-verify a new gate before trusting it: break the rule, watch it
fail, restore (`scripts/probe_*_gate.sh`). Benchmarks are not in the gate;
`tools/benchmark/README.md` says how to run them. Records (changelog,
ADRs, research) follow `docs/guide/contributing.md`.
