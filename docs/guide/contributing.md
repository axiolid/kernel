# Contributing

Axiolid values small, proven boundaries over broad claims. A change is not complete when it compiles; it is complete when its capability, dependency direction, and failure behavior are evident.

## Before you file anything

[Where things go](../contributing/where-things-go.md) decides the surface: new
functionality and optimizations start as **discussions**, papercuts and bugs are
**issues**, and every decision is recorded rather than silently dropped.

## Licence of contributions

Axiolid is `MPL-2.0 OR LicenseRef-Commercial`. Submitting a contribution means you agree to the [contributor terms](./licensing.md#contributor-terms): you keep your copyright, and the copyright holder may also license your contribution commercially.

## Local checks

```bash
cargo fmt --all -- --check
cargo build --workspace
cargo test --workspace
cargo clippy --workspace --all-targets --all-features -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --workspace --no-deps
scripts/geometry-feature-matrix.sh
scripts/probe_layering_gate.sh
```

The feature matrix protects minimal builds. The probe mutates the declared layering boundary and proves the gate fails, then verifies byte-accurate restoration.

### The gate: tiers, slots and the gate's own probes

`scripts/gate.sh` has two tiers ([kernel#443](https://github.com/axiolid/kernel/issues/443)):

```bash
scripts/gate.sh --branch   # while a branch is in flight
scripts/gate.sh            # full gate: before landing workspace-wide changes
GATE_META=1 scripts/gate.sh   # full gate plus the gate's own mutation probes; the release flow
```

- **`--branch`** runs the formatting, `Cargo.lock`, probe-anchor, architecture, closure, context, docs, roadmap, capability and ledger checks, then build, test (default and all features), clippy (both) and doc for the crates changed since the merge base with `origin/main` and every workspace crate that depends on them, directly or transitively (`scripts/gate-affected.py`, from `cargo metadata`). Uncommitted and untracked files count as changed.
- It **escalates to the full gate by itself** when gate infrastructure changed: `scripts/`, `tools/xtask/`, any `Cargo.toml` or `Cargo.lock`, `rust-toolchain.toml`, `.cargo/`, `crates/facade/axiolid-capi/`, `native/`, `tests/consumers/`, `tests/downstream/`, `.github/` and `docs/architecture/*.toml`. A change it cannot determine (no `origin/main`, a failing `git`) counts as infrastructure.
- The **full gate** is unchanged, except that the gate's own mutation probes (`probe_*_gate*`, the semver, minimal-versions, isolated-build, capi-soname and integration-contract probes) run only when gate infrastructure changed since the merge base, or when `GATE_META=1`. Everything else, including every crate's own `probe_*_mutants.py`, is not a gate step to begin with.
- **Slots.** A gate takes one of `GATE_SLOTS` (default 3) machine-wide slots with `flock` before building, and prints that it waits when all are taken. The lock is released on exit, failure and signal. Gates from different worktrees therefore queue instead of running six or seven at once.
- **Build directory.** Each slot has its own persistent `CARGO_TARGET_DIR`, `<repository>/gate-target/slot-N`, so consecutive gates reuse compiled dependencies and a worktree does not hold a target directory of its own. The gate sets `CARGO_TARGET_DIR` itself; every script under `scripts/` honours it. Slot directories are not cleaned automatically: delete `gate-target/` to reclaim the disk.
- **Log.** The full output also goes to `gate-target/gate-logs/<worktree>-<time>.log`, or to `GATE_LOG`; its path is the first line printed.
- `GATE_DRY=1` replaces every step by a sleep (`GATE_DRY_SLEEP` seconds), for testing the tier and slot logic: `scripts/test_gate_tiers.py` and `scripts/probe_gate_tiers.py` do. `GATE_LOCK_DIR` and `GATE_ROOT` relocate the locks and the gate root.

The `scripts/probe_*_mutants.py` mutation probes share one runner, `scripts/probe_lib.py`: each takes `--list`, `--only NAME`, `--timeout SECS` and `--check-anchors`. A mutant whose anchor text is missing or not unique stops the probe before any test runs, and `scripts/check-probe-anchors.sh` checks every probe's anchors in the gate, so a refactor that moves probed code must re-anchor the probe (or remove the mutant with a stated reason) in the same change.

## Releases

Crates version independently ([ADR 0067](../adr/0067-crates-version-independently.md)): `scripts/prepare-crate-release.py --crate <name> --release <version>` bumps one crate and rolls its changelog, and the release workflow publishes every crate whose version is not on crates.io yet. Before that, `scripts/verify-packages.py` (a gate step) packages every publishable crate and checks the to-be-released ones at their **minimal internal dependency versions**.

A crate that is not re-released keeps the requirements it was published with, so a stale floor passes every workspace build and fails only a consumer whose lockfile holds an older version. `axiolid-mesh-compile` 0.3.28 shipped that way: it used `axiolid-evaluate` 0.3.17 APIs through `axiolid-reference`, whose published 0.3.8 admitted `axiolid-evaluate` 0.3.12 (kernel#337). The check builds each to-be-released crate, with its default features and then all of them, in its own scratch workspace against crates.io plus the other to-be-released archives through `[patch.crates-io]`. It lowers every `axiolid-*` crate in the graph, transitively, to the oldest version the published requirements admit (`cargo +nightly-2026-09-25 -Z minimal-versions update -p ...`); external crates resolve normally. The build uses the workspace toolchain.

A failure names the crate that does not build and the minimal internal versions it was built against. Raise its requirement, or that of the crate re-exporting the missing API, and release that crate in the same release. Run the check alone, or for any crate:

```bash
rustup toolchain install nightly-2026-09-25 --profile minimal   # once
python3 scripts/verify-packages.py --minimal-versions-only --allow-dirty
python3 scripts/verify-packages.py --minimal-versions-only --allow-dirty --crates axiolid-mesh-compile
```

It needs the crates.io index, and it stops rather than guess when the index cannot be read. `scripts/probe_minimal_versions_gate.py` lowers `axiolid-mesh-compile`'s floors back to the 0.3.28 situation and proves the check rejects it, also when the check itself is mutated.

Every publishable crate ships its licence texts ([ADR 0089](../adr/0089-dual-licence-with-a-commercial-option.md), kernel#358). A new crate directory needs relative symlinks `LICENSE` and, unless its `license` is `MPL-2.0` only (as `axiolid-mesh-boolean-boolmesh`, [ADR 0047](../adr/0047-absorb-mesh-boolean.md)), `LICENSE-COMMERCIAL`, both pointing at the repository root files (`ln -s ../../../LICENSE .`, with as many `../` as the depth). Cargo packages the target's contents. Crates with `publish = false` get neither. `scripts/verify-packages.py` fails when an archive lacks a file, ships `LICENSE-COMMERCIAL` without the commercial licence, or holds contents that differ from the root files, for example a symlink a Windows checkout turned into a text file.

## Design rules

- Keep Axiolid format-agnostic: source semantics belong in adapter projects.
- Do not add a concrete backend dependency to a representation or format boundary.
- Treat scalar paths as correctness oracles; benchmark and differentially test a faster path before claiming a performance win.
- Keep operation capability tied to an executable provider implementation.
- Record irreversible architecture choices under [`docs/adr/`](../adr/README.md).
- Put open work in a GitHub issue, not in a checked-in plan. A code marker is written `TODO(#N)`. A crate's `README.md` says what the crate is for and makes no coverage claim that the code and tests do not back ([ADR 0078](../adr/0078-context-lives-beside-the-code.md)).

## Documentation

The site is VitePress and deploys from `main` through GitHub Pages. Preview it locally:

```bash
npm --prefix docs ci
npm --prefix docs run docs:dev
```

### Records

- `CHANGELOG.md` files follow [Keep a Changelog](https://keepachangelog.com/).
- Write an ADR from `docs/adr/_template.md` before an irreversible boundary or dependency change.
- A research record in `docs/research/` keeps measured evidence apart from proposals.
- `cargo xtask docs` generates every page derived from the crates: the [crate reference](/reference/) (from each crate's manifest, `README.md` and `CHANGELOG.md`), the [per-crate changelog](/reference/changelog), and `docs/architecture/crate-map.md` and `dependency-graph.md`. Do not edit them by hand; change the source and regenerate. `cargo xtask docs --check` runs in the gate and before every deploy. Package metadata does not prove a capability: the implementation and its conformance tests do.
- A published crate's `README.md` is its crates.io page. It ends with its docs.rs, reference-page and source links, and states no version or publication status; `cargo xtask docs --check` enforces both.

### Diagrams, geometry, and equations

- Use fenced `mermaid` for architecture, state, and data-flow diagrams. Include
  `accTitle` and `accDescr`, and explain the same conclusion in nearby prose.
- Use fenced ASCII `stl` only for small, reviewable geometry examples. GitHub
  renders the source natively; the VitePress theme parses the same source locally
  with an event-driven viewer. Keep capability evidence in tests, not in the
  picture.
- Use `$...$` and `$$...$$` for equations, define every symbol, and state the
  assumptions under which the equation applies.
- Keep code, benchmark tables, file trees, and protocol records as code blocks;
  Mermaid is not a replacement for literal data.

See GitHub's primary documentation for
[Mermaid diagrams](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/creating-diagrams),
[STL models](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/creating-diagrams#creating-stl-3d-models),
and [mathematical expressions](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
