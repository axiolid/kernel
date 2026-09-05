# Pointcloud capability (#59 umbrella, #60–#65)

Status: **COMPLETE** — all of #60–#65 landed 2026-09-05.

Architecture check: 52 packages. Layering probe: MUTATION MATRIX PASSED.
Workspace: 1,250 tests / 0 failed with `--all-features`.

## Outcome

- [x] #60 `axiolid-pointcloud` — representation value, 10 tests
- [x] #61 `PointIndex` KNN/radius in `axiolid-spatial` — 13 tests, brute-force verified
- [x] #62 reconstruction contract — request/evidence/refusal + exported conformance suite
- [x] #63 `axiolid-pointcloud-reconstruction-sdf` reference provider — 10 tests
- [x] #64 dispatch registry — 9 tests incl. the dropped-provider mutation
- [x] #65 facade features + gate + measurement table
- [x] ADR 0044 promoted Proposed -> Accepted with as-built notes

## Decisions made during implementation

1. **The provider is not an adopted dependency.** ADR 0044 anticipated
   vendoring a reconstruction library. Instead the provider composes
   `PointIndex` (#61) with `axiolid-levelset` (#87) — both already in-tree
   and audited. No new licensing surface, and the contract is verifiable now
   rather than after a vendoring decision.
2. **Refusal is not a fallback trigger.** A provider saying "this data cannot
   support that request" is an answer about the data, not a provider failure.
   Falling through would search for a provider willing to guess.
3. **Positions-only reconstruction is weaker and says so.** With no normals
   there is no way to determine an inside, so the result is a shrink-wrap
   offset by the sample spacing, reported through `used_normals`.
4. **`matches_device` was promoted to a shared module** rather than copied
   into the new registry, so routing semantics cannot drift per operation.

## Goal

Make the kernel *optionally* capable of pointcloud work: represent point
samples, query them, and reconstruct meshes from them — without admitting
any source-format (LAS/LAZ/E57/PCD/COPC) type into `crates/`.

## Order (dependency-forced)

1. **#60** `axiolid-pointcloud` — representation only, deps = `axiolid-core`.
2. **#61** KNN + radius queries in `axiolid-spatial`
   (`crates/algorithms/query/spatial`), callback-based.
3. **#62** `axiolid-pointcloud-reconstruction-contract` under
   `crates/contracts/operations/pointcloud-reconstruction/`.
4. **#63** provider under `crates/providers/pointcloud/<backend>/`.
5. **#64** dispatch feature in `axiolid-dispatch`.
6. **#65** facade features + architecture gate + docs + size measurement.
7. **#59** close umbrella once 60–65 are closed. Needs an **ADR** recording
   the ingestion boundary.

## Hard constraints

- `axiolid-pointcloud` is representation-only: no topology, no algorithms,
  no source-format types, no deps beyond `axiolid-core`.
- Typed refusal everywhere; never a silent empty result.
- Queries must not allocate on the hot path — callback-based, mirroring the
  existing spatial queries.
- Broad-phase candidates must never be reported as exact adjacency.
- Facade features additive only; **no change to default features**.
- Layering: `algorithms` may NOT depend on `providers`. Cross that seam
  through a contract, as `decompose::split::Splitter` does.

## Verification per issue

- `cargo xtask architecture check`
- `scripts/probe_layering_gate.sh`
- `bash scripts/gate.sh`
- mutation evidence per new capability
- #65 additionally: `default-features = false, features = ["pointcloud"]`
  size measurement recorded in `current-target-crate-map.md`.

## Progress

- [ ] #60 representation
- [ ] #61 spatial queries
- [ ] #62 contract
- [ ] #63 provider
- [ ] #64 dispatch
- [ ] #65 facade + gate
- [ ] ADR: ingestion boundary
- [ ] #59 umbrella closed
