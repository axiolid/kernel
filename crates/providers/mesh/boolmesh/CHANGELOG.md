# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- Unions of overlapping axis-aligned boxes no longer refuse inside the
  solve with an odd edge-point count (#203). The winding number's xy
  broad phase rejected queries beyond `min + cell * dim`, a rounded
  product that can fall one ulp short of the operand's true bounding box,
  so a vertex lying exactly on the other operand's extreme plane lost a
  face from its winding number. The grid now rejects against the exact
  bounding box. Grid unions at pitches 0.6 to 0.8 with k = 5 and 6
  complete, sequentially, through `union_many` and on the fast winding
  path, with the exact volume and a closed, consistently wound result
  (`tests/overlapping_grid.rs`).

### Changed

- `tests/solve_failure.rs` is replaced: no admissible input is known to
  reach a refusal inside the solve now, so the mapping to
  `BackendContractViolation` is pinned by a unit test where
  `compute_boolean`'s error lands.

## [0.3.2] - 2026-09-28

### Changed

- **Behaviour change:** a refusal inside the solve (an odd edge-point
  count in `pair_up`, #101) is reported as
  `GeomError::BackendContractViolation` naming `boolmesh`, not
  `Degenerate`. The operands passed every input gate, so the failure is
  this provider's defect and no longer reads as the caller's.
  `tests/solve_failure.rs` pins it on a grid union that still reaches the
  refusal (#203).
- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Fixed

- The scratch probe discards one warmup boolean before measuring, and
  measures peaks above the bytes already live, so the first operation is
  no longer charged for process startup (#110). A `scratch_bound` test
  fails if any measured peak exceeds the declared 4 KiB per triangle.

## [0.3.0] - 2026-09-23

### Added

- The pairwise boolean carries attribute channels (#116). Each result triangle's source triangle is tracked through the CSG core, so a corner that is a source corner copies its value and a corner on a cut is derived in its source triangle under the channel's `Blend`. Output is corner-indexed; faces from a tool without the channel are `UNMAPPED`. Fates: `Preserved` when nothing was derived, `Interpolated` otherwise, `Dropped(NotBlendable)` for a `Blend::None` channel a cut needs. Every path carries channels (#116, completed): the analytic box path (sources recovered by plane lookup), grouped `subtract_many` (fused cutters keep theirs), tree `union_many` (each solid conformed to the first's channel set) and the empty result.

### Fixed

- A result corner was sampled in its RECORDED source triangle, but simplification merges coplanar faces, so a result face can span several source triangles: corners were extrapolated (weights down to -0.65), wrong for piecewise data such as atlas UVs or per-face ids. Each corner is now located in the coplanar source region by a point just inside its face, which also picks the right side of a seam.
- `dedupe_edge` pushed per-face data indexed by a vertex id when pinching a vertex, desynchronising face normals and provenance from the faces (inherited from upstream; Manifold pushes only per-vertex data there). The `ProviderLimitation` comment claiming the boolean returns positions only, stale since ADR 0047, is corrected.
