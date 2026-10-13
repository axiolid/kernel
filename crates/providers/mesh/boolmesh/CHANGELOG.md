# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.5] - 2026-10-13

### Fixed

- A memory budget holds for `union_many` and `subtract_many` (#259). The
  declared bound covers one boolean in its operands' triangles; a batch has
  no such bound in its input, because union intermediates outgrow their
  operands and a cut subject grows with every group, and with
  `parallel-batch` a worker waiting inside one boolean's join steals other
  pairs, so every pair of a level can be in flight on fewer workers (5 on
  4 workers was observed). Each step is now charged when its operands
  exist: the declared bound for every boolean that can be live (all pairs
  of a `parallel-batch` level at once, one otherwise, the per-worker term
  once), plus what the batch holds (its copies, the previous level's
  results, the cut subject, a fused tool, its node and index lists). The
  first step that does not fit is refused with `BudgetExceeded` before it
  starts. Before, `union_many` was charged once for its input triangles
  and `subtract_many` not at all beyond each boolean's own check, which
  only ran with `parallel`.
- A boolean given a budget directly, without dispatch, is checked in every
  build, not only with `parallel`.

### Changed

- `union_many` borrows the solids that need no channel conformance and
  moves an odd trailing solid to the next level instead of copying it, and
  `subtract_many` borrows the subject until the first group cuts it: each
  copy was scratch held through the batch.
- The per-worker term no longer includes the 64 KiB base with
  `parallel-batch`. That was the batch's concurrency charged per worker,
  which stolen work can exceed; the batch now charges it per pair. With
  `parallel-batch` alone the term is the 16 KiB of a worker's bookkeeping.
- The `scratch_bound` test and the `scratch_probe` binary measure both
  batch paths at 14 workloads that grow (overlapping boxes, a slab cut by
  overlapping tools, at several sizes, with and without an attribute
  channel) on 1 to 64 workers under the forced schedule: each run at the
  smallest budget it admits stays within it (at most 73% of it), the run a
  byte below is refused within it, and a 1 KiB budget is refused before
  any copy is made. A window now counts only its own allocations: memory
  of an earlier pool freed inside it (rayon retires deque buffers through
  an epoch collector) hid up to 7 KB of a peak.
- `tests/batch_budget.rs` pins what each step is charged, and
  `scripts/probe_scratch_bound_mutants.py` kills a charge that drops any
  of its terms.

## [0.3.4] - 2026-10-04

### Fixed

- The declared scratch is an upper bound on a pool as well (#226). The old
  `PerElement { bytes_per_element: 4096 }` covered large inputs but not
  small ones on rayon: each worker allocates its own bookkeeping when it
  starts, and on a loaded machine that lands inside whichever boolean is
  running. 24 triangles peaked at 113,608 bytes against the declared
  98,304 in one gate run. The provider now declares
  `ScratchRequirement::Affine` with a 64 KiB base, 1,536 bytes per input
  triangle of all operands, and 16 KiB per worker with `parallel` (plus
  the base per worker with `parallel-batch`, where each worker can run a
  boolean of its own; zero without either feature). Each term is 1.7x to
  2.9x its measurement; the doc comment of `scratch_requirement` has the
  table.
- With a memory budget, a boolean (and a `parallel-batch` `union_many`)
  re-checks that bound against the width of the rayon pool it runs in and
  refuses with `BudgetExceeded` when the pool is wider than the budget
  check before dispatch assumed.

### Changed

- `tests/scratch_bound.rs` and the `scratch_probe` binary measure under a
  forced worst-case schedule instead of whatever the machine's load gives:
  with rayon each boolean runs in a pool of 1, 2, 4, 16 or 64 workers
  started inside the measured window, with allocations stalled so join
  halves are stolen. The test checks each peak against the bound for its
  worker count, and checks that worker start-up is really inside the
  window. Under the old declaration it now fails every run instead of
  occasionally under load.
- `scripts/probe_scratch_bound_mutants.py` kills a dropped base term, a
  dropped per-worker term, a worker count ignored in the contract, the
  dispatcher or the provider, a boolean budgeted for its subject alone,
  and the old declaration.

## [0.3.3] - 2026-09-30

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
