# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- The `MESH_HIDDEN_LINES` capability identifier
  (`org.axiolid.geometry.mesh-hidden-lines.v1`) and the
  `Operation::MeshHiddenLines` diagnostic variant, for the hidden-line
  drawing contract (ADR 0090, #138).

## [0.3.4] - 2026-10-28

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.3] - 2026-10-04

### Added

- `ScratchRequirement::Affine { base_bytes, bytes_per_element,
  bytes_per_worker }` (#226): at most `base + per_element * elements +
  per_worker * workers`. A purely per-element bound cannot be an upper
  bound for an operation with fixed setup cost or one that runs on a
  thread pool, whose workers each cost memory however small the input is.
  Additive: the enum is `#[non_exhaustive]`.
- `ScratchRequirement::upper_bound_bytes_on(elements, workers)` and
  `ScratchRequirement::fits_budget_on(options, elements, workers)`, for a
  caller that knows the width the operation runs on.
- `Parallelism::worker_bound`: the most workers a preference allows
  (`Serial` 1, `Threads(n)` n, `Auto` the available parallelism).

### Changed

- `ScratchRequirement::fits_budget` charges a per-worker term for
  `Parallelism::worker_bound` of the options. `upper_bound_bytes` reports
  `None` for an affine requirement with a per-worker term, which has no
  bound without a worker count. Existing variants are unchanged.

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-25

### Added

- `ExecutionOptions::with_chord_error` and `ExecutionOptions::chord_error`
  (#165): an explicit bound on how far a provider's straight chords may sit
  from the curve they replace, separate from the linear tolerance. The
  tolerance is a coincidence test; used as a chord budget it leaves a 5 mm
  arc a few chords at `Tolerance::MILLIMETRE`, and small profiles mesh
  percent-level off. `None` (the default) keeps each provider's previous
  behaviour. A non-finite or non-positive budget is refused (`None`).
