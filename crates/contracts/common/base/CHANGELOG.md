# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

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
