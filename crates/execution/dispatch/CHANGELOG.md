# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.4] - 2026-10-14

### Fixed

- The mesh-section registry uses the shared device matching instead of its
  own copy, so a build with only the `mesh-section` feature no longer warns
  that the shared one is unused, and the routing rules cannot drift apart.

## [0.3.3] - 2026-10-13

### Fixed

- A one-thread context from `with_execution` is charged one worker and
  runs the provider on a pool of one (#259). It installed no pool, so the
  call ran on whatever pool the caller was in, and the budget fell back to
  the options' `Parallelism`.

## [0.3.2] - 2026-10-04

### Fixed

- `MeshBooleanRegistry::boolean` budgets both operands (#226). It counted
  the subject's triangles alone, so a provider's declared scratch was
  checked at up to half the input it would really be given; `subtract_many`
  and `union_many` already counted every operand.
- A per-worker scratch term (`ScratchRequirement::Affine`) is charged for
  the width of the pool configured with `with_execution` when it has more
  than one thread, since every provider call runs inside it, and otherwise
  for the options' `Parallelism`.

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

