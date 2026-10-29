# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `MeshPlaneSectionRegistry::section_request` accepts open meshes for
  `SectionRole::Sheet`: boundary edges are allowed, degenerate, non-manifold
  and inconsistently wound input is refused before dispatch, the evidence must
  echo the role and the audited `boundary_edges`, and a solid answer must not
  end a chain on the mesh boundary, nor may any answer for a mesh without
  boundary edges (ADR 0091, #359).

## [0.3.7] - 2026-10-28

### Added

- `MeshPlaneSectionRegistry::section_request`: the source checks of `section`,
  fallback only on unsupported or unavailable providers, and outcome checks of
  the frame, chain and ring invariants, canonical order, budgets (in-plane
  regions, region vertices, in-plane pairs) and evidence counts (ADR 0091,
  #359).
- `MeshPlaneSectionRegistry::register_request_conformant` also runs the
  region suite and rejects a claimed failure (#359).

## [0.3.6] - 2026-10-28

### Fixed

- `MeshPlaneSectionRegistry` refuses a mesh whose index count is not a
  multiple of three with `InvalidInput`; it used to drop the partial
  triangle and could accept the rest as a closed solid.

### Added

- `MeshHiddenLinesRegistry` behind the `mesh-hidden-lines` feature: shape,
  budget and topology checks before dispatch, outcome validation, fallback
  only on unsupported or unavailable providers, and refusals that name the
  mesh at fault (ADR 0090, #138). `axiolid-hidden-lines` is the first
  provider behind it. Adjacent pieces of one source edge must differ in
  visibility when hidden pieces are requested.

## [0.3.5] - 2026-10-28

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
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

