# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `MeshHiddenLinesRegistry` and `MeshVisibleFacesRegistry` route by
  `supports_occlusion`: a request with a see-through triangle skips every
  provider that does not claim `Occlusion::SeeThrough` and is `Unsupported`
  when none does, before the topology audit. A `PerTriangle` occlusion slice
  whose length differs from the triangle count is `InvalidInput` naming the
  mesh, before routing in `hidden_lines` and `visible_faces` and as a listed
  mesh in `audit`, and a visible-face
  group that mixes see-through and opaque triangles is `InvalidInput` naming
  the mesh (ADR 0111, #496).
- `register_conformant` of both registries also runs the see-through suite
  when the provider claims `Occlusion::SeeThrough`, and reports a failure as a
  `ReturnedError` of the case `see-through claim`, so an unverified claim is
  refused (ADR 0111, #496).
- Both registries validate the see-through part of an outcome (ADR 0111,
  #496). A hidden-line segment counts see-through triangles only when it is
  visible and the request has one, neighbouring pieces of one edge differ in
  visibility or count, and `seen_through_segments` and `see_through_triangles`
  must match. A visible-face `behind` list may only follow a see-through
  region, names request meshes, starts at layer 0 or 1 and steps by at most
  one, ascends by mesh position and source within a layer, never repeats the
  region's own source at layer 0, keeps an occluding surface (a cap, or a mesh
  without see-through triangles) in its last layer and a cap to a solid cut
  by the near plane; regions equal in mesh, source and outer ring are ordered
  by `behind`; the entries are counted in `behind_surfaces` and over
  `max_behind_surfaces` refused as `BudgetExceeded` (`visible-face behind
  surfaces`).

## [0.3.13] - 2026-11-13

### Added

- `MeshHiddenLinesRegistry` and `MeshVisibleFacesRegistry` route by
  `supports_projection`: a provider that does not support the view's
  projection is skipped, and a request no provider supports is `Unsupported`
  before the topology audit. A perspective view needs a `near` below zero
  (`InvalidInput` otherwise). `register_perspective_conformant` registers a
  provider after the base suite and the perspective suite. `register_conformant`
  of both registries also runs the perspective suite when the provider claims
  perspective, so an unverified claim is refused (ADR 0107, #387).
- `MeshVisibleFacesRegistry` accepts `FaceGroups::Planar`: it checks that a
  region's group id is a triangle of the mesh and no later than any triangle
  the region lists, and does not recompute the groups (ADR 0106, #384).
- `MeshHiddenLinesRegistry` accepts `EdgeKind::GroupBoundary` segments only
  when `EdgeSelection::group_boundary` is set, and then refuses a
  `PerTriangle` group slice of the wrong length naming the mesh (ADR 0106,
  #385).

## [0.3.12] - 2026-11-03

### Added

- `MeshVisibleFacesRegistry` (feature `mesh-hidden-lines`): registration,
  conformance-gated registration, fallback on `Unsupported`/`Unavailable`
  only, memory admission, and an outcome check (mesh provenance, rings,
  ascending triangles of the region's group, caps only for a solid with a
  near plane, order, limits and evidence). A group slice whose length is not
  the triangle count is refused naming the mesh (ADR 0094, #383).

### Changed

- The shape, topology and memory checks that the hidden-line registry ran
  moved into helpers shared with the visible-face registry; the refusals
  are unchanged.

## [0.3.11] - 2026-11-03

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).
- `MeshHiddenLinesRegistry::audit` lists every input mesh that
  `hidden_lines` would refuse, in request order, without searching (#370).
  It shares the per-mesh audit with `hidden_lines`; refusals of the whole
  request stay `Err`.

## [0.3.10] - 2026-11-03

### Changed

- `MeshPlaneSectionRegistry::section_request` no longer refuses cut regions
  returned beside a `Both` region; it checks that each `Both` region equals
  one cut region, and carries `PlaneSectionEvidence::in_plane_crossings` into
  the evidence it derives from the answer (ADR 0093, #366).

## [0.3.9] - 2026-10-29

### Added

- `MeshPlaneSectionRegistry::section_request` checks cut regions: a provider
  must answer them exactly when asked, with canonical rings, budgets and
  evidence, and no `Both` region; a sheet request with cut regions is refused
  before dispatch (ADR 0092, #372). A provider that ignores the request and
  answers with `cut_regions: None` is treated as `Unsupported`, so the next
  provider is tried and a registry without one refuses.

## [0.3.8] - 2026-10-29

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

