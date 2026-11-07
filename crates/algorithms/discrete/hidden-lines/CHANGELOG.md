# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- `ScalarVisibleFaces` no longer refuses a scene as `Degenerate` when exact
  ring vertices round to the same binary64 value (#410, ADR 0100). The
  regions are the exact regions snap-rounded to the binary64 grid: an edge
  that passes through the rounding cell of another vertex is routed through
  it, and a region or hole whose area vanishes is dropped and counted. The
  evidence reports `snap_displacement`, `rerouted_edges`, `collapsed_regions`
  and `collapsed_holes`; where nothing is rerouted the vertices are still the
  exact points rounded once. The two seeded oracle scenes that were refused
  (seed 1 oblique, seed 9 isometric, each with a near plane) are answered and
  pinned, with frozen copies in the conformance suite.

### Changed

- Hidden-line drawing does less exact work for the same output (#369): an
  edge stops looking at further candidates, and skips the near-plane caps,
  once one hidden span covers all of its kept part; a sign that a repeated
  vertex id makes exactly zero (`orient2d` with a shared vertex, a point on
  its own triangle's plane, a zero coordinate difference) is not evaluated;
  and a triangle whose plane is certainly not in front of the edge at either
  end, by the interval filter alone, is skipped before its three image
  clips. Visible face regions share this resolver, so their
  `candidate_pairs` and `exact_escalations` evidence and their candidate-pair
  budget fall the same way; the earlier visible-face timing figures predate
  this change. The drawing's segment hash is unchanged on every
  `scenario/hidden-lines` case. Examined candidate pairs and exact
  escalations fall, as the early exit stops counting pairs: isometric 20x20
  125,577 to 79,064 pairs and 187,483 to 43,402 escalations; cut plan 20x20
  23,680 to 14,240 pairs and 113,440 to 51,200 escalations; isometric 50x50
  829,935 to 511,016 pairs and 1,215,313 to 287,386 escalations; cut plan
  50x50 148,000 to 89,000 pairs and 709,000 to 320,000 escalations. A loaded
  machine, three runs each (the BVH provider, ranges of the three
  estimates): isometric 20x20 0.73-0.85 s to 0.35-0.47 s, cut plan 20x20
  0.17-0.48 s to 0.08-0.11 s, isometric 50x50 3.9-6.2 s to 1.9-2.6 s, cut
  plan 50x50 1.0-2.3 s to 0.51-0.79 s. Single runs vary by more than the
  gain on some cases, so read the work counts as the firm result. Because
  fewer pairs are examined, a `max_candidate_pairs` limit that the old
  count reached may no longer be reached.

## [0.1.3] - 2026-11-03

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).
- `ScalarVisibleFaces`, the exact provider for `MeshVisibleFaces` (ADR 0094,
  #383): the visible face regions of caller-grouped triangles in the view of
  the hidden-line provider, decided in the exact planar arrangement of
  `axiolid-exact-arrangement`, with region vertices that are exact points
  rounded once. Checked by the faces conformance suite, an independent
  binary64 ray-casting oracle and a mutation probe. On
  `scenario/visible-faces` (a loaded machine, ranges over three runs) the BVH
  beats a full scan on a cut plan by 1.4-2.1x at 400 boxes and 4.1-6.6x at
  2,500 boxes, and the depth prune of hidden edges saves 2.7-3.2x in an
  isometric view at 400 boxes; the BVH is not consistently faster than the
  scan in an isometric view or on interpenetrating boxes.

### Changed

- Request validation, the mesh audit and the building of view-space
  vertices, triangles and near-plane caps moved into one `prepare` step that
  the visible-face provider will share (ADR 0094, #383). The drawings and
  the refusals are unchanged, except that a crease cosine outside `[-1, 1]`
  is now refused before the other request checks, not after the depth range.
- The rounding of an exact ratio to binary64 moved into
  `axiolid-exact-arrangement` (`round_ratio`), so that the endpoints of
  hidden lines and the vertices of visible-face regions round alike. The
  drawings are unchanged.

## [0.1.2] - 2026-10-29

### Changed

- Candidate triangles for an edge, and the cap triangles of a cut solid, come from an `axiolid-spatial` BVH over view-space boxes instead of a linear scan. The predicate, the candidate set, its ascending order and `candidate_pairs` are unchanged (#138). Wall clock from `scenario/hidden-lines` (`cargo bench -p axiolid-benchmark --bench scenario -- hidden-lines`, which times the BVH and the linear scan, `ScanHiddenLines`, on each case; a loaded shared machine, so read these as ranges), linear scan / BVH: isometric 20x20 boxes 412 / 389 ms (a reviewer's run: 459 / 426), cut plan 20x20 174 / 140 ms (117 / 100; very noisy), isometric 50x50 3.34 / 2.52 s (3.62 / 2.88), cut plan 50x50 1.59 / 0.76 s (1.67 / 0.82). The gain is small at 400 boxes, because exact evaluation of the candidates dominates there, and 25-30% (isometric) to about 50% (cut plan) at 2,500.

## [0.1.1] - 2026-10-28

### Changed

- The near-plane caps of closed solids occlude (#138, ADR 0090): a
  `Solid` that the near plane cuts no longer fails as `UnsupportedInput`.
  A point with `W < near` (farther than the plane) is hidden when the point
  of the plane above it lies in the solid, by exact ray parity with a
  symbolic perturbation of the image point. A face of a cut solid in the
  near plane is drawn, not refused.

## [0.1.0] - 2026-10-28

### Added

- `ScalarHiddenLines`, the exact scalar provider for `MeshHiddenLines`
  (ADR 0090, #138): boundary, silhouette and crease edges of triangle
  meshes split into visible and hidden pieces by exact interval and dyadic
  decisions, with endpoints that are the exact points rounded once. A near
  clipping plane that cuts a `Solid` mesh is refused as `UnsupportedInput`
  naming it.
