# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.1.5] - 2026-11-14

### Added

- `ScalarHiddenLines` and `ScanHiddenLines` draw perspective views and claim
  them through `supports_projection`, replacing the by-name refusal of 0.1.4:
  a pinhole at the frame origin, image `(U / -W, V / -W)`, `near` required
  and below zero. Every visibility predicate is an exact determinant of
  vectors from the eye, so triangles behind the eye need no case, endpoints
  are the exact rationals rounded once, and the near-plane cap is decided by
  ray parity along the eye ray. The perspective conformance fixtures and a
  binary64 eye-ray oracle (seeded scenes, eyes outside and inside) check it;
  the provider refuses a perspective view whose `near` is missing, zero,
  positive or not finite with `InvalidInput`. The visible-face providers
  still refuse a perspective view by name (ADR 0107, #387).
- The broad phase scales the depth coordinate of its hierarchy to the size of
  the image in a perspective view, where the image is in units of one over the
  depth and the longest-extent split would otherwise cut by depth alone. On
  the 50x50 perspective grid of `scenario/hidden-lines` the tree then returns
  4.2M boxes instead of 44.3M before the exact test (540,676 pairs either
  way; the isometric grid returns 3.5M). The drawing is unchanged and
  orthographic views are untouched. A vertex at or behind the eye plane
  (`W >= 0`) has an unbounded image enclosure, so a triangle beside or behind
  the eye sends the whole scene to the linear scan (same pairs, slower
  candidate search); this is documented in the crate docs (#387).
- An edge-on perspective triangle whose plane contains the eye may contain the
  eye itself, strictly inside, on a side or at a vertex. The eye in the closed
  triangle is decided exactly (no two of the three signs of `v_i × v_j` along
  the normal strictly opposite) and every ray from the eye then starts in the
  triangle, so its fibre reaches from the eye to the side the ray leaves by:
  the span of each non-skipped side counts as in front without the nearer-than
  test, while the near plane still has to be reached and the point still has
  to be strictly beyond it (ADR 0107, #387).
- `ScalarHiddenLines` and `ScanHiddenLines` honour see-through triangles and
  claim them through `supports_occlusion` (ADR 0111, #496). A see-through
  triangle hides nothing; a visible piece counts in `seen_through` the
  see-through triangles of nonzero image area whose plane is strictly in
  front of it and whose image holds the point shifted by `(ε, ε²)`, decided
  with the perturbed image-side clips of the caps, in orthographic and
  perspective views. A visible gap is split where that count changes and a
  hidden piece is not split. Without a see-through triangle no step of this
  runs: the drawing, the evidence and the exact work are unchanged. A
  `PerTriangle` occlusion slice of the wrong length, or an unknown occlusion
  mode, is `InvalidInput` naming the mesh.
- `ScalarVisibleFaces`, `ScanVisibleFaces` and `UnprunedVisibleFaces` honour
  see-through triangles and claim them (ADR 0111, #496). A side whose front
  holds a see-through triangle gets its stack: the other candidates grouped
  into exact depth layers down to the first one holding an opaque triangle or
  a cap, found by one linear pass for the occluder and a binary insertion for
  the layers between; each of these extra depth comparisons counts against
  `max_candidate_pairs`. A region is a maximal region of one source and one
  stack, `VisibleRegion::behind` lists the stack without the region's own
  layer-0 entry, the output breaks a tie of mesh, source and outer ring by
  `behind`, and `tied_sides` also counts a side with a tie behind its front.
  The entries are bounded by `max_behind_surfaces` (`BudgetExceeded`,
  `visible-face behind surfaces`). The prune keeps every edge behind a
  see-through triangle, since it hides nothing. A face group that mixes
  see-through and opaque triangles is `InvalidInput` naming the mesh, and
  `FaceGroups::Planar` joins two triangles only when their occlusion is
  equal. Without a see-through triangle no extra comparison is made and the
  regions and evidence are unchanged.
- `scenario/hidden-lines` and `scenario/visible-faces` gain a
  `glazed-isometric` grid, the isometric one with every third box wholly
  see-through, and print a drawing hash that leaves out an empty count or
  stack, so the opaque scenes print what they printed before (#496); they
  did, hash, escalations and pairs alike. Three runs each on a loaded machine
  (median of each run, BVH provider, against the opaque isometric grid of the
  same size): hidden lines 400 boxes 0.50-0.52 s against 0.27-0.39 s (79,064
  to 94,509 pairs), 2,500 boxes 2.33-3.29 s against 2.50-2.51 s (511,016 to
  617,263 pairs); visible faces 100 boxes 0.34-0.36 s against 0.23-0.30 s
  (1,266 to 2,001 pieces, 142 surfaces behind), 400 boxes 1.18-1.92 s against
  1.25-1.30 s (4,938 to 11,097 pieces, 383,574 to 789,043 pairs, 1,931
  surfaces behind). Edges behind glass lose the early-out of full cover and
  are no longer pruned from the arrangement, which the counts show; the run
  to run spread on this machine is as large as the cost.

## [0.1.4] - 2026-11-13

### Added

- `ScalarVisibleFaces` resolves `FaceGroups::Planar` exactly before labelling
  (new `groups` module: union-find over edges of exactly two triangles with
  parallel, same-facing normals, polled for cancellation, memory-budgeted;
  ADR 0106, #384).
- `ScalarHiddenLines` draws `EdgeKind::GroupBoundary` when
  `EdgeSelection::group_boundary` is set, resolving the groups only then;
  with it off the output is unchanged (ADR 0106, #385).

- The hidden-line and visible-face providers refuse a perspective view by
  name (`UnsupportedInput`, input `perspective view`) when called directly,
  and do not claim it through `supports_projection` (ADR 0107, #387).

### Changed

- The providers' view-space vertices carry exact homogeneous coordinates
  `[U, V, W, h]` with `h = 1` in an orthographic view, preparing perspective
  (ADR 0107, #387). The drawings are bit for bit the same: the segment hash
  of every `scenario/hidden-lines` case is unchanged, and the wall clock
  there is unchanged (isometric and cut plan at 20x20 and 50x50 agree with
  the previous commit within run-to-run noise).

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

- Request validation, the mesh audit and the building of view-space
  vertices, triangles and near-plane caps moved into one `prepare` step that
  the visible-face provider will share (ADR 0094, #383). The drawings and
  the refusals are unchanged, except that a crease cosine outside `[-1, 1]`
  is now refused before the other request checks, not after the depth range.
- The rounding of an exact ratio to binary64 moved into
  `axiolid-exact-arrangement` (`round_ratio`), so that the endpoints of
  hidden lines and the vertices of visible-face regions round alike. The
  drawings are unchanged.

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
