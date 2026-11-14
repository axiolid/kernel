# 0109 — Free faces close gaps with reported bridges on the exact arrangement

- **Status:** Accepted
- **Date:** 2026-11-14
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #490, from rheonda/cad#568: a floor plan is cut from a mesh into wall
pieces whose ends miss each other by float noise, about 1e-5 m, and the rooms
between them must still close. The caller gives obstacle polygons (walls) and
barriers (open segments such as a door's line) and wants the bounded free
faces, each ring edge with the inputs it lies on, and a point-in-face query.

The issue proposed snapping boundaries closer than a width `gap` together. That
moves input vertices by up to the width, which at 1e-4 m is not an ulp. A
cluster of close vertices is transitive and depends on processing order, a
moved edge can create crossings that were not there, and the result is no
longer a function of the inputs as given.

`ArcArrangement` (ADR 0104's machinery) takes closed rings only, reads
membership by ring winding, has no open segments and no point location, and
`exact_arc/arrangement.rs` is being rewritten in flight (#344). The editable
`axiolid-arrangement` is built from one polygon. `axiolid-exact-arrangement`
has segments, crossings, tracing and snap rounding, which the visible faces
already use (ADR 0094, ADR 0100), but no query for the faces of the complement.

## Decision

`free_faces` runs one pipeline on `axiolid-exact-arrangement`:

1. Arrange the obstacle edges and the barriers exactly.
2. Close gaps by adding **bridges**: exact straight segments from a vertex to
   the nearest point of a piece it is not on. The pair is bridged when that
   point is strictly nearer than `gap` and the open segment meets no piece.
   The minimum distance of two segments that do not cross is reached at an end
   of one of them, so every narrowest place of a gap is bridged. Every
   qualifying pair is bridged, so no choice is made and nothing is ambiguous.
   Input geometry is never moved. Two edges that meet at a vertex are never
   bridged to each other there: a sharp corner stays sharp.
3. Arrange inputs and bridges, number the faces of the complement (`faces`), and
   read each face as obstacle or free by ring membership (ADR 0104). A bridge
   inside an obstacle is discarded.
4. A bounded free face that fits in a strip narrower than `gap` is a **closed
   gap**: for some edge of the convex hull of its exact ring, every hull vertex
   is strictly nearer than `gap` to that edge's line. It is reported, not
   returned as a face.
5. Rings are traced, snap-rounded and certified as ADR 0100 does for visible
   faces, so the evidence bounds how far the binary64 ring moved.
6. Every ring edge carries the list of inputs and bridges it lies on, since
   coincident inputs are legal.

Details that define the contract:

- `gap` is strict (`<`), `0` is valid and closes nothing, and there is no upper
  bound: the memory budget (ADR 0088) and cancellation bound the work, and a
  guessed constant would be a guessed default. The memory budget bounds what
  grows with crossings, T-junctions and bridges: the cuts of both
  arrangements on top of 256 bytes per input segment or bridge, each exact
  foot and each bridge the search keeps (256 bytes each), and the cells and
  fragments of snap rounding; beyond it the operation is refused with
  `BudgetExceeded` for `memory`. It does not bound the arrangement's vertices
  and pieces, the faces, the rings or the output, which are linear in the size
  of the arrangement, nor the candidate pairs the search examines, which are
  quadratic in the worst case (a comb of staggered parallel barriers) and are
  bounded only through the feet and bridges they keep, and by cancellation.
- A ring on one line is refused (`ZeroArea`), although ADR 0104 leaves it out
  of a region: a wall that vanished silently would make a room leak. The caller
  can pass the line as barriers.
- A piece with the same face on both sides, a dangling barrier or bridge, is in
  no ring. A point on one is `Inside` the face, since the piece bounds nothing.
- Each bridge is a `Bridge` record with its rounded ends, the inputs it joins
  and its length, so every closure is visible and a consumer can drop the edges
  whose only source is a bridge.
- Output is in the repository's order (`point_order`, `canonical_start`,
  `sequence_order`: `x`, then `y`) and depends on the input order only through
  the indices inside `PieceSource`.
- `locate` is exact against the returned binary64 rings and the input
  obstacles: on a ring edge, inside a face, in a closed gap, in an obstacle,
  else unbounded. A point within the snap displacement of a boundary may be
  answered as either neighbour.
- The new edge `axiolid-overlay -> axiolid-exact-arrangement` (which brings
  `axiolid-predicates`). Every closure profile already lists both crates, so no
  `expected_internal` changes. `faces` and the distance predicates live in
  `axiolid-exact-arrangement`; `free_faces` and `locate` in `axiolid-overlay`.
- No C ABI in this change, and no performance claim: the only growth evidence is
  a count of candidate pairs in the tests, not a wall-clock measurement.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Snap boundaries closer than `gap` together (the issue's wording) | Moves input vertices by up to `gap`. Clusters are transitive and depend on processing order, a moved edge can create new crossings, and outputs are no longer bit-identical to inputs: substitute geometry. |
| Morphological closing (dilate then erode by `gap/2`) | Rounds corners, produces arcs and irrational points, and loses the per-edge provenance. |
| Snap rounding on a coarse grid of pitch `gap` | Moves every vertex, does not guarantee closing (two points 1e-5 apart can sit in different cells), and was rejected in ADR 0100. |
| Build on overlay's `exact_arc` core or `ArcArrangement` | Closed rings only, membership by ring winding, no open segments, no point location, and #344 is rewriting it. |
| `axiolid-arrangement` (editable DCEL) | Built from one polygon; no exact multi-segment build, no holes, no snap rounding. |
| Bridge only the nearest pair per vertex | A choice that can leak (the second-nearest piece may bound another neck) and creates ties. Bridging every qualifying pair needs no choice. |
| A new crate for free faces | A topology change (ADR 0035) for one function whose dependencies every closure profile already carries. |

## Consequences

**Positive**

- Inputs are never moved, and the output is independent of the input order.
- Every closure is a record with its ends, the inputs it joins and its length.
- A room that leaks, its walls farther apart than `gap`, shows as the absent
  face and is never guessed closed.
- `faces` and the distance predicates are available to other users of the exact
  arrangement.

**Negative / costs**

- A room ring has a short bridge edge at each open joint instead of a clean
  corner, so its area differs from the nominal room by about the gap times the
  joints' length. A consumer that wants the clean rectangle drops the edges whose
  only source is a bridge.
- A redundant collinear vertex within `gap` of an adjacent edge cuts an acute
  tip off the face and reports it as a closed gap.
- The strip rule measures the width of a cell, not the clearance of the wall
  that made it. At a corner of angle `theta` under 90 degrees that a wall wraps
  at clearance `w` under `gap`, the bridges from the wall's corner to the two
  sides cut off a kite of width `2 w cos^2(theta / 2)`, which is returned as a
  face once it exceeds `gap`, although the wall is everywhere closer than `gap`
  (an equilateral hole of inradius 1 round an island of inradius `1 - w`: at
  `w = 0.6 gap` six gaps and no face, at `w = 0.8 gap` three faces and three
  gaps). By the definition of width that cell is a face, and a rule by inscribed
  disc would keep it as well; the tests pin both cases.
- Component nesting in `faces` is quadratic in the number of components in the
  worst case, and the closed-gap scan is quadratic in the hull in the worst
  case.
- No performance claim.

**Follow-ups / risks to watch**

- Reporting gaps just wider than the closing width, so a caller whose room
  leaks can see where.
- Nesting components in `faces` without the quadratic scan, which needs a
  committed benchmark before any claim.
- Exposing free faces and `locate` through the C ABI, once rheonda/cad has used
  the Rust API.
- At the release, the workspace floor of `axiolid-exact-arrangement` in the root
  `Cargo.toml` must rise to the version that adds `faces`, or the
  minimal-versions check will fail for overlay.

## Relation to existing code

- `crates/algorithms/planar/exact-arrangement/src/faces.rs`: `faces`,
  `Faces`, `FaceId`.
- `crates/algorithms/planar/exact-arrangement/src/distance.rs`: `Nearest` and
  the exact distance predicates on `ExactPoints`.
- `crates/algorithms/planar/overlay/src/free_faces.rs` and `free_faces/`:
  `free_faces`, `free_faces_within`, `FreeFaces::locate` and their types; the
  bridges, the classification, the rings and `locate` each have a module.
- `scripts/probe_free_faces_mutants.py` lists the faults the tests in
  `tests/free_faces.rs`, `tests/faces.rs` and `tests/distance.rs` must catch.
- ADR 0094, ADR 0100 (visible faces: arrangement, tracing, snap rounding),
  ADR 0104 (ring membership), ADR 0088 (memory budget).
