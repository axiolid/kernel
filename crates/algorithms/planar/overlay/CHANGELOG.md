# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `segment_intersections` reports every intersection among many segments
  with a Bentley-Ottmann sweep (#146): each point where two or more
  segments meet, with every segment through it and where on it (start,
  end, interior, or a zero-length segment), and each collinear overlap
  with the segments covering it, in `O((n + k) log n)` time for `k`
  reported incidences. Every decision is exact (interval filter, then
  dyadic arithmetic), so shared endpoints, T-junctions, verticals, many
  segments through one point, overlaps and zero-length segments are
  handled rather than assumed away. Crossings are `ExactPoint2` rationals,
  rounded once and correctly on request; input endpoints come back bit for
  bit. Non-finite input is refused with the segment's index. 64,000 short
  grid segments take about 0.2-0.4 s (`cargo bench -p axiolid-overlay
  --bench segment_sweep`).

## [0.3.8] - 2026-09-30

### Changed

- Exact points share their coefficients, so copying one costs no
  arithmetic, and an input vertex builds its exact form only when a
  question about it gets past the interval filter, which most never do.
  Only the vertices of the result are rounded, not every vertex of the
  subdivision, and settling checks a hole against an outer ring's box
  before its edges. Results are unchanged; a soup of 2,000 overlapping
  triangles that share no edges takes about a fifth fewer instructions.

## [0.3.7] - 2026-09-30

### Changed

- Straight-edge booleans on soups are much faster, with results bit for
  bit the same. Each operand's rings are first reduced as a chain: an edge
  given once each way between the same two exact points cancels, which
  changes no winding number, so a mesh given as a soup of its triangles
  shrinks to its outline before the exact subdivision. What is left is
  split into cycles through distinct vertices and used when every cycle is
  simple, decided exactly; otherwise the rings are kept as given. The
  union of two overlapping 60 x 60 triangle meshes (14,400 triangles) now
  takes about 10 ms, from about 1 s in 0.3.5 (the grid backend before it:
  about 2.5 ms); 160,000 triangles take about 85 ms. Also: the
  subdivision interns its vertices in a tree instead of a sorted vector,
  seeds a point's parameter on a segment in closed form, and settling finds
  repeated vertices and touching edges by hashing and a sweep instead of
  asking every pair. Soups of triangles that share no edges gain less
  (2,000 random triangles: about 370 ms, from about 430 ms).

## [0.3.6] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.5] - 2026-09-28

### Changed

- Straight-edge booleans are exact (#173). `overlay`, `union_soup` and
  `Region` no longer go through `i_overlay`'s integer grid, which snapped
  every output coordinate, including untouched input vertices, to a step
  of about 1.5e-8 of the operands' extent. They now run on the exact
  subdivision the arc path uses: every ring of both operands cut at once,
  each piece classified by exact signs, kept by the operand's winding
  number under the fill rule. An input vertex the operation does not move
  comes back bit-identical (so `[0,4]x[0,0.2]` clipped by a box around it
  has area `0.8`, not `0.800000011920929`), and a crossing of two segments
  is the double nearest to the exact crossing. Fill rules and ring
  orientation keep their meaning (tested against the old backend). A vertex
  where the boundary runs straight on is dropped only when exactly
  straight.
- Arc and arrangement output: rational vertices (every segment crossing)
  are now correctly rounded instead of rounded to about 50 bits.
- The exact subdivision indexes edges and rings in box trees, so building
  it over thousands of rings (a projected mesh's triangles) is no longer
  quadratic in the ring count.

## [0.3.4] - 2026-09-27

### Added

- `minimum_enclosing_circle` (#118): the least circle holding a point set,
  by Welzl's algorithm in a fixed visiting order. Every in/out decision is
  exact, so the support points are the exact minimum circle's; the centre
  is enclosed from exact dyadic values and the radius rounded up, so the
  returned circle holds every point, and `CircleEvidence::error` bounds
  both the centre's offset and the radius's excess. Refuses empty and
  non-finite input (`CircleError`).

## [0.3.3] - 2026-09-27

### Fixed

- Unions of many small cells are operands again (#191): a hole touching
  its outer ring at a vertex was called outside whenever that vertex was
  the hole's first and lay where the boundary-exclusive ray test counts it
  out, so the next operation refused the region with `HoleOutsideOuter`.
  A hole is now outside only if a vertex of it lies strictly outside.
- `union_soup` settles its output like every other operation, so its
  polygons are accepted by `Region::new` (#191).

## [0.3.2] - 2026-09-27

### Added

- `minimum_area_rectangle` (#182): the least-area rectangle enclosing a
  point set, by rotating calipers over the exact convex hull. Caliper
  steps and the area comparison (`W H / |d|^2`, cross-multiplied) are
  exact; ties are broken by the least angle of the first axis, turned into
  `[0, 90)` degrees, so the result does not depend on input order. Returns
  an `OrientedRectangle` (centre, unit axes, half extents) and
  `RectangleEvidence` (hull size, count of tied orientations, and a bound
  on the rounding of the output -- for an axis-aligned rectangle the
  rounding actually done, measured exactly, so zero for a box with
  representable coordinates). Collinear input gives an exactly zero-width
  rectangle; empty or non-finite input is a `RectangleError`.
- `Region::visibility_polygon` (#184): the part of a region with holes in
  sight of a point inside it. An angular sweep round the viewpoint with
  every decision exact -- the order of vertex directions, which edges a
  wedge's middle ray meets and which is nearest -- so walls and holes cast
  exact shadows; only the shadow ends are rounded, once, and the result is
  presented like the other region operations. A viewpoint on the boundary
  or outside is refused (`VisibilityError::NotInside`).

### Fixed

- Every region an operation returns is now one `Region::new` accepts and
  the next operation takes: unions, intersections, differences,
  morphology, Minkowski sums and visibility polygons of shapes that share
  collinear edges, touch at vertices or pinch holes against their outer
  rings could come back with edges shorter than the tolerance or rings
  touching themselves, and were then refused as `RepeatedVertex` or
  `SelfIntersection`. Outputs are settled by the same tests validation
  applies -- short edges merged, a vertex touching another part of its
  ring put on it, rings split where they pass a point twice, holes given
  back to their outer rings -- each move within the tolerance.

## [0.3.1] - 2026-09-27

### Added

- `Region::minkowski_sum` and `Region::minkowski_erosion` (#145): the
  Minkowski sum and erosion of a straight-edged region -- non-convex, with
  holes -- by a convex polygon, at any rotation. Built in one
  `ArcArrangement` with exact decisions: the region translated by a vertex
  of the polygon united with the convex hull of the polygon at both ends of
  every boundary edge, and erosion as the region minus the sum of its
  complement (within a box) with the reflected polygon. Vertex sums are
  rounded once, output vertices once. A non-convex structuring polygon is
  refused as `MinkowskiError::NotConvex`.
- `Region::dilate_inner`, `dilate_outer`, `erode_inner` and `erode_outer`
  (#163): disc morphology on a stated side of the exact result, by
  inscribed and circumscribed 64-gons moved a margin further out or in
  (past vertex and output rounding and the presentation's dropping of
  edges shorter than the tolerance). `Region::bound` gives the side
  (`BoundSide::Inner` or `Outer`) and the greatest deviation from the exact
  disc morphology (`MorphologyBound`). An inner erosion proves reachability
  of a route found in it; an outer erosion proves unreachability where none
  is. `dilate` and `erode` keep their round joins, side unstated. A
  200-edge region takes about a second.

- `ArcArrangement` (#120): the plane cut by several arc rings at once.
  Every crossing, shared boundary piece and ring membership is decided
  exactly (the same predicates as `arc_overlay`); crossing points are
  rounded once, into one vertex table. Each piece records the input edges
  it came from and which rings contain the region on either side, and
  `regions(predicate)` links the pieces bounding any membership set into
  outer rings and holes. Faces built from one arrangement therefore share
  vertices by index, which is what a stepped or stacked solid needs.

### Changed

- `ArcArrangement` skips rings whose box holds neither a piece's sample nor
  an edge: winding numbers and crossings are computed only near each ring.
  Many-ring arrangements are several times faster.
- `arc_overlay` is exact (ADR 0070, #155). Crossings, their order along
  each edge, inside/outside/on classification, linking into rings and
  hole ownership are all exact sign decisions on the given input, via
  `axiolid-exact`; none reads the tolerance. Crossing points of two curves
  are rounded to `f64` once, in the output, and edges shorter than the
  tolerance after rounding are merged. The public API is unchanged.
- `cavalier_contours` is no longer a dependency. The arc path costs about
  70 to 100 us per boolean on typical sections instead of about 1 us
  (`benches/arc_overlay.rs`).
- `arc_overlay` skips edge pairs whose padded bounding boxes are apart,
  and links pieces through a sorted index instead of a scan. Decisions are
  unchanged (still exact); cost now grows close to linearly with edge
  count: two overlapping 256-edge rings went from 779 ms to 7 ms, a
  4096-edge outline against a small disc from 53 ms to 17 ms.

### Fixed

- `overlay` no longer rejects a U-shape or comb as `SelfIntersection`.
  The ring check treated an endpoint on the infinite line through another
  edge as touching it, so two collinear edges that share a line without
  meeting (the two ends of a U) were refused. A touching endpoint must now
  lie on the edge itself. Rings that genuinely touch are still refused.
- `arc_overlay` results no longer depend on drawing units. The arc
  backend's thresholds are fixed in drawing units, so a 5 um gap survived
  a union drawn in millimetres but vanished in metres. The drawing is now
  scaled by a power of two so those thresholds sit at the caller's linear
  tolerance, capped so coordinates stay within what f64 resolves
  (ADR 0069). Superseded by the exact core above, which needs no scaling.
