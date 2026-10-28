# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.5] - 2026-10-28

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.4] - 2026-10-12

### Added

- `MergedPoint` and `Triangulation::merged_points` (#300): the input
  points merged into an exactly coincident one, each as its input index and
  the kept index that stands for it (the lowest input index at that
  position), in ascending input order.
- `TriangulationError::CoincidentConstraint { a, b }` (#300): a constraint
  whose two endpoints are distinct input points at exactly the same
  position. Merged, they are one vertex, and the constraint would be an
  edge from it to itself.

### Fixed

- `triangulate` merges exactly coincident input points instead of
  splitting each into the mesh at its twin's position (#300), which left
  zero-area and overlapping triangles that used one directed edge twice --
  in a 3x3 grid with its first row repeated, edge (3, 10) twice -- and
  reported nothing. A point found at a vertex already inserted is recorded
  as merged into it; the triangles and `Triangulation::constraints` use the
  kept indices, `points` still holds every input point, and
  `merged_points` reports the merge, so it is lossless. Duplicates in free
  points, outlines, holes and constraints triangulate validly; a constraint
  between two coincident points is refused by name. A constraint from an
  index to itself is unchanged.
- The debug adjacency check no longer skips a triangle list with a
  repeated directed edge; it panics on one (#300). Only duplicate points
  produced such lists, and they are merged now.

- `triangulate` covers the exact convex hull of its input (#299). The
  build bracketed the points with a finite super triangle (32 times their
  span) and dropped every triangle touching it, which lost the thin hull
  triangles whose circumcircles reached its corners -- 3,977 instead of
  3,979 triangles for 2,000 random points, with nothing reported. The hull
  is now closed by one symbolic vertex at infinity, decided by the
  predicates without coordinates, so the output always has `2n - 2 - h`
  triangles for `n` points with `h` on the hull boundary, collinear ones
  included, and its boundary is the exact hull. Constrained triangulations
  recover their constraints in that complete triangulation.
- Where the Delaunay triangulation is unique (no four points cocircular)
  and the old build kept every hull triangle, the output is the same set of
  triangles as before; their order in the list changes, because the build
  starts from the first three points that span a triangle instead of from
  the super triangle. With cocircular points either build may choose a
  different, equally Delaunay triangulation.

## [0.3.3] - 2026-10-10

### Changed

- `triangulate` keeps its adjacency current locally (#268): a split of a
  triangle into three, of an edge and its two triangles into four, and
  each Lawson flip -- after insertion, in constraint recovery and in the
  pass that restores Delaunay -- re-link only the twin pointers of the
  triangles they rewrite and of those triangles' outside neighbours,
  including the old neighbours of a reused triangle slot. The build used
  to recompute every twin pointer after each insertion and flip, so n
  points cost O(n^2) or worse. Debug builds compare the maintained
  adjacency with a full recomputation after every step (up to 1024
  triangles outside the crate's own tests) and panic on the first
  difference.
- Point location starts its walk next to the last vertex inserted in the
  point's cell of a coarse grid, instead of from the last triangle, whose
  walk crossed O(sqrt n) triangles per point. A point on an edge or a
  vertex is located by the old walk, so the output is unchanged: the same
  triangles in the same order.
- Together: 16,000 uniform random points took 72 s and take 18 ms; from
  1,000 to 16,000 points the time grows 29x (n log n: 22x, the old build:
  over 300x). Measured by `tools/benchmark`'s `triangulate_points`.

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Fixed

- `triangulate`: after recovering constraint edges, unconstrained edges are
  flipped back to locally Delaunay (Lawson's flips). Recovery used to leave
  long triangles whose circumcircles held points no constraint hid, so the
  result was not constrained Delaunay as documented -- with finely sampled
  walls, triangles spanned whole rooms (#139).
- `triangulate` no longer loops for ever on outlines like a square turned
  45 degrees (#190): legalisation after inserting a point checked the new
  diagonal instead of the two edges across from the point, so the real
  edges were never checked and thin quadrilaterals flipped back and forth.
  It now checks those edges, flips only strictly convex quadrilaterals,
  splits the edge (and both triangles beside it) when a point lands on
  one, and carries a flip bound.
- Constraint recovery no longer gives up at the first crossing edge it
  cannot flip (#190): it follows Anglada's queue, retrying edges that
  cannot flip yet and requeuing new diagonals that still cross, and
  reports `CrossingConstraints` only when a full pass flips nothing.
