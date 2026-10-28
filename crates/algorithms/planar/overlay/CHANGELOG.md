# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.19] - 2026-10-20

### Added

- `overlay_within` and `union_soup_within` take `OverlayOptions` and the
  caller's `ExecutionOptions` (#338, ADR 0088): the tolerance is the
  options', and with `ExecutionOptions::with_memory_budget` the
  subdivision charges what it is about to hold before allocating it and is
  refused with `OverlayRefusal::BudgetExceeded { needed, budget }` as soon
  as the charges pass the budget. `OverlayRefusal::Invalid` carries every
  refusal `overlay_with` and `union_soup_with` make. Without a budget
  nothing new is refused, and an admitted call gives the unbudgeted result
  bit for bit. `tests/grown_soup.rs` measures with a counting allocator
  that a call never holds more heap than its budget, admitted or refused,
  and that four times the unbudgeted high-water mark is admitted. The
  crate now depends on `axiolid-contracts`.

### Changed

- Ring validation no longer compares every pair of a ring's edges (#341).
  A ring of more than 32 edges asks the crate's box tree for the pairs
  whose boxes, grown by the tolerance, overlap, and the unchanged exact
  `segments_intersect` decides each, in the old order: the refusals, their
  order and the first pair found are as before, which a test checks
  against the all-pairs search on rings valid, crossing, touching and
  within the tolerance. And a `Region` carries the proof that its rings
  are valid under the tolerance it was built with: `Region::new` sets it,
  and a boolean on regions and the disc morphology set it from their
  settled output, so `Region::union`, `intersection` and `difference`
  under the same tolerance validate no ring again (holes are still checked
  against their outer rings). A region translated, or used under another
  tolerance, is validated as before; the proof is not part of `Region`
  equality. Measured with `cargo bench -p axiolid-overlay --bench
  ring_validation` (release), on a 40 m by 30 m floor with toothed walls:
  `Region::new` at 10,000, 50,000 and 100,000 corners took 0.17 s, 4.2 s
  and 17 s and now takes 0.006 s, 0.028 s and 0.059 s; a
  `Region::difference` with a 1 m square took 0.20 s, 4.4 s and 18 s per
  call and now takes 0.029 s, 0.18 s and 0.38 s.
  `scripts/probe_ring_validation_mutants.py` lists the faults the tests
  must catch.

### Fixed

- A soup of grown triangles whose union is one polygon no longer needs
  gigabytes (#338). axioval/engine's soup -- the strips of a faceted bar,
  each corner grown by a circumscribed 16-gon, hulled and fanned -- took
  4.1 GB and 26 s for `union_soup` at 128 facets and 1.9 GB and 10 s for
  `overlay(Intersection)`, and at 256 facets, and for the soup's union
  with itself at 128, failed to allocate. Three causes, each fixed with
  results unchanged:
  - The chain reduction walked the whole soup once and gave up on every
    ring when one cycle crossed itself: two overlapping hulls sharing
    vertices let the walk pass from one into the other. It now reduces
    such a soup group by group, the rings joined by edges given each way,
    so every fan becomes its hull's outline.
  - Edges with the same two ends (hulls' flat sides lined up on one line
    and split there by the snap) were paired with one another for every
    copy. Only the first copy keeps their pieces, so the later copies are
    no longer paired among themselves.
  - Where many edges meet at one vertex, every pair recorded the vertex
    for both edges. A point at an edge's own end is no longer recorded
    for that edge: it is its first or last stop, and the rings met there
    are never read.

  Measured with `cargo bench -p axiolid-overlay --bench grown_soup`
  (release, peak resident memory): circumscribed `union_soup` at 64, 128
  and 256 facets takes 0.08 s, 0.32 s and 1.6 s at 22 MB, 64 MB and
  245 MB, and the intersection 0.09 s, 0.37 s and 1.8 s at 24 MB, 70 MB
  and 254 MB, against 0.07 s, 0.32 s and 0.94 s at 20 MB, 58 MB and
  217 MB for the inscribed soup. The self-union at 256 facets takes 5.5 s
  at 482 MB. Every output that completed before is the same bit for bit,
  pinned by hash in `tests/grown_soup.rs`.
  `scripts/probe_grown_soup_mutants.py` lists the faults the tests must
  catch.

## [0.3.18] - 2026-10-19

### Added

- A fast mode for the disc morphology (#330, ADR 0087).
  `Region::erode_inner_with`, `erode_outer_with`, `dilate_inner_with` and
  `dilate_outer_with` take a `MorphologyMode`:
  - `MorphologyMode::Exact`, the default, is what `erode_inner` and its
    siblings give, unchanged bit for bit. A CAD kernel takes it.
  - `MorphologyMode::fast(budget)` first fills the summed set's notches with
    chords proved to lie outside it (`src/fill.rs`), so a wall of teeth
    becomes one edge. A chord is used only when the excess it adds stays
    within `budget`, in metres; a refused chord is split into shorter ones.
    The result keeps its inner or outer side. The excess is part of
    `MorphologyBound::deviation`, and `Region::excess` reports it on its
    own. A rule engine takes this mode.

  The excess is a proved bound on how much larger the filled set makes the
  sum: about the square of the spacing of the corners on the chord over
  eight times the radius, and never more than the notch's depth. A budget
  that is not finite and nonnegative is refused with
  `OverlayError::InvalidOffsetDistance`.

  Measured on #330's cases (#323's room with 0.5 mm teeth, as a room and
  as a hole, a 0.4 m disc, a 1 mm budget, release) with
  `cargo bench -p axiolid-benchmark --bench disc_erosion`. At 10,002
  corners, exact mode takes 3–11 s and about 1.2 GB, and grows about with
  the square of the corners (2.9–7.8 s and 1.15–1.19 GB in this run).
  Fast mode takes 0.010–0.013 s, 0.020–0.024 s and 0.036–0.045 s at 2,402,
  5,002 and 10,002 corners, at 4–7 MB, so it grows about linearly. Its
  excess is 35 µm, 8 µm and 2 µm at those sizes, and 0.9 mm at 10,002
  corners in the four cases where a chord within the budget cuts a corner
  of the room. Exact mode's hashes are those of 0.3.17 in every case. Tests check fast results pointwise against the true disc
  morphology on the certified side, within the stated deviation, for rooms
  and holes, inner and outer, erosion and dilation.
  `scripts/probe_fast_morphology_mutants.py` lists the faults they must
  catch.

## [0.3.17] - 2026-10-16

### Fixed

- Minkowski sums and erosions, and the disc morphology built on them, no
  longer refuse with `OverlayError::ZeroArea` footprints that 0.3.15
  handled (#331, a regression of #323 in 0.3.16). Where a wall runs
  parallel to a side of the convex piece -- a tessellated round body
  against the 64-gon of `dilate_outer`, say -- with a corner an ulp off
  the chord across it, trimming an edge hull by the chord's hull left a
  sliver of three nearly collinear vertices, with an area under the
  tolerance squared, and the arrangement refused it. A chord whose hull or
  any trimmed hull would be such a ring is now not used: its hull and its
  trimmed hulls are dropped and the edges it spans keep their whole hulls.
  The consumer's footprints that failed now give 0.3.15's results bit for
  bit; `tests/minkowski.rs` pins synthetic walls of that kind, 15 of whose
  36 results 0.3.16 refused, against 0.3.15 by hash. Results 0.3.16 gave
  are unchanged, and so is the cost. `scripts/probe_minkowski_degenerate_mutants.py`
  lists the faults the tests must catch.

## [0.3.16] - 2026-10-16

### Changed

- Minkowski sums and erosions, and the disc morphology built on them
  (`Region::erode_inner`, `erode_outer`, `dilate_inner`, `dilate_outer`),
  no longer grow with the square of the corners along walls of many short
  edges (#323). On the 4.5 m by 5 m room whose long walls carry 0.5 mm
  teeth, eroded by a 0.75 m disc, 102, 602 and 1,202 corners took
  0.045 s, 1.15 s and 4.4 s at 27 MB, 416 MB and 1.7 GB, and now take
  9 ms, 68 ms and 0.23 s at 6 MB, 25 MB and 98 MB. With the corners
  exactly collinear they took 0.050 s, 1.4 s and 5.7 s and now take
  2 ms, 8 ms and 12 ms, the rectangle's 2 ms plus reading the corners.
  Measured with `cargo bench -p axiolid-benchmark --bench disc_erosion`,
  interleaved runs against 0.3.15 on a shared machine. Each boundary edge
  gave one hull, each overlapping hundreds of its neighbours. Now a vertex
  exactly between its neighbours on one line (an exact orientation, no
  tolerance) is merged away before the hulls are built, and a chord
  between two vertices of a ring, proved exactly to lie in the summed set,
  gives one hull that covers nearly all of each edge hull it spans; those
  keep only what it does not hold, a sliver about the depth of the teeth.
  The union of the terms is the same set, so the results are too: bit for
  bit on the toothed room and on every sum and erosion of the crate's
  tests but two, whose results lose vertices lying exactly on the line
  between their neighbours, as the merged edges no longer split there.
  `tests/minkowski.rs` pins the toothed room's erosions against 0.3.15
  by hash, a unit test compares every route with the untrimmed hulls
  bit for bit, and `scripts/probe_minkowski_cover_mutants.py` lists the
  faults they must catch.

## [0.3.15] - 2026-10-14

### Changed

- Exact straight-edge booleans of soups of overlapping slender parts,
  the plan footprint of an assembly, take time in proportion to their
  crossings, and an intersection with a small operand in proportion to the
  parts near it, with results bit for bit the same (#315). On `k`
  slender boxes (1 m by 4 cm by 2 cm) over a square metre, the union of
  the soup with itself took 0.92 s, 6.3 s and 39 s at `k` = 100, 200
  and 400 and now takes 0.23 s, 0.94 s and 3.9 s, growing as the
  crossings do (about `k^2`), and its intersection with a 5 cm square
  took 0.32 s, 2.1 s and 13 s and now 5 ms, 19 ms and 62 ms. Measured with `cargo bench -p axiolid-benchmark
  --bench overlay_slender`, interleaved runs against 0.3.14 on a shared
  machine. The time had grown about `k^2.65`, because the turned edges'
  boxes hold much of the soup and every piece of the subdivision searched
  `O(k)` of them for its carriers and counted `O(k)` rings' windings. A
  piece's carriers are now searched only among the edges sharing two or
  more points with its edge, and across any stop only the rings met there
  are counted again, the others carried over; the points met along an
  edge are merged by one sort, and the crossings of a segment given more
  than once are computed once per pair and shared. An intersection leaves
  out, before the subdivision and after the snap, every ring whose box
  misses the other operand's: the module docs of `exact_overlay` give the
  argument that this changes no output. Debug builds check the carriers
  and the rings carried across a stop against the full search and count.
  `tests/soup_snapshots.rs` pins the outputs on the consumer's soups, a
  torus, a cylinder, random triangles and an arc subdivision against
  0.3.14 by hash, and `scripts/probe_overlay_slender_mutants.py` lists
  the faults it and the suite must catch.

### Fixed

- The crate's internal exact booleans refuse a ring with a repeated
  vertex as `RepeatedVertex` instead of hanging on it (#309). Every
  public entry validates first, so no caller could reach it: a piece of
  no length had no point strictly inside, and the exact bisection for
  its sample never ended. The sample now refuses such a piece by name,
  and the straight-edge boolean checks for equal consecutive vertices
  before it starts. Tests call both with a ten-second limit, and
  `scripts/probe_repeated_vertex_mutants.py` lists the faults they must
  catch.

## [0.3.14] - 2026-10-13

### Added

- An explicit opt-in to drop degenerate input rings instead of refusing
  the whole overlay, each drop reported (#256). `overlay_with` and
  `union_soup_with` take `OverlayOptions`; their default is `overlay` and
  `union_soup` exactly, a named refusal.
  `OverlayOptions::new().dropping_degenerate_rings()` leaves out a ring
  that has fewer than three points, an edge no longer than the tolerance
  (`RepeatedVertex`), or an area of at most `tolerance.linear()^2` under
  the exact `ZeroArea` rule, and lists it in `ReportedOverlay::dropped` /
  `ReportedUnion::dropped` as a `DroppedRing`: operand, polygon, outer or
  hole index, `DropReason`, and the magnitude of its signed area, summed
  exactly and rounded once, with a proven upper bound. A hole of a dropped
  outer is dropped and reported with it (`OuterDropped`); a degenerate
  hole of a kept outer is dropped alone. A self-intersecting ring is never
  dropped, nor a `RepeatedVertex` ring that still crosses itself once its
  short edges are merged, since its signed area is not what leaving it out
  removes; non-finite points, holes outside a kept outer and invalid
  frames stay refusals too. The overlay of the kept rings is exact as
  usual, and its evidence counts the kept rings.
  `scripts/probe_degenerate_drop_mutants.py` lists the faults
  `tests/degenerate_drop.rs` must catch.

### Changed

- Exact straight-edge booleans are faster on soups whose rings share no
  edges, with results bit for bit the same (#218). The union of 2,000
  random overlapping triangles takes about 94 ms, from about 290 ms (the
  grid backend: about 26 ms), and two crossing 1000-gons about 8.1 ms,
  from 12.3 ms, two thirds of it now operand validation; a union of two
  60 x 60 meshes given as triangle soups takes about 19.5 ms, from
  144 ms. Measured with `cargo bench -p axiolid-benchmark --bench
  overlay_soup`, five interleaved rounds on a shared machine. The
  subdivision takes a segment no other boundary runs along as carried by
  its own ring only, and carries the rings containing a piece across a
  stop where exactly one other segment crosses, instead of searching for
  carriers and counting windings from a sample for every piece; any
  other stop counts afresh, and debug builds check both shortcuts
  against the search and the count. Samples on segments, crossings of two
  segments and input vertices keep their exact coefficients unbuilt until
  a question gets past their boxes; copies of a point share an identity
  that settles equality and interning at once; a ring's edges share their
  vertices, and each segment its direction. Rounding a crossing compares
  the offset from its guess in a tight interval before any exact
  midpoint, the subdivision's queries reuse their buffers, and the snap
  to the tolerance (#222) finds the vertices near an edge through a grid,
  not its whole band of `x`, which made mesh soups slow. Arc booleans
  share the changes and are unchanged or faster in `benches/arc_overlay`.
  `tests/soup_membership.rs` checks soups by point membership, and
  `scripts/probe_overlay_soup_mutants.py` lists the faults it and the
  suite must catch.

## [0.3.13] - 2026-10-10

### Fixed

- `overlay` accepts its own settled output as an operand when that output
  holds two collinear edges that do not touch (#272). Operand validation
  asks whether non-adjacent edges meet, and its crossing test read the
  side of a rounded cross product with `0.0` counted as a side: two
  disjoint edges on one line whose cross products rounded to `0.0` and to
  a tiny positive residue counted as crossing, and the union of a turned
  U fed back was refused as `SelfIntersection` (6 of 20,000 turned copies
  in the report). Every contact decision is now exact, in the new
  `contact` module that validation and settling share: a proper crossing
  needs each segment's ends strictly on both sides of the other's line,
  an end touches another edge when its cross product is within the
  tolerance and it lies in the edge's box widened by the tolerance, both
  decided exactly, and collinear edges therefore meet only when their
  extents do. The test whether a hole's corner lies inside its outer ring
  is an exact orientation sign too, and settling puts a proper crossing at
  the double nearest the exact one instead of at a point interpolated from
  rounded cross products, which edges crossing at a hair's angle put
  anywhere along the edge. A ring decides the same at any distance from
  the origin. `scripts/probe_edge_contact_mutants.py` lists the faults the
  tests must catch.

## [0.3.12] - 2026-10-10

### Fixed

- Minkowski sums and erosions of large regions use memory linear in their
  subdivision, not quadratic in the region's corners (#292). Every piece
  of an `ArcArrangement` kept a flag for every ring on each side, and the
  sum reads one hull per boundary edge, so a floor's walls of 3,400
  corners and 58 rooms summed with a 0.9 m square needed more than 12 GB.
  A piece now keeps only the rings that hold it, as the subdivision
  computes them, and the sums and erosions read each face from that list:
  a polygon holds a side when its outer ring is listed and none of its
  holes, the union when any piece or translated copy does. Results are
  bit for bit the same. On a synthetic floor (`cargo bench -p
  axiolid-benchmark --bench minkowski_plan`, interleaved runs against
  0.3.11 on a shared machine), 3,370 corners: sum 0.89 s and 531 MB peak
  before, 0.61 s and 124 MB after; with ribbed walls, 3,522 corners: sum
  4.1 s and 2.2 GB before, 3.0 s and 469 MB after, and erosion of the
  rooms 4.7 s and 2.5 GB before, 3.5 s and 528 MB after. Time and memory
  now grow about linearly with the corners. `ArcArrangement::regions`
  fills one buffer of flags per side instead of storing them.

## [0.3.11] - 2026-10-09

### Fixed

- Small rings far from the origin keep their area and orientation (#274).
  A ring's orientation and its `ZeroArea` check came from a shoelace over
  the coordinates as given, whose rounding at georeferenced positions
  (6e5, 5.6e6 m) swamps the area of a small ring: a 1e-4 m² triangle was
  refused as `ZeroArea`, and a clockwise sliver of a chorded round corner
  read as counter-clockwise made the union of two triangles sharing an
  edge fail with `SelfIntersection`. Both are now exact signs of the area
  taken as a fan of triangles from the ring's first vertex, at every site
  that decides with them: validation (`overlay`, `union_soup`,
  `Region::new`, and `validate_arc_ring` for a ring without arcs), the
  exact boolean's winding weights, canonical output orientation,
  settling, the outer-or-hole decision when linking result rings, and the
  orientation `arc_overlay` and `ArcArrangement` give their operands. A
  ring now decides the same at any distance from the origin.
  `arc_ring_area`, `ring_area`, `polygon_area` and `Region::area` sum from
  the first vertex too, so their rounding scales with the ring's extent.
- A ring whose vertices all lie exactly on one line is a valid operand
  wherever it lies; a boolean leaves it out, as it encloses nothing (#219).
  Before, its rounded area decided whether validation refused it as
  `ZeroArea`, so the plan shadow of a vertical face passed at `x = 4.2`
  and was refused at `x = 4`. A ring with arcs, or one passed to
  `validate_arc_ring`, is still refused, as the arc path has no rule
  leaving it out.

## [0.3.10] - 2026-10-02

### Fixed

- Features closer than the tolerance touch (#222). Before the exact
  boolean, a vertex within `tolerance.linear()` of an earlier vertex moves
  onto it, and one within it of another ring's edge is inserted into that
  edge. A door opening computed as `3.9999999999999996..4.199999999999999`
  between rooms ending at 4 and starting at 4.2 used to leave a gap one ulp
  wide, so the union came back in three pieces and coverage could not pass
  through the opening; the grid backend before 0.3.5 had closed it. Inputs
  with nothing that near another feature come back bit for bit as before.

### Added

- Minkowski sums and erosions with a non-convex polygon (#145):
  `Region::minkowski_sum` and `Region::minkowski_erosion` now take any
  simple polygon, where a non-convex one was refused with
  `MinkowskiError::NotConvex`. New `Region::minkowski_sum_region` and
  `Region::minkowski_erosion_region` take a whole region as the other
  operand -- non-convex, with holes, several components; the region sum is
  symmetric and cuts the operand with fewer vertices. A non-convex shape is
  cut into convex pieces of its own vertices (holes bridged in, ear
  clipping, then Hertel-Mehlhorn merging across diagonals while convex),
  every decision an exact orientation, and the cut is certified to tile the
  shape (the triangles' boundaries sum to the shape's); the pieces' sums
  with every boundary edge are united in the same exact arrangement as the
  convex case, so vertex sums are rounded once and output vertices once.
  Should no certified cut be found, the sum is taken from edge-pair
  parallelograms instead. Convex polygons take the path they took before,
  and give the same results bit for bit. A comb-shaped room of 152
  vertices with a hole per two teeth is summed with a rectangle in about
  30 ms, an L-shape (two pieces) in about 55 ms and a plus (three pieces)
  in about 180 ms; erosion takes about 1.3 times as long. The cost grows
  with the number of pieces times the region's edges.
- `MinkowskiError::EmptyStructuring`: an erosion by an empty region, under
  which every point would qualify.

### Fixed

- `Region::minkowski_erosion` by a polygon that does not hold the origin
  was cut to the region itself, so a translate of the polygon lying
  inside the region from a point outside it was missed (eroding
  `[0, 1] x [0, 10]` by `[5, 5.5] x [0, 1]` came out empty instead of
  `[-5, -4.5] x [0, 9]`). The erosion is now anchored on the region moved
  by a point of the polygon. Polygons holding the origin, and the disc
  erosions, are unchanged.

## [0.3.9] - 2026-09-30

### Fixed

- A ring whose vertices all lie on one line -- the plan shadow of a
  vertical face -- is left out of a boolean, as it encloses nothing (#219).
  Its rounded area need not be zero, so validation let it through, and
  since 0.3.5 the exact subdivision refused it as self-intersecting; the
  grid backend before had answered as if it were absent. A union of such a
  ring with itself is empty again, and one with a polygon is that polygon.

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
