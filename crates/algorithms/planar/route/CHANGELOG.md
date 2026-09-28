# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Seeded weighted maps (#198): `weighted_distance_map_seeded` and
  `weighted_distance_map_seeded_within` take `(point, weight)` targets,
  each starting at its own non-negative cost, as `distance_map_weighted`
  does for plain maps (#197); both bounds start there, and
  `WeightedReach::cost` counts the weight. With every weight zero the map
  and its answers are `weighted_distance_map`'s exactly.
- `weighted_forced_walk` and `weighted_forced_walk_within` (#198): the
  cheapest walk that enters a polygon over two weighted maps with the same
  costs, bracketed with a witness as `forced_walk` is. The cell bound falls
  at twice the steepest factor meeting the cell, and the pair bound runs
  over both maps' vertices and cost-edge intervals with their lower
  bounds, at the cell's own factor when no cost edge meets it.
  `WeightedForcedWalk::shortest` brackets the cheapest walk overall; the
  result is never narrower than the maps' brackets at the witness, and the
  search stops at the tolerance plus those.

### Changed

- A cost region touching the region's boundary up to rounding is taken as
  touching it (#198), no longer refused with `MapError::CostCrossing`: a
  cost vertex within 2^-24 of the region's extent (and a few ulps) of a
  region edge, on its free side, is moved just beyond it, and a cost edge
  may cross a region edge that near one of either edge's ends. An edge so
  left on or beyond a wall is wall-borne, as one exactly along it. So a
  footprint clipped to the free region, its corners rounded by the
  overlay, builds a map, and leaves no sliver along a wall costing 1.

### Fixed

- Weighted maps on cost edges at an angle (#198). The points cutting such
  an edge are interpolated and lie off its line by rounding, which made
  three decisions go wrong: an interval's sides were taken of its rounded
  ends, so could both read the factor of one side; an edge could appear to
  cross every hop to its own intervals and block them -- both could put
  the lower bound above the distance; and every hop along the edge was
  charged the greatest factor above. Sides are now taken of the edge's
  exact ends, an edge never blocks hops from its own line, and a hop along
  an edge is costed above as the walk along the edge itself, which lies
  that close to it.

## [0.3.4] - 2026-09-28

### Added

- Weighted targets (#197): `distance_map_weighted` and
  `distance_map_within_weighted` take `(point, weight)` targets, each
  starting at its own non-negative distance -- a stair landing carrying the
  rest of a walk beyond it. `nearest` answers the least route length plus
  weight, in the new `Reach::distance`, and names the target the route
  ends at, which may be a lighter target the heavier one is reached from;
  `route.length` is the route's own. `farthest_point` and `forced_walk`
  keep their brackets (the weighted distance is still 1-Lipschitz), and
  `ForcedWalk::shortest` counts the origins' weights. With every weight
  zero the map is `distance_map`'s. `MapError::InvalidWeight` refuses a
  negative or non-finite weight.

- `forced_walk` (#196): the shortest walk from an origin to a target that
  enters a polygon, `min over p of d_origin(p) + d_targets(p)`, bracketed
  to a tolerance from two distance maps over the same free space. Branch
  and bound as for `farthest_point`, with a second lower bound from the
  maps' own vertices (`D(u) + D'(v) + max(|u - v|, dist(T, u) + dist(T,
  v))` over the vertices a cell may see), which is exact on every cell the
  walk runs straight through. `ForcedWalk::shortest` is the shortest walk
  overall; a lower end above it proves no shortest walk enters the
  polygon. A polygon no walk reaches is infinitely far, not an error.
- `FarthestError::MismatchedMaps`, for two maps over different region or
  barriers.
- Weighted distance maps (#195): `weighted_distance_map` takes cost
  regions (`CostRegion`: a polygon and a factor of at least 1; overlaps
  take the greatest, and a walk along a cost edge pays its cheaper free
  side) and points along cost edges `spacing` apart. `WeightedMap::nearest`
  brackets the weighted distance to the nearest target
  (`WeightedReach::cost`) and gives the walk whose cost is the upper end.
  The upper bound is the exact visibility graph with those points added,
  each edge costing its weighted length; the lower bound is a graph over
  vertices and intervals of cost edges whose hops cost no more than any
  piece of an optimal walk between them, with the states an optimal walk
  cannot take -- turning back at a cost edge, two pieces in a row along
  one line -- barred. Square crossings are exact; otherwise the gap is
  first order in the spacing. `weighted_farthest_point` brackets the
  farthest point of a subregion by the same branch and bound, each
  triangle's slope its greatest factor.
- `MapError::InvalidFactor`, `InvalidSpacing` and `CostCrossing`: cost
  regions may nest, touch, share edges and run along the region's
  boundary; one crossing an obstacle or another's edge, or running along a
  barrier, is refused.

### Changed

- Segment crossing tests reject disjoint bounding boxes before any
  predicate, and a segment's pass-through check visits only vertices with
  obstacle rays: `distance_map` builds several times faster.

## [0.3.3] - 2026-09-27

### Added

- `skeleton` (#139): the corridors of a region with holes as a graph --
  path ends, junctions and paths -- for circulation checks. Nodes are
  Voronoi vertices of the walls sampled at most `spacing` apart, over a
  constrained Delaunay triangulation, pruned of spurs into corners as in
  the lambda-medial axis (`prune`; 1.5 drops spurs into right-angled
  corners). Every node is
  decided inside the region exactly and carries a `clearance` interval
  proven to contain its distance to the nearest wall; each path end names
  the wall it runs into (`ahead`). The skeleton's position approximates
  the medial axis and is not certified.

## [0.3.2] - 2026-09-27

### Added

- `distance_map` and `distance_map_within` (#186): shortest-path distances
  to the nearest of several targets, from one multi-source Dijkstra over
  the visibility graph. `DistanceMap::nearest` gives, for any point of the
  region, the nearest target's index and the route there (`Reach`).
  Refusals are typed: `MapError::NoTargets`, `MapError::TargetOutside`.
- `farthest_point` and `farthest_point_within` (#186): the greatest
  distance to the nearest target over a polygon subregion, as a
  `LengthInterval` that contains the true value, with a witness point
  (`Farthest`). Branch and bound over a constrained triangulation of the
  free space, bounding each cell by the 1-Lipschitz property of the
  distance from anchors proven inside it; widened for the rounding of
  lengths and midpoints. Part of the subregion that no target reaches is
  refused with the triangle as evidence (`FarthestError::Unreachable`);
  crossing barriers are refused (`FarthestError::CrossingObstacles`).

### Fixed

- Routes no longer squeeze through a point where obstacles meet (#189):
  a barrier's foot on a wall, two holes touching at a corner, a barrier
  bent or joined at a vertex. Each graph vertex is split into the free
  sectors between the obstacle rays leaving it, ordered exactly; a route
  moves over (vertex, sector) states, each edge is taken on one side of
  travel and attached to the sector beside it at each end, and passes
  through a vertex only on a side with no obstacle ray and free space.
  Two polygons sharing an edge are walkable along it on either side.
  `shortest_path` and `DistanceMap` share the graph.

## [0.3.1] - 2026-09-27

### Fixed

- A route could run along one wall, through a vertex and on across a gap
  outside the region where two walls line up (#187): visibility tested only
  proper crossings and the segment's midpoint. A segment is now cut at
  every obstacle vertex lying on it, decided exactly, and each stretch is
  either along an obstacle edge or has its midpoint in the region. Two
  rooms whose corridor is cut are `DisconnectedComponents` again.
