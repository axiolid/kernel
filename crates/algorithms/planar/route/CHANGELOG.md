# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

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
