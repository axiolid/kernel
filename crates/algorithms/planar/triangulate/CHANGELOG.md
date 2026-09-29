# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

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
