# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Operands that touch (#167, ADR 0075 stage 2): faces on one surface share
  their overlap (each face's edges are imprinted on the other, and a region
  on the other solid's boundary is kept once by normal agreement); sections
  along an existing edge split only the other face; tangent contact adds no
  section; pieces leaving a vertex in one direction are ordered by
  curvature; solids meeting along an edge are paired radially around it so
  each stays manifold. Cavities go to the smallest solid around them, in
  results of several solids too.

- `section_edges` (#167, ADR 0075 stage 1): the exact intersection curves of
  two exact B-reps' faces, each trimmed to where it lies inside both faces.
  Crossings with a boundary edge are found against the adjacent face's
  surface, or across a seam against the plane through the ruling.
- `boolean(a, b, operator, tolerance)` (#167, ADR 0075 stage 1): the exact
  union, intersection and difference of two exact solids whose faces lie on
  planes and cylinders and meet in lines, circles and ellipses -- not only
  vertical columns. Regions are classified by exact ray parity with
  certified face membership and sewn into shells; cavities become voids.
  Every result audits clean and measures exactly.
- `split_face` (#167): a plane or cylinder face cut along its section edges
  into regions, traced in the face's parameters with exact pcurves (lines,
  conics, rulings, circles about the axis, `Sinusoid2` for oblique cuts).
