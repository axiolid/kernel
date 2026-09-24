# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `ArcArrangement` (#120): the plane cut by several arc rings at once.
  Every crossing, shared boundary piece and ring membership is decided
  exactly (the same predicates as `arc_overlay`); crossing points are
  rounded once, into one vertex table. Each piece records the input edges
  it came from and which rings contain the region on either side, and
  `regions(predicate)` links the pieces bounding any membership set into
  outer rings and holes. Faces built from one arrangement therefore share
  vertices by index, which is what a stepped or stacked solid needs.

### Changed

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
