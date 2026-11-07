# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `snap_round`, `Snapped`, `ClassBoundary`, `HOT_POINT_SCRATCH_BYTES` and
  `FRAGMENT_SCRATCH_BYTES` (ADR 0100, #410): snap rounding of straight,
  pairwise interior-disjoint edges to the binary64 tensor grid. Each edge is
  routed through the centre of every hot cell it meets, in order, decided
  exactly (half-open cells, ties to even); the fragments are arranged and a
  crossing left among them is refused. `Snapped` carries the displacement
  bound in the maximum norm and the number of rerouted edges, and
  `Snapped::class_boundary` nets a class's chain per fragment so that parts
  of zero width cancel.
- `check_winding`, which refuses a traced class that winds more than once
  at some point, and `ExactPoints::frame_sign`.

## [0.1.0] - 2026-11-03

### Added

- The exact point store, segment arrangement and ring tracing that the plane
  section of `axiolid-reference` used (ADR 0093, #366), as a crate of their
  own, so that the visible faces of `axiolid-hidden-lines` decide on the same
  exact signs (ADR 0094, #383).
