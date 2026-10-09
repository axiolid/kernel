# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Station seam sides (#263, ADR 0082 amendment): `OrientedCurveStation`
  gains `seam: SeamSide` (re-exported from `axiolid-curve`), the piece a
  station on a seam of its basis reads; `OrientedCurveStation::new` sets
  `SeamSide::Outgoing`, the stated default and the evaluators' behaviour,
  and `with_seam_side` changes it. `CurveStation` has public fields, so it
  keeps reading the outgoing piece; `CurveStation::with_seam_side` returns
  the oriented station in its base frame reading the given side, and
  `From<CurveStation> for OrientedCurveStation` the unturned outgoing one.
  The `station` module documents the seam rule and how runs of sections
  and offsets are mitred across a seam instead of taking a side.

## [0.3.5] - 2026-10-03

### Added

- Explicit station orientation (#246, ADR 0082 amendment):
  `StationOrientation { axis, ref_direction }`, both optional and given as
  components in the station's base frame `(tangent, lateral, up)` -- the
  reading of a linear placement's axes relative to its curve -- with the
  axis the oriented up and the reference direction orthonormalised
  against it (Gram-Schmidt, axis primary); `unit_axes` fills the defaults
  `(0, 0, 1)` and `(1, 0, 0)`. `ORIENTATION_TOLERANCE` bounds the sine
  below which the two count as parallel. Offsets stay in the base frame.
- `GeometryNode::OrientedCurveStation(OrientedCurveStation)`: a curve
  station with an orientation.
- `SectionAtStation` (built by `new`, `with_tags`, `with_orientation`; it
  is `#[non_exhaustive]`, so fields can follow additively) and its two
  relations, `SolidOperation::SectionsAtStations` (closed profiles, the
  general form of `StationedSpine`, matched by tag when tagged) and
  `SurfaceRelation::OpenSectionsAtStations` (open sections, the general
  form of `SectionedSurface`). Both variants are appended to their enums.
  `From<StationedSection>` and `From<StationedOpenSection>` convert the
  #241 sections.
- `GraphError::InvalidStation` also names a zero, non-finite or parallel
  orientation and a run of sections that mixes tagged and untagged ones.

- Stations (#241, ADR 0082): `Station` (a distance and `StationOffsets`
  lateral/vertical/longitudinal) along a basis curve, measured in that
  curve's convention -- plan distance on an elevated or banked curve, arc
  length on any other -- with offsets in its section frame
  (`StationFrame::Section`, or the upright `StationFrame::Plan`).
  `GeometryNode::CurveStation(CurveStation)` is a point and frame at a
  station; `CurveRelation::OffsetByStations` a 3D curve through offsets
  at stations, interpolated linearly in distance;
  `SolidOperation::StationedSpine` closed profiles standing at stations
  (`StationedSection`), matched by ring and vertex index;
  `SurfaceRelation::SectionedSurface` open sections at stations
  (`StationedOpenSection`) joined by tag.
- `GraphError::InvalidStation` names a malformed station when the node is
  pushed: a non-finite or negative distance, a non-finite offset, fewer
  than two stations in a run, distances that do not increase strictly,
  sections whose tags differ or repeat. A distance beyond the basis
  curve's length is refused when the station is resolved.

### Changed

- Tags are matched as sets (#246): every section of a run carries the
  same tags, none repeated; an open section's tags may run in the first
  section's order or in reverse, so a `SectionedSurface` whose sections
  are authored in opposite directions is now accepted and joined
  reversed. Any other order is still refused by name.

## [0.3.4] - 2026-10-03

### Added

- `TrimSelector::ArcLength(Scalar)` (#239): selects the basis parameter
  where the arc length along the basis, measured from its parameter `0` in
  the trim's sense, equals the value. Stored exactly; an evaluator resolves
  it by quadrature and a root find to a stated tolerance. It is a
  parameter-kind selector: it satisfies `TrimmingPreference::Parameter`,
  a non-finite value is refused by validation, and equal arc lengths are an
  empty trim.

## [0.3.3] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.2] - 2026-09-27

### Fixed

- A `SolidOperation::BoundedHalfSpace` boundary must be a 2D curve (#162).
  The graph accepted a 3D curve, which the compiler refuses, so such a
  graph validated and then could never compile; it is now refused when the
  graph is built.

## [0.3.1] - 2026-09-27

### Added

- A `Curve2::QuadraticGraph` or `Curve2::AngleGraph` with finite
  coefficients is accepted as a trim basis (#119).

- A `Curve2::Sinusoid` is a valid trim basis when its three coefficients
  are finite (ADR 0071).
