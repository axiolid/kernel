# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- JSON wire payloads (ADR 0085) now read every real back bit for bit, as the
  ADR promises. `serde_json` was used without its `float_roundtrip` feature,
  whose fast parser lands one ulp off for about a quarter of ordinary
  doubles (12.838242182209141 read back as its lower neighbour), in every
  format version. CBOR was unaffected. A test now sweeps 6,000 doubles across
  twenty magnitudes through both encodings. The fix is in parsing only, so
  payloads already written are unchanged and now read back exactly.

## [0.3.12] - 2026-10-13

### Added

- Run stations at their basis's start or end, measured by the kernel
  (#303, ADR 0082 amendment): `RunStation` gains `basis_end:
  Option<BasisEnd>` (a new `#[non_exhaustive]` enum, `Start` or `End`),
  set by the constructors `RunStation::at_basis_start(offsets)` and
  `RunStation::at_basis_end(offsets)`; `None`, the default, is the stated
  `station.distance` as before. Such a station states the distance `0`,
  which is not read; a resolver reads the start as `0` and the end as the
  basis's length by the kernel's own measure. The graph refuses, as
  `GraphError::InvalidStation` when the node is pushed, a start station
  that is not the run's first, an end station that is not its last, and
  either stating a distance other than `0`; the stated distances still
  increase strictly, and whether the end comes after them is the
  resolver's check.
- Wire format 1.3 (ADR 0085 note): `RunStation.basis_end`, written only
  when it is set. A graph using it is labelled 1.3; one not using it is
  written exactly as before, as 1.0, 1.1 or 1.2. A 1.2 payload holding it
  is refused as `UnknownField` `basis_end`, by the same test the writer
  labels by. Golden payload `tests/wire/v1.3/every-kind.json`.

## [0.3.11] - 2026-10-12

### Added

- Seam-snapping windows and seam sides on run stations (#301, ADR 0082
  amendment): `SectionAtStation` gains `seam` and `seam_window`, set by
  `with_seam_side` and `with_seam_window`, for `SectionsAtStations` and
  `OpenSectionsAtStations`; a new `RunStation` (`#[non_exhaustive]`, built by
  `new` or `From<Station>` and the same two methods) carries them in a new
  variant, `CurveRelation::OffsetsAtStations`, the general form of
  `OffsetByStations`. Defaults are `Outgoing` and `0`, today's reading;
  `Station`, `StationedSection`, `StationedOpenSection` and the #241
  variants are unchanged. A negative or non-finite window on any run
  station is refused as `GraphError::InvalidStation` when the node is
  pushed.
- Wire format 1.2 (ADR 0085 note): the two fields are written only when
  they are not `Outgoing` and `+0.0`, and `OffsetsAtStations` is a 1.2
  variant. A graph using either is labelled 1.2; one using neither is
  written exactly as before, as 1.0 or 1.1. A 1.0 or 1.1 payload holding
  them is refused by name (`UnknownField` `seam` or `seam_window`,
  `UnknownKind` `OffsetsAtStations`), by the same test the writer labels
  by. Golden payload `tests/wire/v1.2/every-kind.json`.

## [0.3.10] - 2026-10-10

### Fixed

- The wire writer labels each payload with the lowest format version that
  represents it (#297, ADR 0085), not always with `FORMAT_VERSION`: a graph
  with no seam-snapping window is written as `"version": "1.0"`, byte for
  byte what 0.3.8 wrote, so a 0.3.8 reader reads it again; 0.3.9 labelled
  every payload `1.1`, which such a reader refuses. A graph holding a
  window anywhere (an `OrientedCurveStation` node or the station of an
  `InstanceAtStation`, a window of `-0.0` included, since it is written)
  is written as `1.1`, as before. JSON and CBOR alike.

## [0.3.9] - 2026-10-10

### Added

- A seam-snapping window on an oriented station (#294, ADR 0082
  amendment): `OrientedCurveStation::seam_window`, set by
  `with_seam_window`, a non-negative distance in the basis's station
  measure (default `0`, the arc-length tolerance alone, so every existing
  station reads as before). A resolver reads a station within the window
  of one seam whose frame may jump on that seam, from its `SeamSide`, and
  refuses a window holding two seams or a seam whose certified position
  straddles its edge. It also applies to the station of an
  `InstanceAtStation`. A negative or non-finite window is refused as
  `GraphError::InvalidStation` when the node is pushed.
- Wire format 1.1 (ADR 0085): the window is written as `seam_window` and
  omitted at `+0.0`, so a graph without windows writes the 1.0 graph; an
  absent window reads as `0`, so every 1.0 payload still reads to the same
  graph. A payload claiming 1.0 that holds a window is refused as
  `WireError::UnknownField`. `FORMAT_VERSION` is `1.1`; the 1.0 golden
  payloads stay and a 1.1 one is added.

## [0.3.8] - 2026-10-10

### Added

- A versioned wire format for `GeometryGraph` behind an optional `serde`
  feature, off by default (#267, ADR 0085): JSON (`serde_json`) and CBOR
  (`ciborium`) encodings of one envelope,
  `{"format": "axiolid-geometry-graph", "version": "1.0", "graph":
  {"nodes": [...], "roots": [...]}}`, with node references as indices and
  enums externally tagged by their variant names. `GeometryGraph::to_json`,
  `from_json`, `to_cbor` and `from_cbor`, and the `wire` module with
  `FORMAT_NAME`, `FORMAT_VERSION`, `FormatVersion` and the typed
  `WireError`. A reader checks the version before the graph and refuses a
  newer minor, another major, an unknown kind, variant or field, a
  non-finite number and malformed input by name; a decoded graph is
  rebuilt through `GeometryGraphBuilder::push`, so it passes the same
  validation as any other. A writer refuses a graph holding a non-finite
  number. Golden payloads under `tests/wire/` pin format 1.0.
- Both readers accept an integer where a real is expected when it is
  exactly a double (`i as f64 as i128 == i`, so every `|i| <= 2^53`), as
  JavaScript and Python CBOR encoders write whole-number doubles as
  integers. A CBOR integer that is not exactly a double is refused as
  `WireError::Malformed` naming its path, wherever it stands. Integer zero
  has no sign and reads as `+0.0`. The CBOR reader unwraps the
  self-describe tag and refuses other tags, byte strings, simple values,
  non-text and repeated map keys.
- With the `serde` feature, every node payload type derives `Serialize`
  and `Deserialize`. A `NodeId` deserialises only inside a graph payload.

## [0.3.7] - 2026-10-09

### Changed

- A `SolidOperation::BoundedHalfSpace` boundary may be a `Profile` node
  as well as a 2D curve (#277, ADR 0084): a contour of lines and exact
  circular arcs states a boundary a polyline cannot. Graph validation
  accepts either and refuses anything else as `"curve2 or profile"` (a
  3D curve was refused as `"curve2"`). No type, variant or field changed.

- The `station` module documents curve relations as station bases
  (#285): a composite, a trim and a curve placed at a station, nested,
  measured end to end in their pieces' common convention, each joint a
  seam. No API change.

## [0.3.6] - 2026-10-09

### Added

- Nodes placed at stations (#264, ADR 0082 amendment):
  `InstanceAtStation { source, station: OrientedCurveStation }`
  (`#[non_exhaustive]`, built by `InstanceAtStation::new`) and the node
  variant `GeometryNode::InstanceAtStation`, appended last. It reuses a
  curve, solid or surface in the frame of an oriented station, the
  station kept symbolic with its frame, orientation and seam side; the
  source's local `x`, `y`, `z` map onto the oriented tangent, left
  lateral and up, its origin onto the station's point (a linear
  placement's reading), as the `station` module documents. A placed
  curve is a 3D curve whatever its source (accepted where a curve or 3D
  curve is, refused in 2D slots), a placed solid or surface keeps its
  family, a placed profile is no longer a profile. The station is
  validated as an oriented station's is when the node is pushed; the
  node references its source, then the station's basis.

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
