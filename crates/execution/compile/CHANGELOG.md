# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- A solid's void shells are tessellated (#120). They were dropped as
  "boolean intent", which silently filled every cavity of an authored
  `IfcFacetedBrepWithVoids`-style B-rep and forced the exact booleans to
  refuse cavities. Each void is emitted facing into the cavity, so the mesh
  encloses the outer volume less every cavity; a void authored facing out
  of it (the STEP convention, reversed on use) is turned round, since a
  cavity can only remove material. A void shell that is not closed is
  refused rather than leaving the mesh open.

## [0.3.3] - 2026-09-25

### Added

- The reference compiler honours `ExecutionOptions::with_chord_error`
  (#165) everywhere it flattens a curve: profile arcs, circles and ellipses,
  sweep directrices, curved B-rep faces and edges, and CSG primitives. An
  instance scales the budget with its transform, like the tolerance, so it
  stays a world-space distance. Without a budget the chord error is the
  linear tolerance, exactly as before. Measured on a 5 mm disc extruded 1 m
  at `Tolerance::MILLIMETRE`: 10 % short by default, 2.6 % at a 0.1 mm
  budget, 0.16 % at 10 um and 0.01 % at 1 um.

### Fixed

- A surface model with a face whose outer bound encloses no area
  tessellates (#171): that face covers nothing, so it is skipped instead of
  refusing the whole model with `planar face bound has zero or non-finite
  area`. Real Nova MEP exports write pipe-fitting end caps as bowtie quads
  through the pipe axis (an annulus with a negative inner radius), each with
  signed area exactly 0; 42 fittings, pumps and valves in two models were
  refused over them. A declared solid still refuses such a face, and a
  zero-area hole or a non-finite bound is still refused everywhere.
- Closed authored meshes stay closed (#170). Planar faces of a
  `PolygonMesh` and of a B-rep were triangulated with earcut, which drops
  corners on a straight run and runs diagonals and hole bridges over corners
  of the same face; the neighbouring face still split that edge at the
  corner, so the mesh cracked (T-junctions). Every edge earcut invents is now
  split at each face corner on it, within a band of a thousandth of the
  linear tolerance (1 um at `Tolerance::MILLIMETRE`), so export noise on
  shared corners (1e-8 to 5e-8 m on real files) is judged the same on both
  sides. Authored ring edges are never split, and a thin triangle whose long
  side is authored is kept. Closure of an authored mesh is now read from its
  index connectivity instead of `audit_mesh`, which dropped real faces below
  its area threshold before counting edges. On two real ArchiCAD models this
  turns 505 authored-closed `IfcPolygonalFaceSet` products from `Surface`
  into `Solid`; no product that compiled before is refused.

## [0.3.2] - 2026-09-25

### Fixed

- A directrix trimmed from a circle or ellipse ACROSS its seam sweeps the arc
  the trim names (#168). The trimmed curve runs from `start` the way
  `sense_agreement` says, wrapping past the seam if it must; the directrix
  path used to sort the two trims and sample the complementary arc. Standalone
  that was silently wrong geometry (a `315 -> 45` degree bend swept the 270
  degree arc); inside a composite the ends no longer met and the sweep was
  refused as `composite directrix has a N unit gap`. In a real Revit rebar
  model that was 1,494 bent bars, all writing a bend as `270 -> 45`,
  `270 -> 15` or `270 -> 360` rounded just past the seam. A full turn rounded
  past the seam stays a full turn. A sweep `parameter_range` on such a trim is
  read in the trim's unwrapped interval, so `(330, 30)` and `(330, 390)`
  degrees name the same sub-arc of a `315 -> 45` trim; an end off the arc is
  refused as before. Profiles already honoured this.

## [0.3.1] - 2026-09-24

### Added

- `PolygonMesh` faces that are not plain triangles compile (#160): n-gons,
  concave faces and faces with holes (IFC4 `IfcIndexedPolygonalFaceWithVoids`)
  are triangulated in their own plane, keeping the authored positions and
  winding. Plain triangles keep their exact corner order as before. A face
  whose corners leave its plane by more than the linear tolerance, that has
  no area, or whose rings cross is refused with an error naming its index.
- B-reps with shells but no solid tessellate (#161): every shell is
  tessellated as authored and the result is reported as
  `MeshClosure::Surface` through `compile_mesh_reported`, even when the
  shell is closed. Collections are `Solid` only if every member is, and a
  boolean with a surface operand is refused. Authored meshes report
  `Solid` exactly when they are closed, consistently wound two-manifolds.

### Changed

- A `PolygonMesh` with non-triangular faces used to fail with
  `GeomError::Unsupported`; it now compiles. A B-rep with no solid and no
  shell is refused as "neither a solid nor a shell" instead of "has no
  solid".

## [0.3.0] - 2026-09-23

### Fixed

- `Instance` and `Collection` nodes no longer drop attribute channels and normals (#115). Both used to rebuild the mesh from positions and indices only, so a textured item lost its `uv` channel as soon as a product had a second item or was instanced — silently.
  - `Instance` carries channels through unchanged and now transforms normals by the inverse transpose instead of dropping them. Under a mirroring transform, corner-indexed channels and normals swap corners with the triangle.
  - `Collection` merges channels by name. A channel on only some members becomes corner-indexed, with the other members' triangles `UNMAPPED`; a channel on every member as per-vertex stays per-vertex. Members defining one name with a different width or blend drop it as `DropReason::IncompatibleChannels`.
  - Booleans keep the provider's channel fates, composed onto what each operand already went through, instead of discarding the evidence.

### Added

- `ReferenceMeshCompiler` implements `MeshCompiler::compile_mesh_reported`: the compiled mesh plus each channel's fate on its way to the root (worst over parallel members, sequential through booleans).
