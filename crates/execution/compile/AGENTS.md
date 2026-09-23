# axiolid-mesh-compile instructions

Purpose: the reference mesh and exact compilers. `ReferenceMeshCompiler` orchestrates
`GeometryGraph` to `TriMesh`; `ReferenceExactCompiler` owns a separate per-batch
`NodeId -> ExactBRep` cache and never accepts discrete values. Graph traversal,
transform composition, and operation dispatch live here; construction algorithms
and their local invariants remain owned by `axiolid-construct`.

## Invariants

Extrusion output must be **closed, edge-manifold, and outward-oriented**,
because that is exactly `axiolid-mesh-boolean-boolmesh`'s input precondition. Volume alone does
NOT verify this: a cap lying in the z = 0 plane contributes nothing to the
divergence integral, so a flipped bottom cap is invisible to a volume check.
Use the directed-edge parity gate in
`crates/algorithms/construction/construct/tests/extrusion.rs` — every directed
edge exactly once, every edge with exactly one opposing half-edge.

Unsupported profile and solid families return `GeomError::UnsupportedInput` naming
the family, never a silent approximation. Exact compilation currently accepts only
supported extrusion roots and must never delegate to mesh compilation.

No default tolerance or chord budget. The caller supplies both, because
acceptable error depends on source units and downstream use.

**Channels ride with the geometry (#115).** The cache holds
`channels::Built` (mesh + per-channel fates), not a bare `TriMesh`, so every
node kind must say what it did to each channel. `Instance` and `Collection`
never create surface points and so never derive a value; do not rebuild a
mesh from positions and indices on those paths -- that is the #115 bug. New
graph paths go through `channels::{transform, merge, after_boolean}` or wrap
a freshly made mesh in `Built::leaf`. Gate: `tests/graph_channels.rs`.

Curve flattening is **not owned here**. `segment_points`, `circle_rings`, and
`ellipse_rings` all delegate to `axiolid_reference::curve::flatten2` (ADR 0018),
which subdivides adaptively on measured sagitta. The old private
`circle_segments`/`circle_ring` pair is gone — do not reintroduce a
closed-form segment count, it only models circles and cannot express an
ellipse or a rational spline.

`crates/algorithms/construction/construct/tests/extrusion_volume.rs` pins the
identity `volume == area * depth` for every supported profile family and asserts
the chord budget actually bounds the volume error (measured: error is O(chord),
constant under 5). Volume and area come from `axiolid-measure`, never a local
divergence sum: that crate audits closed-two-manifold first, so a hand-rolled
integral would silently measure a torn shell.

**Tolerance must scale with the chord budget.** `audit_mesh` calls a triangle
degenerate when `2A <= tolerance.linear()^2`. A cylinder flattened at chord
`c` has side quads about `sqrt(8*r*c)` wide and cap slivers far smaller, so a
fixed `Tolerance::MILLIMETRE` rejects perfectly correct geometry as soon as a
caller asks for sub-millimetre accuracy. Use `tolerance_for(chord)` in tests;
in production pass a tolerance derived from the same budget that drove
flattening. This is not a test artefact -- it is a real API contract.

## Adopted dependencies

`earcut` (ADR 0015) is owned by
`crates/algorithms/construction/construct/src/profile.rs` and is not re-exported.
`axiolid_reference::triangulate_simple` audits it differentially on hole-free
polygons in `crates/algorithms/construction/construct/tests/oracle.rs` — the
adopted crate is verified, not trusted. This crate's graph-level integration
coverage lives in `tests/pipeline.rs`, `tests/generation_boolean.rs`,
`tests/brep_tessellation.rs`, and the other current `tests/*.rs` targets.

## Layer

L3, an implementation crate alongside `axiolid-backend-cpu` and `axiolid-mesh-boolean-boolmesh`.
It may depend on representation and contract crates; nothing in L0–L2 may
depend on it.
