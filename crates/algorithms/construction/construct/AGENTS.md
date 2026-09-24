# `axiolid-construct`

Scalar geometry-construction algorithms over neutral representations. This is **L2**:
it accepts exact profiles, curves, primitives, certified NURBS traces, and explicit
policy; it creates the broad `TriMesh` reference path plus focused exact extrusion and
analytic B-rep arrangements; it owns no DAG, cache, execution context, or operation-provider dispatch.

## Entry points

- `profile`: lower profile values to sampled rings and triangulate caps.
- `extrude`: mesh extrusion plus exact extrusion of every `Profile` variant --
  rectangle (sharp/hollow), circle, ellipse, contour (arcs and holes), section,
  centre-line, derived and composite. Exact output owns every 3D support,
  pcurve, and native span.
- `revolve`: exact full-turn revolution of any profile that lowers to a contour,
  sweeping cylinders, cones, planar annuli and tori; `sweep`, `loft` place/stitch
  station rings into discrete solids.
- `center_line`: turn constant-width centre-line profiles into rings.
- `half_space`: construct a finite clipping proxy for an unbounded half-space.
- `trimmed_intersection`: promote one certified affine trace into two closed trimmed
  faces on its boundary-owned patch and an explicit embedded pcurve on the containing
  unsplit face; see ADR 0029.

`BACKEND_ID` is `scalar-generate`. Use it for every diagnostic raised here; do
not report `scalar-compile` after this split.

## Invariants

- Every tolerance-sensitive entry point receives an explicit `Tolerance` and/or
  chord budget. Never introduce a global default epsilon.
- Profile mismatch, unbounded construction, invalid frames, insufficient rings,
  and unmet subdivision budgets refuse with `GeomError`; never invent a mesh.
- Shared loft/stitching logic owns winding and cap pairing. Do not duplicate it
  in individual sweep families.
- This crate must not depend on `axiolid-model`, an execution/backend crate, or
  any L3 crate. `cargo xtask architecture check` enforces the declared internal
  dependency allowlist and production role-DAG edge; `scripts/probe_layering_gate.sh`
  mutation-verifies that enforcement.
- Discrete sweeps remain the broad reference path for `sweep` and `loft`. Exact
  solid output covers every profile variant for extrusion, and full-turn
  revolution for any profile that lowers to a contour; the certified affine
  trimmed-intersection arrangement is a separate exact surface slice.
  What still refuses is geometry the kernel cannot represent exactly rather
  than unwritten work -- partial-turn revolution, a profile straddling the
  revolution axis, oblique circle/ellipse extrusion, non-conformal derived
  transforms, disjoint composite members, and Boolean families beyond
  vertical columns -- and must refuse rather than tessellate; see ADR 0020,
  ADR 0023, ADR 0024, ADR 0029, and ADR 0053-0059.
- Coaxial booleans and plane cuts that are not one prism (stepped spans,
  a cut crossing a cap) go through `column.rs` over one exact
  `ArcArrangement` (ADR 0072); `boolean_column.rs` adapts the entry points.
  Every face of a column solid takes its vertices from that one
  arrangement -- do not build bands separately and glue them, the rounded
  crossings will not agree. Cavities are refused because tessellation
  reads only the outer shell.

## Tests

Unit-like generation tests live in `tests/` here. Tests that verify a generated
mesh is accepted by an L3 Boolean provider stay in `axiolid-mesh-compile/tests/`:
letting this L2 crate dev-depend on that provider violates the tier boundary.
