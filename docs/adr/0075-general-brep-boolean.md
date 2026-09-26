# 0075 — General B-rep boolean: our own general-fuse pipeline

- **Status:** Accepted
- **Date:** 2026-09-26
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

Decides the approach for axiolid/kernel#167 (ledger row C9); implementation
follows #119.

## Context

Exact booleans today cover planar-faced polyhedra (`boolean_polyhedra_exact`)
and vertical columns: coaxial arc prisms, stepped spans, sloped cuts and
cavities (ADR 0070–0072, #120). Everything else in C9 is refused by name:

- curved solids that are not coaxial: two cylinders at an angle, a round
  column meeting a sloped wall, a pipe tee;
- solids bounded by cones, spheres, tori or B-spline faces: domes,
  revolved profiles, fillets.

None of these reduce to columns. They need the general pipeline OCCT
implements in `BOPAlgo`/`BRepAlgoAPI`: intersect every face pair, split each
face by the curves it meets, classify the pieces against the other solid,
select by operator, and sew the result.

What exists to build on:

- **Exact intersection curves** (#119): lines, circles and ellipses where
  planes meet planes, cylinders, spheres and cones (circle/ellipse cases),
  where spheres meet, and for every coaxial pair of surfaces of revolution.
  Hits, tangency and containment are decided exactly on `axiolid-exact`.
  General quadric/quadric curves (quartics: a pipe tee) are refused.
- **Exact curve/curve and curve/surface intersection** for lines, circles
  and ellipses (#119), which is what splitting a face along several section
  curves and classifying a point by ray parity need.
- **Oracles that did not exist a month ago:** exact volume, centroid and
  moments of any curved solid (ADR 0073) and a certified boundary distance
  (ADR 0074). A boolean result can now be checked by the identity
  `|A u B| + |A n B| = |A| + |B|` exactly, not through a mesh.
- **Void shells and multi-solid results** are carried end to end (#120,
  #111).

## Decision

We build **our own general-fuse boolean** in a new crate,
`axiolid-brep-boolean` (layer `algorithms`, role `algorithm.construction`),
in stages gated by which intersection curves #119 can construct exactly. We
do not adopt or wrap another kernel.

### Pipeline

1. **Broad phase.** Face pairs whose bounding boxes overlap (the certified
   patch spheres of ADR 0074 serve for curved faces).
2. **Section curves.** `exact_surface_intersection` on each pair's support
   surfaces; each resulting line, circle or ellipse is trimmed to where it
   lies inside BOTH faces' domains, using the certified domain classifier
   (ADR 0074) and exact curve/pcurve-boundary intersections. Tangent and
   coincident faces are classified exactly (`ExactIntersectionCurve`
   already distinguishes them) and handled as shared regions, never as
   near-misses.
3. **Face splitting in 3D, pcurves derived.** All curves on one face are
   conics in 3D, so the face's arrangement is built from exact 3D
   curve/curve intersections, with parameters along each curve; the face's
   own loops join the arrangement the same way. Each new edge's pcurve is
   the inverse image of its 3D curve on the face's surface. Where that image
   is an existing `Curve2` family (line, circle, ellipse, sinusoid on a
   cylinder, latitude on a sphere) it is stored as such; otherwise a new
   `Curve2::OnSurface` variant carries the 3D curve and the surface and
   evaluates the inverse exactly on demand. That variant needs its own ADR
   before it lands.
4. **Classification.** Each split piece is classified in/out/on the other
   solid by ray parity against the other solid's exact faces
   (`exact_curve_surface_intersection` on a line), retrying from a fixed
   set of directions when a ray meets an edge or a tangency, and refusing
   by name when every direction is degenerate — the discipline
   `boolean_polyhedra_exact` already uses.
5. **Selection and assembly.** Pieces are selected by operator, coincident
   pieces kept once by normal agreement, and sewn along shared edges into
   shells; shells with negative enclosed volume become voids of the solid
   that contains them, and disconnected pieces separate solids (#111).
6. **Verification in every test:** geometric audit, closed manifold, and the
   exact volume identity via `exact_properties`.

### Stages

- **Stage 1 (after #119's refusals are named, which they are):** planes
  with cylinders, spheres, cones and tori where the section curves are
  lines, circles or ellipses — a round column meeting a sloped wall, a dome
  cut by a plane, coaxial revolved solids, spheres with spheres. Every other
  pair is refused by the name `exact_surface_intersection` gives it.
- **Stage 2 (after #119 constructs quadric/quadric quartics):** two
  cylinders at an angle, a pipe tee, cone/cylinder junctions. Needs a
  certified representation of the quartic section curve, decided in #119.
- **Stage 3:** B-spline faces, on the certified bounded intersection tier
  that already exists for NURBS pairs.

`boolean_arc_prisms_exact` and the column builder stay as the fast exact
path for vertical columns; the general pipeline must agree with them where
both apply, which becomes a differential test.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Adopt `truck-shapeops` (pure Rust, Apache-2.0) | Its B-rep is NURBS-only: analytic cylinders, spheres and tori would be converted to splines and back, which is exactly the loss of exactness Axiolid exists to avoid. Its crates.io releases are two years behind `master`, and its whole type vocabulary sits on `cgmath`, unmaintained since 2021 (`docs/research/geometry-kernel-landscape.md`). |
| Wrap OCCT (`opencascade-rs`) | C++ and LGPL; the workspace has no C++ dependency path. |
| Mesh boolean (`boolmesh`) behind the exact API | Returns a tessellation; the exact path must never substitute a mesh for exact output. It stays the discrete tier. |
| Extend the column builder | Only vertical columns reduce to columns (ADR 0072); none of the #167 cases do. |
| Wait for a general quartic representation before starting | Stage 1 covers the common building cases (round columns against sloped walls, domes, revolved elements) with curves #119 already constructs exactly. |

## Consequences

**Positive**

- C9 moves towards `implemented` stage by stage, each stage with its
  refused subset named by the intersection it lacks.
- Every result is checkable exactly (ADR 0073), not only against a mesh.

**Negative / costs**

- OCCT-scale work: `BOPAlgo` is about 24k lines. The staged plan keeps each
  landing reviewable, but the pipeline is the largest single construction
  in the kernel.
- `Curve2::OnSurface` adds a pcurve family every consumer of `Curve2` must
  handle (it is `#[non_exhaustive]`, so the change is additive).
- Classification by ray parity needs robust retries; degenerate inputs
  (every direction grazing) are refused, not guessed.

**Follow-ups / risks to watch**

- #119: quadric/quadric quartic curves (gates stage 2).
- ADR for `Curve2::OnSurface` before stage 1 needs it.
- Tolerant (near-coincident) inputs from real IFC models: faces within
  tolerance but not exactly coincident need a merge policy; stage 1 treats
  them as distinct and refuses slivers below tolerance by name.

## Relation to existing code

- `crates/algorithms/parametric/nurbs/src/exact_surface_intersection.rs`,
  `exact_curve_surface_intersection`, `exact_curve_curve_intersection3`: the
  intersections the pipeline consumes.
- `crates/algorithms/query/measure/src/exact_domain.rs`: certified face
  domains for trimming section curves.
- `crates/algorithms/query/measure/src/exact.rs`: the volume identity.
- `crates/algorithms/construction/construct/src/polyhedron.rs`: the planar
  general boolean whose classification discipline this generalises.
- `crates/algorithms/construction/construct/src/column.rs`: the column path
  kept alongside.
