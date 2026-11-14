# 0103 — B-rep faces on swept and ruled surface relations

- **Status:** Accepted
- **Date:** 2026-11-07
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amended by:** [0110](./0110-edges-without-pcurves-are-inverted-onto-surface-relations.md): an edge use without a pcurve is inverted onto the relation (D6).
- **Amends:** [0097](./0097-ruled-surfaces-are-a-graph-relation-between-two-curve-sides.md) (a B-rep face on a ruled surface)

## Context

Issue #411 (consumer: openbimrs/revit#108). Revit's surfaces of revolution
and ruled surfaces arrive as graph relations (`SurfaceRelation::Revolution`,
`LinearExtrusion`, `Ruled`), and their faces carry exact iso-line edges. The
B-rep tessellator was typed on the atomic `axiolid_surface::Surface` and
refused any face over a relation (ADR 0097: "a face on one needs its own
tessellation and a trim"), so the 1,906 revolved faces of one Revit model
could not be meshed.

`axiolid-model` did not document how `Revolution` and `LinearExtrusion`
parameterise their points, and nothing in the kernel evaluated them. The
consumers agree: the IFC adapter (`IfcSurfaceOfRevolution`,
`IfcSurfaceOfLinearExtrusion`), Revit's `SurfRev` and OCCT's
`Geom_SurfaceOfRevolution` and `Geom_SurfaceOfLinearExtrusion` all take the
ISO 10303-42 form.

## Decision

**D1. Parameterisation.** Written into the model's documentation:

- `Revolution`: `S(u, v) = O + R_a(u) (c(v) - O)`, `c` the swept curve, `O`
  the axis origin and `a` the axis direction normalised (its length is not
  read). `R_a(u) w = (w.a) a + cos u (w - (w.a) a) + sin u (a x w)` turns
  right-handed about `a`; `u` is in radians, `u = 0` is the swept curve
  itself, `v` is the curve's own parameter.
- `LinearExtrusion`: `S(u, v) = c(u) + v d`, `u` the curve's parameter, `v`
  in units of the direction's length (it is not normalised).
- `Ruled`: ADR 0097, unchanged.

The turn is about the axis as given, so a producer must emit the axis that
matches its own sense of `u`: Revit's `Surface::Revolved` turns about `x × y`
of its frame, so a lowering must emit `axis_direction = x × y` (a mirrored
frame reverses `u`). A reflecting placement on the IFC side reverses `u`
relative to the pcurves; the adapter must negate the axis or refuse.

**D2. A crate-private support view.** `axiolid-mesh-compile` reads a face's
surface as a `Support`: an atomic surface, or a `RelationSurface` that
evaluates the three relations. The tessellator asks it only for a point, the
parameter periods, a family name, identity, and, through a `SupportOracle`,
the `SurfaceBounds` it already reads from an atomic surface. The curved-face
machinery (the per-triangle Taylor bound, the pcurve lenses and gaps, the
refinement) is unchanged. No dependency edge, no public API and no closure
changes; ADR 0097's `pub(crate)` evaluator stands (no consumer outside the
crate needs one yet).

**D3. Bounds, in closed form from the swept curve.** Rotations preserve
lengths, so over a box `v0..v1`, with `d` the certified derivative bounds of
the swept curve:

- *Revolution:* `|S_u| = |S_uu| = dist(c(v), axis)`, `|S_v| = |c'|`,
  `|S_uv| = |c'_perp|` (the part of `c'` normal to the axis; a line's is
  exact, else `d.first`) and `|S_vv| = |c''|`. The distance from the axis is
  read at the box's middle and widened by half the box times the normal
  speed.
- *Extrusion:* `|S_u| = |c'|`, `|S_v| = |d|`, `|S_uu| = |c''|`, and no mixed
  or `v` second derivative.
- *Ruled:* per side, `a1 = k sup|c'|` and `a2 = k^2 sup|c''|` over the box,
  `k = |t1 - t0|`; `|S_uv| <= |B' - A'|(um) + h/2 (a2_A + a2_B)` at the
  middle `um` of a box of width `h`, `|S_v| <= |B - A|(um) + h/2 |S_uv|`,
  `|S_u|` and `|S_uu|` a blend of the sides' bounds, largest at an end of the
  `v` range, and `S_vv = 0`.

`smooth` is false where a tangent break of the swept curve (`continuity_
breaks3`; a side's, mapped to `u`) lies strictly inside the box, so the
tessellator falls back to its first-order bound there. Every value is
finite-checked and widened by `1 + 1e-12`.

**D4. Domain.** A B-spline swept curve is read only in its domain, within a
slack of `1e-12 max(1, span)` that is clamped; a read beyond is
`InvalidInput`. A box reaching past the domain is bounded over its part
inside, since evaluation refuses the rest. A ruled surface is read on its
unit square the same way. A conic is periodic and is read anywhere.

**D5. Periods.** A revolution has `u` period `2 pi`, and `v` period `2 pi`
when its swept curve is a circle or ellipse; an extrusion has `u` period `2
pi` for a conic, `v` none; a ruled surface has none.

**D6. Refusals.** The swept curve is a line, circle, ellipse or B-spline (the
set ADR 0097 reads for a side). A swept curve that is not a `Curve3` node is
`UnsupportedInput` ("a swept curve that is not an atomic 3D curve"), as is any
other `Curve3` ("a swept curve other than a line, circle, ellipse or
B-spline"). A zero or non-finite axis or direction is `InvalidInput`. A face
on a relation needs a pcurve on every edge use: an edge use without one is
`UnsupportedInput` ("an edge use without a pcurve on a face over a surface
relation"), raised before any meshing, because inverting an edge onto a
relation is not done (#454). Amended by
[0110](./0110-edges-without-pcurves-are-inverted-onto-surface-relations.md):
such a use is inverted onto the relation along its swept parameter, and the
refusal is gone. Any other relation (curve-bounded, rectangular
trimmed, offset, sectioned) stays `UnsupportedInput` ("a B-rep face on a
surface relation other than a revolution, a linear extrusion or a ruled
surface"; rectangular trims are #455). Exact compilation is unchanged: a
B-rep with relation faces fails the polyhedral-warp check of the referee and
the exact compiler refuses a bare relation, so the mesh stands.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Lower a relation to an atomic `Surface` where it coincides (cylinder, cone, sphere, torus, plane) | Only special placements; the general case is needed anyway, and two evaluations of one relation could disagree. |
| Convert to a NURBS surface | A conic's angle is not rational (ADR 0097), so pcurves would need reparameterising. |
| A public evaluator and oracle in `axiolid-evaluate` | ADR 0097 kept the ruled evaluator `pub(crate)` for want of a consumer; revisit when measure, the exact compiler or hidden lines need relation surfaces. |
| A trait object in the tessellator | No gain over a closed `Copy` enum. |
| Include the inversion of edges without pcurves | Roughly doubles the work and adds certified inversion code; Revit's iso-line edges carry exact pcurves. Filed as #454. |

## Consequences

**Positive**

- Revit's revolved and ruled faces mesh with a certified bound, checked
  against closed-form oracles and against the atomic cylinder, cone, sphere
  and torus they coincide with.
- The model keeps a single input seam; the refusals name what is missing.

**Negative / costs**

- A face over a relation needs pcurves; IFC advanced faces on swept surfaces
  carry none until #454 lands.
- A Hermite profile's piecewise-Bezier nodes are tangent breaks, so
  triangles straddling a node use the first-order bound: sound, but it may
  report over the budget, as an atomic B-spline surface does today.
- A hole on a curved face is bounded above the 1e-3 budget on the fixture
  used (an atomic cylinder with a circular pcurve hole gives the same
  1.47e-3); the relation matches the atomic bound exactly.

**Follow-ups / risks to watch**

- #454: invert edges without pcurves onto relation faces.
- #455: rectangular trims over swept surfaces.
- Curve-relation and elevated or banked swept curves.

## Relation to existing code

`crates/execution/compile/src/surface_relation.rs` (`RelationSurface`,
`RelationOracle`), `src/brep/support.rs`, `src/ruled.rs` (`Side`),
`tests/relation_faces.rs`, `tests/ruled_surface.rs` and
`scripts/probe_relation_faces_mutants.py`. The model's
`SurfaceRelation::Revolution` and `LinearExtrusion` document the
parameterisation. The ledger rows B5 and C8 record the level.
