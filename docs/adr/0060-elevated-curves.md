# Elevated curves: composing a planar layout with an elevation law

Status: accepted
Date: 2026-09-14

Closes axiolid/kernel#105.

## Context

`Curve2::Intrinsic` carries a transition spiral exactly: a `CurvatureLaw`
anchored to a start frame over an arc length, covering the clothoid, Helmert
and sine-corrected families. `Curve3` had no counterpart, and `CurveRelation`
composed existing nodes without any way to pair a planar curve with an
independent elevation law.

A road or rail centreline is authored as exactly that pair: a horizontal
layout carrying the spirals, and a vertical profile giving height as a
function of distance along that layout. With no 3D form, a consumer holding
an exact plan and an exact profile had to approximate the spiral as a
`BSpline3`, keep the halves apart and make every consumer recompose them, or
refuse. The kernel refused.

## Decision

Add `Curve3::Elevated(Elevated3)`, holding a boxed `Curve2` plan and an
`ElevationLaw`. This is option 2 from the issue: composition, not re-encoding.

Two things fell out of reading the existing code that changed the shape of
the work.

**`Intrinsic2` had no evaluator at all.** It exposes `total_turning` and
documents that *position does not integrate in closed form, which is why this
method exists and an `evaluate` does not*. The `evaluate2` dispatch has no
`Intrinsic` arm and refuses by name. So the issue's third criterion -- evaluate
a point and frame at a distance -- was unmet in 2D as well, and a 3D
composition alone would not have delivered it. `arc_length.rs` now supplies
`intrinsic_point`/`intrinsic_tangent`, and the 3D path is a thin layer over
them.

**Height is a function of PLAN distance, not 3D arc length.** The vertical
profile a surveyor writes is indexed by chainage on the horizontal layout.
The two differ whenever grade is non-zero, since `ds3 = sqrt(1 + g^2) ds_plan`.
Measured on a 120 m curve at 2% entry grade: 3D length exceeds plan length by
12.5 mm, which reads the profile 0.35 mm off. Small, but systematic and
compounding along a chain of segments, so the convention is named in the type's
documentation rather than left to the caller.

### Position is quadratured; the value stays exact

Heading is the integral of curvature and is exact in closed form for every law
in the family -- `heading_at` reuses the existing `turning_over`. Position is
the integral of `(cos th, sin th)` and has no elementary antiderivative; for
the clothoid it is the Fresnel integral. Gauss-Legendre 8-point discharges it,
subdivided by total turning so a tight spiral gets more panels, bounded so a
malformed law refuses rather than hangs.

This is the same bargain as evaluating `sin`: the stored curve is not
approximated, its evaluation is computed to tolerance. Against a Fresnel series
computed independently in the test, a 120 m clothoid into R=300 agrees to
better than 1e-9.

### Only arc-length parameterisations may carry an elevation law

`Line`, `Circle` and `Intrinsic` qualify. A B-spline's parameter is not arc
length, so pairing one would silently mean something other than what the law
says; it refuses. This keeps the composition honest rather than accepting
anything shaped like a curve.

## Consequences

- An exact spiral plan and an exact vertical profile compose into one
  `Curve3` with neither approximated, and either half is recoverable unchanged.
- The composition is an ordinary `GeometryNode::Curve3`, so it round trips
  through `GeometryGraphBuilder::finish` and validates as a 3D curve with no
  special case in the graph.
- `ElevationLaw::Piecewise` mirrors `CurvatureLaw::Piecewise`: each piece is
  written in its own distance, restarting at zero at its seam, so moving a
  piece never rewrites its coefficients. A seam belongs to the piece that
  starts there.
- Mismatched piece lists report `None` rather than guessing, matching the rest
  of the curve crate.
- `IfcGradientCurve` and `IfcSegmentedReferenceCurve` now have a target
  representation. Lowering them stays with the consumer; this ADR adds the
  value they lower into, not road semantics.
- Torsion is still not representable. An alignment does not state torsion, so
  `Curve3::Intrinsic` with curvature and torsion laws (option 1) remains open
  for a future need rather than being built speculatively.

## Verification

- 9 evaluation tests, 1 graph round-trip test, 1586 workspace tests passing.
- Clothoid position pinned against an independently computed Fresnel series,
  not against another run of the kernel.
- A constant-curvature law is checked against the elementary quarter-arc, and
  a straight law against the line.
- Mutation testing, 5/5 killed: grade dropped from the tangent; `sqrt(1+g^2)`
  normalisation removed; spiral plan refused; Horner degenerated to a sum;
  piecewise piece not rebased to its own zero.
