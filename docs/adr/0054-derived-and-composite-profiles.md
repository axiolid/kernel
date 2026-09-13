# Derived and composite profile extrusion

Status: accepted
Date: 2026-09-12

## Context

After ADR 0053, `extrude_profile_exact` accepted `Rectangle`, `Circle` and
`Contour`. `Derived` and `Composite` still refused, and they are the two
variants an imported model produces most: `Derived` is how the profile crate
records PLACEMENT (its `AGENTS.md` states "Profile placement uses Derived, not
baked tessellation"), and `Composite` is the `IfcCompositeProfileDef` case.

## Decision

Lower both onto the concrete builders rather than teaching each extruder about
transforms or member lists.

`Derived` pushes its transform down onto the basis geometry and re-enters
`extrude_profile_exact` with the lowered profile. Nested derivations COMPOSE
their transforms instead of recursing on already-lowered output, so a
repeatedly re-placed profile costs one pass.

`Composite` unions its members through `union_soup` and extrudes the single
resulting polygon, outer ring plus holes.

## What the transform may and may not do

Measured before implementing:

| transform        | det   | conformal | circle survives |
|------------------|-------|-----------|-----------------|
| rotate 30 deg    | +1.00 | yes       | yes             |
| uniform scale 2  | +4.00 | yes       | yes             |
| non-uniform 2x1  | +2.00 | no        | NO -> ellipse   |
| shear            | +1.00 | no        | NO -> ellipse   |
| mirror x         | -1.00 | yes       | yes             |

Two consequences, both load-bearing:

- A circle survives only a CONFORMAL linear part. The test is `M^T M` being a
  positive multiple of the identity, NOT the determinant: the shear above has
  determinant exactly 1 and still maps the unit circle to radii spanning
  `[0.781, 1.281]`. A determinant check would admit it and report a wrong
  radius.
- A negative determinant MIRRORS, reversing ring orientation. The extruders
  assume counter-clockwise outer rings, so mirrored rings are reversed during
  lowering -- both the segment order and each segment's own sense -- rather
  than handed over inside-out.

## Why composite members are unioned, not concatenated

Members of one section may touch or overlap. Passing them as separate rings
to the multi-ring extruder would treat member 1 as the outer boundary and the
rest as holes, which is simply a different shape. Measured union behaviour:

- overlapping members -> 1 polygon, 0 holes
- edge-touching members -> 1 polygon
- four bars as a picture frame -> 1 polygon, 1 hole
- disjoint members -> 2 polygons

The last case is REFUSED. `Solid` holds one outer shell plus voids, so two
separate bodies have nowhere honest to live; returning either alone would
silently discard the other.

## Consequences

- A transformed rectangle lowers to a `Contour`, not a `Rectangle`: under a
  general affine map it is a parallelogram, which `RectangleProfile` cannot
  express. The contour path handles it exactly.
- Contour segments are transformed as GEOMETRY, not sampled. A line stays a
  line and an arc stays an arc; taking segment endpoints would silently turn
  every arc into a chord.
- A translated circle is refused: the circle extruder places its cylinder at
  the origin, so an off-origin circle has nowhere to record its centre.
- Composite members carrying arcs are refused rather than chord-sampled, since
  the union operates on polygons.
- `Section`, `Ellipse` and `CenterLine` still refuse. `CenterLine` needs an
  offset of an open path, which is a genuinely different operation.

## Verification

Volumes are checked against closed forms computed independently: the
overlapping pair reads `6 * depth` (not `8 * depth`, which double-counting
would give) and the picture frame reads `32 * depth` (not `36 * depth`, which
a dropped hole would give). Every built solid is also checked with the ADR
0052 geometric audit.

Mutation, 5/5 killed: conformality weakened to a determinant check; disjoint
composite silently taking the first polygon; composite holes dropped; nested
derived applying only the outer transform; circle scale ignored.
