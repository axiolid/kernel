# Variable-radius and multi-corner fillets

Status: accepted
Date: 2026-09-12

## Context

`fillet_polygon_corner` blends ONE corner at a CONSTANT radius. Two gaps
remain from the thesis feedback: filleting several corners in one call,
and letting the radius vary along the extrusion.

Multi-corner is a generalisation of existing machinery. Variable-radius is
not: it changes the SURFACE TYPE of the blend, so it was measured before
any code was written.

## What a tapered blend actually is

With radius varying linearly, `r(z) = r0 + k z`, at each height the blend
section is a circular arc whose centre sits at `r(z)/sin(theta/2)` along
the corner bisector. Measured, not assumed:

- The surface IS a cone: every ruling passes through a single apex at
  `z = -r0/k`, to `2.3e-16`.
- The surface is NOT a right circular cone. The minimum achievable spread
  of ruling angles about ANY axis is `0.045 rad` at `theta = 90`, found by
  direct optimisation over the axis rather than by plane fitting.
- Obliqueness scales with taper and vanishes only as `k -> 0`, which is the
  constant-radius cylinder already handled.

`Surface::Cone` is right circular (`frame`, `radius`, `semi_angle`), so it
CANNOT represent a tapered blend. Storing one would be a silent lie.

## Decision

Represent the tapered blend as a RATIONAL B-spline surface: degree 2 in
`u` (the standard 3-point rational-quadratic arc) times degree 1 in `v`
(a linear loft between the bottom and top arcs). Both the centre and the
radius are affine in `z`, so the loft is exact rather than approximate.

Measured radial error of the lofted surface against the true fillet:
`1.1e-16`. The weights are the conic weights `cos(sweep/2)` and are
constant along `v`, which is what makes the linear loft exact.

The wall boundary stays straight: setback is `r(z)/tan(theta/2)`, linear
in `z`, so the adjacent walls remain planar. Straightness residual
`5.6e-16`. This is why the tapered fillet does not force curved walls.

## Multi-corner

`fillet_polygon_corners` takes `(corner, radius)` pairs. Corners are
solved independently, then all replacements are applied in one pass.
Independence is only valid when no two blends overlap, so the setbacks of
adjacent filleted corners are checked against the shared edge length and
the call is refused when they collide -- rather than emitting a
self-intersecting ring.

## Consequences

- A tapered fillet yields a `Surface::BSpline`, a constant one yields
  `Surface::Cylinder`. Callers matching on surface type see the difference,
  which is honest: the geometry genuinely differs.
- Exactness claim: the NURBS blend is exact in the same sense as the
  arc path, agreement to machine precision, not integer-predicate exact.

## Findings

### The tapered blend is not a cone

The obvious guess is that a linearly tapered fillet sweeps a cone, so it could
reuse `Surface::Cone`. Measured, it does not. The rulings do all meet a single
apex -- residual `2.3e-16`, so it is a cone in the general sense -- but it is
an OBLIQUE cone: an independent optimisation over all candidate axes found a
minimum half-angle spread of `0.045 rad`, far above numerical noise.
`Surface::Cone` is right-circular (`frame`, `radius`, `semi_angle`) and cannot
represent it.

The obliqueness scales with the taper and vanishes only as `k -> 0`, which is
the constant-radius cylinder already handled. Measured spread by interior
angle at `k = 0.2`: 60 deg `0.0696`, 90 deg `0.0451`, 120 deg `0.0212`,
150 deg `0.0054`. So there is no useful sub-case where a right circular cone
would do; a fixed shortcut would be wrong everywhere except the limit.

### A rational quadratic loft is exact

Every horizontal section of the blend is a full circular arc, and a circular
arc is exactly a rational quadratic Bezier with shoulder weight `cos(sweep/2)`.
Lofting two such arcs linearly in `v` reproduces the closed-form fillet to
`1.1e-16`. This is exact representation, not approximation: no tolerance
parameter and no tessellation appear anywhere in the construction.

The surface is therefore `Surface::BSpline` with `u_degree = 2` rational,
`v_degree = 1`, a 3x2 control net and weights `[1, cos(sweep/2), 1]` repeated
along `v`.

### Setback grows with radius, and the wall boundary stays straight

Because setback is `r/tan(theta/2)` and `r` varies with height, the seam
between blend and wall is not vertical -- it is a slanted line. Worth checking
that it is still a straight line, since a curved seam could not be a planar
wall boundary. Measured residual from straightness: `<= 5.6e-16` at 60, 90 and
120 degrees. So the adjacent walls remain planar quadrilaterals and no wall
needs to become a spline.

### Multi-corner needs a pairwise check, not a per-corner one

The single-corner path checks each setback against the whole adjacent edge.
That is insufficient once two filleted corners share an edge: their setbacks
must sum to less than the edge length. A concrete case that passes per-corner
checks and still collides: edge length 2.0 with radii 1.1 and 0.8, whose
setbacks sum to 1.9 -- fine -- versus radii 0.5 and 1.6, summing to 2.1, which
overruns. Corners that are not being filleted contribute zero setback, so the
check is over the requested set, not all corners.
