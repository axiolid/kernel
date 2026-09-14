# Multi-ring arc caps and tapered sections

Status: accepted
Date: 2026-09-12

## Context

Two gaps were recorded when ADR 0053 and ADR 0057 landed.

1. The arc extruder built one ring only. A cap face carrying several
   bounds was not written, so `Profile::Contour` with holes was refused
   and arc sections could not have through-passages. The polygon path
   already supported holes, so the two paths disagreed on what a profile
   could express.

2. Parameterised sections refused a declared flange, web or leg taper.
   The stated reason was that a taper moves fillet tangency onto an
   inclined face.

## Decision

### One cap face, several bounds

`extrude_arc_rings` replaces `extrude_arc_ring`. Ring 0 is the outer
boundary and the rest are through-holes, matching the polygon path's
existing convention exactly: `outer: index == 0` on each `FaceBound`,
and one cap face per end carrying every ring's loop.

Winding is enforced rather than demanded. An imported contour may hand
holes in either direction, so `orient_arc_ring` reverses a ring whose
signed area has the wrong sign -- reversing the vertex order, rotating
the bulge assignment onto the correct segment, and negating each bulge.
A hole handed counter-clockwise builds the same solid as one handed
clockwise, instead of an inside-out passage.

Signed area for a bulge ring is the polygon shoelace plus each arc's
circular-segment correction, so the sign is the winding for curved rings
too -- a polygon-only shoelace would misjudge a ring whose arcs bulge
outward.

### Taper needs a corner list, not new geometry

The claim in ADR 0057 was wrong, and measurement showed it. The corner
rounding built there places the arc centre along the angle bisector at
`r / cos(turn/2)`, which is tangent to BOTH adjacent edges at any angle.
Probing it at 5, 8 and 14 degrees of taper gave tangency distances of
exactly 0.021000000 against both the inclined flange face and the
vertical web, for a stated radius of 0.021.

So a taper changes WHERE the corners sit, not how they are rounded.
Each tapered face pivots about the MID-POINT of its run, which keeps the
declared thickness the MEAN thickness -- the value section tables state.
Pivoting about the tip or the web face instead would silently change the
declared thickness and with it the area.

Implemented for I, AsymmetricI, T, U and L. A T's web and flange tapers
are independent and are carried in a `TaperT` struct rather than as two
adjacent `Scalar` arguments, because same-typed neighbours in an argument
list are exactly the shape of a silent swap.

## Consequences

- `Profile::Contour` with holes extrudes, with arcs, as does any arc
  section with through-passages: one fix, two features as predicted.
- Declared tapers build instead of refusing. `None` and `Some(0.0)` both
  mean a parallel face and give identical geometry.
- A slope at or beyond 0.9 of a quarter turn is refused: the inner face
  would approach parallel with the web and leave no flange.
- Two stale refusal tests were flipped to capability assertions rather
  than deleted, so the record shows what changed.

## Evidence

- Plate with a round hole: net section area 16 - pi to within 1e-9,
  four exact cylindrical walls at the stated radius.
- Both hole windings produce identical topology and wall counts.
- Tapered flange preserves the parallel-case area to 1e-12, confirming
  the mid-point pivot keeps the mean thickness.
- Tapered root fillet sits one radius from each ADJACENT face to 1e-9.
- Mutation testing: 7 of 7 mutants killed -- holes dropped, hole ring
  left counter-clockwise, outer winding not enforced, every ring marked
  outer, taper ignored, taper pivoted about the tip, steep-slope check
  removed.

## Notes

A tangency assertion must be made against the segments the arc actually
joins. The first version compared each fillet centre to the nearest line
anywhere in the outline and failed at 0.0136 on geometry that was exactly
correct, because an unrelated face sat closer than the radius.
