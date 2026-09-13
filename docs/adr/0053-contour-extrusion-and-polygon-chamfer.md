# Contour extrusion and general-polygon chamfer

Status: accepted
Date: 2026-09-12

## Context

`extrude_profile_exact` accepted 2 of the 8 `Profile` variants: `Rectangle`
and `Circle`. `Contour` -- the variant carrying an arbitrary exact boundary --
was refused outright, which is the one a boolean result or an imported model
most naturally produces.

Separately, `chamfer_extruded_profile` was rectangle-only, the same
restriction the fillet carried before ADR 0051.

## Decision

Lower `Profile::Contour` onto the existing arc-ring extruder rather than
writing a new builder. `Line2` segments become planar walls and `Circle2`
segments become cylindrical ones, which is exactly what the arc path already
produces. Contours of only straight segments route to the polygon path
instead, which carries face naming and hole support the arc path lacks.

Curve kinds that cannot be carried exactly -- ellipse, spline, polyline,
intrinsic -- are REFUSED. Sampling them into short chords would produce a
solid that closes, validates, audits clean and is not the requested shape.

Add `chamfer_polygon_corners` mirroring `fillet_polygon_corners`: same reflex
and pairwise-edge constraints, straight cut instead of an arc.

## Consequences

- A contour with holes is refused: the arc cap-loop builder emits a single
  bound, so a hole would vanish from the caps while still appearing in the
  walls. That is a solid that closes and is wrong, which is worse than a
  refusal.
- A circular segment sweeping half a turn or more is refused. `bulge` is
  `tan(sweep/4)`, and beyond a half turn the chord no longer determines the
  arc. Callers split such an arc into two segments, as the stadium test does.
- Contour extrusion results are checked by the ADR 0052 geometric audit in
  tests, not only by closure.

## Findings

### The adapter exposed a sign bug in shipped arc code

Deriving the `Curve2::Circle` to bulge conversion meant checking it against
the shipped `arc_geometry`. Positive sweeps agreed to `1e-16`; negative
sweeps were wrong by order 1.

The cause, confirmed by an independent property test: `arc_geometry` placed
the centre with an UNSIGNED apothem,

    centre = mid + normal * (radius * cos(half))

and `cos` is even, so `bulge = +0.5` and `bulge = -0.5` produced the
IDENTICAL centre. Every arc bulged the same way regardless of its stated
direction. Verified concretely: both produced centre `(1.0, +0.75)` for the
same chord.

Every existing arc test used positive bulges only -- `ArcRing::circle` and the
hand-written square fixtures -- so nothing caught it. Fixed by carrying the
sign:

    centre = mid + normal * (radius * cos(half)) * sweep.signum()

Verified at `6.1e-16` across both sweep signs AND both frame handednesses. The
regression test was confirmed to FAIL against the old code before being kept.

### Frame handedness is part of the conversion, not a detail

`Circle2` stores `frame.x` and `frame.y` independently, so the frame can be
left-handed. When it is, the parameter runs clockwise in world orientation and
the world sweep is the negation of the parameter sweep. Dropping
`handedness.signum()` puts the arc on the wrong side of its chord while the
ring still closes -- invisible to any topological check.

### A mutation harness can lie

Three mutants initially "survived". Two were genuine test gaps (handedness,
half-turn). The third was a harness error: the chamfer pairwise check and the
fillet pairwise check share identical `if` text, so a
`replace(old, new, 1)` mutated the FILLET path and left the chamfer check
intact. Line-targeted mutation killed it immediately.

The lesson generalises: when a mutant survives, confirm the mutation actually
landed where intended before concluding the test is weak.

Mutation, all killed: handedness sign ignored; half-turn no longer refused;
contour closure not checked; chamfer pairwise check removed; chamfer reflex
check removed; arc centre sign dropped.
