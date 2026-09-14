# Exact centre-line offsetting

Status: accepted
Date: 2026-09-12

## Context

`Profile::CenterLine` was the last refusing variant reachable by adapting
existing machinery. It denotes the set of points within `half_width` of an
open path, so resolving it is an offsetting problem, not a lowering one.

An offsetter already existed: `center_line.rs` offsets the FLATTENED path
with miter joins and butt caps, feeding the tessellating pipeline. That is
correct there and wrong for the exact extruder, because flattening bakes a
chord budget into a result that is supposed to be exact. The output would
close, audit clean, and not be the requested shape.

## Decision

Add `center_line_exact.rs`, which offsets the path AS CURVES and produces a
`ContourProfile`. The contour path (ADR 0053) then extrudes it exactly,
including genuine cylindrical walls for curved paths.

The existing flattening offsetter is untouched and still serves the
tessellating pipeline. The two coexist deliberately: one is parameterised by
a chord budget, the other refuses anything it cannot do exactly.

## What can be offset exactly

Measured, not assumed (`offset.py`):

| segment | offset | exact |
|---------|--------|-------|
| line | parallel line | yes |
| circular arc (ccw) | concentric arc, `r - d` | yes |
| circular arc (cw) | concentric arc, `r + d` | yes |
| ellipse | NOT an ellipse | no -- refused |

The circle result was verified by offsetting 200 sampled points and measuring
the distance to the original centre: constant to `6.7e-16`. The ellipse case
was verified by fitting the best ellipse to the offset curve: residual
`0.074` for a 3:1 ellipse at distance `0.3`, so it is simply a different kind
of curve. Splines are refused for the same reason.

An offset arc keeps its centre, frame and parameter domain, changing only the
radius. That is what makes the two offset endpoints correspond at equal
parameters, so the caps join the sides without a search.

## Consequences

- Butt caps, matching the flattening offsetter. The source states a width and
  an extent, not an end treatment; a round or square cap would add material
  the author never declared.
- A CLOSED path is refused: it denotes an annulus, which needs a hole rather
  than a single outer ring. Welding the ends would silently fill it.
- A disconnected path is refused rather than bridged.
- A tangent reversal is refused: the two offsets would cross and the enclosed
  area would depend on where they self-intersect.
- A half-width that drives an inner offset radius to zero or below is refused
  by name; that is not a curve.
- `Profile::Section` remains the only refusing variant.

## Findings

### Traversal order decides orientation, and only a SIGNED check sees it

The first implementation walked the LEFT offset forward and the RIGHT offset
back. That traces the boundary CLOCKWISE, so the solid came out inside-out.

The volume was `-4` where the closed form gives `4`: correct in magnitude,
wrong in sign. A test comparing magnitudes -- or an `abs()` anywhere in the
chain -- would have passed. The shoelace formula over the four corners
confirmed the winding independently before the fix.

Correct order is RIGHT side forward, then LEFT side back.

### Frame handedness decides which side grows

`Circle2` stores `x` and `y` independently, so an arc frame may be
left-handed, and then the parameter runs clockwise in world orientation. The
offset side follows the WORLD turn, so dropping `handedness.signum()` swaps
which side grows while the ring still closes.

That mutant initially SURVIVED: every arc fixture used a right-handed frame.
A left-handed fixture was added, asserting the two walls still straddle the
path radius at 1.8 and 2.2, and it now kills the mutant.

## Mutation evidence

6/6 killed:

- offset side ignores turn direction
- handedness dropped from the world sweep
- collapsed arc no longer refused
- closed path silently welded
- disconnected path accepted
- elliptical segment silently accepted
