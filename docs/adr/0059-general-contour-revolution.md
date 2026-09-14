# General contour revolution

Status: accepted
Date: 2026-09-12

## Context

`revolve_profile_exact` refused `Section`, `Contour`, `CenterLine`, `Derived`
and `Composite`. Those refusals predated ADR 0053-0058: every one of those
variants now lowers to a contour, so the refusals described a missing module
rather than geometry the kernel cannot represent. Left in place they would
read as deliberate limits.

`revolve_rectangle` builds one shape from four known corners and cannot
generalise: it hard-codes two cylindrical walls and two annular caps.

## Decision

Add `revolve_contour`, which revolves an arbitrary closed section a full turn
about the profile's local y axis. Each profile SEGMENT sweeps a surface
determined by its own geometry, verified numerically before implementation:

| segment | swept surface | check |
|---|---|---|
| parallel to the axis | cylinder | radius constant, spread 0 |
| perpendicular | planar annulus | height constant, spread 0 |
| oblique | cone | radius affine in height, residual 2e-13 |
| arc | torus | implicit residual 1.05e-15 |

Cylinder and cone share one wall builder: both parameterise as
`(angle, height)`, so only the surface differs. The torus needs its own
because its second parameter wraps around the tube.

`Section`, `Contour`, `CenterLine` and `Derived` route through the shared
`profile_to_contour`, so extrusion and revolution cannot drift apart on what
a profile means.

## Consequences

- Volumes are checked against Pappus (`V = 2*pi*R_centroid*A`), computed from
  the section independently of the kernel. A rectangular contour is also
  checked against `revolve_rectangle`: two independent paths must agree.
- `Composite` still refuses. `lower_composite` returns point rings, not
  contours -- it lowers past the contour stage -- so routing it through would
  mean re-deriving arcs from points. Refused by name rather than approximated.
- `Circle` and `Ellipse` still refuse: revolving a closed circle about an
  external axis is a torus the existing path does not build, and an ellipse
  sweeps a surface the kernel has no type for.
- Partial turns, profiles crossing the axis, and sections with holes remain
  refused. The first is a different topology; the last would need nested cap
  bounds on a revolved annulus.

## Findings

Two bugs that only the geometric audit caught, both invisible to a volume
check because both produce a closed solid of the right size:

1. **Reversed cap loops need reversed pcurve intervals.** An inner hole loop
   traversed `Reversed` with a forward interval put the pcurve start
   diametrically opposite its 3D edge. Reported error was exactly `2*r`.

2. **A torus seam is an arc, not a ruling.** Cylinders and cones have a
   straight v-direction, so their seam is a line. On a torus the seam runs
   around the tube: a straight chord sags by the sagitta
   `r*(1 - cos(sweep/2))`, which the audit reported to the digit
   (0.14644660940672635 against 0.14644660940672621).

A third finding concerns testing, not geometry: the arc end angle must come
from the arc's own sweep, never from `atan2` of the endpoint. The two agree
for every arc that avoids the tube's branch cut at `pi`, so the fillet
fixtures could not see the difference -- a mutant survived until a fixture
whose arc crosses the axis-facing side of the tube was added, where the span
differs by a full turn (+90 vs -270 degrees) at the same endpoint.

## Stale-claim sweep

Landing this exposed documentation that no longer matched behaviour:

- `construct/src/lib.rs` claimed exact generation was limited to
  "sharp rectangular and axial circular extrusions".
- `construct/AGENTS.md` listed revolution among families that "must refuse",
  and described `revolve` as discrete-only.
- ADR 0053's hole refusal and ADR 0057's taper refusal are marked superseded
  in place rather than rewritten: an ADR records what was decided when.

`PLAN.md`'s two open entries are Boolean issues, untouched by this work and
correctly still open.
