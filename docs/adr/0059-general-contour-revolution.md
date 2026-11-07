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

## Amendment 2026-11: a profile touching the axis (#379)

A full turn may touch the axis. The rules below replace the refusal of a
profile that touches it, in `revolve_arc_ring` and `revolve_profile_exact`.

**Contact is exact.** "On the axis" is the exact `side_of_axis` sign, never
a tolerance. A vertex on the axis is admitted when one of its two edges lies
on the axis (an axis edge, which sweeps no face) and the other is another
axis edge, a disc cap (a straight edge of exactly equal heights, which sweeps
a plane) or any other straight edge (which sweeps a cone up to its apex). A
vertex with no axis edge (a pinch point, a non-manifold vertex) is refused by
name; the mesh path refuses it too (#363). A straight edge ending on the axis
within the linear tolerance of perpendicular to it but not exactly is refused
by name: a cone of rise 1e-12 builds but fails the geometric audit by 6e-5.

**An apex has no degenerate edge.** The B-rep has no degenerate edges. An
apex is a vertex, and the wall's loop runs
`[rim circle] -> seam forward -> [rim circle] -> seam reversed`: in the
parameters (u, v) it leaves a gap along the pole line, where every u names
the same point. The partial-turn apex and brep-boolean's whole-sphere fixture
already used this convention, and `builder.finish`, `audit_brep`,
`geometric_audit` and `exact_properties` accept it. Consumers close the gap
instead of the builder modelling it: brep-boolean's seam planner reads a step
over a pole as no turn round the surface (it used to snap the jump to whole
periods and add a second seam, so every such face was refused with
`TangentSplit` or `Assembly`), and the deviation oracle treats a loop that
leaves a pole at the height it reached it as closed and ends a cone's box at
its apex. The B-rep tessellator has no pole row of samples yet and refuses
such a face by name (#463).

**The apex frame.** A cone with an end on the axis is framed at its apex with
radius zero, as the partial-turn path does, so the apex is at v = 0 exactly.
A cone between two circles keeps its frame at its start, bit for bit.

**Nesting.** A piece's range along the axis, which decides which shell is the
void of which, takes the heights of both stations of each of its faces, not
only the stations off the axis: a cavity inside a spinning top lies inside
the range of its two apexes, which no station off the axis spans.

**Poles of a sphere, and the quadrant snap.** An arc with exactly one end on
the axis is a sphere's meridian when the centre of its bulge circle lies
within `CENTRE_ROUNDING = 64 * EPSILON` times the scale of its coordinates
of the axis; otherwise it is a spindle or horn torus and is refused by name, as
are an arc with both ends on the axis and an arc whose sweep disagrees with
its latitude span by a half turn or more (it runs across the axis). The sphere
is placed exactly on the axis at the height equidistant from both arc ends,
`c = h_pole + (r^2 + d^2) / (2 d)`, with its frame like the cone's and `v` the
latitude, `+-pi/2` at the poles. The gate is a rounding gate, not a
tolerance: no bulge encodes a quarter arc exactly (`tan(pi / 8)` is
irrational), so the centre computed from a true sphere's arc is about
`1e-16 R` off the axis; measured at most a sixteenth of the gate in the worst
case tried (R 0.37, axis at 12.5). The only alternative to the gate is a
refusal, never a different solid.

`lower_arc` evaluates a quarter-turn parameter exactly (within
`QUADRANT_ROUNDING = 2 * EPSILON * max(1, |t|)`, `sin` and `cos` of 0 and
`+-1`): at a quadrant `|cos t|` is `0.18 * t * EPSILON` for every form a
producer writes, and it decides whether the arc's end is on an axis through
its circle's centre. This moves a point by less than the error of the point
itself, the same doctrine as the joint rule of #250. Neither gate is
reported: each chooses between a sphere and a refusal, and moves points below
their own evaluation error.

Rejected: carrying the source circles through `ArcRing` (a cross-crate type
change; composites lose their sources anyway; and it still could not make an
f64 pole exact), and a rounding band on the on-axis test itself (it breaks
the rule that on-axis tests are exact). Spindle and horn tori and pinch
points stay refused: supporting them needs a new torus patch domain or
non-manifold vertices.

**The axis frame (#450).** The exact revolution builds about the axis as
given. An axis exactly along +y or -y builds in place, as before. Any other
axis lying exactly in the profile plane (`z` exactly zero in its direction and
origin) gets the right-handed frame `(n, d, z)` with `d` its unit direction
and `n = (d.y, -d.x)`, origin on the axis: the section is mapped rigidly into
it (`p' = ((p - o) . n, (p - o) . d)`, the bulges kept), revolved about +y
there with the angle as given, and the solid placed back with
`ExactBRep::transformed`. The mapping keeps the exact decisions of the
revolution: a vertex exactly on the axis (the exact sign against the axis as
given) maps to `x' = 0`, an off-axis vertex keeps the sign of its side or is
refused as closer to the axis than the frame's rounding, an edge exactly
perpendicular to the axis with one end on it copies its off-axis end's height
(an exact `along_axis` sign), and a perpendicular-looking edge that is not
exactly perpendicular is refused. Any out-of-plane component is refused by
name: the edges would sweep hyperboloids, so there is no tolerance band to
offer.

Follow-up: exact booleans of a wall minus a revolution about a slightly
tilted axis are not sewn or classified yet (#470).

Rejected: refusing everything not exactly +-y (it throws away correct
in-plane axes, such as the profile's own x axis, which the mesh path
supports), and snapping the axis within the tolerance and reporting a
reading (a tilt of 1.4 mrad over a 10 m extent moves geometry by 14 mm, which
is not within tolerance, and construct has no reading channel). This agrees
with the rule that an axis above rounding keeps its exact coordinates;
re-expressing an axis within rounding, with a reading, waits on a reading
kind for it.
