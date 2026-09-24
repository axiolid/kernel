# 0071 — Sinusoid pcurves for plane cuts across cylinders

- **Status:** Accepted
- **Date:** 2026-09-24
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

Part of axiolid/kernel#120.

## Context

An exact B-rep stores every edge twice: as a 3D curve, and as a 2D pcurve in
the parameters of each face that uses it. `ExactBRep` refuses a face whose
edge uses lack a pcurve, and `axiolid-brep-audit` checks that every pcurve,
lifted through its surface, lands on its 3D edge.

The next #120 case is a curved prism cut by a sloped plane: a round column
under a sloped roof. #119 already gives the 3D edge exactly: a plane cuts a
cylinder in an `Ellipse3`. The missing piece is the pcurve on the cylinder.

A cylinder is parameterised by angle `u` and height `v`. The plane
`z = h + g . (x, y)` meets it where `x` and `y` are `r cos u` and `r sin u`
(plus the centre), so the cut, read in `(u, v)`, is

```text
v(u) = mean + a cos(u) + b sin(u)
```

None of the existing `Curve2` variants holds that exactly. It is not a conic in
`(u, v)`. A B-spline can only approximate it, because `u` is an angle and the
wave is transcendental in it. An `Intrinsic2` is defined by curvature against
arc length, and the wave's arc length is an elliptic integral.

## Decision

Add `Curve2::Sinusoid(Sinusoid2 { mean, cosine, sine })`: the graph
`t -> (t, mean + cosine cos t + sine sin t)`.

- **The parameter is the first coordinate.** A pcurve use states the angle
  span it covers directly as its interval, with no reparameterisation.
  Inversion is exact: the parameter of a point is its `u`, then its height
  is checked.
- **Evaluation is closed form**, with first and second derivatives and a
  one-turn conventional domain, like a circle. Graph validation accepts it as
  a trim basis when its three coefficients are finite.
- **`Curve2` is `#[non_exhaustive]`**, so adding a variant is additive for
  downstream matches.

The first consumer is `clip_arc_prism_exact(prism, half_space, tolerance)`.
When the plane passes cleanly between the prism's caps over the whole
section, the result is the prism with one cap replaced by the cut. Each
cylindrical wall stays a `Cylinder`, its sloped rim an `Ellipse3` taken from
`exact_surface_intersection` with a `Sinusoid2` pcurve. Each planar wall gets
a sloped straight edge, and the new cap is a plane face in the cutting plane.

Whether the plane stays clear of the caps is decided on the true range of
the plane's height over the section. An affine function over an arc edge
takes its extremes at the endpoints or where the arc is tangent to the
plane's contour lines, so those candidates are checked exactly rather than
sampled.

## Alternatives considered

| Option | Why not |
| --- | --- |
| B-spline pcurve within a stated error bound (what OCCT does) | Stores an approximation. The kernel's exact tier would carry an edge whose two representations disagree by a tolerance, and the audit would have to accept that tolerance. |
| New ruled-surface type for the wall | Keeps both trims as lines, but the wall stops being a `Cylinder`. IFC export and other downstream users read "this face is a cylinder" to recover intent, and that would be lost. |

## Consequences

**Positive**

- A plane cut through a cylinder or elliptical cylinder is exact end to end:
  the 3D edge comes from #119, the trim is the wave, and the audit verifies
  both against each other.
- The same pcurve serves every later plane cut through a cylindrical wall:
  sloped roofs, oblique openings, mitred round members.

**Negative / costs**

- One more `Curve2` variant for every consumer to handle or refuse. Because
  the enum is non-exhaustive, the ones that don't know it already refuse it
  by name.
- Arc-length queries (`arc_length.rs`) refuse the wave: its length is an
  elliptic integral, and nothing needs it yet.

**Follow-ups / risks to watch**

- Refused by name for now: a plane that crosses a cap inside the section,
  which needs a face with both an original and a cut cap (built since
  ADR 0072); a plane parallel to the axis, which is a plan cut and not a
  cap cut; and cuts of cones, spheres and tori.
- Elliptical-cylinder walls use the same wave. The builder only produces
  circular cylinders today, because arc sections have circular arcs.

## Relation to existing code

- `crates/representations/analytic/curve/src/sinusoid.rs`: the value type.
- `crates/algorithms/parametric/evaluate/src/curve.rs`: evaluation,
  derivatives, domain, inversion.
- `crates/representations/modeling/graph/src/validation.rs`: trim-basis
  validation.
- `crates/algorithms/construction/construct/src/extrude_arc.rs`: `Level`
  and the sloped-cap builder.
- `crates/algorithms/construction/construct/src/boolean_exact.rs`:
  `clip_arc_prism_exact`.

## Verification

- 14 clip tests. Volumes are checked against a closed form (section area
  times the plane's mean height over the section centroid), including an
  off-axis column, a section mixing straight and curved walls, and a
  concave arc (a notch), whose cut edge runs backwards along its ellipse. The cut
  edges are sampled and checked to lie on both the cylinder and the plane.
  Every result passes the geometric audit, and every refusal is covered.
- 6 wave tests. A lifted wave lands on its plane, on a circular and an
  elliptical cylinder. Derivatives match finite differences. Inversion
  refuses off-curve points.
- One graph-validation test covering non-finite coefficients.
- An 11-fault mutation probe, `scripts/probe_sloped_cap_mutants.py`, catches
  all 11. Its first run caught 8 of 12: the in-arc range extremes, the
  cap naming and the crossing-cap refusal had no test that could fail, and
  those tests were added. The twelfth fault forced the ellipse direction
  check to always return +1. It is equivalent today, because both frames
  point up, so it was dropped from the probe. Forcing the direction to -1
  instead fails 9 of the 14 clip tests, so a reversed edge is caught.
