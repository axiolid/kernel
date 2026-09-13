# Elliptical cylinder surface

Status: accepted
Date: 2026-09-12

## Context

`Profile::Ellipse` refused exact extrusion because the wall it sweeps has no
surface type. The kernel had `Cylinder { frame, radius }`, which is circular
by construction, so an ellipse could only be approximated by a spline or a
polygon -- both of which discard the exact identity of the surface.

## Decision

Add `Surface::EllipticalCylinder { frame, semi_axis_x, semi_axis_y }` as a
distinct variant rather than generalising `Cylinder`, and extrude
`Profile::Ellipse` onto it exactly.

### Why not generalise `Cylinder`

`Cylinder`'s outward normal is its radial direction, and every consumer is
entitled to rely on that. For an ellipse the two coincide only at the four
axis points; elsewhere they diverge, by up to **53.1301 deg** for a 3:1
ellipse (measured, `ellcyl.py`). Widening `Cylinder` would silently break
that guarantee for existing code. A separate type forces every consumer to
answer the question.

The normal is therefore computed as `S_u x S_v`, not copied from the circular
case. Verified perpendicular to both partials to `2.2e-16`.

## Consequences

- `Profile::Ellipse` extrudes to an exact solid with a genuine
  `EllipticalCylinder` wall; the wall satisfies `(x/a)^2 + (y/b)^2 = 1` to
  `<1e-15` over 320 sampled points.
- Cap pcurves are `Ellipse2` and cap curves `Ellipse3`, sharing the surface's
  `(semi_axis_x, semi_axis_y)` parameterisation, so the cap boundary and the
  wall agree by construction. Confirmed by `geometric_audit`.
- `surface_periods` reports a u-period of `TAU`. This is load-bearing: the
  variant would otherwise fall into a catch-all arm reporting `(None, None)`,
  leaving the seam open instead of wrapping.
- `exact_properties` stays planar-only and now refuses the new wall **by
  name** (`NonPlanarFace("elliptical-cylindrical")`), exactly as it already
  refuses circular cylinders.
- Equal semi-axes are accepted and produce a geometrically circular
  elliptical cylinder. It is NOT rewritten to `Surface::Cylinder`: silently
  changing the surface type would make the output type depend on the input
  values, which callers matching on the variant cannot predict.

## Findings

### A new variant on a `#[non_exhaustive]` enum fails silently

`Surface` is `#[non_exhaustive]`, so consumers already carry wildcard arms and
adding a variant **compiles cleanly**. It does not fail loudly; it falls into
whatever the wildcard does. Two such sites were found by audit before writing
the surface:

- `surface_periods` (`compile/src/brep.rs`) returned `(None, None)`, which
  would have left the wall seam open.
- `family` (`measure/src/exact.rs`) returned the generic `"non-planar"`,
  which would have made the refusal message useless.

Both were handled explicitly. The lesson generalises: on a `#[non_exhaustive]`
enum, adding a variant demands an audit of every wildcard arm, because the
compiler will not do it for you.

### The normal shortcut is invisible at the obvious sample points

The disagreement between the true normal and the radial direction is exactly
zero at `u = 0, pi/2, pi, 3pi/2` -- the four points a test is most likely to
sample. A test checking only those would pass with a completely wrong normal.
The test therefore samples off-axis and asserts the peak disagreement matches
the independently computed 53.1301 deg.

## Mutation evidence

4/4 killed:

- normal copied from the circular (radial) case
- semi-axes swapped in the u-partial
- wall built with equal semi-axes (a circle)
- elliptical cylinder loses its u-period

## Alternatives considered

- **Generalise `Cylinder` to carry two semi-axes.** Rejected: silently breaks
  the radial-normal guarantee existing consumers rely on.
- **Represent the wall as a rational B-spline.** An exact ellipse needs
  rational quadratics, so this is representable -- but it discards the
  surface's identity, blocks any future analytic intersection, and costs more
  to evaluate for a shape the kernel can name exactly.
