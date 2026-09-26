# 0076 — Ruled quadric sections: exact quartic intersection curves

- **Status:** Accepted
- **Date:** 2026-09-26
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

Part of axiolid/kernel#119 (ledger row B10); gates stage 2 of ADR 0075.

## Context

`exact_surface_intersection` derived only the pairs whose intersection is a
line, circle or ellipse. Everything else between elementary surfaces was
refused by name: two cylinders of different radii crossing (a pipe tee), a
sphere and a cylinder off each other's axis, a plane cutting a cone off its
axis (an oblique ellipse, a parabola or a hyperbola), elliptical cylinders
against anything curved. These are space quartics or open conics; no
`Curve3` variant held them, and fitting a spline would be an approximation
in the one module that promises none.

These are the pairs a general boolean (ADR 0075) meets first after the
conic cases: pipes joining pipes, columns piercing domes, cones cut
obliquely.

## Decision

Represent the section of a **ruled carrier** by a quadric as a root branch
of a quadratic over the carrier's angle.

A cylinder, elliptical cylinder or cone is linear in `v` along each ruling,
`P(u, v) = A0(u) + v A1(u)`. Substituted into a quadric `Q`, it gives
`a(u) v^2 + b(u) v + c(u) = 0` with `a`, `b`, `c` trigonometric polynomials
of degree at most two. The curve is the graph of one root.

- **`Curve2::QuadraticGraph(QuadraticGraph2)`**: `a`, `b`, `c` as `Trig2`
  coefficients and a `Branch` (`Plus`, `Minus`). The parameter is the first
  coordinate, as for `Sinusoid2` (ADR 0071), which is its `a = 0` case.
  Evaluated as `v = 2c / (-b - sign sqrt(D))`, finite where `a = 0` and free
  of cancellation; derivatives by implicit differentiation.
- **`Curve3::RuledSection(RuledSection3)`**: the carrier (`RuledCarrier`:
  frame, two radii, slope) evaluated along the graph. The same point, never
  a separate approximation, so pcurve and edge agree by construction.
- **Spans are decided exactly.** Coefficients are built in dyadic
  arithmetic from the operands' doubles. Under `t = tan(u/2)` the
  discriminant times `(1 + t^2)^4` is an integer polynomial of degree 8;
  its roots are isolated with Sturm sequences (`axiolid-exact`), and the
  sign on each span is decided at a dyadic point. A root at `u = pi`
  (`t = infinity`) shows as a vanishing value there and is added as a
  break. Whether the surfaces meet in one loop, two loops, none, or only
  touch is therefore not a rounding question.
- **Result shape.** `ExactIntersectionCurve` gains `spans`, aligned with
  `branches`: `None` for curves on their whole natural domain, `Some` for a
  ruled piece. Over a bounded span, a `Plus` and a `Minus` piece join at
  both ends into one loop; a whole-turn span gives two closed curves.
  `Derivation::RuledQuadricSection` names the identity.
- **Scope.** Carriers: cylinder and elliptical cylinder against plane,
  sphere, cylinder, elliptical cylinder and cone; a cone against plane,
  sphere and cone. The closed forms still win where they apply; the ruled
  path runs only when they refuse.
- **Tori.** A torus is not ruled, but at a fixed tube angle `v` its
  points form a circle about the axis, which a plane or sphere meets where
  `A(v) cos u + B(v) sin u = C(v)` (degree one in `v`). The section is `u`
  as a function of `v`: `Curve2::AngleGraph(AngleGraph2)` (the parameter is
  the second coordinate) and `Curve3::TorusSection(TorusSection3)`. It
  exists where `E = A^2 + B^2 - C^2 >= 0`; the returned angle wraps at
  `u = pi`, so spans are also split at the roots of `B^2 C^2 - A^2 E`,
  keeping each pcurve piece continuous. Both are decided by the same exact
  root isolation. `Derivation::TorusAngleSection`.
- **Nappes.** A cone's implicit equation holds both nappes, and a cone
  carrier's parameterisation reaches the other nappe past its apex. Each
  cone taking part adds a condition `h0(u) + h1(u) v >= 0` (the carrier's
  `r + s v`, the other's `r' + s' Z' . (P - O')`). Along a root branch it
  has the sign of `a (P + sign Q sqrt D)` with `P = 2 a h0 - b h1`,
  `Q = h1`, decided exactly at each span's dyadic sample with
  `axiolid_exact::sign_root`; its sign changes are roots of
  `P^2 - Q^2 D`, which join the span breaks.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Fit a B-spline to traced points | An approximation with an error bound, in a module whose contract is exactness. |
| Levin's pencil parameterisation (a ruled quadric in the pencil of the two) | Also exact, but needs choosing a ruled member of the pencil and radicals of radicals; the carrier is already ruled when one operand is a cylinder or cone, which covers every pair in scope. |
| A generic "inverse image on the surface" curve type | Evaluates by inversion, so exactness rests on an iterative solver; the quadratic is closed form. |
| Parabola and hyperbola `Curve3` variants | Cover only plane/cone; the graph family covers those and every quartic in scope with one representation. |

## Consequences

**Positive**

- Pipe tees, off-axis sphere/cylinder junctions, oblique cone cuts and
  elliptical cylinder junctions have exact intersection curves, with the
  branch structure decided exactly.
- One representation for all of them, with the plane cut as its degenerate
  member.

**Negative / costs**

- Two new `Curve2`/`Curve3` variants (both enums are `#[non_exhaustive]`);
  consumers that do not handle them refuse by name, as for any unknown
  family.
- A piece's parameter `u` is singular at its branch ends (vertical
  tangent): positions are defined there, derivatives are not. Consumers
  stop at the span.
- Pieces going to infinity (a hyperbola, a parabola) are returned over open
  spans; a consumer must trim them.

**Follow-ups / risks to watch**

- A torus against a cylinder, cone or torus is quartic in `u` at each `v`
  and stays refused unless coaxial.
- B-spline pairs stay on the certified bounded tier.

## Relation to existing code

- `crates/representations/analytic/curve/src/quadric_section.rs`: the
  representation.
- `crates/algorithms/parametric/evaluate/src/curve.rs`: evaluation arms.
- `crates/algorithms/parametric/nurbs/src/ruled_section.rs`: the exact
  derivation; `exact_surface_intersection.rs`: dispatch and `spans`.
- `crates/algorithms/parametric/nurbs/tests/ruled_section.rs`,
  `crates/algorithms/parametric/evaluate/tests/quadratic_graph.rs`:
  oracles; `scripts/probe_ruled_section_mutants.py`: 16 faults, all caught.
