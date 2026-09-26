# 0077 — Implicit section curves: certified tracing of analytic sections

- **Status:** Accepted
- **Date:** 2026-09-26
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

Part of axiolid/kernel#119 (ledger rows B8–B10); gates stage 2 of ADR 0075
(#167).

## Context

After ADR 0076, `exact_surface_intersection` built every analytic pair
with a ruled carrier (a cylinder, elliptical cylinder or cone) and a torus
against a plane or sphere. Left over were a torus against a cylinder,
elliptical cylinder, cone or another torus off its axis. These sections are
quartic or octic in space, with no closed form of either kind used so far:
not a root of a quadratic along a ruling, nor an angle solving
`A cos u + B sin u = C`.

The general boolean (ADR 0075) needs more than that. Every section edge
on a sphere, cone or torus face needs an exact pcurve on that face. Most of
these cannot be expressed in the existing pcurve families, for example an
oblique plane's ellipse on a cone or a cylinder's cut across a sphere.

Both needs have the same shape. Read in one surface's parameters `(u, v)`,
the section is where the other surface's implicit equation, composed with
the first surface's point, vanishes. For planes, quadrics and tori that
composition is always a finite sum of products of powers of a linear
parameter and harmonics of an angle, with coefficients computed from the
two surfaces' own numbers.

## Decision

Carry such sections as **stretches of a field's zero set, in monotone
cells**. Find all of them in a parameter window by **certified subdivision**,
never by marching.

- **`Field2`** (`axiolid-curve`) is `sum c[i][j] B_i(u) B_j(v)` with a
  `Basis` per parameter: `Power` (`x^k`) or `Fourier` (`1, cos kx, sin kx`).
  It evaluates its value, gradient and Hessian in closed form.
- **`Curve2::Implicit(ImplicitCurve2)`** is a field and a chain of
  `ImplicitCell`s. In each cell, one parameter (`axis`) runs linearly from
  `from` to `to`, and the field is *strictly monotone* in the other parameter
  over the bracket `[low, high]`. So at every free value the curve is the
  field's *unique* zero in the bracket. The parameter `t` runs over
  `[0, cells.len()]`, one unit per cell.
  - A point is that root, found to full precision by a Newton iteration
    safeguarded by bisection. It cannot leave the bracket or reach a
    different branch.
  - Derivatives come from the implicit function theorem.
  - Inversion locates the cell whose box holds the point.
- **`Curve3::ImplicitSection(ImplicitSection3)`** is the same curve on its
  **`Carrier`**: a copy of a plane, ruled surface (reusing `RuledCarrier`),
  sphere or torus (reusing `TorusCarrier`) in its `axiolid_surface`
  parameterisation. The point at `t` is the carrier at the pcurve's point.
- **`section_field`** (`axiolid-nurbs`, `field.rs`) builds the field:
  1. The carrier's point becomes three fields.
  2. The other surface's implicit equation is written in its own local
     coordinates. That is `l_z` for a plane, or `x² + y² − r²`,
     `b²x² + a²y² − a²b²`, `x² + y² − (r + z tan α)²`, `|p|² − r²`, or
     `(|p|² + R² − r²)² − 4R²(x² + y²)`.
  3. Products stay in the basis by the product-to-sum rules.
- **`trace`** (`implicit_trace.rs`) finds every component in a window:
  1. Start from an 8 × 8 grid, with its lines offset from round fractions.
  2. For each box, bound the field by interval arithmetic over the terms'
     exact ranges, tightened by the mean-value form and widened by a rounding
     margin. If the bound excludes zero, the box is dropped.
  3. If one partial's bound excludes zero, the box is **regular**: the curve
     in it is a graph over the other parameter.
  4. Otherwise the box is split into four.
  5. In a regular box the graph exists exactly where the field has opposite
     signs on the two sides across it. Those sides' roots are isolated in
     1D with a certificate: monotone and changing sign, or excluded. The
     stretches between them are the cells.
  6. Cells of neighbouring boxes meet where the curve crosses their common
     side, so chaining by end points gives the components: closed loops,
     loops that wind around a periodic parameter (unwrapped by whole turns),
     or chains leaving the window.
  7. Completeness needs no step size. The regular boxes cover every point of
     the zero set, and each regular box's cells are exact.
- **Singular points are refused by name.** A box where the field and both
  partials may all vanish, at the smallest size allowed, may hold a point
  where the surfaces touch and branches cross. The trace answers
  `NotRegularCurve`, as every other tangency in #119 does.
- **Where it is used:**
  - `exact_surface_intersection` falls back to it after the closed forms and
    ADR 0076. It is used for analytic pairs with a compact surface to carry
    the trace, a torus first (`Derivation::ImplicitTrace`, one span per
    curve).
  - `implicit_surface_intersection` is public, with an explicit window for
    any analytic pair.
  - The general boolean uses it for exact pcurves of section edges on
    sphere, cone and torus faces.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Marching (predictor–corrector) as OCCT's `IntWalk` does | Step size is a guess: loops smaller than a step are missed, branches can be jumped near tangencies, and nothing certifies that every component was found. |
| Fit B-splines to the traced points | Turns an exact definition into an approximation with a tolerance, in the module that promises none; pcurves and edges would disagree by the fit error. |
| Resultants and root isolation of the full algebraic curve | Exact topology, but the degrees are large: a torus against a torus is degree 8 in each of two variables before the resultant. Isolating the discriminant's roots with big integers per face pair is far slower than subdivision, and the curve still needs evaluating afterwards. |
| A guide polyline with normal projection | Needs a step size to build and a certificate that the normal lines meet the curve once. Monotone cells give both by construction. |

## Consequences

**Positive**

- Every analytic pair's section now exists exactly, whether or not a closed
  form does. The only exception is the singular case (touching surfaces),
  which is refused by name.
- One pcurve family serves every section on every analytic face, so the
  boolean needs no per-pair pcurve derivations beyond the cheap closed forms
  it already has.
- The pcurve and the edge are the same definition, one read in parameters
  and one in space. They agree by construction, not up to a fit.

**Negative / costs**

- Evaluation solves a 1D root per point. It is bracketed and converges in a
  few Newton steps, but it is not closed form.
- A curve has one cell per regular box it crosses, typically tens to a few
  hundred, and the parameter is uniform per cell, not in arc length.
- Coefficients are computed in `f64` from the operands' doubles. As for
  ADR 0076's curve coordinates, the field is the rounded composition. The
  certificates hold for that field, and the rounding margin covers
  evaluating it.

**Follow-ups / risks to watch**

- Singular points (branches crossing where surfaces touch): isolate them as
  vertices of the section graph instead of refusing.
- Very thin sections or nearly tangent pairs can exhaust the trace's box
  budget (400 000 boxes). That is also refused as `NotRegularCurve`.

## Relation to existing code

- `crates/representations/analytic/curve/src/implicit.rs`: `Field2`,
  `ImplicitCurve2`, `ImplicitSection3`, `Carrier`.
- `crates/algorithms/parametric/nurbs/src/field.rs`: field construction,
  interval bounds, partials.
- `crates/algorithms/parametric/nurbs/src/implicit_trace.rs`: the certified
  trace.
- `crates/algorithms/parametric/nurbs/src/implicit_section.rs`:
  `implicit_surface_intersection`, the fallback of
  `exact_surface_intersection`.
- `crates/algorithms/parametric/evaluate/src/curve.rs`: evaluation and
  inversion of the new families. Inversion now also covers the ADR 0076
  graphs.
- `crates/algorithms/parametric/nurbs/tests/implicit_section.rs`: points on
  both surfaces, continuity, and completeness against a dense scan of tube
  circles.
