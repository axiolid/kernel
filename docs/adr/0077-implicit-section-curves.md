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

### B-spline surfaces (stage 3 of ADR 0075)

- **`Field2` is an enum.** `Series(SeriesField2)` holds powers and
  harmonics, on analytic carriers. `Patches(PatchField2)` holds a
  tensor-product Bernstein polynomial per cell of a grid, on B-spline
  carriers.
- **Building the patch field.**
  1. The spline is split into its rational Bezier patches by knot insertion.
     On each, the homogeneous point `(X, W)` is a Bernstein polynomial.
  2. The analytic surface's implicit equation is homogenised to degree `d`
     in `(X, W)`: `W^d Q(X / W)`, a Bernstein polynomial of degree
     `(d p, d q)`, since Bernstein products stay Bernstein. Its zero set is
     the section wherever `W > 0`, which a B-spline's weights guarantee.
  3. The coefficients bound the field over any box by the convex hull
     property, after de Casteljau restriction to the box. That gives
     tighter bounds than the series' interval arithmetic, and the same
     certified trace runs unchanged.
  4. The edge cells' polynomials continue past the grid, so a window may
     reach slightly beyond the spline's domain. Traced curves are clipped
     to the domain afterwards.
- **`Carrier::Spline`.** `BSplineSurface` moved to `axiolid-curve`, and
  `axiolid-surface` re-exports it unchanged. A spline carrier evaluates its
  point and partials to second order in closed form (rational tensor
  product). Its inverse is iterative: `evaluate::surface::locate`.
- **Pcurves on the analytic face: `Curve2::Lifted`.** A B-spline has no
  implicit equation to read in an analytic surface's parameters. So the
  section's pcurve on the analytic face is the section's own space curve
  read back through the analytic surface's closed-form inverse, sharing the
  edge's parameter. A guide of unwrapped parameters picks the whole turn at
  each point. Its turning points come from a dense scan of the derivative's
  signs; that step is not certified.

### Two B-spline surfaces

Neither surface has an equation the other can be read in, so the section
is carried on both: `Curve3::PairSection` (`PairSection3`).

- **Representation.** A chain of nodes, each a point on both surfaces with
  its parameters on each, corrected to the last bits. Between two nodes the
  curve is defined, not interpolated. At local parameter `s` it is where
  both surfaces meet on the plane across the chord at `P0 + s (P1 - P0)`:
  four equations (`S1(a) = S2(b)` and the plane) in the four parameters,
  solved by Newton from the nodes' parameters. Its rates come from the same
  system. The parameter runs one unit per chord.
- **Finding every component** (`spline_pair_intersection`).
  1. Pairs of rational Bezier sub-patches are split until their control
     hulls' boxes are apart (no section there), or their normal cones are
     apart. The cones hold the Bernstein coefficient vectors of the
     normal's own polynomial, so they are certain.
  2. A pair with apart normal cones has no two parallel normals, so it
     holds no closed loop of the section: every piece of the section in it
     crosses an edge of one of the sub-patches (Sederberg and Meyers).
  3. Each sub-patch edge is intersected with the other sub-patch (hull
     pruning, then Newton on three unknowns). These crossings seed every
     component.
  4. From each seed not already on a traced curve, the curve is followed
     both ways with steps held to a few degrees of turning, each node
     corrected onto both surfaces. It ends where it closes on itself or
     where it leaves either face's window. The last node is solved exactly
     on that window's edge.
  5. Pairs that never separate within the depth limit are refused as
     `NotRegularCurve`: the surfaces touch there, or come closer than the
     search resolves.
- **What is certified.** The pruning and the loop-free test are certain.
  The edge crossings and the following are Newton with step control, not
  certified. Every node lies on both surfaces to rounding, whichever way
  the steps went.
- **Pcurves.** On either spline face the pcurve is `Curve2::Lifted` on a
  spline carrier. When the lifted curve is a `PairSection` on one of its
  own surfaces, evaluation reads that surface's parameters straight from
  the solve.
- **Curve against surface.** A B-spline curve against a B-spline surface
  uses the same hull pruning over the Bezier pieces of both, then Newton.
  The boolean uses it to cut a pair section where it crosses a spline
  face's boundary edge.

### The section families against each other and against distance (B8, B7)

- **Plane curves** (`section_curve_curve_intersection2`). Every family with
  a defining field (line, circle, ellipse, `Sinusoid2`, `QuadraticGraph2`,
  `AngleGraph2`, `ImplicitCurve2`) is traced over its span. The other
  curve's field has its roots isolated along the cells, and a root is kept
  when it lies on both pieces, so a graph's other branch is filtered out.
- **Space curves** (`section_curve_curve_intersection3`). The first curve
  is intersected, with certification, with each surface the second lies on
  (its carrier, a conic's plane, a line's two planes). A point is kept
  where it lies on the second curve within the tolerance. In space the
  question is only meaningful up to a tolerance.
- **Extrema** (`extrema::minimum_distance`). Branch and bound over points,
  curve spans and surface patches. The lower bounds are boxes certain to
  hold each piece's image:
  - affine images for lines and planes;
  - exact harmonic ranges for conics and for cylinders, cones, spheres and
    tori;
  - control hulls for B-splines;
  - certified cells mapped through the carrier for traced sections.

  The upper bounds are evaluated witness points. The result brackets the
  minimum within the accuracy asked for.

### Singular points

- An isolated point where the surfaces touch without crossing is a
  non-degenerate extremum of the field on its zero set (definite Hessian).
  The trace finds it by Newton on the gradient and drops the boxes around
  it that rounding cannot decide: no curve runs through it. A trace that
  finds only such points answers `NotRegularCurve`, as the closed forms do.
- Two branches crossing at a saddle (indefinite Hessian) stay refused by
  name. Within about `sqrt(rounding / curvature)` of the crossing, which
  is micrometres at metre scale, the field is below its own rounding, so
  where the branches run there is not decidable in doubles. Ending them at
  that distance would leave gaps larger than a micrometre tolerance.
- Windows a whole turn wide start an irrational fraction of a radian past
  `-pi` (and past a face's own seam angle), so a symmetric section's
  special points never sit on the window's edge.

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
  form does, and so does every section involving B-splines. The only
  exception is the singular case (touching surfaces), which is refused by
  name.
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
- Two B-splines: seeding and following are not certified (see above). A
  certified version would isolate the edge crossings with interval Newton
  and bound each step's corridor.

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
- `crates/representations/analytic/curve/src/pair_section.rs`:
  `PairSection3`, the section of two B-splines.
- `crates/algorithms/parametric/nurbs/src/pair_trace.rs`:
  `spline_pair_intersection` and B-spline curve/surface crossings.
- `crates/algorithms/parametric/evaluate/src/curve.rs`: evaluation and
  inversion of the new families. Inversion now also covers the ADR 0076
  graphs.
- `crates/algorithms/parametric/nurbs/tests/implicit_section.rs`: points on
  both surfaces, continuity, and completeness against a dense scan of tube
  circles.
- `crates/algorithms/parametric/nurbs/tests/pair_section.rs`: two sheets
  crossing (completeness against a dense scan), a closed loop inside one
  patch pair, and a spline curve crossing a spline surface.
