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
- **Singular points become vertices.** A box where the field and both
  partials may all vanish, at the smallest size allowed, holds a point
  where the surfaces touch. A touching point with no curve through it is
  dropped; branches crossing or touching there end at it (see *Singular
  points*). A trace that exhausts its budget, or meets a singular point
  whose branch ends do not match the field's sign changes about it, is
  `Undecided`; one that finds only touching points is `NotRegularCurve`.
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
  0. Each Bezier patch is first cut to the faces' windows, so the windows'
     edges are sub-patch edges too.
  1. Pairs of rational Bezier sub-patches are split until their control
     hulls' boxes are apart (no section there), or their normal cones are
     apart. The cones hold the Bernstein coefficient vectors of the
     normal's own polynomial, so they are certain.
  2. A pair with apart normal cones has no two parallel normals, so it
     holds no closed loop of the section: every piece of the section in it
     crosses an edge of one of the sub-patches (Sederberg and Meyers).
  3. Each sub-patch edge is intersected with the other sub-patch by
     certified root isolation: boxes of (edge parameter, other surface's
     parameters) are halved until Krawczyk's test proves each holds no
     crossing or exactly one. These crossings seed every component.
  4. From each seed not already on a traced curve, the curve is followed
     both ways with steps held to a few degrees of turning, each node
     corrected onto both surfaces. It ends where it closes on itself or
     where it leaves either face's window. The last node is solved exactly
     on that window's edge.
  5. Every chord is then proven: Krawczyk's test on `S1(a) = S2(b)` and the
     chord's plane, with the plane's level an interval over the whole
     chord, shows a box in both surfaces' parameters holding exactly one
     point of the section at every level. The chord's two nodes are that
     point at its first and last level, so they lie on one arc, and the
     arc is the curve `PairSection3` defines there. A chord that cannot be
     proven is halved at its solved middle; a curve that still cannot be
     is followed again with shorter steps, and after that refused.
  6. A seed is on a curve already traced exactly when it lies in one of
     its chords' boxes at a level of that chord.
  7. Pairs that never separate within the depth limit are refused as
     `NotRegularCurve`: the surfaces touch there, or come closer than the
     search resolves.
  8. An edge crossing that no box can prove, where the edge touches the
     other surface (the Jacobian is singular there, so no proof can exist),
     is found by damped Newton and kept as a seed; a regular crossing that
     cannot be proven is refused.
- **Enclosures** (`pair_certify`). Each surface's point and first partials
  over a parameter box are bounded by Bernstein coefficients: the point by
  its rational control points (the weights are positive), each partial by
  `(W X_u - W_u X) / W^2`, numerator and `W^2` bounded by their own
  coefficients. Past the domain the edge cells' polynomials continue, and
  the certificates evaluate points from those same polynomials. Every
  bound is widened for the rounding of its coefficients.
- **What is certified.** Every component is found, and every chord is
  proven to follow one arc. The following itself is still a march; the
  proof is what makes its result trustworthy, and a march that strayed
  fails its proof instead of passing a wrong curve on.
- **Pcurves.** On either spline face the pcurve is `Curve2::Lifted` on a
  spline carrier. When the lifted curve is a `PairSection` on one of its
  own surfaces, evaluation reads that surface's parameters straight from
  the solve.
- **Curve against surface.** A B-spline curve against a B-spline surface
  uses the same certified root isolation, over the curve's own Bezier
  pieces.
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
- Two branches crossing at a saddle (indefinite Hessian) meet at a vertex
  of the section graph.
  1. Newton on the gradient finds the crossing and its Hessian. Rounding is
     measured on the terms that make up the value there
     (`Field2::scale_at`), not the coefficients: a power series read far
     from its origin cancels terms much larger than its coefficients.
  2. The pieces around the crossing that cannot be certified form a
     square. It holds at least the `sqrt(rounding / curvature)` where the
     field is below its own rounding, and it grows to hold every piece
     that fails certification there, up to a thousandth of the window.
     Only pieces wholly inside it are dropped.
  3. Exactly four branch ends must lie at the square, two along each
     direction where the Hessian's form vanishes, on opposite sides.
     Otherwise the point is treated as degenerate (below).
  4. Each end is joined to the crossing by a *bridge* cell: the cubic that
     matches the branch's value and slope at its certified end and its
     tangent at the crossing (the Hessian's direction). It is smooth, its
     Bezier control values bound it, and it leaves the branch by about its
     length to the fourth power, far below rounding for these lengths.
  5. Every branch then ends at the crossing, so the chains stop there: the
     crossing is a vertex, as the Steinmetz ellipses' crossings are.
     A boolean cuts faces there like any other section vertex.
- Where the surfaces touch to higher order (a singular Hessian), any
  even number of branches may end at the point, or none.
  1. Newton on the gradient is damped (Levenberg-Marquardt), so it still
     converges, if only linearly.
  2. The branch ends at the square are checked against the field itself:
     there must be as many as its sign changes round a square three times
     the size. None means the surfaces only touch there.
  3. A tacnode's branches (x ~ c z^2) all arrive along the Hessian's null
     direction, and the cubic bridge reproduces them exactly to second
     order; with no direction left (a higher crossing), each arrives along
     its chord.
  4. Along such a contact the surfaces agree to rounding over a stretch of
     about `rounding^(1/4)`: where along it the branches meet cannot be
     decided in doubles. Every point, the vertex too, lies on both
     surfaces to rounding.
- Where the surfaces are tangent along a whole curve, every piece along
  it is singular; more than a few dozen singular points in one trace hand
  over to a search for such curves.
  1. Along the curve the field vanishes as `c n^m` in the distance `n`
     across it (`m` even: touching; odd: crossing). So every derivative of
     order `m - 1` vanishes exactly on it, and generically one crosses
     zero there regularly: its trace, certified as usual, holds the
     curve. Orders are tried from the first up; a derivative whose own
     trace meets a degenerate point is passed over.
  2. Each curve of that trace is a *line of contact* only where the field
     vanishes along the whole of it (to rounding), its gradient far below
     the gradient just beside it, and changes sign across it (crossing) or
     not (touching). Curves of the derivative that are regular zeros of
     the field, or not zeros at all, are dropped.
  3. The field is traced again with tubes left out about every line of
     contact a thousandth of the window wide: every field has a certified
     tier (below), so no tube needs widening to where rounding stops
     hiding the line. Two branches meeting a line within rounding of each
     other meet it at one vertex. A branch of it ending at a tube
     runs into the line: it is bridged to the line's nearest point (an end
     the line already has, within rounding, is that point), and the line
     cut there, a vertex as at a crossing.
  4. Crossing lines are sections, their points the derivative's zeros
     (the field's own to rounding); touching lines add none.
  5. Where no derivative yields a line (a pole, a curve of touching
     points on the window's edge), the singular points are taken one by
     one as before.
- Finding such points early matters. Where two branches nearly touch,
  certified cells shrink as the square of the distance to the contact. So
  a small piece where both partials may vanish is searched for a singular
  point at once, and its square spares that work. A side's roots are
  bounded by the tighter of the direct and the mean-value bound, and a
  side that needs a deep search makes its piece split instead.
- Windows a whole turn wide start an irrational fraction of a radian past
  `-pi` (and past a face's own seam angle), so a symmetric section's
  special points never sit on the window's edge.

### Certified multiplicity of contacts (#397, #427)

A second field that only touches a traced curve has a root of even
multiplicity, and one that osculates it a root of order four or more.
Bisecting a stretch on the sign of the field alone read f64 noise at the
touch as a pair of spurious sign changes about 2.4e-8 either side of it
(#397), and at a contact of order four it did not finish: the second
derivative's enclosure straddles zero out to the square root of the
stretch's width, so the unsettled stretches multiplied with the depth until
the depth floor, and every leaf the floor reached was labelled a certified
touch (#427).

The roots of the second field `H` along the curve `F = 0` are now isolated
by a Rolle chain of order up to eight, where the enclosures resolve it, on the derivative tower
`H_0 = H`, `H_{j+1} = det(grad F, grad H_j) = F_u ∂_v H_j - F_v ∂_u H_j`
(`curve_tower.rs`, `root_chain.rs`).

- **The tower.** Along the curve, `d(G|C)/ds = ±det(grad F, grad G) / |grad F|`
  with a sign constant along a regular stretch, so `H_{j+1}|C` is `H_j|C`'s
  derivative up to a factor that never vanishes: a zero of `H_j|C` of
  multiplicity `m` is one of `H_{j+1}|C` of multiplicity `m - 1`, and
  `H_{j+1}|C != 0` on a stretch makes `H_j|C` strictly monotone there. The
  values do not depend on which parameter a cell leaves free, so a stretch
  can cross a cell join. (Derivatives of order two and up with respect to
  the free parameter do change at a join; a chain on those could not.) In a
  cell the curve is `(f, w(f))` and `H_{j+1}|C = sigma F_w|C d/df (H_j|C)`,
  `sigma = -1` for `Axis::U` and `+1` for `Axis::V`, so the tower is built
  from univariate Taylor coefficients alone: those of `w` from the implicit
  equation, then of `H|C` and `F_w|C`, each a polynomial of the partials up
  to the order wanted and of the earlier coefficients. Every partial is
  enclosed once over a box certain to hold the stretch's points (the old
  hull, plus the error of the solved value itself, which may reach past
  the end of its bracket), the interval Taylor
  coefficients over the hull being the order-0 model.
- **The chain.** `isolate` splits the parameter range until a leaf has an
  order `k` (at most `MAX_CONTACT_ORDER = 8`) whose `H_k` excludes zero, and
  `chain` climbs from `H_k` down to `H_0`. Between two zones of `H_{j+1}`,
  `H_j` is strictly monotone: a sign change brackets exactly one simple
  zero, which is bisected. A run of zones and ends whose enclosures of `H_j`
  straddle zero is one cluster of multiplicity `1 + sum of the members'`
  (a run of two or more critical points refuses, below).
  The chain allows at most `k - j` zeros of `H_j`, and a cluster between two
  signed neighbours must change sign exactly when its multiplicity is odd;
  either failing refuses.
- **Generic splits.** A stretch an order settles is never split. An
  unsettled one is cut at a cell join, or else at one of 15 fixed fractions
  (the middle first), only where `H_0`, `H_1` and `H_2` are signed, so a cut
  never sits on a root of order two or less; when no fraction qualifies (the
  stretch is inside a flat region), a cut where `H_0` alone is signed. A cut
  that straddles at a higher order makes the chain refuse the leaf, and a
  leaf the chain refuses is split further instead of refusing the whole
  isolation: no zone touches an inner end, and each zero belongs to exactly
  one leaf.
- **Gaps.** Where the enclosures of `H_j` straddle zero over a run of zones
  (from overestimation or a true zero), the gaps between the run and its
  signed neighbours hold no zero of `H_{j+1}`, so `H_j` is monotone there.
  Each gap is searched from the signed neighbour for the last point still
  signed like it, and for a point of the opposite sign, which brackets a
  simple zero. The cluster's extent starts at the last signed point, so every
  zero lies in a bracket or in a cluster. A zero found in a gap counts
  against the cluster's multiplicity. At `H_0` a cluster's extent is halved
  while the enclosure straddles, to the part that does; if the value at the
  cluster is signed, the leaf is refused, not reported (an exact root is never
  signed at the point reported for it, to the tower's rounding).
  On a piece where `H_1` is signed `H_0` is monotone, so the halving
  decides it from the values at its ends. A run of two or more zones of
  `H_{j+1}` whose enclosures of `H_j` all straddle zero (distinct critical
  points the enclosures cannot tell from a zero) refuses, instead of being
  read as one root of the summed order. An inner end that
  straddles zero at an order the chain climbs refuses its leaf before the
  chain starts.
- **Curve ends.** The curve's own ends report. A root there carries the
  multiplicity the chain sees at the end, which depends on how far, at
  rounding level, the end is from the contact. A closed curve's seam is one
  point and is reported once, with the larger of its two multiplicities.
- **Budget and cap.** One call asks at most `ROOT_BUDGET` enclosures (each node costs one more
  than the enclosures it asked), one leaf's chain at most `CHAIN_EVALS` before
  the leaf is split instead, and
  depth is at most 128 with leaves at least `2^-44` of a cell wide.
  Exceeding either, a contact above order eight, or one rounding hides
  refuses as `ExactCurveRefusal::Undecided`: it says nothing of whether the
  curve meets the surface there. There is no degree bound to derive the cap
  from, because the fields are trigonometric series and the curve is
  implicit, so the cap is a named constant.

Soundness:

1. *Enclosure.* For every point of a stretch, `(f, w(f))` lies in the hull,
   so the partials there lie in their enclosures, and each Taylor
   coefficient and each `H_j|C` is a polynomial of them (with a division by
   `F_w`, bounded away from zero over the hull) evaluated in outward-rounded
   intervals. The result contains the true value at every point of the
   stretch at once. The fields are the rounded compositions the module
   already accepts.
2. *Rolle.* If `H_k|C` has no zero on a stretch, `H_j|C` has at most `k - j`
   there, counting multiplicity.
3. *Completeness,* by downward induction: if the zones of `H_{j+1}` cover
   its zeros, then `H_j` is strictly monotone between consecutive zones, a
   zero in an open piece forces rigorously opposite signs at its ends, and
   any other zero lies in a zone whose enclosure straddles. No zero is
   missed.
4. *Nothing invented beyond resolution.* A bisected zone holds a zero by the
   intermediate value theorem. The one claim that is not rigorous is that a
   straddling enclosure is a zero of the multiplicity the chain allows; the
   true count is at most that, with its parity whenever both neighbours are
   signed (checked), so crossing versus touching is always right and two
   zeros closer than the rounding read as one touch.
5. *Termination.* Depth, width and the budget bound the leaves; a chain
   has at most eight levels and at most 100 bisection steps per zone.
6. *Decided when possible.* At a contact of order `mu <= 8` with `H_mu`
   nonzero there, the enclosures converge (they are Lipschitz inclusion
   functions), so a certified stretch exists at finite depth as long as it
   is wider than the contact's flat region, the width below which rounding
   hides it. Otherwise the result is `Undecided`.

A flat stretch of the second field next to crossings, below what the
enclosures resolve, is refused as `Undecided`, never reported as a root
(the box enclosure of `H` over a wide extent overestimates by the extent's
width, so a displacement of a high-order contact by 1e-6 is already beyond
it). Orders 7 and 8 on trigonometric fields and circles are likewise often
refused, and the tests only show polynomial fields decided there.

The multiplicity of a root at the curve's own end is the order seen at the
end, at rounding distance from a contact the curve may stop short of; it is
documented as such on the hit types, not as a crossing/touch class.

The bridge cells of a crossing keep their 64-sample sign scan, which is not
certified (#441). The caller's merge of close cuts (#381) stays, as the
tolerance reading for cuts that are not one root; the roots themselves are
no longer merged within 1e-9, since zones are disjoint by construction.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Marching (predictor–corrector) as OCCT's `IntWalk` does | Step size is a guess: loops smaller than a step are missed, branches can be jumped near tangencies, and nothing certifies that every component was found. |
| Fit B-splines to the traced points | Turns an exact definition into an approximation with a tolerance, in the module that promises none; pcurves and edges would disagree by the fit error. |
| Resultants and root isolation of the full algebraic curve | Exact topology, but the degrees are large: a torus against a torus is degree 8 in each of two variables before the resultant. Isolating the discriminant's roots with big integers per face pair is far slower than subdivision, and the curve still needs evaluating afterwards. |
| A guide polyline with normal projection | Needs a step size to build and a certificate that the normal lines meet the curve once. Monotone cells give both by construction. |
| A node budget alone for contacts of higher order (#427) | It terminates, but it refuses plausible geometry: a cylinder osculating a section, or a plane tangent at a section's extreme. It stays as the backstop. |
| A Rolle chain per cell on `h^(k)` with respect to the free parameter (the direct generalisation of #397) | Derivatives of order two and up change at a cell join where the free axis changes, so a chain cannot cross a join; and a stretch end inside the contact's flat region produces spurious zeros at that end. |
| Exact classification of `h` at the zeros of `h'` | The touch point is an irrational root of the traced field, and the exact series tier only decides at `f64` points; the algebraic route needs resultants, rejected above. |
| Taylor models | They cut the subdivision depth but meet the same rounding floor and need a polynomial-plus-remainder composition facility. The interval Taylor coefficients over the hull are the order-0 model; a later optimisation. |

## Consequences

**Positive**

- Every analytic pair's section now exists exactly, whether or not a closed
  form does, and so does every section involving B-splines, including
  branches crossing or touching each other where the surfaces touch.
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

- A value's sign is trusted only clear of its rounding, on a piece's
  sides and across the piece. Below it:
  - **B-spline fields decide it exactly** (`exact_field`). On a cell
    `[a, b] x [c, d]` the field times `(b - a)^p (d - c)^q` is
    `sum c_kl C(p,k) C(q,l) (u - a)^k (b - u)^(p-k) (v - c)^l (d - v)^(q-l)`,
    which needs no division: its value at any `f64` point is exact in
    dyadic arithmetic. Over a box, de Casteljau subdivision with
    homogeneous weights (`(b - x) P_i + (x - a) P_(i+1)`) gives its
    Bernstein coefficients there, again without division; all of one
    strict sign prove that sign over the box. Derivatives' signs come
    from the coefficients' differences. Each question is asked in
    outward-rounded intervals first and in big integers only where they
    cannot decide. So a B-spline field's trace has no precision limit:
    pieces are excluded, sides found monotone and roots placed on exact
    signs, lines of contact of every order the field's degree allows are
    found (`x^21` against a plane in the tests), and their tubes need no
    widening.
  - **Series fields decide it with certified arithmetic**
    (`exact_series`, #181). Harmonics have no exact value, but a sign
    does not need one:
    - *Zero at a point.* At `f64` parameters every harmonic angle is an
      integer multiple of one dyadic angle `w`, so with `z = e^(iw)` the
      value is a Laurent polynomial in `z` with dyadic complex
      coefficients. `z` is transcendental (Lindemann), so the value is
      zero exactly when every collected coefficient is: decided in dyadic
      arithmetic.
    - *Sign at a point.* Otherwise the value is nonzero, and fixed-point
      intervals (`axiolid_exact::FixedInterval`: big-integer bounds
      rounded outward, `sin` and `cos` from their Taylor series with the
      Lagrange remainder added) at rising precision show its sign.
    - *Sign over a box.* A Taylor form about a corner: the value and
      derivatives there, certified, make a polynomial, and the next
      order's derivative bounds over the box bound the remainder. Its
      Bernstein coefficients all clearing the remainder with one sign
      prove that sign. The order rises until they do or cannot.
    The trace consults this tier only where the field is *flat*, its
    gradient below `1e-4` of its terms' magnitudes: along a line of
    contact or at a singular point, where `f64` bounds fail at every
    size. Elsewhere `f64` and subdivision decide as before, so ordinary
    sections are traced exactly as they were, at their old cost. Lines
    of contact that rounding hides over more than a twentieth of the
    window (`sin(u)^13` in the tests) are traced, no longer
    `Undecided`. Such contact needs inflections quadrics and tori do not
    have, so no known input needed it; the tier is there for the one
    that does.
- A pcurve on a face whose trace is too costly (a tacnode read far from
  the face's parameter origin) falls back to the section's own space
  curve read on the face (`Curve2::Lifted`).
- Very thin sections or nearly tangent pairs can exhaust the trace's box
  budget (400 000 boxes). That is refused as `Undecided`.

## Relation to existing code

- `crates/representations/analytic/curve/src/implicit.rs`: `Field2`,
  `ImplicitCurve2`, `ImplicitSection3`, `Carrier`.
- `crates/algorithms/parametric/nurbs/src/field.rs`: field construction,
  interval bounds, partials.
- `crates/algorithms/parametric/nurbs/src/implicit_trace.rs`: the certified
  trace.
- `crates/algorithms/parametric/nurbs/src/exact_field.rs` and
  `exact_series.rs`: the certified tiers of B-spline and series fields.
- `crates/algorithms/exact/src/fixed.rs`: `FixedInterval`, certified
  `sin` and `cos`.
- `crates/algorithms/parametric/nurbs/src/implicit_section.rs`:
  `implicit_surface_intersection`, the fallback of
  `exact_surface_intersection`.
- `crates/representations/analytic/curve/src/pair_section.rs`:
  `PairSection3`, the section of two B-splines.
- `crates/algorithms/parametric/nurbs/src/pair_trace.rs`:
  `spline_pair_intersection` and B-spline curve/surface crossings.
- `crates/algorithms/parametric/nurbs/src/pair_certify.rs`: enclosures,
  chord proofs and certified crossing isolation.
- `crates/algorithms/parametric/nurbs/src/interval.rs`: the outward-rounded
  interval of those certificates.
- `crates/algorithms/parametric/nurbs/src/root_chain.rs` and
  `curve_tower.rs`: the Rolle chain and the derivative tower along a traced
  curve (`implicit_ops::roots_along`).
- `crates/algorithms/parametric/evaluate/src/curve.rs`: evaluation and
  inversion of the new families. Inversion now also covers the ADR 0076
  graphs.
- `crates/algorithms/parametric/nurbs/tests/implicit_section.rs`: points on
  both surfaces, continuity, and completeness against a dense scan of tube
  circles.
- `crates/algorithms/parametric/nurbs/tests/pair_section.rs`: two sheets
  crossing (completeness against a dense scan), a closed loop inside one
  patch pair, and a spline curve crossing a spline surface.

## Amendment 2026-11-07: certified chord bounds of lifted pcurves (#414)

`bound::chord_bound2` bounded implicit pcurves (#249) and, with #414, the
graphs of ADR 0076, but not the pcurve a B-spline section leaves on an
analytic face (`Curve2::Lifted`), nor a section's pcurve on a B-spline
(`PairSection3` read on one of its surfaces). The compilers refused those
trims as unbounded. No boundary or dependency changes: everything is in
`axiolid-evaluate`.

- **One continuation.** A lifted pcurve is the solution of `S(p(t)) =
  c(t)`. On a piece `[a, b]` of the span, with the evaluated ends `p_a`,
  `p_b` and the box `P` of their hull widened by `rho`, a certified
  separation `sigma` of the surface over `P` gives `|J w| >= sigma |D w|`
  at every point of `P` and `|S(x) - S(y)| >= sigma |D (x - y)|` on it
  (`D` scales the parameters by the columns' norms at the box's centre).
  Then `|D p'| <= sup |c'| / sigma`, and differentiating twice `J p'' = c'' -
  H[p', p']` bounds `|D p''|`. The lift of an exact end moves at most `|p'|
  h / 2` to the middle of the piece; if that and the end's error are below
  `rho`, the lifts from both ends stay in `P` and injectivity makes them
  one, so the evaluated ends lie on one continuous branch. `rho` grows
  until that holds, or the piece is refused. The chord of the pieces'
  bounds is `min(h sup|p'| / 2, h^2 sup|p''| / 8)` over the best of the
  splits into 1 to 64 equal pieces, the graph bounds' scheme.
- **Separation by Weyl's inequality.** The separation is `sqrt(1 - |cos
  theta|)` at the box's centre, `theta` the angle between the columns, less
  the Frobenius norm of the columns' scaled spread over the box (the second
  partials' suprema times the half widths) and a padded jet: the same for
  every carrier, from `SurfaceBoundOracle::for_carrier` (closed forms for
  planes, ruled surfaces, spheres and tori; the derivative nets of a
  B-spline). It is `None` where it cannot be shown positive, so a pole, an
  apex, a fold, a box too wide and a B-spline crease (`!smooth`) refuse.
  The interval Gram matrix and the closed-form first fundamental form were
  rejected: the first needs enclosures of every family, the second is not a
  hull bound and cannot prove injectivity.
- **Ends are certified, not assumed.** Each evaluated end is off its exact
  lift by at most `|S(p_i) - c_i| / (sigma min(|S_u|, |S_v|))`, the residual
  enclosed in outward-rounded arithmetic (closed forms for the analytic
  carriers, the rational tensor sum of an interval basis for a B-spline),
  and the larger of the two is added to the bound, because the chord joins
  the evaluated ends. `chord_bound2`'s contract is now that: the chord of
  what `evaluate2` returns. The space curve is read as the evaluator
  computes it and lies on the carrier. For an intrinsic curve the same
  contract adds the distance of each evaluated (quadrature) end from a
  certified `IntrinsicPositionTable` point, with that point's error.
- **The space curve's derivative bounds** are those of
  `curve_derivative_bounds3`; for the sections of an analytic carrier they
  come through the carrier from the plane curve `q(t)` read through the
  surface (`|c'| <= S_u |q_u'| + S_v |q_v'|` and the second order
  analogue): a graph's root (ADR 0076) or one cell of an implicit curve,
  whose rate along its free parameter is constant and whose bridge cell is
  a cubic bounded by its Bezier control values.
- **Pair sections.** A pair section is the same continuation on `G(s, x)
  = (S1(a) - S2(b), n . (S1(a) - P0) - s |d|)` in the four parameters of
  both surfaces, over the product of two boxes. The smallest singular value
  `mu` of the scaled matrix is `(1 - |I - R M'|) / |R| - |dM'|` (`R` an
  approximate inverse at the centres, certified by the residual, `dM'` by
  both surfaces' spreads, the first counted twice); `|y'| <= |d| / mu`, and
  the end error is the 4D residual over `mu`. A side's pcurve is `y`
  unscaled through its own smaller scale factor.
- **Breaks.** The rate of an implicit section changes where one cell
  joins the next, and a pair section's chord plane at each interior node,
  so neither is `C^2` in its parameter there: `continuity_breaks3` names
  them, a lifted pcurve inherits its space curve's, and a span across one
  is refused as across a B-spline knot.
- **Preconditions and refusals.** The boxes on a B-spline are clamped to
  its knot domain, which the section lies in over the span (the assumption
  of the face bound's lens); an evaluated end outside it is refused. A
  guide that could change the turn of an angle inside a piece, a space
  curve without derivative bounds (`Intrinsic3`, `Polyline3`, a pair section
  read on a third surface) and a B-spline carrier outside a pair section
  (which cannot be evaluated) return `None`.
- **Tightness.** The pcurve's second derivative is the difference of `c''`
  and the carrier's curvature along `p'`, nearly equal for a curve that
  hugs the surface; the bound adds their magnitudes. It is within a few
  times the sampled deviation on a flat carrier and sound but loose on a
  curved one, most near an inflection. A tighter bound needs the tangential
  parts of `c''` and of the carrier's curvature separately, which the
  derivative suprema do not give.

Not done here: the mesh compiler samples pcurve edges at uniform spans and
ignores their breaks, so it still reports such trims unbounded (#436);
`chord_bound3` and derivative bounds of the sections in space (#437); the
chain bound's evaluated ends (#438).
