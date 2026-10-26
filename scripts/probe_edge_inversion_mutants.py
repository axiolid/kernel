"""Mutation probe for curved faces whose edges carry only 3D curves (#310).

Each mutant weakens one decision of `crates/execution/compile/src/brep/
inversion.rs` or its caller in `brep.rs`: the refusal at a cone apex or
sphere pole, the unwrapping and its quarter-turn guard, the strict pass of
the certified B-spline inversion, the polish inside its enclosure, the
domain clamp, the trim of an edge curve at its vertices (#319), and
the terms the deviation report adds per span (the parameter chord's image,
the lens and the containment that sizes it).
`tests/edge_inversion.rs` or the module's unit tests must fail on it. A
build error also counts as killed.

Not listed, as equivalent at test resolution:
- the parameter chord's image left out of a span's term: the report adds
  the face's triangle bound to its trim terms, and that bound already
  covers the surface over every triangle edge, the chord included, so no
  sample can tell the image from it.
- the lens term left out: it is of second order in a width the edge's
  sampling keeps far below the budget, so dropping it moves no bound
  past a measured sample.
  The containment check that sizes it is listed.
- `align` and the B-spline parameter check: both act only on a periodic
  B-spline face, which no fixture has.
- the #265 weld pre-pass skipping planar faces bounded by curved edges:
  such a face takes the curved path, which does not take recorded
  corners, so the skip keeps it from recording any. Its corners touch one
  of its own chords only in a cap whose vertex polygon folds onto itself,
  which no fixture builds; a recorded corner would leave the shell
  `OpenSolid`, never a wrong solid.

First run: 8 / 13. The survivors were the two mutants above, the
quarter-turn guard, the polish and the span terms; the guard and the
polish got unit tests, and the flat patch trimmed by a 3D circle (whose
bound is its trim terms alone) was added for the span terms. 11 / 11
killed after. The planar caps (closed cylinders, #310 follow-up) added two
mutants: 13 / 13. Trimming edge curves at their vertices (#319) replaced
the whole-domain check with seven mutants: 19 / 19. The edge pass (#336),
with the D-shaped half cylinder in every face order, added four: 23 / 23.
B-spline and polyline trims (#333) moved the #319 mutants to `trim.rs` and
added eleven: 33 / 34 on the first run, the seam cut surviving (every
fixture's seam span is matched by an equal one elsewhere); a unit test on
a span across a square's seam corner kills it: 34 / 34.

The early off-surface exit and the overrun retry (#320) added four mutants: 38 / 38 killed.

Not listed for #320, as equivalent: the model pass, run after a strict
overrun, returning an off-surface refusal that is passed on as the
refusal. Dropping the arm turns it into the strict `BudgetExceeded`, and
only a point the strict pass cannot decide yet the model pass proves off
the surface distinguishes the two, which no fixture builds.

Not listed for the trims, as equivalent:
- a circle's or ellipse's span wrapped at its turn: the families evaluate
  periodically, so wrapping names the same points.

Not listed for the edge pass, as equivalent at test resolution:
- a planar face's bound without its samples' distance off their edge: the
  fixtures' straight edges are sampled on their segments, so the term is
  rounding.
- the pcurve use the samples are taken through: every pcurve use must
  agree with them within the budget, and no fixture shares an edge between
  two pcurve uses that differ.
"""
import probe_lib

I = "crates/execution/compile/src/brep/inversion.rs"
B = "crates/execution/compile/src/brep.rs"
E = "crates/execution/compile/src/brep/edges.rs"
T = "crates/execution/compile/src/brep/trim.rs"
N = "crates/algorithms/parametric/nurbs/src/certified_surface_projection.rs"
NURBS_INVERSION = ["-p", "axiolid-nurbs", "--test", "certified_surface_inversion"]
NURBS_PROJECTION = ["-p", "axiolid-nurbs", "--test", "certified_surface_projection"]
EDGE_INVERSION = ["-p", "axiolid-mesh-compile", "--test", "edge_inversion"]
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "edge_inversion"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::inversion"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::trim"],
]

MUTANTS = [
    ("edges without a pcurve refused again", B,
     "            let Some(pcurve) = use_.pcurve else {\n",
     "            let Some(pcurve) = use_.pcurve else {\n"
     "                return Err(GeomError::Unsupported {\n"
     "                    backend: crate::BACKEND_ID,\n"
     "                    operation: axiolid_contracts::Operation::Tessellation,\n"
     "                });\n"
     "                #[allow(unreachable_code)]\n"),
    ("no refusal at an apex or pole", I,
     "        if radial.is_nan() || radial <= self.tolerance.linear() {",
     "        if radial.is_nan() {"),
    ("edge parameters not unwrapped", I,
     "            localized_parameter(params[i - 1].x, params[i].x, pu),",
     "            params[i].x + 0.0 * localized_parameter(params[i - 1].x, params[i].x, pu),"),
    ("quarter-turn guard reads half a period", I,
     "period.is_none_or(|p| d.abs() < 0.25 * p)",
     "period.is_none_or(|p| d.abs() < 0.5 * p)"),
    ("B-spline inverted at the model tolerance only", I,
     "        let outcome = match run(search.strict) {",
     "        let outcome = match run(search.model) {"),
    ("a strict budget overrun not accepted from the model pass", I,
     "                Ok(Ok(certificate)) => Ok(certificate),\n",
     ""),
    # The early off-surface exit of the certified surface search (#320).
    ("no separation exit", N,
     [("lower > threshold", "lower > threshold && false")],
     [NURBS_INVERSION, EDGE_INVERSION]),
    ("separation exit taken on the surface", N,
     [("lower > threshold", "lower >= 0.0 * threshold")],
     [NURBS_INVERSION]),
    ("projection given a separation threshold", N,
     "        surface, target, options, false, false, None,\n",
     "        surface, target, options, false, false, Some(options.distance_tolerance().linear()),\n",
     [NURBS_PROJECTION]),
    ("certified B-spline parameters not polished", I,
     "        for _ in 0..POLISH_STEPS {",
     "        for _ in 0..0 {"),
    ("polish not clamped to the domain", I,
     "        let (u_lo, u_hi) = (u_lo.max(u_domain.0), u_hi.min(u_domain.1));",
     "        let (u_lo, u_hi) = (u_lo, u_hi);"),
    # Edge curves trimmed at their vertices (#319), now in `trim.rs`.
    ("edge curve taken whole, not trimmed at its vertices", T,
     "    Ok(CurveSpan { span, period })",
     "    let _ = span;\n"
     "    Ok(CurveSpan { span: axiolid_reference::curve::domain3(curve), period })"),
    ("a one-vertex edge is half a turn", T,
     "            let sweep = if start_id == end_id {\n                turn\n",
     "            let sweep = if start_id == end_id {\n                0.5 * turn\n"),
    ("arc runs backwards from its start vertex", T,
     "                (t1 - t0).rem_euclid(turn)",
     "                (t0 - t1).rem_euclid(turn)"),
    ("two vertices at one point of a closed curve not refused", T,
     "            } else if (end - start).length() <= linear {",
     "            } else if false {"),
    ("vertices located at a loose tolerance", T,
     "                axiolid_reference::curve::invert3(curve, vertex, tolerance).map_err(|error| {",
     "                axiolid_reference::curve::invert3(curve, vertex, Tolerance::new(1.0, 1e-9).unwrap())\n"
     "                    .map_err(|error| {"),
    ("line clamped to its unit domain", T,
     "        None => Interval { start: t0, end: t1 },",
     "        None => Interval { start: t0, end: t1.min(1.0) },"),
    ("other families' whole domain need not end at the vertices", T,
     "            if bounded && gap(domain.start, start)? <= linear && gap(domain.end, end)? <= linear {",
     "            if bounded {"),
    # B-spline and polyline edge curves trimmed at their vertices (#333).
    ("an off-curve B-spline vertex accepted", T,
     [("    if certificate.distance_lower_bound > linear {", "    if false {"),
      ("    if residual.is_nan() || residual > linear {", "    if false {")],
     TARGETS),
    ("the B-spline cover judged unrefined", T,
     "    let intervals = minimiser_cover(curve, vertex, cells, threshold)?;",
     "    let intervals = cells;\n    let _ = (minimiser_cover, threshold);"),
    ("a separated cover read as unique", T,
     "        _ => false,\n    };\n    if !connected {",
     "        _ => true,\n    };\n    if !connected {"),
    ("a cover across the seam read as separated", T,
     "        ([first, last], Some(domain)) => first.0 <= domain.start && last.1 >= domain.end,",
     "        ([first, last], Some(domain)) => first.0 <= domain.start && last.1 >= domain.end && false,"),
    ("a closed B-spline read as open", T,
     "                    let period = axiolid_nurbs::PeriodicCurve3::try_new(spline, tolerance)\n"
     "                        .ok()\n",
     "                    let period = axiolid_nurbs::PeriodicCurve3::try_new(spline, tolerance)\n"
     "                        .ok()\n                        .filter(|_| false)\n"),
    ("a closed polyline read as open", T,
     "                    let period = polyline.closed.then_some(domain);",
     "                    let period = polyline.closed.then_some(domain).filter(|_| false);"),
    ("B-spline vertex parameter not polished", T,
     "    for _ in 0..POLISH_STEPS {",
     "    for _ in 0..0 {"),
    ("polyline tolerance band ignored", T,
     "        if reach >= 0.0 {",
     "        if false {"),
    ("spans not split at breaks", T,
     "        cuts.retain(|&t| lo < t && t < hi);",
     "        cuts.clear();"),
    ("the seam not a cut", T,
     "            cuts.push(domain.end);",
     "            let _ = domain.end;"),
    ("split points' distance off the chord left out", T,
     "        let chord = chord + corner;",
     "        let chord = chord + 0.0 * corner;"),
    ("span terms left out of the report", I,
     "            self.trim = self.trim.max(near + lens);",
     "            self.trim = self.trim.max(lens);"),
    ("lens containment not checked", I,
     "            if on_box >= sigma {",
     "            if true {"),
    # Planar caps bounded by curved edges.
    ("caps bounded by curved edges keep the planar path", B,
     "    !surface_is_planar(surface) || curved_plane(brep, graph, face)\n",
     "    !surface_is_planar(surface)\n"),
    ("every edge curve reads as straight", B,
     "            Some(axiolid_model::GeometryNode::Curve3(curve)) => !inversion::straight_curve(curve),",
     "            Some(axiolid_model::GeometryNode::Curve3(curve)) => !inversion::straight_curve(curve) && false,"),
    ("sampling ignores the certified terms", I,
     "            if certified.is_none_or(|worst| worst <= budget) {",
     "            if certified.is_none_or(|worst| worst <= budget) || true {"),
    # The edge pass (#336).
    ("no edge pass: the first face to reach an edge samples it", E,
     "    let plans = plan_edges(ctx, shells)?;",
     "    let plans = plan_edges(ctx, &shells[..0])?;"),
    ("planar faces take only their edges' vertices", B,
     "    let mut rings = planar_rings(brep, face, Some((cache, &mesh.positions)))?;",
     "    let mut rings = planar_rings(brep, face, None)?;"),
    ("the last face's requirement wins, not the strictest", E,
     "                plan.segments = plan.segments.max(segments);",
     "                plan.segments = segments;"),
    ("a planar face takes a reversed use's samples forwards", B,
     "            let (indices, _) = shared.in_direction(use_.orientation);",
     "            let (indices, _) = shared.in_direction(Orientation::Forward);"),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
