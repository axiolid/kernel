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
the whole-domain check with seven mutants: 19 / 19.
"""
import probe_lib

I = "crates/execution/compile/src/brep/inversion.rs"
B = "crates/execution/compile/src/brep.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "edge_inversion"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::inversion"],
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
     "        let outcome = match run(search.strict)? {",
     "        let outcome = match run(search.model)? {"),
    ("certified B-spline parameters not polished", I,
     "        for _ in 0..POLISH_STEPS {",
     "        for _ in 0..0 {"),
    ("polish not clamped to the domain", I,
     "        let (u_lo, u_hi) = (u_lo.max(u_domain.0), u_hi.min(u_domain.1));",
     "        let (u_lo, u_hi) = (u_lo, u_hi);"),
    # Edge curves trimmed at their vertices (#319).
    ("edge curve taken whole, not trimmed at its vertices", I,
     "    let domain = vertex_span(curve, (edge.start, start), (edge.end, end), ctx.tolerance)?;",
     "    let domain = vertex_span(curve, (edge.start, start), (edge.end, end), ctx.tolerance)\n"
     "        .map(|_| axiolid_reference::curve::domain3(curve))?;"),
    ("a one-vertex edge is half a turn", I,
     "            let sweep = if start_id == end_id {\n                turn",
     "            let sweep = if start_id == end_id {\n                0.5 * turn"),
    ("arc runs backwards from its start vertex", I,
     "                (t1 - t0).rem_euclid(turn)",
     "                (t0 - t1).rem_euclid(turn)"),
    ("two vertices at one point of a circle not refused", I,
     "            } else if (end - start).length() <= linear {",
     "            } else if false {"),
    ("vertices located at a loose tolerance", I,
     "        axiolid_reference::curve::invert3(curve, vertex, tolerance).map_err(|error| {",
     "        axiolid_reference::curve::invert3(curve, vertex, Tolerance::new(1.0, 1e-9).unwrap())\n"
     "            .map_err(|error| {"),
    ("line clamped to its unit domain", I,
     "        Curve3::Line(_) => axiolid_core::Interval { start: t0, end: t1 },",
     "        Curve3::Line(_) => axiolid_core::Interval { start: t0, end: t1.min(1.0) },"),
    ("other families' whole domain need not end at the vertices", I,
     "        if bounded && gap(domain.start, start)? <= linear && gap(domain.end, end)? <= linear {",
     "        if bounded {"),
    ("span terms left out of the report", I,
     "            self.trim = self.trim.max(near + lens);",
     "            self.trim = self.trim.max(lens);"),
    ("lens containment not checked", I,
     "            if on_box >= sigma {",
     "            if true {"),
    # Planar caps bounded by curved edges.
    ("caps bounded by curved edges keep the planar path", B,
     "        if !surface_is_planar(surface) || curved_plane(brep, graph, face) {",
     "        if !surface_is_planar(surface) {"),
    ("every edge curve reads as straight", B,
     "            Some(axiolid_model::GeometryNode::Curve3(curve)) => !inversion::straight_curve(curve),",
     "            Some(axiolid_model::GeometryNode::Curve3(curve)) => !inversion::straight_curve(curve) && false,"),
    ("sampling ignores the certified terms", I,
     "            if certified.is_none_or(|worst| worst <= budget) {",
     "            if certified.is_none_or(|worst| worst <= budget) || true {"),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
