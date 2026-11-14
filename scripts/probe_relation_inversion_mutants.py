"""Mutation probe for edges without pcurves on surface relations (#454, ADR 0110).

Each mutant weakens one decision of `crates/execution/compile/src/
surface_relation/invert.rs` (the cover, its local speed and its work cap,
uniqueness, the spread rule for folds and turning generators, the closed
seam's span check, the pole and input refusals, the closed-form other
parameter and its wrap, the revolved line's bracket), the tangents and open
domain a relation hands the lens, the inverter's periods, or the span's
local arc bound in `brep/trim.rs`. The module's unit tests,
`tests/relation_edge_inversion.rs`, `tests/relation_faces.rs` or the
tessellator's own tests must fail on it. A build error also counts as
killed.

Not listed, as not observable (no fixture, not equivalent):
- R8, the coarse cell without a certified speed bound, and the same refusal
  in the spread: the four supported swept families always have one.
- the `InvalidInput` refusals of non-finite arithmetic (a distance, the
  revolved line's bracket, the extruded plane's parameters, the residual,
  the spread): finite, valid curves give finite values throughout, and no
  fixture builds one that does not. The non-finite sample's own refusal is
  listed.
- the revolved line's rho bracket without its `m^2` term (equivalent): it
  only widens the bracket to a superset, whose extra cells lie more than
  `tau` from the sample and are dropped; the answer is the same, at more
  work.
- `COARSE_CELLS` (8 instead of 64 survives): the cells are only where the
  bisection starts and where the speed bound is taken; fewer of them cost
  more splits and looser bounds, not another answer.
- the `ROUNDING` inflation of the speed bounds: a relative 1e-12 widening.
"""
import probe_lib

V = "crates/execution/compile/src/surface_relation/invert.rs"
S = "crates/execution/compile/src/surface_relation.rs"
P = "crates/execution/compile/src/brep/support.rs"
I = "crates/execution/compile/src/brep/inversion.rs"
T = "crates/execution/compile/src/brep/trim.rs"
# The span's arc bound has its own unit test; the faces it sizes also mesh
# without it, only dearer.
TRIM = [["-p", "axiolid-mesh-compile", "--lib", "brep::trim"]]
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--lib", "surface_relation"],
    ["-p", "axiolid-mesh-compile", "--test", "relation_edge_inversion"],
    ["-p", "axiolid-mesh-compile", "--test", "relation_faces"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::"],
]

MUTANTS = [
    # The cover.
    ("the cover never drops a cell", V,
     "            if cell.lower() > self.tau.min(best + slack) {\n",
     "            if false {\n"),
    ("the Lipschitz term is zero", V,
     "        0.5 * self.speed * (self.b - self.a)\n",
     "        0.0 * self.speed * (self.b - self.a)\n"),
    ("the drop bound loses tau", V,
     [("            if cell.lower() > self.tau.min(best + slack) {\n",
       "            if cell.lower() > best + slack {\n"),
      ("        let threshold = self.tau.min(best + slack);\n",
       "        let threshold = best + slack;\n")],
     TARGETS),
    ("the minimiser slack is ten times wider", V,
     "const MINIMISER_SLACK: Scalar = 1e-3;\n",
     "const MINIMISER_SLACK: Scalar = 1e-2;\n"),
    ("an off-surface floor is not sharpened", V,
     "            let settled = cell.lower() >= FLOOR_SHARPNESS * best\n",
     "            let settled = true\n"),
    # Uniqueness.
    ("two separated parts are accepted", V,
     "            Some((start, (last.0 - period, first.1)))\n"
     "        }\n"
     "        _ => None,\n",
     "            Some((start, (last.0 - period, first.1)))\n"
     "        }\n"
     "        _ => Some((nearest, parts[0])),\n"),
    ("the periodic seam join is removed", V,
     "        [first, last] if bracket.periodic && first.0 <= bracket.start && last.1 >= bracket.end => {\n",
     "        [first, last] if false && first.0 <= bracket.start && last.1 >= bracket.end => {\n"),
    ("a fold's spread generators are accepted", V,
     "        if moving + turning > self.tau {\n",
     "        if false {\n"),
    ("the spread drops its first-order term", V,
     "        let (moving, turning) = (r * speed * up, 0.5 * r * r * second * up);\n",
     "        let (moving, turning) = (0.0, 0.5 * r * r * second * up);\n"),
    ("the spread drops its second-order term", V,
     "        let (moving, turning) = (r * speed * up, 0.5 * r * r * second * up);\n",
     "        let (moving, turning) = (r * speed * up, 0.0);\n"),
    ("generators turning about the sample are named a fold", V,
     "            return Err(if turning > moving {\n",
     "            return Err(if false {\n"),
    ("a child keeps its coarse cell's speed bound", V,
     "                cell.speed.min(local)\n",
     "                cell.speed\n"),
    ("a child's local speed drops its acceleration term", V,
     "            (self.generator_speed(middle)? + 0.5 * cell.second * (b - a)) * (1.0 + ROUNDING);\n",
     "            (self.generator_speed(middle)? + 0.0 * cell.second * (b - a)) * (1.0 + ROUNDING);\n"),
    ("R6, the work cap, never stops the cover", V,
     "            if *work > RELATION_COVER_WORK {\n",
     "            if false {\n"),
    ("a span across a closed spline's seam is accepted", V,
     "            .any(|w| (swept(&w[1]) - swept(&w[0])).abs() > half)\n",
     "            .any(|w| (swept(&w[1]) - swept(&w[0])).abs() > 2.0 * half)\n"),
    ("a closed swept spline is read as open", V,
     "                Some(domain) if closes(curve.point(domain.start)?, curve.point(domain.end)?) => {\n",
     "                Some(domain) if false && closes(curve.point(domain.start)?, curve.point(domain.end)?) => {\n"),
    ("a periodic parameter is not wrapped into its turn", V,
     "            t.rem_euclid(bracket.end - bracket.start) + bracket.start\n",
     "            t\n"),
    ("a non-finite sample is not refused up front", V,
     "        if !point.is_finite() {\n",
     "        if false {\n"),
    ("a swept line without a direction is not refused", V,
     "    if !(direction.is_finite() && direction.length() > 0.0) {\n",
     "    if false {\n"),
    # Poles and the residual.
    ("R1, a sample on the axis, is not refused", V,
     "                if probe.rho.is_nan() || probe.rho <= self.tau {\n",
     "                if false {\n"),
    ("R2, a sample at a ruled apex, is not refused", V,
     "                        if near.is_nan() || near <= self.tau {\n",
     "                        if false {\n"),
    ("R3, a sample on a collapsed ruling, is not refused", V,
     "                if length.is_nan() || length <= self.tau {\n",
     "                if false {\n"),
    ("the residual check accepts everything", V,
     "        if off > self.tau {\n"
     "            return Err(self.off_surface(floor));\n",
     "        if false {\n"
     "            return Err(self.off_surface(floor));\n"),
    ("R7, a line extruded along itself, is not refused", V,
     "                if normal.length() <= PARALLEL * e.length() {\n",
     "                if false {\n"),
    # The other parameter.
    ("the revolution's atan2 arguments are swapped", V,
     "                let u = axis\n"
     "                    .cross(q_perp)\n"
     "                    .dot(probe.perp)\n"
     "                    .atan2(q_perp.dot(probe.perp));\n",
     "                let u = q_perp.dot(probe.perp).atan2(axis.cross(q_perp).dot(probe.perp));\n"),
    ("the revolution's sine changes sign", V,
     "                let u = axis\n"
     "                    .cross(q_perp)\n"
     "                    .dot(probe.perp)\n"
     "                    .atan2(q_perp.dot(probe.perp));\n",
     "                let u = (-axis.cross(q_perp).dot(probe.perp)).atan2(q_perp.dot(probe.perp));\n"),
    ("the extrusion's v divides by |d|, not d.d", V,
     "                    (probe.point - c).dot(*direction) / direction.dot(*direction),\n",
     "                    (probe.point - c).dot(*direction) / direction.length(),\n"),
    ("the extruded plane's v divides by |d|, not d.d", V,
     "        let v = (probe.point - c).dot(*direction) / direction.dot(*direction);\n",
     "        let v = (probe.point - c).dot(*direction) / direction.length();\n"),
    ("the extruded plane's t changes sign", V,
     "        let t = -w_perp.dot(e_perp) / e_perp.dot(e_perp);\n",
     "        let t = w_perp.dot(e_perp) / e_perp.dot(e_perp);\n"),
    ("the ruled v is not clamped", V,
     "                let v = ((probe.point - a).dot(d) / (length * length)).clamp(0.0, 1.0);\n",
     "                let v = (probe.point - a).dot(d) / (length * length);\n"),
    ("the Gauss-Newton polish is disabled", V,
     "        for _ in 0..POLISH_STEPS {\n",
     "        for _ in 0..0 {\n"),
    # The revolved line's bracket.
    ("the h bracket's centre drops the line's height", V,
     "            let (centre, half) = ((probe.height - h0) / e_a, tau / e_a.abs());\n",
     "            let (centre, half) = (probe.height / e_a, tau / e_a.abs());\n"),
    # What the inverter and the lens read of a relation.
    ("a relation's inverter has no periods", I,
     "                periods: surface.periods(),\n",
     "                periods: (None, None),\n"),
    ("a relation's support has no tangents", P,
     "            Self::Relation { surface, .. } => surface.tangents(u, v),\n",
     "            Self::Relation { .. } => Err(axiolid_contracts::GeomError::Degenerate(\n"
     "                \"no tangents\".to_owned(),\n"
     "            )),\n"),
    ("a ruled surface's open domain is unbounded", S,
     "            Self::Ruled { .. } => ((0.0, 1.0), (0.0, 1.0)),\n",
     "            Self::Ruled { .. } => (unbounded, unbounded),\n"),
    ("the revolution's S_u drops the rotation", S,
     "                    axis.cross(rotate(*axis, u, w)),\n",
     "                    axis.cross(w),\n"),
    ("the revolution's S_v drops the rotation", S,
     "                    rotate(*axis, u, curve.tangent(v)?),\n",
     "                    curve.tangent(v)?,\n"),
    ("the extrusion's S_v is the unit direction", S,
     "            Self::Extrusion { curve, direction } => Ok((curve.tangent(u)?, *direction)),\n",
     "            Self::Extrusion { curve, direction } => Ok((curve.tangent(u)?, direction.normalize())),\n"),
    ("the ruled S_u reads the first side only", S,
     "                    a.derivative(u)? * (1.0 - v) + b.derivative(u)? * v,\n",
     "                    a.derivative(u)?,\n"),
    ("the ruled S_v drops the first side", S,
     "                    b.at(u)? - a.at(u)?,\n",
     "                    b.at(u)?,\n"),
    ("a revolved spline's open domain is unbounded", S,
     "            Self::Revolution { curve, .. } => (unbounded, swept(curve)),\n",
     "            Self::Revolution { .. } => (unbounded, unbounded),\n"),
    ("an extruded spline's open domain is unbounded", S,
     "            Self::Extrusion { curve, .. } => (swept(curve), unbounded),\n",
     "            Self::Extrusion { .. } => (unbounded, unbounded),\n"),
    # The span's arc bound (`brep/trim.rs`), which the lens reads.
    ("a span's arc keeps the curve's fastest speed", T,
     "    Some(local.map_or(global, |local| local.min(global)))\n",
     "    Some(global)\n", TRIM),
    ("a span's local arc drops its acceleration term", T,
     "            (speed.length() * h + 0.25 * bounds.second * h * h)\n",
     "            (speed.length() * h)\n", TRIM),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
