"""Mutation probe for edges without pcurves on surface relations (#454, ADR 0110).

Each mutant weakens one decision of `crates/execution/compile/src/
surface_relation/invert.rs` (the cover and its certificate, uniqueness, the
fold and pole refusals, the closed-form other parameter, the revolved line's
bracket), the tangents and open domain a relation hands the lens, or the
inverter's periods. The module's unit tests, `tests/relation_edge_inversion.rs`,
`tests/relation_faces.rs` or the tessellator's own tests must fail on it. A
build error also counts as killed.

Not listed, as not observable (no fixture, not equivalent):
- R8, the coarse cell without a certified speed bound: the four supported
  swept families always have one.
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
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--lib", "surface_relation"],
    ["-p", "axiolid-mesh-compile", "--test", "relation_edge_inversion"],
    ["-p", "axiolid-mesh-compile", "--test", "relation_faces"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::"],
]

MUTANTS = [
    # The cover.
    ("the cover never drops a cell", V,
     "            if cell.lower().is_nan() || cell.lower() > self.tau.min(best + slack) {\n",
     "            if false {\n"),
    ("the Lipschitz term is zero", V,
     "        0.5 * self.speed * (self.b - self.a)\n",
     "        0.0 * self.speed * (self.b - self.a)\n"),
    ("the drop bound loses tau", V,
     [("            if cell.lower().is_nan() || cell.lower() > self.tau.min(best + slack) {\n",
       "            if cell.lower().is_nan() || cell.lower() > best + slack {\n"),
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
    ("a fold's wide cover is accepted", V,
     "        if (piece.1 - piece.0) * speed > self.tau {\n",
     "        if false {\n"),
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
     "        let off = (self.surface.evaluate(uv.x, uv.y)? - point).length();\n"
     "        if off.is_nan() || off > self.tau {\n",
     "        let off = (self.surface.evaluate(uv.x, uv.y)? - point).length();\n"
     "        if false && off > self.tau {\n"),
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
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
