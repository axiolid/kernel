"""Mutation probe for certified lengths of offsets by distances along an
ellipse or a B-spline, and adaptive certified panels (#307, ADR 0082
amendment).

Each mutant breaks one step the new certificate rests on. The law read
through the base's certified length: the varying law refused again, read
as constant, read from the wrong end on a reversed base, read through a
coarse table (the panel's ends, not a sub-panel to the point) or a table
far too coarse for the tolerance. The deficit bound: the third-derivative
bound ignored (the order drops back to `h^3`), its lower bound overstated,
its upper bound understated, a conic's fourth derivative taken as zero
and a rational spline's fourth derivative read without its weights. The
adaptivity: the narrowest panel halved first, refinement stopped at once,
and the midpoint rule's panels never refined. Each must turn
`evaluate/tests/station_certify_varying.rs`, `station_certify.rs`, the
crate's own unit tests (`bound::tests`) or
`compile/tests/station_window_varying_offset.rs` red.

Not probed, because no test can tell them apart: a point's error from the
table dropped from the chain's rounding and the velocity's (they are
orders of magnitude below the tolerances the tests ask; the containment
holds without them), the unit vector's third-derivative constant (loose
by design: the supremum is not attained), `G`'s own error term, and the
cuts at third-derivative breaks (no test base has a knot whose third
derivative jumps inside an offset piece). The panel budget itself is not
probed: halving it only refuses earlier.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
BESIDE = "crates/algorithms/parametric/evaluate/src/station/certify/beside.rs"
PANELS = "crates/algorithms/parametric/evaluate/src/station/certify/panels.rs"
CERTIFY = "crates/algorithms/parametric/evaluate/src/station/certify.rs"
BOUND = "crates/algorithms/parametric/evaluate/src/bound.rs"

VARYING = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_certify_varying"]
CERT = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_certify"]
LIB = ["-p", "axiolid-evaluate", "--all-features", "--lib"]
DECK = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_window_varying_offset"]

MUTANTS = [
    # The law read through the base's certified length.
    ("an offset by distances along an ellipse or a B-spline is refused", BESIDE,
     "            let table = BaseTable::build(base_curve, t0, t1, &breaks, target)?;",
     "            let table = BaseTable::build(base_curve, t0, t1, &breaks, target)?;\n"
     "            return Err(uncertified());", [VARYING, DECK]),
    ("the varying law is read as its start", BESIDE,
     "                let middle = 0.5 * (ratio.0 + ratio.1);",
     "                let middle = 0.0 * (ratio.0 + ratio.1);", [VARYING]),
    ("a reversed base reads its law from the wrong end", BESIDE,
     "                    vb - span * middle\n",
     "                    va + span * middle\n", [VARYING]),
    ("the base's length to a point is read from its panel's ends", BESIDE,
     "        let part = base.panel(panel.a, base.node(t)?)?;",
     "        let part = base.panel(panel.a, base.node(t)?)?;\n"
     "        let part = Panel { chord: Scalar::NAN, ..part };", [VARYING]),
    ("the base table is far too coarse for the tolerance", BESIDE,
     "let target = (PANEL_SHARE * tolerance).min(tolerance / (16.0 * rsum * TURNING));",
     "let target = 1e4 * tolerance;", [VARYING]),
    # The deficit bound.
    ("the third-derivative bound is ignored", PANELS,
     "            if let Some(third) = third.filter(|m| m.is_finite()) {",
     "            if let Some(third) = third.filter(|m| m.is_finite() && false) {", [VARYING]),
    ("the deficit's lower bound is overstated", PANELS,
     "                        lower = moment * moment / h3;",
     "                        lower = 4.0 * moment * moment / h3;", [VARYING, CERT]),
    ("the deficit's upper bound is understated", PANELS,
     "                        alpha * alpha * h3",
     "                        0.25 * alpha * alpha * h3", [VARYING, CERT]),
    ("a conic's fourth derivative is taken as zero", BOUND,
     "        Curve2::Ellipse(e) => Some(inflate(conic_second(\n"
     "            e.semi_axis_x,\n"
     "            e.semi_axis_y,\n"
     "            lift2(e.frame.x),\n"
     "            lift2(e.frame.y),\n"
     "            a,\n"
     "            b,\n"
     "        ))),\n"
     "        Curve2::BSpline(spline) => {\n"
     "            bspline_bounds(spline, |p| [p.x, p.y], a, b).map(|(.., fourth)| fourth)",
     "        Curve2::Ellipse(_) => Some(0.0),\n"
     "        Curve2::BSpline(spline) => {\n"
     "            bspline_bounds(spline, |p| [p.x, p.y], a, b).map(|(.., fourth)| fourth)",
     [LIB, VARYING]),
    ("a rational spline's fourth derivative ignores its weights", BOUND,
     "    [first, second, third, fourth]\n}",
     "    [first, second, third, a[4] / w_min + 0.0 * fourth]\n}", [LIB]),
    # The adaptivity.
    ("the narrowest panel is halved first", PANELS,
     "        self.0.total_cmp(&other.0).then(other.1.cmp(&self.1))",
     "        other.0.total_cmp(&self.0).then(other.1.cmp(&self.1))", [VARYING, CERT]),
    ("the panels are never refined", PANELS,
     "        if count + 2 > MAX_CERTIFIED_PANELS {",
     "        if count + 2 > 0 {", [VARYING, CERT]),
    ("the midpoint rule's panels are never refined", CERTIFY,
     "            |p| p.3 - p.2,",
     "            |p| 0.0 * (p.3 - p.2),", [CERT]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False)
