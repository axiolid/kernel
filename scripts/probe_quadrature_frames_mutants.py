"""Mutation probe for certified frame-law bounds of frames read by
quadrature (#316, ADR 0082 amendment).

Each mutant breaks one step the composed bounds rest on. The frame law in
the exact curve's arc length: the curvature bound short a power of the
speed floor, the tangent's second derivative without the unit vector's, a
3D curve's or a vertically varying offset's frame taken as planar. The
reads: their certified true distance dropped from the measure, a reversed
piece or an offset beside a reversed base measured from the wrong end, the
corner rule never applied, the read point's and the read axes' computed
distances from the exact frame dropped. The composition: the joint's
certified error not added on a composite, a piece read by quadrature
handed to the closed-form bounds, and the run's frames certified to the
whole budget rather than a thousandth of it. Each must turn
`evaluate/tests/station_frame_law.rs`, the crate's own unit tests
(`reading::tests`) or `compile/tests/station_bound_quadrature.rs` red.

Not probed, because no test can tell them apart: the read errors' share of
a cell's slide in `station/bound.rs` and the turned axes' `sqrt 3` lever
(the errors are about `1e-9`, five orders below any budget the tests ask),
the tangent's radius from the velocity's slip and the plan reading's
widening of the axis error (the same scale), and the speed floor's span
term (the bounds stay above the differences without it on every test
span; it is what keeps them sound on a long one).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
READING = "crates/algorithms/parametric/evaluate/src/station/certify/reading.rs"
FB = "crates/algorithms/parametric/evaluate/src/station/frame_bound.rs"
SB = "crates/execution/compile/src/station/bound.rs"

LAW = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_frame_law"]
LIB = ["-p", "axiolid-evaluate", "--all-features", "--lib", "reading"]
QUAD = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_bound_quadrature"]

MUTANTS = [
    # The frame law in the exact curve's arc length.
    ("the curvature bound is short a power of the speed", READING,
     "            let k1 = second / (sigma * sigma);",
     "            let k1 = second / sigma;", [LAW]),
    ("the tangent's second derivative drops the unit vector's", READING,
     "            let k2 = unit_second / (sigma * sigma) + x * x / (sigma * sigma);",
     "            let k2 = 0.0 * unit_second / (sigma * sigma) + x * x / (sigma * sigma);",
     [LAW]),
    ("every curve read by quadrature is framed as planar", READING,
     "                    planar: matches!(curve, StationCurve::Two(_)),",
     "                    planar: true,", [LAW]),
    ("an offset whose vertical law varies is framed as planar", READING,
     "                    OffsetLaw::Linear { start, end, .. } => start.vertical == end.vertical,",
     "                    OffsetLaw::Linear { .. } => true,", [LAW]),
    # The reads.
    ("the measure drops the certified true distance", READING,
     "        let measure = (u - into_lo).max(into_hi - u).max(0.0);",
     "        let measure = 0.0 * (u - into_lo).max(into_hi - u).max(0.0);", [LAW]),
    ("a reversed piece is measured from its far end", READING,
     "                    increasing: !piece.reversed,",
     "                    increasing: true,", [LAW]),
    ("an offset beside a reversed base is measured from its far end", READING,
     "                    increasing: base.reversed == piece.reversed,",
     "                    increasing: !piece.reversed,", [LAW]),
    ("the corner rule never applies", READING,
     "        let corner =\n            self.corners.contains(&t)",
     "        let corner =\n            false && self.corners.contains(&t)", [LAW]),
    ("the read point's distance from the exact one is dropped", READING,
     "        let point_error = (frame.point - point).length() + error + POINT_ROUNDING * scale;",
     "        let point_error = error + POINT_ROUNDING * scale;", [LIB]),
    ("the read tangent's distance from the exact one is dropped", READING,
     "            ((frame.tangent - exact.tangent).length() + r_t)",
     "            (r_t)", [LIB]),
    # The composition.
    ("a joint's certified error is not added", FB,
     "                        bounds.measure += self.joint_error(composite, k)?;",
     "                        bounds.measure += 0.0 * self.joint_error(composite, k)?;", [LAW]),
    ("a piece read by quadrature takes the closed-form bounds", FB,
     "            if !by_quadrature(piece.curve) {",
     "            if true || !by_quadrature(piece.curve) {", [LAW, QUAD]),
    ("the run's frames are certified to the whole budget", SB,
     "const CERTIFICATION_SHARE: Scalar = 1e-3;",
     "const CERTIFICATION_SHARE: Scalar = 1e3;", [QUAD]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False)
