"""Mutation probe for the exact fold decision of turned sections and frames
read by quadrature (#406, ADR 0082 amendment).

Each mutant breaks one step the decision rests on.

The frame law's rates (`FrameLaw::rates`): the twist's reference-up term
dropped, the roll's rate dropped from the twist (and, under the vertical-rise
convention, its `hypot(G, 1)` factor), the plan reading's rise dropped, the
section reading given a rise, the twist's Lipschitz term dropped, a span
across a curvature break of a frame read by quadrature given rates, and
`has_bend` blind to frames read by quadrature.

The frames read by quadrature (`PieceReader::jet`, `Measured::acceleration`):
the curvature divided by the speed instead of its square, a reversed base's
direction of travel dropped, the native acceleration read as the velocity, an
ellipse's third derivative with the wrong sign, a B-spline's third derivative
read as its second, an offset's `2 c' a1` term and `unit_accel`'s `2 u1 r1`
dropped.

The decision (`SectionLaw::decide`, `floors`): the orientation's up turn and
lateral turn dropped, the profile point dropped from the vertex, the rise term
dropped, the turned components not widened over the span, the turned frame
replaced by the base's, the depth-limit fallback reporting a span unbounded
without trying the conservative bound, and the bound of the turned tangent's
turn dropped.

Each must turn `evaluate/tests/station_frame_rates.rs`, a unit test,
`compile/tests/station_fold_turned.rs` or `station_fold_quadrature.rs` red.
A mutant that does not compile counts as BROKEN, not killed.

Not probed, because no test can tell them apart:

- the witness window's widening by the measure (about `1e-8`, below any
  distance the tests put between a vertex and a fold);
- the `64 eps` rounding allowances, the acceleration's and the jet's error
  terms, the measure widening of the Lipschitz span and the knot term of the
  native acceleration's error (all far below the distances the tests put
  between a vertex and a fold, as `probe_station_fold_mutants.py` states of
  its own rounding terms).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
EVAL = "crates/algorithms/parametric/evaluate"
RATES = f"{EVAL}/src/station/frame_bound/rates.rs"
BEND = f"{EVAL}/src/station/frame_bound/bend.rs"
BANKED = f"{EVAL}/src/banked.rs"
READING = f"{EVAL}/src/station/certify/reading.rs"
CERTIFY = f"{EVAL}/src/station/certify.rs"
BESIDE = f"{EVAL}/src/station/certify/beside.rs"
PANELS = f"{EVAL}/src/station/certify/panels.rs"
CURVE = f"{EVAL}/src/curve.rs"
COMPILE = "crates/execution/compile"
FOLD = f"{COMPILE}/src/station/bound/fold.rs"
RATE = f"{COMPILE}/src/station/bound/fold/rate.rs"
TURNED = f"{COMPILE}/src/station/bound/fold/turned.rs"
BOUND = f"{COMPILE}/src/station/bound.rs"

FRAME_RATES = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_frame_rates"]
STATION_BEND = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_bend"]
EVAL_LIB = ["-p", "axiolid-evaluate", "--all-features", "--lib"]
FOLD_TURNED = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_fold_turned"]
FOLD_QUADRATURE = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_fold_quadrature"]
TURNED_UNIT = ["-p", "axiolid-mesh-compile", "--all-features", "--lib", "turned"]

MUTANTS = [
    # The frame law's rates.
    ("the twist's reference-up term is dropped", RATES,
     "            let middle = phi + dt.dot(up0.cross(Vec3::Z)) / rho;\n",
     "            let middle = phi + 0.0 * dt.dot(up0.cross(Vec3::Z)) / rho;\n",
     [FRAME_RATES, FOLD_TURNED]),
    ("the roll's rate is dropped from the twist", RATES,
     "            let (phi, phi_error) = roll?;\n",
     "            let (phi, phi_error) = (0.0 * roll?.0, roll?.1);\n",
     [FRAME_RATES, FOLD_TURNED]),
    ("the vertical-rise factor is dropped from the roll's rate", BANKED,
     "            let scale = v.z.hypot(1.0);\n",
     "            let scale = 1.0;\n",
     [EVAL_LIB, FRAME_RATES]),
    ("the plan reading's rise is dropped", RATES,
     "                v.dot(frame.up),\n",
     "                0.0 * v.dot(frame.up),\n",
     [FRAME_RATES, FOLD_TURNED]),
    ("the section reading is given a rise", RATES,
     "        let rise = if plan {\n",
     "        let rise = if true {\n",
     [FRAME_RATES]),
    ("the twist's Lipschitz term is dropped", RATES,
     "            let slope = spanned.up[0] * spanned.lateral[0] + spanned.lateral[1];\n",
     "            let slope = 0.0 * (spanned.up[0] * spanned.lateral[0] + spanned.lateral[1]);\n",
     [FRAME_RATES]),
    ("a span across a curvature break is given rates", RATES,
     "                if reader.bend_breaks().iter().any(|&k| lo < k && k < hi) {\n",
     "                if reader.bend_breaks().iter().any(|&k| lo < k && k < hi && false) {\n",
     [FRAME_RATES]),
    ("has_bend is blind to frames read by quadrature", BEND,
     "            Some(Slot::Read(reader)) => matches!(reader.jet(&at.piece, mid, at.table), Ok(Some(_))),\n",
     "            Some(Slot::Read(_)) => false,\n",
     [STATION_BEND, FOLD_QUADRATURE]),
    # The frames read by quadrature.
    ("the curvature is divided by the speed, not its square", READING,
     "                curvature: (acceleration - unit * acceleration.dot(unit)) / (speed * speed),\n",
     "                curvature: (acceleration - unit * acceleration.dot(unit)) / speed,\n",
     [FRAME_RATES, FOLD_QUADRATURE]),
    ("a reversed base's direction of travel is dropped", READING,
     "                tangent: self.travel() * unit,\n",
     "                tangent: unit,\n",
     [FRAME_RATES]),
    ("the native acceleration is read as the velocity", CERTIFY,
     "            let a = self.second_derivative(t)?;\n            let d = self.bounds(t, t)?;\n",
     "            let a = self.derivative(t)?;\n            let d = self.bounds(t, t)?;\n",
     [EVAL_LIB, FRAME_RATES]),
    ("an ellipse's third derivative has the wrong sign", CURVE,
     "        Curve2::Ellipse(e) => Ok(-conic_tangent2(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),\n",
     "        Curve2::Ellipse(e) => Ok(conic_tangent2(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),\n",
     [EVAL_LIB, FRAME_RATES]),
    ("a B-spline's third derivative is read as its second", CURVE,
     "        (a3[axis] - 3.0 * w1 * c2[axis] - 3.0 * w2 * c1[axis] - w3 * c0[axis]) / w0\n",
     "        c2[axis]\n",
     [EVAL_LIB, FRAME_RATES]),
    ("an offset's 2 c' a1 term is dropped", BESIDE,
     "2.0 * q * rate * a1,",
     "0.0 * q * rate * a1,",
     [FRAME_RATES]),
    ("unit_accel drops the 2 u1 r1 term", PANELS,
     "    let u2 = (w2 - 2.0 * u1 * r1 - u * r2) / r;\n",
     "    let u2 = (w2 - u * r2) / r;\n",
     [EVAL_LIB, FRAME_RATES]),
    # The decision.
    ("the orientation's up turn is dropped", TURNED,
     "    let up_rate = axis_rate.dot(normal) / size;\n",
     "    let up_rate = 0.0 * axis_rate.dot(normal) / size;\n",
     [TURNED_UNIT, FOLD_TURNED]),
    ("the orientation's lateral turn is dropped", TURNED,
     "        up,\n        lateral_rate,\n        up_rate,\n",
     "        up,\n        lateral_rate: 0.0,\n        up_rate,\n",
     [TURNED_UNIT, FOLD_TURNED]),
    ("the profile point is dropped from the vertex", RATE,
     "                    [0, 1, 2].map(|k| q[k] + x * frame.lateral[k] + y * frame.up[k])\n",
     "                    [0, 1, 2].map(|k| q[k] + 0.0 * x * frame.lateral[k] + 0.0 * y * frame.up[k])\n",
     [FOLD_TURNED]),
    ("the rise term is dropped from the rate", RATE,
     "        sum.term(1.0, &[c.rise, n[2]]);\n",
     "",
     [FOLD_TURNED]),
    ("the turned components are not widened over the span", RATE,
     "            let radius = half * rate + ROUNDING;\n",
     "            let radius = 0.0 * half * rate + ROUNDING;\n",
     [FOLD_TURNED]),
    ("the turned frame is replaced by the base's", RATE,
     "            turned_frame(ta, tb, 0.5 * (u0 + u1), delta)?\n",
     "            TurnedFrame::BASE\n",
     [FOLD_TURNED]),
    ("the depth-limit fallback gives up without the conservative bound", FOLD,
     "                    if !matches!(self.advancing(part, reading, margin, identity)?, Ok(true)) {\n"
     "                        verdict = Err(FOLD_NOT_EXCLUDED);\n"
     "                    }\n",
     "                    let _ = part;\n                    verdict = Err(FOLD_NOT_EXCLUDED);\n",
     [FOLD_QUADRATURE]),
    ("the turned tangent's turn bound is dropped", BOUND,
     "        tangent: [x1, x2],\n",
     "        tangent: [0.0, 0.0],\n",
     [FOLD_TURNED]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False, reject_uncompiled=True)
