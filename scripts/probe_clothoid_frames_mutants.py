"""Mutation probe for certified frame laws along clothoid and other
intrinsic plans (#361, ADR 0082 amendment).

Each mutant breaks one step the bounds rest on. The position table: the
Gauss-Legendre remainder dropped, the partial panel's remainder dropped, a
checkpoint that forgets the remainder so far. The reader: the read point's
computed distance from the exact point dropped, a banked curve taken as an
elevated one, the orthonormal start frame and the level plane of a space
curve not required. The frame law: the tangent's second derivative without
the curvature's rate. The soundness rule: `piece_frame_bounds`, whose
`point_error` is zero, handed an intrinsic arm. Each must turn
`evaluate/tests/intrinsic_position_table.rs`,
`evaluate/tests/station_frame_law_intrinsic.rs` or the reader's own unit
tests (`intrinsic::tests`) red.

Not probed, because no test can tell them apart:

- the `3 k'^2` term of the position table's fourth-derivative bound: the
  table's step falls with the fourth root of the bound and the rule's real
  error sits well under its worst case, so dropping one of the terms of a
  slowly turning law leaves every point inside its error;
- the table error's share of a cell's slide in `station/bound.rs`, about
  `1e-7` below any budget the tests ask, as in #316;
- the heading's rounding in the axes' charge (about `1e-13`), and the side
  a kink of the profile is read from (the stations land off the kinks);
- the plan reading's widening of the axis error.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
EVAL = "crates/algorithms/parametric/evaluate"
TABLE = f"{EVAL}/src/bound/position_table.rs"
READER = f"{EVAL}/src/station/certify/intrinsic.rs"
FB = f"{EVAL}/src/station/frame_bound.rs"

POSITIONS = ["-p", "axiolid-evaluate", "--all-features", "--test", "intrinsic_position_table"]
LAW = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_frame_law_intrinsic"]
UNIT = ["-p", "axiolid-evaluate", "--all-features", "--lib", "intrinsic"]

BANKED_ARM = (
    "            StationCurve::Three(Curve3::Banked(banked)) => {\n"
    "                let Curve2::Intrinsic(plan) = banked.base.plan.as_ref() else {\n"
    "                    return None;\n"
    "                };\n"
    "                Some(Self::Elevated(&banked.base, plan))\n"
    "            }\n"
    "            StationCurve::Three(Curve3::Intrinsic(space)) => {\n"
)

INTRINSIC_ARM = (
    "        Curve2::Intrinsic(plan) => {\n"
    "            let b = crate::bound::intrinsic_derivative_bounds2(plan, lo, hi)?;\n"
    "            let mid = crate::curve::derivative2(curve, 0.5 * (lo + hi)).ok()?;\n"
    "            Some(TangentLaw {\n"
    "                point: [b.first, b.second],\n"
    "                tangent: [b.second, b.third],\n"
    "                middle: unit(lift(mid))?,\n"
    "                kind: Kind::Planar,\n"
    "                measure: 0.0,\n"
    "            })\n"
    "        }\n"
    "        _ => None,\n"
    "    }\n"
    "}\n"
    "\n"
    "/// An elevated or banked curve's point law in plan distance"
)

MUTANTS = [
    # The position table.
    ("the Gauss-Legendre remainder is dropped", TABLE,
     "    width.powi(5) / 4320.0 * b4",
     "    0.0 * width.powi(5) / 4320.0 * b4", [POSITIONS]),
    ("the partial panel's remainder is dropped", TABLE,
     "            state.remainder += remainder(rest, span.b4);",
     "            state.remainder += 0.0 * remainder(rest, span.b4);", [POSITIONS]),
    ("a checkpoint forgets the remainder so far", TABLE,
     "                table.checkpoints.push(state);",
     "                table.checkpoints.push(Checkpoint {\n"
     "                    remainder: 0.0,\n"
     "                    ..state\n"
     "                });", [POSITIONS]),
    # The reader.
    ("the read point's distance from the exact one is dropped", READER,
     "        let point_error = (frame.point - exact.point).length() + table + POINT_ROUNDING * scale;",
     "        let point_error = table + POINT_ROUNDING * scale;", [UNIT]),
    ("a banked curve is read as an elevated one", READER,
     "            StationCurve::Three(Curve3::Intrinsic(space)) => {\n",
     BANKED_ARM, [LAW]),
    ("the start frame need not be orthonormal", READER,
     "    (frame.x.length() - 1.0).abs() <= NEAR\n",
     "    true || (frame.x.length() - 1.0).abs() <= NEAR\n", [LAW]),
    ("a space curve need not be level", READER,
     "                let level = start.x.z == 0.0 && start.y.z == 0.0;",
     "                let level = true;", [LAW]),
    # The frame law.
    ("the tangent's second derivative drops the curvature's rate", FB,
     "        tangent: [bounds.second, bounds.third],",
     "        tangent: [bounds.second, bounds.second * bounds.second],", [LAW]),
    # The soundness rule.
    ("the closed-form bounds take an intrinsic arm", FB,
     "        _ => None,\n    }\n}\n\n/// An elevated or banked curve's point law in plan distance",
     INTRINSIC_ARM, [LAW]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False)
