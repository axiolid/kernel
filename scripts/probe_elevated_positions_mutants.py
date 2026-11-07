"""Mutation probe for the certified positions of elevated and banked curves
over intrinsic plans (#393, ADR 0082 amendment).

Each mutant breaks one step the bounds rest on. The table's point: the
height dropped, the pivot dropped. Its uniform error: the remainders
dropped. The chord bound of the evaluated ends: the gap of the ends dropped,
the smaller of the two ends' gaps taken, the zero-width span answered with
zero. The tube: its centre read from the evaluator instead of the table,
the named refusal of a chain plan and of an intrinsic profile removed. Each
must turn `evaluate/tests/elevated_positions.rs`,
`evaluate/tests/intrinsic_position_table.rs`, the tube wall's unit test
(`deviation::paths::tests`) or `compile/tests/elevated_directrix.rs` red.

Not probed, because no test can tell them apart:

- the tube's `+ e`, the table's error added to the certified distance: it
  is a thousandth of the budget, below any sampling of the tube;
- the closed-form height's and pivot's own rounding, which the table does
  not charge (about `1e-14`, as for every closed-form point in the tube
  certificate).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
EVAL = "crates/algorithms/parametric/evaluate"
TABLE = f"{EVAL}/src/bound/position_table.rs"
ELEVATED = f"{EVAL}/src/bound/elevated_positions.rs"
BOUND = f"{EVAL}/src/bound.rs"
PATHS = "crates/execution/compile/src/deviation/paths.rs"

POSITIONS = ["-p", "axiolid-evaluate", "--all-features", "--test", "elevated_positions"]
TABLES = ["-p", "axiolid-evaluate", "--all-features", "--test", "intrinsic_position_table"]
WALL = ["-p", "axiolid-mesh-compile", "--lib", "a_tube_wall_over_a_clothoid"]
NAMED = [
    "--release",
    "-p",
    "axiolid-mesh-compile",
    "--test",
    "elevated_directrix",
    "a_disk_along_a_chain_plan",
]

MUTANTS = [
    # The table's point.
    ("the height is dropped from the point", ELEVATED,
     "        let z = height + pivot;",
     "        let z = pivot;", [POSITIONS]),
    ("the pivot is dropped from the point", ELEVATED,
     "        let z = height + pivot;",
     "        let z = height;", [POSITIONS]),
    # The uniform error.
    ("the uniform error drops the remainders", TABLE,
     "            self.total_remainder + self.end * (self.heading_rounding + 16.0 * Scalar::EPSILON),",
     "            0.0 * self.total_remainder\n"
     "                + self.end * (self.heading_rounding + 16.0 * Scalar::EPSILON),",
     [TABLES]),
    # The chord bound of the evaluated ends.
    ("the chord bound drops the ends' gap", ELEVATED,
     "    let out = inflate(Iv::point(law).add(Iv::point(g)).hi);",
     "    let out = inflate(Iv::point(law).add(Iv::point(0.0 * g)).hi);", [POSITIONS]),
    ("the chord bound takes the smaller end gap", ELEVATED,
     "    let g = gap(lo)?.max(gap(hi)?);",
     "    let g = gap(lo)?.min(gap(hi)?);", [POSITIONS]),
    ("a zero-width span has no gap", BOUND,
     "        if over_intrinsic_plan(curve) {",
     "        if false && over_intrinsic_plan(curve) {", [POSITIONS]),
    # The tube.
    ("the tube reads the evaluator's centre", PATHS,
     "        Some(table) => table.point3(s).map(|(point, _)| point),",
     "        Some(_) => axiolid_reference::curve::evaluate3(curve, s).ok(),", [WALL]),
    ("a chain plan is not refused by name", PATHS,
     "    (matches!(base.plan.as_ref(), Curve2::Chain(_)) || has_intrinsic_profile(&base.elevation))",
     "    (has_intrinsic_profile(&base.elevation))", [NAMED]),
    ("an intrinsic profile is not refused by name", PATHS,
     "    (matches!(base.plan.as_ref(), Curve2::Chain(_)) || has_intrinsic_profile(&base.elevation))",
     "    (matches!(base.plan.as_ref(), Curve2::Chain(_)))", [NAMED]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
