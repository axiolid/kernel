"""Mutation probe for the column builder (#120, ADR 0072).

Covers the arrangement's membership labels, the column builder's wall,
cap, height and shell decisions, and the boolean/clip adapters. Every
mutant runs the column-level test targets that own its claim.

A mutant counts as killed only if some test FAILS; a build error also
counts (the tree must not compile a wrong variant into a passing suite).

Dropped after the first run: "ledge named after the wrong end" mutated the
subject's opposite-cap naming branch, which no boolean can reach (see the
comment on `cap_name` in boolean_column.rs); the branch was removed rather
than kept untested.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
COL = "crates/algorithms/construction/construct/src/column.rs"
BC = "crates/algorithms/construction/construct/src/boolean_column.rs"
ARR = "crates/algorithms/planar/overlay/src/exact_arc/arrangement.rs"
HOLD = "crates/algorithms/planar/overlay/src/arrangement.rs"
MEAS = "crates/algorithms/query/measure/src/exact.rs"
TESTS = [
    "-p", "axiolid-construct",
    "--test", "column_booleans",
    "--test", "boolean_exact",
    "--test", "boolean_arc",
    "--test", "boolean_disconnected",
    "--test", "clip_arc_prism",
    "--test", "boolean_stepped",
    "--test", "exact_measure",
]
OT = ["-p", "axiolid-overlay", "--test", "arrangement"]

MUTANTS = [
    # Builder: walls.
    ("wall where both sides agree", COL,
     "            if l == r {\n                gap += 1;\n                continue;\n            }",
     "            if l == r && gap > 0 {\n                gap += 1;\n                continue;\n            }", TESTS),
    ("wall faces the solid side", COL,
     "        let mut uses = vec![(rim(wall.low), wall.left)];",
     "        let mut uses = vec![(rim(wall.low), !wall.left)];", TESTS),
    ("wall run not split on side change", COL,
     "                && d.solid(&d.left_stack, gap) == l\n                && d.solid(&d.right_stack, gap) == r",
     "                && d.solid(&d.left_stack, gap) != d.solid(&d.right_stack, gap)", TESTS),
    # Builder: caps.
    ("cap on the wrong end of a block", COL,
     ".any(|&(lo, hi)| if up { hi == plane } else { lo == plane })",
     ".any(|&(lo, hi)| if up { lo == plane } else { hi == plane })", TESTS),
    ("cap facing inverted", COL,
     "                let orientation = if *up {\n                    Orientation::Forward\n                } else {\n                    Orientation::Reversed\n                };",
     "                let orientation = if *up {\n                    Orientation::Reversed\n                } else {\n                    Orientation::Forward\n                };", TESTS),
    ("cap holes dropped", COL,
     "                let mut plan = vec![region.outer];\n                plan.extend(region.holes);",
     "                let plan = vec![region.outer];", TESTS),
    # Builder: shells.
    ("cavity taken as an outer shell", COL,
     "        if volume > 0.0 {", "        if volume != 0.0 {", TESTS),
    ("non-manifold edge accepted", COL,
     "            2 => sets.join(faces[0], faces[1]),",
     "            2 | 3 | 4 => sets.join(faces[0], faces[1]),", TESTS),
    # Adapters.
    ("difference keeps the tool", BC,
     "            _ => a && !b,", "            _ => a && b,", TESTS),
    ("operand span open at the top", BC,
     "        self.bottom <= gap && gap < self.top",
     "        self.bottom <= gap && gap <= self.top", TESTS),
    ("stepped difference with covering tool", BC,
     "            !(tool.0 <= subject.0 + tolerance.linear() && tool.1 >= subject.1 - tolerance.linear())",
     "            !(tool.0 <= subject.0 + tolerance.linear() || tool.1 >= subject.1 - tolerance.linear())", TESTS),
    ("wall ordinal ignores earlier rings", BC,
     "            let ordinal = u32::try_from(before + source.edge).ok()?;",
     "            let ordinal = u32::try_from(source.edge).ok()?;", TESTS),
    ("clip keeps the cut side", BC,
     "                vec![(if under_bottom { flat_bottom } else { cut }, flat_top)]",
     "                vec![(if under_bottom { cut } else { flat_bottom }, flat_top)]", TESTS),
    # Arrangement labels.
    # The side a carrier holds is read in `ArrangementEdge::holding` since
    # the sparse labels of #292 replaced `RawEdge::sides`.
    ("shared piece on the wrong side", HOLD,
     "            if on_left == left {",
     "            if on_left != left {", OT),
    ("containment ignored", ARR,
     "                    .filter(|&other| windings[other].winding(&piece.sample, &mut scratch) != 0)",
     "                    .filter(|&other| windings[other].winding(&piece.sample, &mut scratch) != 0 && false)", OT),
    # Measure fix.
    ("measure ignores face orientation", MEAS,
     "            ^ (face.orientation == Orientation::Reversed);",
     "            ^ false;", TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
