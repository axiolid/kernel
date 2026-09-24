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
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
COL = "crates/algorithms/construction/construct/src/column.rs"
BC = "crates/algorithms/construction/construct/src/boolean_column.rs"
ARR = "crates/algorithms/planar/overlay/src/exact_arc/arrangement.rs"
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
    ("shared piece on the wrong side", ARR,
     "                        if same {\n                            left[other] = true;\n                        } else {\n                            right[other] = true;\n                        }",
     "                        if same {\n                            right[other] = true;\n                        } else {\n                            left[other] = true;\n                        }", OT),
    ("containment ignored", ARR,
     "                    } else if winding(&piece.sample, &parts[other]) != 0 {",
     "                    } else if winding(&piece.sample, &parts[other]) != 0 && false {", OT),
    # Measure fix.
    ("measure ignores face orientation", MEAS,
     "            ^ (face.orientation == Orientation::Reversed);",
     "            ^ false;", TESTS),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode


only = set(sys.argv[1:])
survivors = []
for name, rel, old, new, target in MUTANTS:
    if only and name not in only:
        continue
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(target)
        except subprocess.TimeoutExpired:
            code = -1
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)

print(f"{len(MUTANTS) - len(survivors) if not only else '-'} / {len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
