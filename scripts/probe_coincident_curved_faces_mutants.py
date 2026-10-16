"""Mutation probe for a cutter face on a curved face of its host (#324).

A quarter-ring wall less a niche whose inner arc lies on the wall's inner
cylinder, or one unit in the last place inside it. Each mutant breaks one
step that makes it certified -- the other operand's coincident face read
with its winding, the winding itself, the reading reported as a deviation
contribution -- and must turn a test red.

Equivalent mutants, deliberately not listed: dropping the chord that
closes a gap at a pole in `face_against` (it changes the winding's size,
not its sign, for every face these tests build: a pole gap only closes a
loop the seam already leaves open), and dropping the
`(Orientation::Reversed, true)` arm (no operand here has a reversed face
whose loops wind clockwise).
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/construction/brep-boolean/src/classify.rs"
S = "crates/algorithms/construction/brep-boolean/src/split.rs"
D = "crates/execution/compile/src/deviation/boolean.rs"
BREP = ["-p", "axiolid-brep-boolean", "--test", "coincident_cylinders"]
MESH = ["-p", "axiolid-mesh-compile", "--test", "curved_flush_niche"]
TESTS = [BREP, MESH]

MUTANTS = [
    ("a coincident face read without its winding", C,
     "                        self.against[face],",
     "                        false,"),
    ("the winding read the other way", S,
     "    Ok(swept < 0.0)\n}\n\n/// `1/2 int",
     "    Ok(swept > 0.0)\n}\n\n/// `1/2 int"),
    ("the coincidence not reported", D,
     "        .find(|d| d.kind == ToleranceDecisionKind::CoincidentSupports)",
     "        .find(|d| d.kind == ToleranceDecisionKind::MergedPoints)",
     MESH),
    ("the coincidence reported as no move", D,
     "            DeviationBound::Certified(moved),",
     "            DeviationBound::Certified(0.0 * moved),",
     MESH),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
