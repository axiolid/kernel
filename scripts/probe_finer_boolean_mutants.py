"""Mutation probe for curved walls with abutting openings and a box tangent
to a curved wall's inner face (#332, #327).

A curved wall less a row of wedge openings, abutting ones reaching 15 um or
0.6 mm past each other's sides, is cut again at the finer tolerance where
the caller's readings are refused; a box whose back face touches a quarter
ring's inner cylinder along a ruling is classified with the ray's start
read exactly, at the caller's tolerance. Each mutant breaks one step that
makes them exact -- the retry, its tolerance, its fresh report, the report's
flag and its merge, the exact pass and its rounding floor, a split's
vertices welded only within the tolerance on the surface -- and must turn a
test red.

Equivalent mutant, deliberately not listed: keeping only the rounding floor
for the faces read as the point's own support in the exact pass
(`Start::Exact(own) if own.contains(&face)`). A point on a region of such a
face is settled by `on_face` before the ray is cast where it lies in the
other face, and outside it the hit there is outside the face either way.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
L = "crates/algorithms/construction/brep-boolean/src/lib.rs"
R = "crates/algorithms/construction/brep-boolean/src/report.rs"
C = "crates/algorithms/construction/brep-boolean/src/classify.rs"
S = "crates/algorithms/construction/brep-boolean/src/split.rs"
BREP = ["-p", "axiolid-brep-boolean", "--test", "curved_wall_openings"]
FILLET = ["-p", "axiolid-brep-boolean", "--test", "fillet_tangent"]
MESH = ["-p", "axiolid-mesh-compile", "--test", "curved_wall_openings"]
TESTS = [BREP, MESH, FILLET]

MUTANTS = [
    ("no retry at the finer tolerance", L,
     "                && session.retry_at(finer) =>",
     "                && false && session.retry_at(finer) =>"),
    ("the retry at the caller's tolerance", L,
     "Tolerance::new(tolerance.linear() * FINER_FACTOR, tolerance.angular())",
     "Tolerance::new(tolerance.linear(), tolerance.angular())"),
    ("the retry keeps the refused attempt's readings", R,
     "                s.tolerance = tolerance;\n"
     "                s.report = BooleanReport::default();",
     "                s.tolerance = tolerance;"),
    ("a finer result not flagged", R,
     "                s.finer = true;",
     "                s.finer = false;"),
    ("the flag lost when reports merge", R,
     "        self.finer |= other.finer;\n",
     ""),
    ("no exact pass over a region's points", L,
     "        classify::Start::Exact(coincident),",
     "        classify::Start::Tolerance,"),
    ("the exact pass reads the start within the tolerance", C,
     "                                    <= crate::report::rounding()",
     "                                    <= tolerance.linear()"),
    ("a split's vertices welded by parameters alone", S,
     "(*q - p).length() <= slack && (*b - at).length() <= reach",
     "(*q - p).length() <= slack"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
