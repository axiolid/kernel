"""Mutation probe for the seam planner's reading of a pole or apex (#379).

A loop that reaches a pole or apex jumps in u there, and the jump is not a
turn round the surface. Each mutant weakens one decision of
`crates/algorithms/construction/brep-boolean/src/seams.rs` (the joint rule,
the closing rule, the pole test, the whole-period snap) and must turn a test
of `poles.rs`, `quadrics.rs` or `seamless.rs` red. A build error counts as
killed.
"""
import probe_lib

S = "crates/algorithms/construction/brep-boolean/src/seams.rs"
TESTS = [
    "-p", "axiolid-brep-boolean",
    "--test", "poles", "--test", "quadrics", "--test", "seamless",
]

MUTANTS = [
    ("joint pole not recognised", S,
     "        let u = if at_pole(surface, end) {\n            0.0",
     "        let u = if false {\n            0.0"),
    ("closing pole not recognised", S,
     "    let u = if at_pole(surface, last) {\n        0.0",
     "    let u = if false {\n        0.0"),
    ("pole test too loose", S,
     "su.length() <= 1e-9 * (1.0 + sv.length())",
     "su.length() <= 1e9 * (1.0 + sv.length())"),
    ("u snap dropped everywhere", S,
     "            (x / TAU).round() * TAU",
     "            0.0 * (x / TAU).round() * TAU"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1200)
