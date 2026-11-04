"""Mutation probe for the Rolle certificate of a touching root (#397).

A traced section against a plane tangent to it has one root of multiplicity
2; lowered it has two simple roots, raised none. Each mutant breaks one step
of the certificate and must turn a test red: the first the boolean of #381 at
`Tolerance::ZERO` (whose double root only the certificate settles), the others
`double_root`.
"""
import probe_lib

N = "crates/algorithms/parametric/nurbs/src/implicit_ops.rs"
TESTS = ["-p", "axiolid-nurbs", "--test", "double_root"]
BOOLEAN = ["-p", "axiolid-brep-boolean", "--test", "curved", "a_box_less_a_torus"]

MUTANTS = [
    ("the second-derivative exclusion forced false", N,
     "if !second.lo.is_finite() || !second.hi.is_finite() || second.straddles_zero() {",
     "if true || second.straddles_zero() {",
     BOOLEAN),
    ("the critical point sought on the wrong sign", N,
     "            if (d > 0.0) == rising {",
     "            if (d > 0.0) != rising {",
     TESTS),
    ("a touch emitting sign changes", N,
     "            out.push((sc, 2));\n            return true;",
     "            out.push((self.bisect(s0, sc, a), 1));\n"
     "            out.push((self.bisect(sc, s1, -a), 1));\n            return true;",
     TESTS),
    ("an extremum near zero read as a touch", N,
     "        if enclosure.straddles_zero() {",
     "        if enclosure.straddles_zero() || hc.abs() <= 1e-9 {",
     TESTS),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
