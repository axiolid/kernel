"""Mutation probe for exact curve intersection (#119): each mutant must be caught."""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/parametric/nurbs/src/exact_curve_intersection.rs"

MUTANTS = [
    ("multiplicity stops at one",
     "    while !q.is_zero() && root.sign_of(&q) == Sign::Zero {",
     "    while false && root.sign_of(&q) == Sign::Zero {"),
    ("cone nappe condition ignored",
     "        if !allowed(&root) {\n            continue;\n        }",
     "        if false && !allowed(&root) {\n            continue;\n        }"),
    ("antipode never reported",
     "        if lost >= 1 && on_nappe {",
     "        if lost >= 2 && on_nappe {"),
    ("common roots from first equation only",
     "        .fold(polys[0].0.clone(), |acc, (p, _)| acc.gcd(p));",
     "        .fold(polys[0].0.clone(), |acc, (_p, _)| acc);"),
    ("partial overlap accepted",
     "        .any(|root| multiplicity(root, &int) % 2 == 1);",
     "        .any(|root| multiplicity(root, &int) % 2 == 7);"),
    ("theta order ignores the lower half",
     "        let band = if is_line || root.cmp_dyadic(&Dyadic::zero()) != Sign::Negative {",
     "        let band = if true || root.cmp_dyadic(&Dyadic::zero()) != Sign::Negative {"),
]

probe_lib.main(MUTANTS, tests=[["--release", "-p", "axiolid-nurbs",
         "--test", "exact_curve_intersection"]], path=SRC, timeout=900)
