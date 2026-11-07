"""Mutation probe for the Rolle chain of contacts along a section (#397, #427).

A traced section against a surface that touches it has a root of multiplicity
2 (a tangent plane) or 4 (an osculating cylinder); moved a little it has two
simple roots or none. The roots are isolated by a Rolle chain on the
derivative tower along the curve (`root_chain.rs`, `curve_tower.rs`). Each
mutant breaks one step of the chain or the tower and must turn a test red:
`double_root` (the traced half of #381 against planes and cylinders),
`root_chain` (the chain against polynomials) or `curve_tower` (the tower
against closed forms).
"""
import probe_lib

R = "crates/algorithms/parametric/nurbs/src/root_chain.rs"
T = "crates/algorithms/parametric/nurbs/src/curve_tower.rs"
DOUBLE = ["-p", "axiolid-nurbs", "--test", "double_root"]
CHAIN = ["-p", "axiolid-nurbs", "--lib", "root_chain"]
TOWER = ["-p", "axiolid-nurbs", "--lib", "curve_tower"]

MUTANTS = [
    ("the order cap lowered to 2", R,
     "pub(crate) const MAX_CONTACT_ORDER: usize = 8;",
     "pub(crate) const MAX_CONTACT_ORDER: usize = 2;",
     DOUBLE),
    ("a cluster's multiplicity read as 1", R,
     "        let m = run.iter().map(|q| q.m).sum::<usize>() + 1;",
     "        let m = 1;",
     CHAIN),
    ("an inner end allowed to sit on a root", R,
     "                if run.iter().any(|q| q.end == Some(End::Inner)) {",
     "                if false {",
     CHAIN),
    ("the depth floor returning no roots", R,
     "    if depth >= MAX_DEPTH || t1 - t0 <= MIN_WIDTH {\n"
     "        return Err(RootRefusal::Unresolved);",
     "    if depth >= MAX_DEPTH || t1 - t0 <= MIN_WIDTH {\n"
     "        return Ok(());",
     CHAIN),
    ("the node budget never spent", R,
     "    budget.spend()?;",
     "    let _ = &budget;",
     CHAIN),
    ("a split point that need not be generic", R,
     "        if generic(tower, m) {",
     "        if true {",
     CHAIN),
    ("the parity of a cluster unchecked", R,
     "                    if change != (cluster.m % 2 == 1) {",
     "                    if false {",
     CHAIN),
    ("the sign of the tower flipped on the u axis", T,
     "            Axis::U => -1.0,",
     "            Axis::U => 1.0,",
     TOWER),
    ("the derivative factor of a Taylor coefficient dropped", T,
     ".scale((m - i + 1) as Scalar)",
     ".scale((m - i) as Scalar)",
     TOWER),
    ("no order above 1 ever certifying a stretch", R,
     "(2..=MAX_CONTACT_ORDER).find(",
     "(2..2).find(",
     DOUBLE),
    ("a straddle read within a nanometre of zero", R,
     "        let straddles: Vec<bool> = e.iter().map(|x| !x.excludes_zero()).collect();",
     "        let straddles: Vec<bool> = e\n"
     "            .iter()\n"
     "            .map(|x| !x.excludes_zero() || x.lo.abs().min(x.hi.abs()) <= 1e-9)\n"
     "            .collect();",
     DOUBLE),
]

probe_lib.main(MUTANTS, tests=DOUBLE, timeout=1800, reject_uncompiled=True)
