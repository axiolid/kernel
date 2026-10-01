"""Mutation probe for axiolid-numeric: each mutant must fail at least one test.

Run: python3 scripts/probe_numeric_mutants.py   (a few minutes; not in gate.sh)

Each entry weakens one property the crate promises: a refusal, a bracket
or enclosure, an error estimate, or a result. A survivor means the tests do
not pin that property; add a test, never delete the mutant. Anchors must
match exactly once, so a refactor that moves code fails loudly here
instead of silently skipping a mutant.
"""
import pathlib, subprocess, sys, os

K = pathlib.Path(__file__).resolve().parents[1]
SRC = K / "crates/algorithms/numeric/src"

MUTANTS = [
    # Refusals.
    ("negative tolerance accepted", "error.rs",
     "    if value < 0.0 {\n        return Err(NumericError::InvalidArgument {",
     "    if false {\n        return Err(NumericError::InvalidArgument {"),
    ("root: unbracketed start accepted", "root.rs",
     "    if fa.signum() == fb.signum() {",
     "    if false {"),
    ("root: exact zero at the lower end not recognised", "root.rs",
     "    if fa == 0.0 {\n        return Ok(exact(a, 0, evaluations));\n    }\n",
     ""),
    # Brent's bracket and tolerance.
    ("root: bracket not re-established after a same-sign step", "root.rs",
     "        if fb.signum() == fc.signum() {\n            c = a;",
     "        if false {\n            c = a;"),
    ("root: x tolerance applied ten times over", "root.rs",
     "let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * x_tolerance;",
     "let tol1 = 2.0 * f64::EPSILON * b.abs() + 5.0 * x_tolerance;"),
    ("root: interpolation step always accepted", "root.rs",
     "            if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min((e * q).abs()) {",
     "            if true {"),
    # Polynomial roots.
    ("poly: rounding bound shrunk", "poly.rs",
     "let bound = 2.0 * gamma * magnitude",
     "let bound = 0.01 * gamma * magnitude"),
    ("poly: unresolved region undercounts", "poly.rs",
     "                        max_count: crit_count + 1,",
     "                        max_count: crit_count.max(1),"),
    ("poly: sign changes between critical points ignored", "poly.rs",
     "if i + 1 < points.len() && signs[i + 1] != 0 && signs[i] != signs[i + 1] {",
     "if i + 1 < points.len() && signs[i + 1] != 0 && signs[i] == signs[i + 1] {"),
    ("poly: Brent bracket trusted without certified signs", "poly.rs",
     "                let lo = if self.sign(root.lower) != 0 {",
     "                let lo = if true {"),
    ("poly: derivative coefficients off by one", "poly.rs",
     "            .map(|(i, &a)| a * i as f64)",
     "            .map(|(i, &a)| a * (i as f64 + 1.0))"),
    # Quadrature.
    ("quad: no adaptive refinement", "quad.rs",
     "    while total_error > target(total_value) {",
     "    while false && total_error > target(total_value) {"),
    ("quad: Gauss weights misapplied", "quad.rs",
     "            gauss += WG[j / 2] * (f1 + f2);",
     "            gauss += WG[j / 2] * (2.0 * f1);"),
    ("quad: rounding floor of the estimate dropped", "quad.rs",
     "        error = error.max(50.0 * f64::EPSILON * resabs);",
     "        let _ = resabs;"),
    ("quad: reversed interval not negated", "quad.rs",
     "        (upper, lower, -1.0)",
     "        (upper, lower, 1.0)"),
    ("quad: panel budget ignored", "quad.rs",
     "        if heap.len() + settled.len() >= options.max_panels {",
     "        if false {"),
    ("gauss-legendre: weights halved", "quad.rs",
     "            let w = 2.0 / ((1.0 - x * x) * derivative * derivative);",
     "            let w = 1.0 / ((1.0 - x * x) * derivative * derivative);"),
    # Linear algebra.
    ("lu: no row pivoting", "linalg.rs",
     "                .unwrap_or(k);\n            if lu[(p, k)] == 0.0 {",
     "                .map(|_| k)\n                .unwrap_or(k);\n            if lu[(p, k)] == 0.0 {"),
    ("lu: transposed solve unpermuted", "linalg.rs",
     "        for (i, &p) in self.perm.iter().enumerate() {\n            c[p] = w[i];",
     "        for (i, &_p) in self.perm.iter().enumerate() {\n            c[i] = w[i];"),
    ("lu: numerically singular matrices accepted", "linalg.rs",
     "        if !out.condition.is_finite() || out.condition * f64::EPSILON >= 1.0 {\n            return Err(NumericError::Singular {\n                name: \"system matrix\",",
     "        if !out.condition.is_finite() {\n            return Err(NumericError::Singular {\n                name: \"system matrix\","),
    ("solve: residual rounding left out of the backward error", "linalg.rs",
     "    let eta = (r_norm + residual_rounding) / denominator;",
     "    let eta = r_norm / denominator + 0.0 * residual_rounding;"),
    ("cholesky: symmetry not checked", "linalg.rs",
     "                if (a[(i, j)] - a[(j, i)]).abs() > symmetry_tolerance {",
     "                if (a[(i, j)] - a[(j, i)]).abs() > f64::INFINITY {"),
    ("cholesky: negative pivot accepted", "linalg.rs",
     "            if d <= 0.0 || !d.is_finite() {",
     "            if !d.is_finite() {"),
    ("qr: rank threshold dropped", "linalg.rs",
     "            if best <= (m.max(n) as f64) * f64::EPSILON * first_pivot || best == 0.0 {",
     "            if best == 0.0 {"),
    ("constrained: solution left in the rotated basis", "linalg.rs",
     "    let mut x = y;\n    ct.apply_q(&mut x);",
     "    let x = y;"),
    ("constrained: redundant constraints accepted", "linalg.rs",
     "    if ct.rank < p {",
     "    if false {"),
    # Minimisation.
    ("lm: steps that raise the cost accepted", "optimize.rs",
     "        if predicted > 0.0 && actual > 0.0 {",
     "        if predicted > 0.0 && actual > -1e300 {"),
    ("lm: one iteration over budget", "optimize.rs",
     "        if iterations >= options.max_iterations {\n            break Termination::IterationLimit;",
     "        if iterations > options.max_iterations {\n            break Termination::IterationLimit;"),
    ("lm: non-finite start accepted", "optimize.rs",
     "    if r.iter().any(|v| !v.is_finite()) {",
     "    if false {"),
    ("fmin: bracket updated on the wrong side", "optimize.rs",
     "            if u >= x {\n                a = x;",
     "            if u < x {\n                a = x;"),
    ("fmin: stopping test ten times too loose", "optimize.rs",
     "        let tol2 = 2.0 * tol1;",
     "        let tol2 = 20.0 * tol1;"),
]

env = {k: v for k, v in os.environ.items()
       if k not in ("CARGO_HOME", "GIT_DIR", "GIT_WORK_TREE")}
env["PATH"] = str(pathlib.Path.home() / ".cargo/bin") + ":" + env["PATH"]

def tests_pass():
    try:
        r = subprocess.run(["cargo", "test", "-q", "-p", "axiolid-numeric"], cwd=K, env=env,
                           capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        # A mutant that makes a test hang is detected, not survived.
        return False
    return r.returncode == 0

assert tests_pass(), "baseline must pass"
survivors = []
for name, fname, old, new in MUTANTS:
    path = SRC / fname
    original = path.read_text()
    assert original.count(old) == 1, f"anchor not unique/missing: {name}"
    path.write_text(original.replace(old, new, 1))
    try:
        killed = not tests_pass()
    finally:
        path.write_text(original)
    print(f"{'killed  ' if killed else 'SURVIVED'}  {name}", flush=True)
    if not killed:
        survivors.append(name)
assert tests_pass(), "tree restored and passing"
print(f"\n{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
