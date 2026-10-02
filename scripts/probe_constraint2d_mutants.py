"""Mutation probe for Apollonius and tangent-circle constructions (#159,
ledger row B14).

Each mutant weakens one decision in
`crates/algorithms/planar/constraint2d/src/lib.rs` -- an equation coefficient,
a sign, a refusal threshold -- and
`crates/algorithms/planar/constraint2d/tests/tangent_circles.rs` must turn a
test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
T_FILE = "crates/algorithms/planar/constraint2d/src/lib.rs"
T = ["-p", "axiolid-constraint2d", "--test", "tangent_circles"]

MUTANTS = [
    (
        "circle constraint drops the sign choice",
        T_FILE,
        "                d: -2.0 * sign * rho,",
        "                d: -2.0 * rho,",
    ),
    (
        "line constraint drops the sign choice",
        T_FILE,
        "                d: -sign,",
        "                d: -1.0,",
    ),
    (
        "line normal not normalized",
        T_FILE,
        "            let n = line.direction.perp().normalize();",
        "            let n = line.direction.perp();",
    ),
    (
        "point radius term not squared away",
        T_FILE,
        "            f: p.x * p.x + p.y * p.y,",
        "            f: p.x + p.y,",
    ),
    (
        "quadratic-elimination base picks the wrong other index",
        T_FILE,
        "        lines.push(if eqs[i].a != 0.0 {\n            sub(&eqs[base], &eqs[i])\n        } else {\n            eqs[i]\n        });",
        "        lines.push(eqs[i]);",
    ),
    (
        "cx(r) numerator uses the wrong equation's d",
        T_FILE,
        "    let m1 = (l1.d * l0.c - l0.d * l1.c) / det;",
        "    let m1 = (l0.d * l0.c - l1.d * l1.c) / det;",
    ),
    (
        "closure quadratic drops the r^2 term",
        T_FILE,
        "    let a_coef = m1 * m1 + m2 * m2 - 1.0;",
        "    let a_coef = m1 * m1 + m2 * m2;",
    ),
    (
        "radius floor removed: zero/negative roots accepted",
        T_FILE,
        "                if r.is_finite() && cx.is_finite() && cy.is_finite() && r > linear {",
        "                if r.is_finite() && cx.is_finite() && cy.is_finite() {",
    ),
    (
        "non-positive circle radius accepted",
        T_FILE,
        "                if circle.radius <= 0.0 {",
        "                if circle.radius < 0.0 {",
    ),
    (
        "zero-direction line accepted",
        T_FILE,
        "                if line.direction.length() <= linear {",
        "                if line.direction.length() < 0.0 {",
    ),
    (
        "3x3 singular system not refused",
        T_FILE,
        "    if det.abs() <= Scalar::EPSILON * 1e8 * scale * scale * scale {",
        "    if det == 12345.0 {",
    ),
    (
        "deduplication distance check dropped",
        T_FILE,
        "            .any(|&(c, r)| (c - centre).length() <= linear && (r - radius).abs() <= linear);",
        "            .any(|&(c, r)| (r - radius).abs() <= linear);",
    ),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(T)
        except subprocess.TimeoutExpired:
            code = -1
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
