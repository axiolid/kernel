"""Mutation probe for exact ruled sections (#119, ADR 0076).

Each mutant weakens one step -- the coefficients of the substituted
quadratic, the exact span decisions, the branch evaluation, or a refusal --
and must turn a test red.

Equivalent mutant, deliberately not listed: testing the discriminant with
`!= Negative` instead of `== Positive` at a span's sample point. Samples lie
strictly between isolated roots of the same square-free polynomial, so the
discriminant is never zero there; touching surfaces are pinned instead by
`apart_and_touching_pairs_are_named` through the empty-span path.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/parametric/nurbs/src/ruled_section.rs"
Q = "crates/representations/analytic/curve/src/quadric_section.rs"
RULED = ["-p", "axiolid-nurbs", "--test", "ruled_section"]
PAIRS = ["-p", "axiolid-nurbs", "--test", "exact_intersection"]
EVAL = ["-p", "axiolid-evaluate", "--test", "quadratic_graph"]

MUTANTS = [
    ('trig derivative sign', Q, '        -self.cos * s + self.sin * c - 2.0 * self.cos2 * s2 + 2.0 * self.sin2 * c2', '        self.cos * s + self.sin * c - 2.0 * self.cos2 * s2 + 2.0 * self.sin2 * c2', EVAL),
    ('branches swapped', Q, '        let denominator = -b - self.branch.sign() * root;', '        let denominator = -b + self.branch.sign() * root;', EVAL),
    ('span-end rounding refused', Q, '        if d < -1e-12 * scale {', '        if d < 0.0 {', RULED),
    ('squares not halved', R, '        f(&x.e0, &y.e0).add(&cc.add(&ss).mul(&half())),', '        f(&x.e0, &y.e0).add(&cc.add(&ss)),', RULED),
    ('discriminant 2ac', R, '        let discriminant = psub(&pmul(&tb, &tb), &pmul(&trig_scale_poly(&ta, 4), &tc));', '        let discriminant = psub(&pmul(&tb, &tb), &pmul(&trig_scale_poly(&ta, 2), &tc));', RULED),
    ('nappe side inverted', R, '            Some(sn == Sign::Zero || (sn == sb))', '            Some(sn == Sign::Zero || (sn != sb))', RULED),
    ('root at pi ignored', R, '        for (start, end) in spans(&breaks, &good, at_pi_ok, vanishes) {', '        for (start, end) in spans(&breaks, &good, at_pi_ok, false) {', RULED),
    ('wrap span dropped', R, '        out.push((last.2, first.2 + core::f64::consts::TAU));', '', RULED),
    ('apex plane accepted', R, '        if curve_carrier.slope != 0.0 && is_zero(&nappe_trig) {\n            return Ok(None);\n        }', '', RULED),
    ('degenerate frame accepted', R, '    if dot(&cross, &z).sign() == Some(Sign::Zero) {\n        return Err(refuse());\n    }', '', PAIRS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode

survivors = []
for name, rel, old, new, target in MUTANTS:
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
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
