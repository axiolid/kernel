"""Mutation probe for the certified deviation report (#232).

Each mutant loosens a step a reported bound rests on: the flattener's
certified acceptance, a derivative or chord bound, the interpolation lemma
over a B-rep triangle, the trim lens, a vertex's offset, the branch and
bound's Taylor remainder or its corner test, a tube's frame derivatives, a
profile's or an instance's stretch, the revolution's turn share. Each must
turn a test red: `deviation_report.rs` samples the exact surfaces against
the compiled meshes and checks every reported bound covers them, and the
`bound`, `certify` and `deviation` unit tests sample the derivative,
chord and remainder bounds directly.

Not probed: dropping the weight term from the rational chord and
interpolation bounds. The proof needs it, but no counterexample turned up
(20000 random rational quadratics stayed under the numerator term alone),
so no test can tell the mutant apart.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/parametric/evaluate/src/bound.rs"
F = "crates/algorithms/parametric/evaluate/src/curve.rs"
B = "crates/execution/compile/src/brep.rs"
K = "crates/execution/compile/src/certify.rs"
P = "crates/execution/compile/src/deviation/paths.rs"
D = "crates/execution/compile/src/deviation.rs"
R = "crates/algorithms/construction/construct/src/profile.rs"
REPORT = ["-p", "axiolid-mesh-compile", "--test", "deviation_report"]
BOUND = ["-p", "axiolid-evaluate", "--lib", "bound"]
CERTIFY = ["-p", "axiolid-mesh-compile", "--lib"]

MUTANTS = [
    ('flattener accepts on the midpoint sagitta alone', F,
     'if sagitta(pa, pb, pm) <= tol && bound(a, b).is_none_or(|certified| certified <= tol) {',
     'if sagitta(pa, pb, pm) <= tol {', [BOUND, REPORT]),
    ('ellipse second derivative from the minor axis only', E,
     '    let main = beta2 + (alpha2 - beta2) * cos2;',
     '    let main = alpha2.min(beta2);', [BOUND, REPORT]),
    ('arc sagitta quartered', E,
     '(orthonormal && h <= core::f64::consts::PI).then(|| radius.abs() * (1.0 - (0.5 * h).cos()))',
     '(orthonormal && h <= core::f64::consts::PI).then(|| 0.25 * radius.abs() * (1.0 - (0.5 * h).cos()))',
     [BOUND, REPORT]),
    ('Taylor chord factor halved', E,
     '    h * h * 0.125 * second',
     '    h * h * 0.0625 * second', [BOUND, REPORT]),
    ('cylinder curvature halved', E,
     '                    duu: s * c.radius.abs(),',
     '                    duu: 0.5 * s * c.radius.abs(),', [BOUND, REPORT]),
    ('torus twist term dropped', E,
     '                    duv: s * r * max_abs_sin(v0, v1),',
     '                    duv: 0.0,', [BOUND, REPORT]),
    ('B-rep interpolation lemma halved', B,
     '                best = best.min(0.5 * radius * radius);',
     '                best = best.min(0.25 * radius * radius);', [REPORT]),
    ('B-rep enclosing circle shrunk', B,
     '    ((ab * bc * ca).sqrt() / (2.0 * twice_area)).max(0.5 * longest.sqrt())',
     '    0.5 * longest.sqrt()', [CERTIFY, REPORT]),
    ('B-rep trim lens ignored', B,
     '            self.trim = self.trim.max(g * d);',
     '            self.trim = self.trim.max(0.0 * g * d);', [REPORT]),
    ('B-rep vertex offset ignored', B,
     '            self.interior = self.interior.max(interp + offset);',
     '            self.interior = self.interior.max(interp + 0.0 * offset);', [CERTIFY, REPORT]),
    ('branch and bound Taylor remainder halved', K,
     '    let e2 = 0.5 * (a * hx * hx + 2.0 * b * hx * hy + c * hy * hy);',
     '    let e2 = 0.25 * (a * hx * hx + 2.0 * b * hx * hy + c * hy * hy);', [CERTIFY, REPORT]),
    ('branch and bound reads the centre, not the corners', K,
     '                        .map(|&c| lattice[c][k])\n                        .fold(0.0, Scalar::max)',
     '                        .map(|&c| lattice[c][k])\n                        .fold(Scalar::INFINITY, Scalar::min)',
     [CERTIFY, REPORT]),
    ('tube frame curvature dropped', P,
     '            d2 + r * (m2 * m2 + b2 * b2).sqrt(),',
     '            d2,', [CERTIFY, REPORT]),
    ('end disk ring curvature dropped', P,
     '        Some([0.0, 1.0, rho.0.abs().max(rho.1.abs())])',
     '        Some([0.0, 1.0, 0.0])', [CERTIFY, REPORT]),
    ('revolution turn share dropped', P,
     '            proven(profile(*id, half)?, half),',
     '            proven(profile(*id, half)?, 0.0),', [REPORT]),
    ('derived profile stretch ignored', R,
     '.scaled(stretch2(transform))',
     '.scaled(1.0)', [REPORT]),
    ('instance stretch ignored', D,
     '        let factor = stretch(transform);',
     '        let factor = 1.0;', [REPORT]),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode

survivors = []
for name, rel, old, new, targets in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    code = 0
    try:
        for target in targets:
            try:
                code = run(target)
            except subprocess.TimeoutExpired:
                code = -1
            if code != 0:
                break
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
