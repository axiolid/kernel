"""Mutation probe for sections of two B-spline surfaces (#119, ADR 0077).

Each mutant weakens one step -- the loop-free test, the seed deduplication,
the curve's rates, the solve on a window's edge, the pcurve read from the
solve, the domain's breaks at nodes, or B-spline curve/surface Newton --
and must turn a test red.

Equivalent mutants, deliberately not listed: pinning the other surface's
parameter in the Jacobian of the edge solve (the pin is reapplied after
every step, so Newton still lands on the edge), and scaling a curve
piece's rate by a constant in curve/surface Newton (it still converges,
linearly, well within its iterations). Nor is sending two B-splines through
`exact_surface_intersection` in the boolean as well. It traces the whole
domains, finds no line or conic, and falls through to the windowed trace:
only slower.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/parametric/nurbs/src/pair_trace.rs"
S = "crates/representations/analytic/curve/src/pair_section.rs"
E = "crates/algorithms/parametric/evaluate/src/curve.rs"
M = "crates/algorithms/query/measure/src/exact_domain.rs"
PAIR = ["-p", "axiolid-nurbs", "--test", "pair_section"]
BOOL = ["-p", "axiolid-brep-boolean", "--test", "splines"]

MUTANTS = [
    ('normal cones never compared', P, 'if angle - ha - hb > 1e-9 && (core::f64::consts::PI - angle) - ha - hb > 1e-9 {', 'if true {', PAIR),
    ('seeds on traced curves traced again', P, '        if covered(&out, &seed) {', '        if false {', PAIR),
    ('rates not per chord', S, '        let x = solve4(matrix(&j1, &j2, d), [0.0, 0.0, 0.0, d.dot(d)])?;', '        let x = solve4(matrix(&j1, &j2, d), [0.0, 0.0, 0.0, 1.0])?;', BOOL),
    ('window edges never crossed', P, '                let outside = if value == w.0[axis] {\n                    x1 < value', '                let outside = if false {\n                    x1 < value', BOOL),
    ('pcurve reads the other surface', E, '                return finite2(if first { a } else { b }, "curve point");', '                return finite2(if first { b } else { a }, "curve point");', BOOL),
    ('no domain breaks at nodes', M, 'axiolid_curve::Curve3::ImplicitSection(_) | axiolid_curve::Curve3::PairSection(_)', 'axiolid_curve::Curve3::ImplicitSection(_)', BOOL),
    ('hit read in the piece, not the curve', P, '                        let t = k0 + (k1 - k0) * local.clamp(0.0, 1.0);', '                        let t = local.clamp(0.0, 1.0);', PAIR),
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
