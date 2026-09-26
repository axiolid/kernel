"""Mutation probe for branches crossing where surfaces touch (#119,
ADR 0077).

Each mutant weakens one step -- bridging the branch ends to the crossing,
the bridge's tangent at the crossing, its reversal, degenerate contact and
its tangent, or telling touching from undecided -- and must turn a test
red.

Equivalent mutants, deliberately not listed:
- loosening the checks that the four ends lie two by two along the
  Hessian's null directions: they refuse only sections that are not a
  clean crossing, which the tests do not build;
- swapping a bridge's two end slopes, or its two tangents: over a bridge
  the branch turns by its curvature times the bridge's length, a few
  millionths, so the cubic hardly changes;
- reading the rounding from the coefficients alone
  (`Field2::scale_at`): the square about the crossing then starts too
  small, but grows to hold every piece that cannot be certified anyway.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
T = "crates/algorithms/parametric/nurbs/src/implicit_trace.rs"
Q = "crates/representations/analytic/curve/src/implicit.rs"
S = "crates/algorithms/parametric/nurbs/src/implicit_section.rs"
CURVE = ["-p", "axiolid-nurbs", "--test", "implicit_section"]
SPLINE = ["-p", "axiolid-nurbs", "--test", "spline_section"]

MUTANTS = [
    ('degenerate contact always refused', T, '    if changes != ends.len() || changes % 2 != 0 {', '    if true {', CURVE),
    ('tacnode bridges arrive along their chords', T, '        let into = null.unwrap_or(Vec2::new(chord.x, chord.y));', '        let into = Vec2::new(chord.x, chord.y);', CURVE),
    ('touching read as undecided', S, '        TraceRefusal::Touching(_) => ExactIntersectionRefusal::NotRegularCurve,', '        TraceRefusal::Touching(_) => ExactIntersectionRefusal::Undecided,', SPLINE),
    ('crossings not bridged', T, '        bridge(&mut curves, z, periodic)?;\n', '', CURVE),
    ('bridge leaves along the normal', T, '        let leaving = Vec2::new(-g.y, g.x);\n        let cell = ImplicitCell::bridge(p, c, leaving, into);', '        let leaving = Vec2::new(g.x, g.y);\n        let cell = ImplicitCell::bridge(p, c, leaving, into);', CURVE),
    ('reversed bridge keeps its ends', Q, '                low: self.high,\n                high: self.low,\n                bridge: Some((m1, m0)),', '                bridge: Some((m1, m0)),', CURVE),
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
