"""Mutation probe for branches crossing where surfaces touch (#119,
ADR 0077).

Each mutant weakens one step -- bridging the branch ends to the crossing,
the bridge's tangent at the crossing, its reversal, degenerate contact and
its tangent, telling touching from undecided, lines of contact where the
surfaces are tangent and cross, and ordering pieces that touch to third
order at a vertex -- and must turn a test red.

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
B = "crates/algorithms/construction/brep-boolean/src/section.rs"
P = "crates/algorithms/construction/brep-boolean/src/split.rs"
ROOF = ["-p", "axiolid-brep-boolean", "--test", "splines"]

MUTANTS = [
    ('lines of contact never sought', T, '    let radius = 1e-3 * size;\n    // Tangent without crossing', '    let radius = 1e-3 * size;\n    if radius > 0.0 {\n        return None;\n    }\n    // Tangent without crossing', SPLINE),
    ('regular zeros taken for contact', T, '                if field.jet(p).gradient.length() > 0.25 * beside {', '                if false {', SPLINE),
    ('crossing contact taken for touching', T, '                if (a < 0.0) != (b < 0.0) {\n                    crossing += 1;', '                if (a < 0.0) == (b < 0.0) {\n                    crossing += 1;', SPLINE),
    ('branches into a contact line left loose', T, '                if (q - p).length() > 4.0 * radius {', '                if true {', SPLINE),
    ('a traced contact dropped as touching', B, '                    if !traced && touching(sa, sb, mid, tolerance)? {', '                    if touching(sa, sb, mid, tolerance)? {', ROOF),
    ('pieces touching to third order refused', P, '                        let Some(t) = c.into_iter().find(|t| *t > 1e-12 && *t < TAU - 1e-12)', '                        let Some(t) = None::<Scalar>', ROOF),
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
