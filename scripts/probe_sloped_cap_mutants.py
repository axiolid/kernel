"""Mutation probe for the exact sloped-cap cut (#120 / ADR 0071).

Covers the three layers the change spans: the Sinusoid2 value, its
evaluation, the sloped-cap builder, the half-space clip decisions, and graph
validation. Each mutant runs only the test target that owns its claim.

Equivalent mutant, deliberately not listed: forcing the cut-ellipse span
"always forward" (turn = 1.0). Both frames handed to the intersection point
up, so that branch always yields +1 today; it is kept as a guard against a
frame-convention change, and clip_arc_prism's start/end-vertex check pins
the invariant end to end.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
A = "crates/algorithms/construction/construct/src/extrude_arc.rs"
B = "crates/algorithms/construction/construct/src/boolean_exact.rs"
C = "crates/algorithms/parametric/evaluate/src/curve.rs"
V = "crates/representations/modeling/graph/src/validation.rs"
S = "crates/representations/analytic/curve/src/sinusoid.rs"
CT = ["-p", "axiolid-construct", "--test", "clip_arc_prism"]
ET = ["-p", "axiolid-evaluate", "--test", "sinusoid"]
GT = ["-p", "axiolid-model", "--test", "open_profile_atomic_trim_basis"]

MUTANTS = [
    ('wave sine term sign', S, 'self.mean + self.cosine * c + self.sine * s', 'self.mean + self.cosine * c - self.sine * s', ET),
    ('wave derivative sign', C, 'Ok(Vec2::new(1.0, -w.cosine * sin + w.sine * cos))', 'Ok(Vec2::new(1.0, w.cosine * sin + w.sine * cos))', ET),
    ('wave inversion unchecked', C, 'Curve2::Sinusoid(_) => verify2(curve, point.x, point, linear),', 'Curve2::Sinusoid(_) => Ok(point.x),', ET),
    ('rim wave sine from x axis', A, 'sine: arc.radius * level.gradient.dot(y),', 'sine: arc.radius * level.gradient.dot(x),', CT),
    ('rim wave mean ignores anchor', A, 'mean: level.at(arc.centre) - anchor,', 'mean: level.at(arc.centre),', CT),
    ('cap ellipse axes unprojected', A, '                    x: vector(ellipse.frame.x),\n                    y: vector(ellipse.frame.y),', '                    x: Vec2::new(ellipse.frame.x.x, ellipse.frame.x.y),\n                    y: Vec2::new(ellipse.frame.y.x, ellipse.frame.y.y),', CT),
    ('range ignores arc extremes', B, '            if turned <= arc.sweep.abs() {', '            if turned <= arc.sweep.abs() && false {', CT),
    ('kept side inverted', B, 'let keeps_above = (normal.z > 0.0) == half_space.agreement;', 'let keeps_above = (normal.z > 0.0) != half_space.agreement;', CT),
    ('cut cap named anyway', B, 'extrude_arc_rings_between(&rings, level, Level::flat(prism.top), (false, true))?', 'extrude_arc_rings_between(&rings, level, Level::flat(prism.top), (true, true))?', CT),
    ('crossing cap accepted', B, '        } else if low > prism.bottom + linear && high < prism.top - linear {\n            extrude_arc_rings_between(&rings, level, Level::flat(prism.top), (false, true))?', '        } else if high < prism.top - linear {\n            extrude_arc_rings_between(&rings, level, Level::flat(prism.top), (false, true))?', CT),
    ('validation accepts non-finite wave', V, 'Curve2::Sinusoid(wave) => wave.is_finite(),', 'Curve2::Sinusoid(_) => true,', GT),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=900,
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
