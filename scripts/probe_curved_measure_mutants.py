"""Mutation probe for exact mass properties over curved faces (#125, C17).

Each mutant weakens one decision in the Green's-theorem face integral or in
how `exact_properties` routes a face to it, and must turn a closed-form test
red. The unit tests in `exact_face.rs` own the pole, apex, tube-winding and
B-spline cases no constructor builds; `curved_measure.rs` owns the solids the
constructors do build, including the revolutions whose orientation this work
corrected.

Not listed, deliberately: the pole-jump branch in `join`. A loop that runs
along a pole by less than half a turn takes the same straight segment with or
without it; the branch only matters for a jump of half a turn or more, which
no fixture here has.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
F = "crates/algorithms/query/measure/src/exact_face.rs"
E = "crates/algorithms/query/measure/src/exact.rs"
R = "crates/algorithms/construction/construct/src/revolve_contour.rs"
UNIT = ["-p", "axiolid-measure", "--features", "exact", "--lib"]
SOLIDS = ["-p", "axiolid-construct", "--test", "curved_measure"]
PLANAR = ["-p", "axiolid-construct", "--test", "exact_measure"]
BOTH = (UNIT, SOLIDS)

MUTANTS = [
    ('area element halved', F, '        n.length(),\n', '        n.length() * 0.5,\n', SOLIDS),
    ('volume field 1/4 not 1/3', F, '        w / 3.0,\n', '        w / 4.0,\n', SOLIDS),
    ('first moment reads y for x', F, '        p.x * w / 4.0,\n', '        p.y * w / 4.0,\n', SOLIDS),
    ('second moment reads x for z', F, '        p.z * p.z * w / 5.0,\n', '        p.x * p.x * w / 5.0,\n', SOLIDS),
    ('green sign dropped', F, '(-tangent.x, inner_along_v(surface, reference, point, &floor)?)', '(tangent.x, inner_along_v(surface, reference, point, &floor)?)', SOLIDS),
    ('tube form sign dropped', F, '(tangent.y, inner_along_u(surface, reference, point, &floor)?)', '(-tangent.y, inner_along_u(surface, reference, point, &floor)?)', UNIT),
    ('tube winding ignored', F, '    let along_v = boundary.wraps[1];', '    let along_v = false;', UNIT),
    ('loop winding counted backwards', F, '            periods(shift.x, chart.u_period),', '            periods(-shift.x, chart.u_period),', UNIT),
    ('pole taken on the wrong side', F, '.filter(|pole| (*pole > boundary.anchor.y) == upward)', '.filter(|pole| (*pole > boundary.anchor.y) != upward)', UNIT),
    ('pole ignored, reference at the loop', F, '    } else if boundary.winding[0] != 0 {\n        pole_on_domain_side(&chart, &boundary)?', '    } else if false {\n        pole_on_domain_side(&chart, &boundary)?', UNIT),
    ('unbounded winding accepted', F, '            "face boundary winds around a surface with no pole on the domain side",\n        ))', '            "face boundary winds around a surface with no pole on the domain side",\n        )).or(Ok(boundary.anchor.y))', UNIT),
    ('reversed bound read forward', F, '        let reversed = bound.orientation == Orientation::Reversed;', '        let reversed = false;', UNIT),
    ('seam periods not unwrapped', F, "            period.map_or(0.0, |period| (gap / period).round() * period)", "            period.map_or(0.0, |_| 0.0)", BOTH),
    ('face flip ignored on curved faces', E, '        let sign = if flip_face { -1.0 } else { 1.0 };', '        let sign = 1.0;', SOLIDS),
    ('curve-bounded plane fanned', E, '            if !matches!(curve, Curve3::Line(_)) {\n                return Ok(false);', '            if false {\n                return Ok(false);', SOLIDS),
    ('planar hole area added', E, '        *area += (b - a).cross(c - a) * 0.5;', '        *area += ((b - a).cross(c - a) * 0.5).abs();', PLANAR),
    ('quadrature bound loosened', F, 'const RELATIVE: Scalar = 1e-13;', 'const RELATIVE: Scalar = 1e-2;', UNIT),
    ('revolution frame left-handed again', R, '        y: Vec3::NEG_Z,', '        y: Vec3::Z,', SOLIDS),
]

def run(targets):
    """Non-zero when any target fails; a tuple names several."""
    if not isinstance(targets, tuple):
        targets = (targets,)
    for target in targets:
        code = subprocess.run(
            ["cargo", "test", "-q", *target],
            cwd=ROOT, capture_output=True, text=True, timeout=900,
        ).returncode
        if code != 0:
            return code
    return 0

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
