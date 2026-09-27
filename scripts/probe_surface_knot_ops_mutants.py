"""Mutation probe for surface knot removal, degree change and iso-curves
(#141).

Each mutant accepts a lossy change the tolerance forbids, understates the
deviation bound, or builds the wrong surface or curve, and must turn a test
red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/parametric/nurbs/src/surface_ops.rs"
TESTS = ["-p", "axiolid-nurbs", "--test", "surface_knot_ops"]

MUTANTS = [
    ('removal ignores the tolerance', S, '        if bound + step > tolerance {', '        if false && bound + step > tolerance {', TESTS),
    ('rational bound drops the weight term', S, '(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() + radius * (a[3] - b[3]).abs();', '(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();', TESTS),
    ('rational bound ignores the least weight', S, '    let bound = worst / least_weight;', '    let bound = worst;', TESTS),
    ('polynomial bound halved', S, '                .fold(0.0, Scalar::max),\n        );', '                .fold(0.0, Scalar::max)\n                * 0.5,\n        );', TESTS),
    ('removal count not limited by times', S, '    while removed < times && columns[0].knots.contains(&parameter) {', '    while columns[0].knots.contains(&parameter) {', TESTS),
    ('end knots removable', S, '    if index == 0 || index + 1 == surface.u_knots.len() {', '    if false && index == 0 {', TESTS),
    ('reduction ignores the tolerance', S, '    if bound.is_nan() || bound > tolerance {', '    if bound.is_nan() {', TESTS),
    ('reduction accepts rational surfaces', S, '    if surface.weights.is_some() {\n        return Err(GeomError::Unsupported {', '    if false {\n        return Err(GeomError::Unsupported {', TESTS),
    ('iso-curve drops the weights', S, '            let w = surface.weights.as_ref().map_or(1.0, |rows| rows[i][v]);', '            let w = 1.0;', TESTS),
    ('iso-curve basis on the wrong span', S, '            let i = span - degree + k;', '            let i = (span - degree + k + 1).min(count - 1);', TESTS),
    ('iso-curve clamps outside the domain', S, '    if !(parameter >= lo && parameter <= hi) {', '    if parameter.is_nan() {', TESTS),
    ('transpose leaves the net unflipped', S, '        control_points: flip(&surface.control_points),', '        control_points: surface.control_points.clone(),', TESTS),
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
