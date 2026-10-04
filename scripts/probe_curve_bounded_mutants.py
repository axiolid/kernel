"""Mutation probe for curve-bounded planes (#192).

Each mutant fills what must be refused, maps boundaries wrongly, or meshes
the plane facing the wrong way, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/execution/compile/src/bounded.rs"
TESTS = ["-p", "axiolid-mesh-compile", "--test", "curve_bounded"]

MUTANTS = [
    ('implicit outer boundary not refused', B, '    if implicit_outer {', '    if false && implicit_outer {', TESTS),
    ('non-planar basis taken as a plane', B, '        Some(GeometryNode::Surface(axiolid_surface::Surface::Plane(plane))) => plane.frame,', '        Some(GeometryNode::Surface(axiolid_surface::Surface::Plane(plane))) => plane.frame,\n        Some(GeometryNode::Surface(axiolid_surface::Surface::Cylinder(c))) => c.frame,', TESTS),
    ('clockwise outer loop kept', B, '        if index == 0 && signed_area(&ring) < 0.0 {', '        if false && index == 0 && signed_area(&ring) < 0.0 {', TESTS),
    ('parameters mapped with the axes swapped', B, '    let map = |p: Point2| frame.origin + frame.x * p.x + frame.y * p.y;', '    let map = |p: Point2| frame.origin + frame.y * p.x + frame.x * p.y;', TESTS),
    ('open polylines closed silently', B, '        !closed && !meets', '        false && !closed && !meets', TESTS),
    ('3D boundaries off the parameter plane flattened', B, '    if let Some(off) = points.iter().find(|p| p.z.abs() > linear) {', '    if let Some(off) = points.iter().find(|p| p.z.abs() > linear && false) {', TESTS),
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
