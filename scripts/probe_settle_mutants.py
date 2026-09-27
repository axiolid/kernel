"""Mutation probe for settled region outputs (#188 follow-up).

Each mutant lets an operation return rings its own validation rejects --
by not settling an output path, not merging short edges, not splitting a
ring that passes a point twice, or not putting a touching vertex on the
edge it touches -- and must turn the round-trip test red.

Equivalent mutant, deliberately not listed: not merging edges shorter than
the tolerance first. The touching pass then snaps the short edge's end
onto its neighbour, and the split drops the empty piece between them.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
L = "crates/algorithms/planar/overlay/src/lib.rs"
O = "crates/algorithms/planar/overlay/src/offset.rs"
S = "crates/algorithms/planar/overlay/src/settle.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "roundtrip"]

MUTANTS = [
    ('overlay output not settled', L, 'let polygons = settle::settle(shapes_to_polygons(shapes), tolerance);', 'let polygons = shapes_to_polygons(shapes);', TESTS),
    ('offset output not settled', O, '        crate::settle::settle(\n            to_kernel(backend_shape(polygons).outline(&style)),\n            tolerance,\n        )', '        to_kernel(backend_shape(polygons).outline(&style))', TESTS),
    ('rings passing a point twice kept whole', S, '        if let Some((i, j)) = split_point(&ring) {', '        if let Some((i, j)) = split_point(&ring).filter(|_| false) {', TESTS),
    ('touching vertices left off the edge', S, '        match touching(&ring, eps) {', '        match None::<Vec<Point2>>.or_else(|| { let _ = touching; None }) {', TESTS),
    ('pieces keep the ring\'s kind regardless of winding', S, '                if same != hole {', '                if !hole {', TESTS),
    ('holes judged by their first vertex', L, '    hole.points\n        .iter()\n        .any(|&q| !on_boundary(q) && !contains(outer, q))', '    let _ = on_boundary;\n    !contains(outer, hole.points[0])', TESTS),
    ('union_soup output not settled', L, '    Ok(settle::settle(shapes_to_polygons(shapes), tolerance))', '    Ok(shapes_to_polygons(shapes))', TESTS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
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
