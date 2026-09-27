"""Mutation probe for region skeletons (#139) and the constrained
Delaunay repair they rely on.

Each mutant leaves the triangulation non-Delaunay after recovering
constraints, lets legalisation or recovery fail on real room outlines
(#190), rounds the exact incircle's differences, keeps spurs into
corners, points an end at the wrong wall, or narrows a clearance
interval, and must turn a test red (a hang past the timeout counts).

Not listed, deliberately: splitting the one triangle a vertex lands on
into three even when the vertex lies on an edge. Legalisation then flips
the empty triangle away in every case the tests reach; splitting the edge
and both triangles beside it is kept so that correctness does not rest on
that.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/planar/route/src/skeleton.rs"
R = "crates/algorithms/planar/triangulate/src/recover.rs"
B = "crates/algorithms/planar/triangulate/src/build.rs"
P = "crates/algorithms/predicates/src/sphere.rs"
PR = ["-p", "axiolid-predicates", "--test", "incircle_exact"]
SK = ["-p", "axiolid-route", "--test", "skeleton"]
TR = ["-p", "axiolid-triangulate", "--test", "cdt"]

MUTANTS = [
    ('Delaunay not restored after recovery', R, '    restore_delaunay(tri);\n', '', TR),
    ('spurs kept: no spread filter', S, '            spread(nodes[i], c, 0.25 * c, &segments) >= prune * c', '            true', SK),
    ('short leaf branches kept', S, '            if adj[here].len() >= 3 && length < clearance[here].1 {', '            if false && adj[here].len() >= 3 && length < clearance[here].1 {', SK),
    ('the farthest wall ahead', S, '        if best.is_none_or(|(bt, _)| t < bt) {', '        if best.is_none_or(|(bt, _)| t > bt) {', SK),
    ('clearance not rounded outward', S, '    let lo = d2.lo.max(0.0).sqrt().next_down().max(0.0);', '    let lo = d2.hi.max(0.0).sqrt();', SK),
    ('legalize checks the new diagonal again', B, '            stack.push(3 * t + ti);\n            stack.push(3 * ot + (oi + 2) % 3);', '            stack.push(3 * t + (ti + 1) % 3);\n            stack.push(3 * ot + (oi + 1) % 3);', TR),
    ('recovery drops a new diagonal that still crosses', R, '            queue.push_back((x, y));', '            let _ = (x, y);', TR),
    ('recovery drops an edge it cannot flip yet', R, '            queue.push_back((u, v));\n            idle += 1;', '            idle += 1;', TR),
    ('incircle differences rounded', P, '    let (d, err) = two_diff(p, q);', '    let (d, err) = (p - q, 0.0 * two_diff(p, q).1);', PR),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=300,
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
