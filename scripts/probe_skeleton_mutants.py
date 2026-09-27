"""Mutation probe for region skeletons (#139) and the constrained
Delaunay repair they rely on.

Each mutant leaves the triangulation non-Delaunay after recovering
constraints, keeps spurs into corners, points an end at the wrong wall, or
narrows a clearance interval, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/planar/route/src/skeleton.rs"
R = "crates/algorithms/planar/triangulate/src/recover.rs"
SK = ["-p", "axiolid-route", "--test", "skeleton"]
TR = ["-p", "axiolid-triangulate", "--test", "cdt"]

MUTANTS = [
    ('Delaunay not restored after recovery', R, '    restore_delaunay(tri);\n', '', TR),
    ('Delaunay not restored (seen by the skeleton)', R, '    restore_delaunay(tri);\n', '', SK),
    ('spurs kept: no spread filter', S, '            spread(nodes[i], c, 0.25 * c, &segments) >= prune * c', '            true', SK),
    ('short leaf branches kept', S, '            if adj[here].len() >= 3 && length < clearance[here].1 {', '            if false && adj[here].len() >= 3 && length < clearance[here].1 {', SK),
    ('the farthest wall ahead', S, '        if best.is_none_or(|(bt, _)| t < bt) {', '        if best.is_none_or(|(bt, _)| t > bt) {', SK),
    ('clearance not rounded outward', S, '    let lo = d2.lo.max(0.0).sqrt().next_down().max(0.0);', '    let lo = d2.hi.max(0.0).sqrt();', SK),
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
