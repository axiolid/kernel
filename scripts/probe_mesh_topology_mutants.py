"""Mutation probe for mesh topology and homology bases (#144).

Each mutant accepts a non-manifold mesh, misclassifies a surface, or hands
back loops that are not a homology basis, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
T = "crates/algorithms/query/inspect/src/topology.rs"
TESTS = ["-p", "axiolid-inspect", "--test", "topology"]

MUTANTS = [
    ('edges on three triangles accepted', T, 'filter(|uses| uses.len() > 2)', 'filter(|uses| uses.len() > 3)', TESTS),
    ('vertices with two fans accepted', T, '        if fans > 1 {', '        if fans > 2 {', TESTS),
    ('degenerate triangles accepted', T, '        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {', '        if false {', TESTS),
    ('orientation conflicts ignored', T, '                    orientable = false;', '                    orientable = orientable;', TESTS),
    ('inconsistent winding reported consistent', T, '                consistent &= !flips;', '', TESTS),
    ('boundary loops left out of the genus', T, '        let deficit = 2 - characteristic - boundary_loops as i64;', '        let deficit = 2 - characteristic;', TESTS),
    ('boundary loops counted as one', T, '        let boundary_loops = count_loops(&boundary);', '        let boundary_loops = usize::from(!boundary.is_empty());', TESTS),
    ('crosscaps halved like a genus', T, '                crosscaps: u32::try_from(deficit).unwrap_or(0),', '                crosscaps: u32::try_from(deficit / 2).unwrap_or(0),', TESTS),
    ('basis offered with boundary', T, '(orientable && boundary.is_empty()).then(', '(orientable).then(', TESTS),
    ('no dual tree: every non-tree edge a loop', T, '        if find(&mut dual, s) == find(&mut dual, t) {', '        if true {', TESTS),
    ('loop repeats its common ancestor', T, '            up_b.pop();\n', '', TESTS),
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
