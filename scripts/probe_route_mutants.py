"""Mutation probe for route visibility (#187).

Each mutant drops one condition a visible segment must meet -- the cuts at
obstacle vertices lying on it (the collinear case), or the proper crossing
test -- and must turn a test red.

Equivalent mutant, deliberately not listed: testing a stretch along an
obstacle edge by its midpoint instead of accepting it. The midpoint of a
stretch on the boundary is on the boundary, which counts as inside.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/planar/route/src/lib.rs"
TESTS = ["-p", "axiolid-route"]

MUTANTS = [
    ('collinear vertices not cut at', R, '                cuts.push(v);', '', TESTS),
    ('proper crossings ignored', R, '        if crosses(a, b, *p, *q)? {\n            return Ok(false);\n        }', '        if crosses(a, b, *p, *q)? && false {\n            return Ok(false);\n        }', TESTS),
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
