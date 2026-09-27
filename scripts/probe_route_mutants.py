"""Mutation probe for route visibility (#187).

Each mutant drops one condition a visible segment must meet -- the proper
crossing test -- and must turn a test red.

Equivalent mutants, deliberately not listed:

- testing a stretch along an obstacle edge by its midpoint instead of
  accepting it. The midpoint of a stretch on the boundary is on the
  boundary, which counts as inside.
- not cutting the segment at obstacle vertices lying on it. Since #189
  every such vertex is a graph vertex whose sectors are checked on each
  side of travel (scripts/probe_route_sector_mutants.py), and a gap
  outside the region beyond it is a sector that is not free; the cuts
  stay as a second, independent guard.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/planar/route/src/lib.rs"
TESTS = ["-p", "axiolid-route"]

MUTANTS = [
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
