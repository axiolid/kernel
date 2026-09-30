"""Mutation probe for the boolmesh winding broad phase (#203).

Each mutant brings back a way a query on the grid's far edge is dropped,
and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/providers/mesh/boolmesh/src/csg/manifold/collider.rs"
TESTS = [
    ["-p", "axiolid-mesh-boolean-boolmesh", "--lib"],
    ["-p", "axiolid-mesh-boolean-boolmesh", "--test", "overlapping_grid"],
]

MUTANTS = [
    ("far edge from the rounded cell product", C,
     "            if px < self.min.x || py < self.min.y || px > self.max.x || py > self.max.y {",
     "            let (hx, hy) = (self.min.x + self.cell * self.dim as Real, self.min.y + self.cell * self.dim as Real);\n            if px < self.min.x || py < self.min.y || px > hx || py > hy {"),
    ("the far edge itself excluded", C,
     "            if px < self.min.x || py < self.min.y || px > self.max.x || py > self.max.y {",
     "            if px < self.min.x || py < self.min.y || px >= self.max.x || py >= self.max.y {"),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode

survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = max(run(t) for t in TESTS)
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
