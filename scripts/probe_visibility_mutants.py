"""Mutation probe for visibility polygons (#184).

Each mutant skips the occlusion test -- which edge is nearest along a
wedge -- or weakens the sweep around it, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
V = "crates/algorithms/planar/overlay/src/visibility.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "visibility"]

MUTANTS = [
    ('occlusion skipped: first edge met wins', V, '        if best.is_none_or(|b| {', '        if best.is_none() || false && best.is_none_or(|b| {', TESTS),
    ('farthest edge wins', V, '                far: b,\n            }) == Sign::Positive', '                far: b,\n            }) == Sign::Negative', TESTS),
    ('edges behind the viewpoint count', V, '        Some(if sn != Sign::Zero && sn == sd {', '        Some(if true || sn != Sign::Zero && sn == sd {', TESTS),
    ('directions not merged', V, '            if at < around.len() && same_direction(v, around[at], a) {\n                continue;\n            }', '', TESTS),
    ('boundary viewpoint accepted', V, '            if s == Sign::Zero && between {\n                return false;\n            }', '', TESTS),
    ('shadow ends never exact', V, '        if sign(&Orient { a: v, b: w, c: end }) == Sign::Zero {\n            return end;\n        }', '', TESTS),
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
