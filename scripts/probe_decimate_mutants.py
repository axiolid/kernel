"""Mutation probe for the decimator's collapse guards (#201).

Each mutant drops one of the two checks that refuse a damaging edge
collapse, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/discrete/decimate/src/collapse.rs"
TESTS = ["-p", "axiolid-decimate", "--test", "decimation"]

MUTANTS = [
    ('link condition dropped', C, '    if shared != 2 {\n        return None;\n    }', '', TESTS),
    ('normal-inversion guard dropped', C, '            if after.length_squared() == 0.0 || before.dot(after) <= 0.0 {', '            if after.length_squared() == 0.0 {', TESTS),
    ('normal-inversion guard admits a flip', C, '            if after.length_squared() == 0.0 || before.dot(after) <= 0.0 {', '            if after.length_squared() == 0.0 || before.dot(after) < -1.0 {', TESTS),
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
