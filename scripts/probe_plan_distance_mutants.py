"""Mutation probe for plan distance between exact bodies (#217).

Each mutant makes a plan bound unsound -- a lower bound above the plan
distance, an upper bound not between boundary points' projections -- or
certifies an overlap or a gap it has not shown, and must turn a test red.

Equivalent mutants, deliberately not listed: the sphere fallback of a
horizontal projection using the space radius (only looser), and the
overlap margin scaled up tenfold (only more cautious).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/query/measure/src/exact_distance.rs"
TESTS = [
    ["-p", "axiolid-construct", "--test", "plan_distance"],
    ["-p", "axiolid-measure", "--all-features", "--lib"],
]

MUTANTS = [
    ("the plan gap measured in space", D,
     "        (flat(b.centre - a.centre), a.plan_radius + b.plan_radius)\n    } else {",
     "        (b.centre - a.centre, a.plan_radius + b.plan_radius)\n    } else {"),
    ("plan discs of no size", D,
     "        (flat(b.centre - a.centre), a.plan_radius + b.plan_radius)\n    } else {",
     "        (flat(b.centre - a.centre), 0.0)\n    } else {"),
    ("normals read whole in plan", D,
     "        if !plan {\n            return Some(n);\n        }",
     "        if true {\n            return Some(n);\n        }"),
    ("horizontal turning bounded by the larger axis alone", D,
     "        (a * flat_length(frame.x)).hypot(b * flat_length(frame.y))",
     "        (a * flat_length(frame.x)).max(b * flat_length(frame.y))"),
    ("a plan patch's vertical speed dropped", D,
     "            turning(&c.frame, c.radius.abs(), c.radius.abs()),\n            flat_length(c.frame.z),",
     "            turning(&c.frame, c.radius.abs(), c.radius.abs()),\n            0.0,"),
    ("plan witnesses measured in space", D,
     "        Metric::Plan => flat_length(p - q),",
     "        Metric::Plan => (p - q).length(),"),
    ("an overlap shown without a margin", D,
     "            length > 0.0 && cross(e0, e1, at) / length > margin",
     "            length > 0.0 && cross(e0, e1, at) / length >= 0.0"),
    ("apart when merely not overlapping", D,
     "        None if found.bounds.lower > 0.0 => PlanOverlap::Disjoint {",
     "        None if found.bounds.lower >= 0.0 => PlanOverlap::Disjoint {"),
    ("faces pruned by the offset in space", D,
     "        Metric::Plan => (\n            flat(other.centre - face.centre),",
     "        Metric::Plan => (\n            other.centre - face.centre,"),
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
