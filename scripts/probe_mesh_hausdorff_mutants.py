"""Mutation probe for the certified mesh Hausdorff distance (#148).

Each mutant makes a Hausdorff bound unsound -- a lower bound above the
distance, an upper bound below it -- or stops refining short of the
requested accuracy, and must turn a test red.

Equivalent mutants, deliberately not listed: the rounding pads (the
16 eps upper margin, the 4 eps support and box shrinks, the midpoint
drift), which move bounds by a few units in the last place that no
closed form here resolves; box pruning weakened (only slower); and a
patch's members bounded by their largest rather than smallest distance
(only looser).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
H = "crates/algorithms/query/measure/src/mesh_hausdorff.rs"
P = "crates/algorithms/query/measure/src/mesh_hausdorff/patches.rs"
TESTS = [
    ["-p", "axiolid-measure", "--test", "mesh_hausdorff"],
    ["-p", "axiolid-measure", "--lib", "mesh_hausdorff"],
]

MUTANTS = [
    ("a piece bounded by its nearest corner, not its farthest", H,
     "                    .fold(f64::INFINITY, f64::min)\n            })\n            .fold(0.0, f64::max)",
     "                    .fold(f64::INFINITY, f64::min)\n            })\n            .fold(f64::INFINITY, f64::min)"),
    ("the support direction pointing away from the triangle", H,
     "    let w = v - p;\n    let length = w.length();",
     "    let w = p - v;\n    let length = w.length();"),
    ("the support taken at the farthest corner", H,
     "        .map(|&corner| u.dot(v - corner))\n        .fold(f64::INFINITY, f64::min);",
     "        .map(|&corner| u.dot(v - corner))\n        .fold(f64::NEG_INFINITY, f64::max);"),
    ("a sample's lower bound read from its upper bound", H,
     "                        lower = lower.min(lower_at(v, p, facet));",
     "                        lower = lower.min(up);"),
    ("a flat patch's excess dropped", H,
     "best.min(self.worst(corners, &patch.members) + patch.excess)",
     "best.min(self.worst(corners, &patch.members))"),
    ("flat patches never used", P,
     "    patches.retain(|patch| patch.excess <= flatness);",
     "    patches.retain(|_| false);"),
    ("a pair's hull excess halved", P,
     "    2.0 * height * (1.0 + 4.0 * EPS)",
     "    height * (1.0 + 4.0 * EPS)"),
    ("a dart paired as a convex quadrilateral", P,
     "    let across_cd = frame.sign(d - c, a - c)? * frame.sign(d - c, b - c)? == -1;",
     "    let across_cd = true;"),
    ("a fan's reflex link corner accepted", P,
     "        if frame.sign(here - v, next - v)? != 1 || frame.sign(here - prev, next - here)? != 1 {",
     "        if frame.sign(here - v, next - v)? != 1 {"),
    ("a doubly wound fan accepted", P,
     "    if (turning - core::f64::consts::TAU).abs() >= 1.0 {",
     "    if turning.is_nan() {"),
    ("refinement stopped at ten times the accuracy", H,
     "        if upper - best.lower <= accuracy || splits >= MAX_SPLITS {",
     "        if upper - best.lower <= 10.0 * accuracy || splits >= MAX_SPLITS {"),
    ("pieces within the accuracy of the lower bound dropped", H,
     "            if child.upper > best.lower {",
     "            if child.upper > best.lower + 1e-6 {"),
    ("the two-sided upper bound taken from the nearer side", H,
     "    distance.upper = forward.upper.max(backward.upper);",
     "    distance.upper = forward.upper.min(backward.upper);"),
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
