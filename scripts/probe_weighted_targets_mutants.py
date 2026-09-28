"""Mutation probe for distance maps with weighted targets (#197).

Each mutant drops a target's weight somewhere it matters -- the seed, the
route's own length, the forced walk's shortest walk -- or resolves the
target wrongly, or accepts a bad weight, and must turn a test red.

Equivalent mutant, deliberately not listed: seeding a state twice keeping
the later start, not the lesser. Only the lightest target on a point seeds
it, so no state is seeded twice; the lesser is kept for the proof.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
M = "crates/algorithms/planar/route/src/map.rs"
F = "crates/algorithms/planar/route/src/forced.rs"
TESTS = ["-p", "axiolid-route", "--test", "weighted_targets", "--test", "distance_map"]

MUTANTS = [
    ('targets seeded at zero', M, 'sources.push((state, weights[t]));', 'sources.push((state, 0.0 * weights[t]));', TESTS),
    ('a seeded state ends the chain', M, 'while next[at] != usize::MAX && steps <= seed.len() {', 'while seed[at] == usize::MAX && next[at] != usize::MAX && steps <= seed.len() {', TESTS),
    ('route length keeps the weight', M, 'length: length - self.weights[target],', 'length,', TESTS),
    ("forced walk drops the origin's weight", F, 'shortest = shortest.min(weight + d);', 'shortest = shortest.min(d + 0.0 * weight);', TESTS),
    ('bad weights accepted', M, 'if let Some(index) = weights.iter().position(|w| !(w.is_finite() && *w >= 0.0)) {', 'if let Some(index) = weights.iter().position(|w| w.is_nan() && false) {', TESTS),
    ('the heavier of two on one point', M, '.min_by(|&a, &b| weights[a].total_cmp(&weights[b]).then(a.cmp(&b)));', '.max_by(|&a, &b| weights[a].total_cmp(&weights[b]).then(a.cmp(&b)));', TESTS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
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
