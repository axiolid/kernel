"""Mutation probe for #198: seeded weighted maps, cost regions touching
walls, turned cost edges, and forced walks over weighted maps.

Each mutant makes a bound unsound or loose past what the tests allow, or
drops a refusal, and must turn a test red.

Equivalent mutants, deliberately not listed:

- the forced walk's stopping goal without the maps' bracket width: it
  only spends more cells, and the tests' budgets absorb them;
- the wall-hugging lower bound without its side check: once vertices are
  moved across the walls no hugged piece lies on a wall's free side.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
W = "crates/algorithms/planar/route/src/weighted.rs"
F = "crates/algorithms/planar/route/src/forced.rs"
TESTS = [
    "-p", "axiolid-route", "--lib",
    "--test", "weighted", "--test", "weighted_touch",
    "--test", "weighted_seeded", "--test", "weighted_forced_walk",
]

MUTANTS = [
    # Seeded targets.
    ('upper search ignores the weights', W, '                sources.push((state, seeds[t]));', '                sources.push((state, 0.0));\n                let _ = seeds[t];', TESTS),
    ('a seeded state ends its chain', W, '            while previous[at] != usize::MAX && steps <= seed.len() {', '            while seed[at] == usize::MAX && previous[at] != usize::MAX && steps <= seed.len() {', TESTS),
    ('a negative weight accepted', W, '    if let Some(index) = seeds.iter().position(|w| !(w.is_finite() && *w >= 0.0)) {', '    if let Some(index) = seeds.iter().position(|w| !w.is_finite()) {', TESTS),
    # Touching walls.
    ('no vertex moved across a wall', W, '                    *v = moved;', '                    let _ = moved;', TESTS),
    ('a touch refused as a crossing', W, '                    if region_edge && near {', '                    if false && region_edge && near {', TESTS),
    ('any crossing of a region edge accepted', W, '                    if region_edge && near {', '                    if region_edge {', TESTS),
    ('an edge beyond a wall gets intervals', W, '        if weights.beyond_wall(p, q, reach)? {', '        if false {', TESTS),
    ('a stretch along a wall at 1 below', W, '                resolved = known;', '                resolved = false && known;', TESTS),
    # Turned cost edges.
    ('a hop along an edge at the greatest factor above', W, '                let mut all_hug = stretch > tiny;', '                let mut all_hug = false;', TESTS),
    ('an edge blocks its own intervals', W, '            if piece.line != la && piece.line != lb && blocks(piece.p, piece.q)? {', '            if blocks(piece.p, piece.q)? {\n                let _ = (la, lb);', TESTS),
    ('interval sides on the interpolated ends', W, '        let d = q - p;\n        let (a, b) = (p, q);\n        let (mut left, mut right) = (1.0f64, 1.0f64);', '        let d = q - p;\n        let _ = (p, q);\n        let (mut left, mut right) = (1.0f64, 1.0f64);', TESTS),
    # Forced walks over weighted maps.
    ('the pair bound at the steepest factor', F, '    let factor = map.inside_factor(t)?;', '    let factor = map.steepest(t)?;', TESTS),
    ('a cell across a cost edge at its factor', W, '            if meets_triangle(piece.p, piece.q, t)? {\n                return Ok(1.0);\n            }', '            let _ = meets_triangle(piece.p, piece.q, t)?;', TESTS),
    ('the cell bound at slope 2', F, '    f - 2.0 * k * radius\n}', '    let _ = k;\n    f - 2.0 * radius\n}', TESTS),
    ('an origin\'s weight left out of the shortest walk', F, '            shortest.lower = shortest.lower.min(weight + lo);\n            shortest.upper = shortest.upper.min(weight + hi);', '            shortest.lower = shortest.lower.min(lo);\n            shortest.upper = shortest.upper.min(hi);\n            let _ = weight;', TESTS),
    ('maps over different costs matched', W, '        self.region == other.region && self.walls == other.walls && self.costs == other.costs', '        self.region == other.region && self.walls == other.walls', TESTS),
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
