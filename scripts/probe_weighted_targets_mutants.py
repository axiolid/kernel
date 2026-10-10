"""Mutation probe for distance maps with weighted targets (#197).

Each mutant drops a target's weight somewhere it matters -- the seed, the
route's own length, the forced walk's shortest walk -- or resolves the
target wrongly, or accepts a bad weight, and must turn a test red.

Equivalent mutant, deliberately not listed: seeding a state twice keeping
the later start, not the lesser. Only the lightest target on a point seeds
it, so no state is seeded twice; the lesser is kept for the proof.
"""
import pathlib

import probe_lib

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

probe_lib.main(MUTANTS, timeout=1800)
