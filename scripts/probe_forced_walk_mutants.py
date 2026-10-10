"""Mutation probe for the shortest walk forced through a region (#196).

Each mutant weakens the bracket -- the Lipschitz radius, the vertex-pair
bound, the proofs that a vertex is hidden from a cell, the floor of the
shortest walk -- or the handling of mismatched maps and unreached parts,
and must turn a test red.

Equivalent mutants, deliberately not listed:

- dropping the widening for the rounding of lengths, a few ulps of the
  length, far inside the tolerance;
- dropping the floor of the shortest walk overall. It is sound and cheap,
  but the vertex-pair bound reaches the shortest walk on every cell the
  walk passes through, so the floor never binds in a closed-form test.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
F = "crates/algorithms/planar/route/src/forced.rs"
TESTS = ["-p", "axiolid-route", "--test", "forced_walk", "--lib"]

MUTANTS = [
    ('Lipschitz constant 1, not 2', F, '    f - 2.0 * radius\n', '    f - radius\n', TESTS),
    ('inherited anchor without radius', F, '(lipschitz_lower(f, radius(a)), (a, f))', '(f, (a, f))', TESTS),
    ('pair bound adds both terms', F, 'best = best.min(du + dv + (u - v).length().max(gu + gv));', 'best = best.min(du + dv + (u - v).length() + gu + gv);', TESTS),
    ('pair bound unused', F, '(bound.max(pairs) - slack(depth, radius(anchor.0)))', '(bound - slack(depth, radius(anchor.0)))', TESTS),
    ('pair search stops too early', F, '                if s + r >= best {', '                if s + r >= best - 1.0 {', TESTS),
    ('hidden when one corner is blocked', F, '            if crosses(v, t[0], p, q)? && crosses(v, t[1], p, q)? && crosses(v, t[2], p, q)? {', '            if crosses(v, t[0], p, q)? {', TESTS),
    ('no obstacle hides a vertex', F, '            if crosses(v, t[0], p, q)? && crosses(v, t[1], p, q)? && crosses(v, t[2], p, q)? {', '            if false {', TESTS),
    ('solid wedge on the wrong side', F, '(small_is_inside == hole).then_some', '(small_is_inside != hole).then_some', TESTS),
    ('solid wedge from one corner', F, '            if inside(t[0])? && inside(t[1])? && inside(t[2])? {\n                continue \'vertex;', '            if inside(t[0])? {\n                continue \'vertex;', TESTS),
    ('upper bound from points outside the polygon', F, '            if crate::map::in_polygon(self.through, a)? && self.upper', '            if self.upper', TESTS),
    ('mismatched maps accepted', F, ') -> Result<ForcedWalk, FarthestError> {\n    if !(tolerance.is_finite() && tolerance >= 0.0) {\n        return Err(FarthestError::InvalidTolerance);\n    }\n    if !from.same_space(to) {', ') -> Result<ForcedWalk, FarthestError> {\n    if !(tolerance.is_finite() && tolerance >= 0.0) {\n        return Err(FarthestError::InvalidTolerance);\n    }\n    if false {', TESTS),
    ('an unreached part is an error', F, '            // No walk reaches this part of the free space.\n            Anchored::Unreached => {}', '            Anchored::Unreached => return Err(FarthestError::Triangulation),', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
