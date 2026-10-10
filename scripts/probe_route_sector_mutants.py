"""Mutation probe for route sectors (#189).

Each mutant lets a route through a point where obstacles meet -- by
ignoring the rays on a side of travel, the freeness of a sector, which
side of a ring the region lies on, or the side an edge is taken on -- and
must turn a test red.

Equivalent mutant, deliberately not listed: passing a vertex on a side
with no ray there without asking whether that side is free. Along a wall
the side's freeness is the same at the segment's ends, where attaching
the edge to a sector already requires it, and it can only change at a
vertex with a ray on that side.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
G = "crates/algorithms/planar/route/src/graph.rs"
TESTS = ["-p", "axiolid-route"]

MUTANTS = [
    ('rays on the side of travel ignored', G, '            if side(a, b, ray.to)? == blocking {\n                return Ok(false);', '            if false && side(a, b, ray.to)? == blocking {\n                return Ok(false);', TESTS),
    ('cut vertices not checked', G, '        if !within(a, b, *v) || side(a, b, *v)? != Sign::Zero {\n            continue;', '        if true || !within(a, b, *v) || side(a, b, *v)? != Sign::Zero {\n            continue;', TESTS),
    ('every sector free', G, '            n => self.rays[sector].ccw_free && self.rays[(sector + 1) % n].cw_free,', '            _ => { let _ = sector; true }', TESTS),
    ('holes treated as outer rings', G, '    Ok(ccw != hole)', '    Ok(ccw)', TESTS),
    ('shared ring edges free on neither side', G, '                self.ccw_free |= other.ccw_free;\n                self.cw_free |= other.cw_free;', '                self.ccw_free &= other.ccw_free;\n                self.cw_free &= other.cw_free;', TESTS),
    ('left edge attached on the right at its far end', G, '                        Side::Left => (stars[i].ccw_of(nodes[j])?, stars[j].cw_of(nodes[i])?),', '                        Side::Left => (stars[i].ccw_of(nodes[j])?, stars[j].ccw_of(nodes[i])?),', TESTS),
    ('arrival on the wrong side', G, '            Side::Left => star.cw_of(from)?,\n            Side::Right => star.ccw_of(from)?,', '            Side::Left => star.ccw_of(from)?,\n            Side::Right => star.cw_of(from)?,', TESTS),
    ('on a ray counts as the sector before it', G, '            Place::On(i) => (i + n - 1) % n,', '            Place::On(i) => i,', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
