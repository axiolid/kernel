"""Mutation probe for the slender-soup speed-ups and the intersection cull in the exact overlay (#315): each mutant must be caught.

The subdivision searches a piece's carriers only among the edges sharing
two or more points with its edge, counts afresh only the rings met at a
stop, merges the points met along an edge by one sort, and shares the
crossings of repeated segments; an intersection leaves out the rings whose
box misses the other operand's. Each mutant below bends one of these. The
tests run with debug assertions off, so the cross-checks against the full
search and count are off and a kill comes from a wrong result, or from
`tests/soup_snapshots.rs` and the cull's own tests.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/planar/overlay/src"
TESTS = [["--config", "profile.dev.debug-assertions=false", "-p", "axiolid-overlay"]]
ARR = "exact_arc/arrangement.rs"
OVERLAY = "exact_overlay.rs"

MUTANTS = [
    # Carriers: the edge itself and the edges sharing a stretch with it.
    ("an edge not a candidate carrier of its own pieces", ARR,
     "        candidates.push(index);\n",
     ""),
    ("edges sharing a stretch not candidate carriers", ARR,
     "                shared[index].push(near);\n                shared[near].push(index);\n",
     ""),
    # Rings carried across any other stop.
    ("a ring met at the stop keeps its membership", ARR,
     "                            .filter(|other| met.binary_search(other).is_err())",
     "                            .filter(|_| true)"),
    ("the rings met at the stop not counted again", ARR,
     "                                inside.push(other);\n",
     ""),
    ("a carrier of the piece counted as containing it", ARR,
     "                            if !carried_by(other)\n",
     "                            if true\n"),
    ("the rings met at the next stop", ARR,
     "&rings_met[w[0].2.clone()]",
     "&rings_met[w[1].2.clone()]"),
    ("rings carried past a piece another ring keeps", ARR,
     "                before = None;\n                continue;",
     "                continue;"),
    # Merging the points met along an edge.
    ("the rings meeting a copy of a point lost", ARR,
     "                    held.1 = None;\n                    rings_met.extend(other);",
     "                    held.1 = None;"),
    # Not listed: an unstable sort, which may keep another copy of a point
    # first. Every copy is the same exact point, and a straight piece's
    # output is its exact ends correctly rounded, so no straight output
    # can tell; the sort stays stable so a sub-arc's bulge, read off the
    # copies' approximations, is the one before #315.
    # Repeated segments' crossings.
    ("segments named by their start alone", ARR,
     "            let key = [a.x, a.y, b.x, b.y].map(f64::to_bits);",
     "            let key = [a.x, a.y, a.x, a.y].map(f64::to_bits);"),
    ("a pair's crossings shared with another pair", ARR,
     "                match met_before.entry((named[index], named[near])) {",
     "                match met_before.entry((named[index], named[index])) {"),
    # Not listed: repeated segments left unmarked. Their crossings are then
    # computed for every pair, as before #315: slower, the same output.
    # The intersection cull.
    ("every operation culled", OVERLAY,
     "    if steps.cull && operation == OverlayOperation::Intersection {",
     "    if steps.cull {"),
    ("boxes touching at a corner taken as apart", OVERLAY,
     "    a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3]",
     "    a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3]"),
    ("a ring culled against its own operand", OVERLAY,
     "            let other = if count.clip { subject } else { clip };",
     "            let other = if count.clip { clip } else { subject };"),
    ("every ring kept against an empty operand", OVERLAY,
     "            other.is_some_and(|other| boxes_meet(b, &other))",
     "            other.is_none_or(|other| boxes_meet(b, &other))"),
    ("a ring's box short of its right side", OVERLAY,
     "[b[0].min(p.x), b[1].min(p.y), b[2].max(p.x), b[3].max(p.y)]",
     "[b[0].min(p.x), b[1].min(p.y), b[2].min(p.x), b[3].max(p.y)]"),
]

probe_lib.main(MUTANTS, tests=TESTS, root=SRC, timeout=300)
