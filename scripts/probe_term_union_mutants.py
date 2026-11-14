"""Mutation probe for the Minkowski term union (#344, ADR 0108).

Each mutant weakens one step of `term_union::region` or of the exact core's
partial mode -- the certificate that drops a term edge strictly inside a
convex piece (an edge along a side line of the piece, an end outside it, or
nothing at all certified), the stops at term vertices (none, or only at the
pieces' vertices), the reading of a piece's sides (a piece holding the
sample not absorbing, the set rings not wound, a carrier read on the wrong
side, carriers counted again as non-carriers) and the order the kept edges
are listed in -- and must turn a test red.

Equivalent mutants, deliberately not listed:
- closed instead of strict containment of the sample in a piece: a
  non-carrier piece's boundary through the sample can only be the open
  segment of a dropped edge, so the sample is strictly inside its
  certifier anyway;
- allowing a ring to certify its own edge: an edge lies on its own ring's
  side line, so the certificate fails anyway;
- the certifier's box test using `overlaps` instead of the two point
  tests: it changes the cost only;
- a split point recorded as a simple stop in the core with winding off:
  stop tags are read only by the winding carry-over;
- testing a piece's sample against the carriers too: a carrier's boundary
  holds the sample, so the strict test fails for it anyway.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
T = "crates/algorithms/planar/overlay/src/term_union.rs"
A = "crates/algorithms/planar/overlay/src/exact_arc/arrangement.rs"
TESTS = ["-p", "axiolid-overlay", "--lib", "--test", "minkowski"]

MUTANTS = [
    (
        "an edge along a side line of a piece is certified",
        T,
        "        a != Sign::Negative && b != Sign::Negative && (a == Sign::Positive || b == Sign::Positive)",
        "        a != Sign::Negative && b != Sign::Negative",
        TESTS,
    ),
    (
        "an edge with an end outside a piece is certified",
        T,
        "        a != Sign::Negative && b != Sign::Negative && (a == Sign::Positive || b == Sign::Positive)",
        "        a != Sign::Negative && (a == Sign::Positive || b == Sign::Positive)",
        TESTS,
    ),
    (
        "no edge is certified",
        T,
        "                        .is_none()\n",
        "                        .is_none()\n                        || true\n",
        TESTS,
    ),
    (
        "no stops at term vertices",
        T,
        "build_partial(&oriented, &keep, &splits)?",
        "build_partial(&oriented, &keep, &[])?",
        TESTS,
    ),
    (
        "stops only at the vertices of piece rings",
        T,
        "    let splits: Vec<Point2> = oriented\n        .iter()",
        "    let splits: Vec<Point2> = oriented[first..]\n        .iter()",
        TESTS,
    ),
    (
        "a piece holding the sample does not absorb",
        T,
        "            if let Some(j) = sample_at {",
        "            if let Some(j) = sample_at.filter(|_| false) {",
        TESTS,
    ),
    (
        "set rings are not wound",
        T,
        "                windings.holding(s, carried, &mut held);",
        "                held.clear();",
        TESTS,
    ),
    (
        "carriers are read on the wrong side",
        T,
        "            if same {",
        "            if !same {",
        TESTS,
    ),
    (
        "a ring's kept edges are listed in reverse",
        A,
        '            edges.retain(|_| *flags.next().expect("a flag per edge"));',
        '            edges.retain(|_| *flags.next().expect("a flag per edge"));\n            edges.reverse();',
        TESTS,
    ),
    (
        "carriers are counted again as non-carriers",
        T,
        [
            (
                "                |j| !carried(first + j) && inside_convex(s, &piece_exact[j]),",
                "                |j| inside_convex(s, &piece_exact[j]),",
            ),
            (
                "                windings.holding(s, carried, &mut held);",
                "                windings.holding(s, |_| false, &mut held);",
            ),
        ],
        TESTS,
    ),
]

# A mutant may make a test loop forever; past this many seconds it counts as
# killed, by timeout.
TIMEOUT = 1800

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
