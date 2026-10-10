"""Mutation probe for exact edge contact in the overlay (#272).

Operand validation and settling share one notion of when two boundary
edges meet (`crates/algorithms/planar/overlay/src/contact.rs`): an end of
one on the other within the tolerance, or a proper crossing, every
decision exact. Each mutant puts back a rounded sign, a rounded
comparison or a missing case at one site and must turn the tests red:
the unit tests in `contact.rs` and `settle.rs` (an integer oracle, exact
translations to 1e5 and 5e6, ulp steps across a line) and
`tests/edge_contact.rs` (the issue's U, collinear disjoint edges with
rounded-zero cross products, overlaps, touches, crossings, a hole corner a
hair inside its outer ring, and settled soups fed back).

Equivalent mutant, deliberately not listed: an exact zero counted as a
side in `segments_intersect` as a whole. A pair with an end exactly on
the other's line and a proper sign change otherwise has that end on the
other segment, which `on_segment` already reports; the unit tests kill it
through `proper_crossing` alone (listed below).

Runs on `scripts/probe_lib.py`: each mutant in its own process group
under a timeout that kills the group, the source always restored, and a
missing or ambiguous anchor aborts the probe. A mutant that does not
compile counts as BROKEN, not killed.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/planar/overlay/src/"
C, L, S = D + "contact.rs", D + "lib.rs", D + "settle.rs"
UNIT = ["-p", "axiolid-overlay", "--lib"]
EDGE = ["-p", "axiolid-overlay", "--test", "edge_contact"]
BOTH = [UNIT, EDGE]
TIMEOUT = 1800

ROUNDED_ON = "(b - a).perp_dot(p - a).abs() <= epsilon"
# The touching test before #272, rounded throughout.
OLD_ON = ("(({b} - {a}).perp_dot({p} - {a}).abs() <= {e}"
          " && {p}.x >= {a}.x.min({b}.x) - {e} && {p}.x <= {a}.x.max({b}.x) + {e}"
          " && {p}.y >= {a}.y.min({b}.y) - {e} && {p}.y <= {a}.y.max({b}.y) + {e})")

MUTANTS = [
    # The strict test.
    ('proper crossing from rounded cross products (#272)', C,
     '    opposite(side(a, b, c), side(a, b, d)) && opposite(side(c, d, a), side(c, d, b))',
     '    let _ = opposite;\n'
     '    let x = |a: Point2, b: Point2, c: Point2| (b - a).perp_dot(c - a) > 0.0;\n'
     '    x(a, b, c) != x(a, b, d) && x(c, d, a) != x(c, d, b)', BOTH),
    ('an exact zero read as the negative side', C,
     '    opposite(side(a, b, c), side(a, b, d)) && opposite(side(c, d, a), side(c, d, b))',
     '    let _ = opposite;\n'
     '    let x = |a: Point2, b: Point2, c: Point2| side(a, b, c) == Sign::Positive;\n'
     '    x(a, b, c) != x(a, b, d) && x(c, d, a) != x(c, d, b)', BOTH),
    ('an exact zero counted as a side', C,
     '            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)',
     '            (Sign::Positive, Sign::Negative)\n'
     '                | (Sign::Negative, Sign::Positive)\n'
     '                | (Sign::Zero, Sign::Positive | Sign::Negative)\n'
     '                | (Sign::Positive | Sign::Negative, Sign::Zero)', BOTH),
    ('proper crossing tests one segment only', C,
     '    opposite(side(a, b, c), side(a, b, d)) && opposite(side(c, d, a), side(c, d, b))',
     '    opposite(side(a, b, c), side(a, b, d))', BOTH),
    # Touching.
    ('collinear ends touch outside the extent', C,
     '    within_extent(a, b, p, epsilon)\n        && cross_sign(a, b, p, epsilon) != Sign::Positive',
     '    true\n        && cross_sign(a, b, p, epsilon) != Sign::Positive', BOTH),
    ('touching from the rounded cross product', C,
     '        && cross_sign(a, b, p, epsilon) != Sign::Positive\n'
     '        && cross_sign(a, b, p, -epsilon) != Sign::Negative',
     '        && ' + ROUNDED_ON, BOTH),
    ('touching only exactly on the line', C,
     '        && cross_sign(a, b, p, epsilon) != Sign::Positive\n'
     '        && cross_sign(a, b, p, -epsilon) != Sign::Negative',
     '        && side(a, b, p) == Sign::Zero', BOTH),
    ('an end on the other segment not checked from the first', C,
     '    on_segment(c, d, a, epsilon)\n        || on_segment(c, d, b, epsilon)\n',
     '    false\n', BOTH),
    ('the widened box rounded (`p >= min - epsilon`)', C,
     '    at_least(p.x, a.x.min(b.x), -epsilon)\n'
     '        && at_least(a.x.max(b.x), p.x, -epsilon)\n'
     '        && at_least(p.y, a.y.min(b.y), -epsilon)\n'
     '        && at_least(a.y.max(b.y), p.y, -epsilon)',
     '    p.x >= a.x.min(b.x) - epsilon\n'
     '        && p.x <= a.x.max(b.x) + epsilon\n'
     '        && p.y >= a.y.min(b.y) - epsilon\n'
     '        && p.y <= a.y.max(b.y) + epsilon', BOTH),
    ('a difference rounding onto the bound decides by itself', C,
     '    error >= 0.0\n}', '    let _ = error;\n    true\n}', BOTH),
    ('the box prefilter drops the tolerance', C,
     '        !at_least(r.max(s), p.min(q), -epsilon) || !at_least(p.max(q), r.min(s), -epsilon)',
     '        !at_least(r.max(s), p.min(q), 0.0) || !at_least(p.max(q), r.min(s), 0.0)', BOTH),
    # Holes: the ray test.
    ('a point inside decided by the rounded ray crossing', L,
     '        if (a.y > p.y) != (b.y > p.y) {\n'
     '            let towards = if b.y > a.y {',
     '        if (a.y > p.y) != (b.y > p.y)\n'
     '            && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x\n'
     '        {\n'
     '            inside = !inside;\n'
     '        }\n'
     '        if false {\n'
     '            let towards = if b.y > a.y {', BOTH),
    ('a hole corner on the outer boundary from the rounded cross product', L,
     '            on_segment(a, b, q, t.linear())',
     '            ' + OLD_ON.format(a='a', b='b', p='q', e='t.linear()'), BOTH),
    # Settling.
    ('settled touches from the rounded cross product', S,
     '            let on = |p: Point2, q: Point2, v: Point2| on_segment(p, q, v, eps);',
     '            let _ = on_segment;\n'
     '            let on = |p: Point2, q: Point2, v: Point2| '
     + OLD_ON.format(a='p', b='q', p='v', e='eps') + ';', BOTH),
    ('a settled crossing found from rounded cross products', S,
     '                let x = rounded_crossing(a, b, c, d)?;',
     '                let _ = rounded_crossing;\n'
     '                let cr = |a: Point2, b: Point2, c: Point2| (b - a).perp_dot(c - a);\n'
     '                let x = a + (b - a) * (cr(c, d, a) / (cr(c, d, a) - cr(c, d, b)));',
     BOTH),
]


probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
