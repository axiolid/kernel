"""Mutation probe for the minimum-area rectangle (#182).

Each mutant skips a rotating-calipers step, an edge, the exact area
comparison or the tie rule, and must turn a test red.

Equivalent mutant, deliberately not listed: advancing a caliper only on a
strict increase instead of a non-decrease. On a tie the two vertices are
equally extreme, so the spans, and the area, are the same.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/planar/overlay/src/rectangle.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "rectangle"]

MUTANTS = [
    ('far caliper not advanced', R, '            hi = advance(hi, d, false, true);', '', TESTS),
    ('top caliper not advanced', R, '            top = advance(top, d, true, true);', '', TESTS),
    ('near caliper not advanced', R, '            lo = advance(lo, d, false, false);', '', TESTS),
    ('calipers move one step at most', R, '        for _ in 0..n {\n            if !step(', '        for _ in 0..1 {\n            if !step(', TESTS),
    ('half the hull edges tried', R, '    for base in 0..n {', '    for base in 0..n.div_ceil(2) {', TESTS),
    ('larger area preferred', R, '            Some(Sign::Negative) => {\n                best = *c;', '            Some(Sign::Positive) => {\n                best = *c;', TESTS),
    ('ties keep the first found', R, '                if clockwise(canonical(best.d), axis) {', '                if false && clockwise(canonical(best.d), axis) {', TESTS),
    ('tied orientations not counted', R, '                    ties.push(axis);', '', TESTS),
    ('collinear hull vertices kept', R, '                if sign(&Orient { a, b, c }) == Sign::Positive {', '                if sign(&Orient { a, b, c }) != Sign::Negative {', TESTS),
    ('segment width left rounded', R, '            rectangle.half_extents[usize::from(axis == d || axis == -d)] = 0.0;', '', TESTS),
    ('measured error never used', R, '            measured_error(&h, &rectangle).unwrap_or(64.0 * f64::EPSILON * size),', '            64.0 * f64::EPSILON * size,', TESTS),
    ('corners left out of the measured error', R, '        (c0.x, exact(x0)),\n        (c0.y, exact(y0)),\n        (c1.x, exact(x1)),\n        (c1.y, exact(y0)),\n        (c2.x, exact(x1)),\n        (c2.y, exact(y1)),\n        (c3.x, exact(x0)),\n        (c3.y, exact(y1)),\n', '', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
