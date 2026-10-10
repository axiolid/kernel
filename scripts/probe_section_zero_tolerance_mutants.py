"""Mutation probe for contour closure at Tolerance::ZERO (#250).

Each mutant drops the evaluation rounding from a joint check, widens it
past rounding, or unshares a tangent point, and must turn a test red.

Equivalent or insensitive mutants, deliberately not listed:

- leaving a line's origin out of its rounding scale. A line's end is
  `origin + direction`, so `|origin| <= |end| + |direction|`: the
  direction's run and the other side's scale, which is at least the
  joint's own magnitude, already cover it to a factor of two.
- adding the rounding to the tolerance instead of taking the larger of
  the two. At ZERO they agree; above it the rounding is ~1e-15 of a
  tolerance of 1e-6.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/construction/construct/src/contour_lower.rs"
S = "crates/algorithms/construction/construct/src/section_lower.rs"
TESTS = ["-p", "axiolid-construct", "--test", "section_zero_tolerance"]
BEAM = ["-p", "axiolid-mesh-compile", "--test", "exact_section_zero_tolerance"]

MUTANTS = [
    ('joint rounding dropped', C, '            if gap > joint_slack(tolerance, previous_scale, lowered.start_scale) {', '            if gap > tolerance.linear() {', TESTS),
    ('joint rounding dropped (beam)', C, '            if gap > joint_slack(tolerance, previous_scale, lowered.start_scale) {', '            if gap > tolerance.linear() {', BEAM),
    ('closing rounding dropped', C, '        if gap > joint_slack(tolerance, last_scale, first_scale) {', '        if gap > tolerance.linear() {', TESTS),
    ('first segment scale never recorded', C, '            first_scale = lowered.start_scale;', '', TESTS),
    ('rounding bound far past rounding', C, 'const JOINT_ROUNDING: Scalar = 8.0 * Scalar::EPSILON;', 'const JOINT_ROUNDING: Scalar = 1e5 * Scalar::EPSILON;', TESTS),
    ('line end scale ignores its run', C, '                end_scale: origin + direction * to.abs(),', '                end_scale: 0.0,', TESTS),
    ('arc scale ignores its centre', C, '            let scale = circle.frame.origin.length() + circle.radius.abs();', '            let scale = circle.radius.abs();', TESTS),
    ('ring vertex taken from the entering end', C, '        let start = lowered.vertices[0].point;', '        let start = previous_end.map_or(lowered.vertices[0].point, |(p, _)| p);\n        let mut lowered = lowered;\n        lowered.vertices[0].point = start;', TESTS),
    ('line runs to the corner, not the shared tangent point', S, '            Some((_, start, _, _)) => start,', '            Some(_) => corners[next].point,', TESTS),
]

probe_lib.main(MUTANTS, timeout=2400)
