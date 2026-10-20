"""Mutation probe for the rectangle corner-radius reach (#339).

A rectangle's corner radius may pass its limit by rounding alone, at most
`RADIUS_REACH_EPSILONS` relative epsilons of the half-extent the limit is cut
from, and is then snapped onto the limit. Each mutant drops the reach, widens
or narrows it, measures it at the wrong scale, drops it from one ring, or
keeps the unsnapped radius, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/construction/construct/src/section_lower.rs"
TESTS = [
    ["-p", "axiolid-construct", "--lib", "section_lower"],
    ["-p", "axiolid-construct", "--test", "rounded_rectangle"],
]

MUTANTS = [
    ('no reach', S, '|| value > limit + reach {', '|| value > limit {', TESTS),
    ('reach a modelling tolerance', S,
     'let reach = RADIUS_REACH_EPSILONS * Scalar::EPSILON * scale;\n    if',
     'let reach = RADIUS_REACH_EPSILONS * Scalar::EPSILON * scale * 1e10;\n    if', TESTS),
    ('reach twice as wide', S,
     'const RADIUS_REACH_EPSILONS: Scalar = 4.0;', 'const RADIUS_REACH_EPSILONS: Scalar = 8.0;', TESTS),
    ('reach half as wide', S,
     'const RADIUS_REACH_EPSILONS: Scalar = 4.0;', 'const RADIUS_REACH_EPSILONS: Scalar = 2.0;', TESTS),
    ('radius not snapped onto the limit', S, '    Ok(value.min(limit))\n', '    Ok(value)\n', TESTS),
    ('outer ring has no reach', S,
     'corner_radius(rectangle.outer_radius, scale, scale, "outer")',
     'corner_radius(rectangle.outer_radius, scale, 0.0, "outer")', TESTS),
    ('inner ring has no reach', S,
     'corner_radius(rectangle.inner_radius, ix.min(iy), scale, "inner")',
     'corner_radius(rectangle.inner_radius, ix.min(iy), 0.0, "inner")', TESTS),
    ('inner reach at the scale of the core', S,
     'corner_radius(rectangle.inner_radius, ix.min(iy), scale, "inner")',
     'corner_radius(rectangle.inner_radius, ix.min(iy), ix.min(iy), "inner")', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
