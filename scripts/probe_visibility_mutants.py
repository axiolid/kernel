"""Mutation probe for visibility polygons (#184).

Each mutant skips the occlusion test -- which edge is nearest along a
wedge -- or weakens the sweep around it, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
V = "crates/algorithms/planar/overlay/src/visibility.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "visibility"]

MUTANTS = [
    ('occlusion skipped: first edge met wins', V, '        if best.is_none_or(|b| {', '        if best.is_none() || false && best.is_none_or(|b| {', TESTS),
    ('farthest edge wins', V, '                far: b,\n            }) == Sign::Positive', '                far: b,\n            }) == Sign::Negative', TESTS),
    ('edges behind the viewpoint count', V, '        Some(if sn != Sign::Zero && sn == sd {', '        Some(if true || sn != Sign::Zero && sn == sd {', TESTS),
    ('directions not merged', V, '            if at < around.len() && same_direction(v, around[at], a) {\n                continue;\n            }', '', TESTS),
    ('boundary viewpoint accepted', V, '            if s == Sign::Zero && between {\n                return false;\n            }', '', TESTS),
    ('shadow ends never exact', V, '        if sign(&Orient { a: v, b: w, c: end }) == Sign::Zero {\n            return end;\n        }', '', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
