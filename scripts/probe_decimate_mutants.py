"""Mutation probe for the decimator's collapse guards (#201).

Each mutant drops one of the two checks that refuse a damaging edge
collapse, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/discrete/decimate/src/collapse.rs"
TESTS = ["-p", "axiolid-decimate", "--test", "decimation"]

MUTANTS = [
    ('link condition dropped', C, '    if shared != 2 {\n        return None;\n    }', '', TESTS),
    ('normal-inversion guard dropped', C, '            if after.length_squared() == 0.0 || before.dot(after) <= 0.0 {', '            if after.length_squared() == 0.0 {', TESTS),
    ('normal-inversion guard admits a flip', C, '            if after.length_squared() == 0.0 || before.dot(after) <= 0.0 {', '            if after.length_squared() == 0.0 || before.dot(after) < -1.0 {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
