"""Mutation probe for exact faces whose loops close within the boolean's
readings (#334).

The deviation measurement admits a face whose pcurve loops close within the
largest linear reading the exact boolean reported, capped by the caller's
tolerance, beyond the rounding of its parameters; the slack holds the gap
measured. Each mutant breaks one of those steps -- the admission, the
reading it takes, the cap, how much it admits, and the slack -- and must
turn a test red.

Equivalent mutant, deliberately not listed: dividing the reading by `1`
instead of the surface's largest first derivative `g`. The faces a flat
wall has are planes in unit frames, `g = 1`, and no curved face of the
tests closes beyond rounding.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/execution/compile/src/deviation/boolean.rs"
UNIT = ["-p", "axiolid-mesh-compile", "--lib", "loops_close_within_the_booleans_readings"]
WALL = ["-p", "axiolid-mesh-compile", "--test", "coincident_wall_openings"]
PLANAR = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation",
          "an_all_planar_boolean_is_bounded_near_rounding"]
TESTS = [UNIT, WALL]

MUTANTS = [
    ("loops close within rounding only", B,
     "let closing = (1e-9 * scale.max(1.0)).max(reading / g);",
     "let closing = 1e-9 * scale.max(1.0);"),
    ("the caller's tolerance taken for the reading", B,
     "    report.linear().min(options.tolerance().linear())",
     "    options.tolerance().linear()"),
    ("the reading not capped by the caller's tolerance", B,
     "    report.linear().min(options.tolerance().linear())",
     "    report.linear()"),
    ("a gap admitted far beyond the reading", B,
     "let closing = (1e-9 * scale.max(1.0)).max(reading / g);",
     "let closing = (1e-9 * scale.max(1.0)).max(1e4 * reading / g);"),
    ("the slack holds the gap admitted, not the gap measured", B,
     "        Ok(Self::from_rings(&rings, delta, gap, scale, affine))",
     "        Ok(Self::from_rings(&rings, delta, closing, scale, affine))",
     [UNIT, WALL, PLANAR]),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
