"""Mutation probe for a boolean's deviation measured against its exact
result (#235).

Each mutant loosens a step the measured bound rests on: a flat region's
off-plane term or its covering test, and a trimmed face's domain
classification. Each must turn a test red: `boolean_deviation.rs` samples
the exact results of booleans against the compiled meshes and checks the
reported bound covers them (stopping once within the budget), the
`deviation::boolean` unit test checks a bound tightened to the sampled
maximum against the cut, and the `certify` unit tests probe flat regions
directly. The branch and bound's own steps are probed by
`probe_deviation_report_mutants.py`.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
FLAT = "crates/execution/compile/src/certify/flat.rs"
BOOL = "crates/execution/compile/src/deviation/boolean.rs"
MEASURED = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"]
UNIT = ["-p", "axiolid-mesh-compile", "--lib"]
REFEREE = ["-p", "axiolid-mesh-compile", "--test", "exact_union_referee"]
WALLS = ["-p", "axiolid-mesh-compile", "--test", "coincident_wall_openings"]

MUTANTS = [
    ('flat region off-plane term dropped', FLAT,
     '            let bound = off + region.thickness;',
     '            let bound = region.thickness;', [UNIT, MEASURED]),
    ('flat region boundary edges ignored', FLAT,
     '            if inside && !touches(&quad, &region.boundary, region.margin) {',
     '            if inside {', [UNIT, MEASURED]),
    ('flat region joins folded neighbours', FLAT,
     '                    if nu.dot(normal) > 0.5 && level {',
     '                    if level {', [UNIT, MEASURED]),
    ('trimmed face drops its boundary cells', BOOL,
     '        if self.near(lo, hi).next().is_some() {\n            return Coverage::Boundary;',
     '        if self.near(lo, hi).next().is_some() {\n            return Coverage::Outside;',
     [UNIT, MEASURED]),
    ('trimmed face parity inverted', BOOL,
     '        if self.inside(centre) {',
     '        if !self.inside(centre) {', [UNIT, MEASURED]),
    # #236: the detail claims the given operands only when nothing was
    # read within tolerance.
    ('detail claims the given operands whatever the report', BOOL,
     '    } else if report.is_exact() {',
     '    } else if true {', [MEASURED]),
    ('detail never claims the given operands', BOOL,
     '    } else if report.is_exact() {',
     '    } else if false {', [MEASURED]),
    # #421: every kind of reading and the rounding floor are named, and
    # widen the boolean's own bound by their sum.
    ('readings of other kinds than coincident supports dropped', BOOL,
     '    let reading = reading_distance(report);',
     '    let reading = report\n        .decisions()\n        .iter()\n        .filter(|d| d.kind == ToleranceDecisionKind::CoincidentSupports)\n        .map(|d| d.linear)\n        .fold(0.0, Scalar::max);', [MEASURED]),
    ('readings not named', BOOL,
     '        (TOLERANCE_READINGS, reading),',
     '        (TOLERANCE_READINGS, 0.0),', [MEASURED]),
    ('rounding floor not named', BOOL,
     '        (BOOLEAN_ROUNDING_FLOOR, floor),',
     '        (BOOLEAN_ROUNDING_FLOOR, 0.0),', [MEASURED]),
    ('bound not widened by the rounding floor', BOOL,
     '            DeviationBound::Certified(widened(widened(value, floor), reading))',
     '            DeviationBound::Certified(widened(value, reading))', [MEASURED]),
    ('bound not widened by the readings', BOOL,
     '            DeviationBound::Certified(widened(widened(value, floor), reading))',
     '            DeviationBound::Certified(widened(value, floor))', [UNIT, MEASURED, WALLS]),
    ('bound the worst of the measurement and the readings', BOOL,
     '            DeviationBound::Certified(widened(widened(value, floor), reading))',
     '            DeviationBound::Certified(value.max(floor).max(reading))', [UNIT, MEASURED, WALLS]),
    ('angular readings not turned into a distance', BOOL,
     '        .map(|(linear, angular)| widened(linear, angular * lever))',
     '        .map(|(linear, _)| linear)', [UNIT]),
    ('angular readings turned by half the lever', BOOL,
     '        .map(|(linear, angular)| widened(linear, angular * lever))',
     '        .map(|(linear, angular)| widened(linear, angular * 0.5 * lever))', [UNIT]),
    ('a reading that moved and turned counts the larger only', BOOL,
     '        .map(|(linear, angular)| widened(linear, angular * lever))',
     '        .map(|(linear, angular)| linear.max(angular * lever))', [UNIT]),
    ('decisions summed across kinds', BOOL,
     '        .fold(0.0, Scalar::max)\n}',
     '        .sum()\n}', [UNIT]),
    ('a proven bound left unwidened', BOOL,
     '        DeviationBound::Certified(value) | DeviationBound::Proven(value) => {',
     '        DeviationBound::Certified(value) => {', [REFEREE]),
]

probe_lib.main(MUTANTS, timeout=1800)
