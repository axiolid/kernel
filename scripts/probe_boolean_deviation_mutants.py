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
     '        if near {\n            return Coverage::Boundary;',
     '        if near {\n            return Coverage::Outside;', [UNIT, MEASURED]),
    ('trimmed face parity inverted', BOOL,
     '        if self.inside(centre) {',
     '        if !self.inside(centre) {', [UNIT, MEASURED]),
    # #236: the detail claims the given operands only when nothing was
    # read within tolerance.
    ('detail claims the given operands whatever the report', BOOL,
     '    let detail = if report.is_exact() {',
     '    let detail = if true {', [MEASURED]),
    ('detail never claims the given operands', BOOL,
     '    let detail = if report.is_exact() {',
     '    let detail = if false {', [MEASURED]),
]

probe_lib.main(MUTANTS, timeout=1800)
