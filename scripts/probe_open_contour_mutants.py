"""Mutation probe for the open-contour refusal of the mesh profile flattening (#368).

Each mutant drops, weakens or misdirects the check in `snap_joints` that a
contour's consecutive segments, and its last and first, meet within the
linear tolerance, or the exemption of an open centre-line path from the
closing joint. Each must turn a test in `tests/extrude_contour.rs` red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/construction/construct/src/profile.rs"
TESTS = ["-p", "axiolid-construct", "--test", "extrude_contour"]

MUTANTS = [
    ('gap beyond the tolerance bridged', P, '        if !near2(end, start, slack) {\n            let gap', '        if false {\n            let gap', TESTS),
    ('closing joint exempt for a closed contour', P, '        if wraps && !closed {', '        if wraps {', TESTS),
    ('open centre-line path joined to its start', P, '        if wraps && !closed {', '        if false {', TESTS),
    ('centre-line path flattened as closed', P, '    contour_points_measured(contour, chord_error, tolerance, &mut merged, false)\n}', '    contour_points_measured(contour, chord_error, tolerance, &mut merged, true)\n}', TESTS),
    ('interior gap read as the closing joint', P, '        let wraps = j <= i;', '        let wraps = true;', TESTS),
    ('closing joint read as an interior gap', P, '        let wraps = j <= i;', '        let wraps = false;', TESTS),
    ('rounding slack of the joint dropped', P, '        if !near2(end, start, slack) {', '        if !near2(end, start, tolerance.linear()) {', TESTS),
    ('gap refused with the wrong variant', P, '            return Err(GeomError::InvalidInput(if wraps {', '            return Err(GeomError::Degenerate(if wraps {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
