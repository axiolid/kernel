"""Mutation probe for certified joint distances and seam-snapping windows
(#294, ADR 0082 amendment).

Each mutant breaks one step the window or the certificate rests on: in
`axiolid-evaluate`'s window module the snap itself, the window's width,
both named refusals and what counts as inside; in its certify module the
bound on a chord's excess over its arc, the midpoint rule's error term,
reading a numerical length as exact, the grade a leaning frame adds to an
offset's plan speed, the curvature's sign, the tolerance check and the
certified interval a seam is snapped by; in the composite which seams
count; in the provider and the contract the query and its default
refusal; in the conformance suite the containment check; in the compiler
the window applied at all; and in `axiolid-model` the window's validation
and its wire format (omitted at zero, refused in a 1.0 payload). Each must
turn `evaluate/tests/station_certify.rs`, `curve_evaluate_path.rs`,
`compile/tests/station_window.rs` or the model's `station`/`wire` tests
red.

Not probed, because no test can tell them apart: halving a bound the
certificate holds with room to spare (the derivative suprema are loose by
design, so tightening one by a constant still contains the truth), the
compensated summation and the rounding of the chord chain (they move a
bound by a few units in the last place), and the panel budget.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
W = "crates/algorithms/parametric/evaluate/src/station/window.rs"
C = "crates/algorithms/parametric/evaluate/src/station/certify.rs"
K = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
P = "crates/algorithms/parametric/evaluate/src/provider.rs"
T = "crates/contracts/operations/curve-evaluate/src/contract.rs"
F = "crates/contracts/operations/curve-evaluate/src/conformance.rs"
S = "crates/execution/compile/src/station.rs"
V = "crates/representations/modeling/graph/src/validation.rs"
G = "crates/representations/modeling/graph/src/wire.rs"
M = "crates/representations/modeling/graph/src/station.rs"

CERT = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_certify"]
PATH = ["-p", "axiolid-evaluate", "--all-features", "--test", "curve_evaluate_path"]
COMP = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_window"]
STATION = ["-p", "axiolid-model", "--test", "station"]
WIRE = ["-p", "axiolid-model", "--all-features", "--test", "wire"]

MUTANTS = [
    # The window.
    ('a seam within the window is not snapped onto', W,
     '    Ok(inside.unwrap_or(distance))', '    Ok(distance)', [CERT, COMP]),
    ('the window is half as wide', W,
     '    let (lo, hi) = (distance - window, distance + window);',
     '    let (lo, hi) = (distance - 0.5 * window, distance + 0.5 * window);', [CERT, COMP]),
    ('two seams in one window are not refused', W,
     '    if ambiguous {', '    if false {', [CERT, COMP]),
    ('a seam on the window\'s edge is not refused', W,
     '    if straddled {', '    if false {', [CERT, COMP]),
    ('a seam on the window\'s edge counts as inside', W,
     '        if interval.lower >= lo && interval.upper <= hi {',
     '        if interval.upper >= lo && interval.lower <= hi {', [CERT, COMP]),
    ('smooth seams are snapped onto', C,
     '        .filter(|seam| !seam.smooth)\n', '', [CERT]),
    ('a joint is snapped by the kernel\'s estimate, not its certificate', K,
     '                    distance: self.starts[k],\n                    interval: run.interval(),',
     '                    distance: self.starts[k],\n'
     '                    interval: DistanceInterval::point(self.starts[k]),', [CERT, COMP]),
    ('the compiler ignores the window', S,
     '    let distance = basis.snapped_distance(station.distance, *seam_window)?;',
     '    let distance = station.distance + 0.0 * *seam_window;', [COMP]),
    # The certificate.
    ('a chord is taken for its arc', C,
     '                chord * (1.0 + 4.0 * Scalar::EPSILON) * 2.0 * half * half / turn.cos()',
     '                0.0 * half', [CERT]),
    ('the midpoint rule claims no error', C,
     '                h * h * h / 24.0 * bend + VALUE_ROUNDING * h * speed_mid',
     '                0.0 * bend + VALUE_ROUNDING * h * speed_mid', [CERT]),
    ('a numerical offset\'s length is read as exact', C,
     '            if offset.measure_is_exact() {\n                return Ok(DistanceInterval::point(piece.length()));',
     '            if true {\n                return Ok(DistanceInterval::point(piece.length()));',
     [CERT]),
    ('a leaning frame adds no grade to the plan speed', C,
     '        let grade = elevation.filter(|_| leans && (vertical != zero || longitudinal != zero));',
     '        let grade = elevation.filter(|_| false && leans && (vertical != zero || longitudinal != zero));',
     [CERT]),
    ('the curvature\'s sign is flipped', C,
     '                .map(|k| sign * k)', '                .map(|k| -sign * k)', [CERT]),
    ('a zero tolerance is accepted', C,
     '    if tolerance.is_finite() && tolerance > 0.0 {', '    if tolerance.is_finite() {',
     [CERT]),
    # The contract.
    ('the provider does not certify', P,
     '        CompositeBasis::from_path(path)?.joint_intervals(tolerance)',
     '        CompositeBasis::from_path(path)?.joint_intervals(tolerance * 0.0 + 1.0)?;\n'
     '        Err(GeomError::InvalidInput("no".into()))', [CERT, PATH]),
    ('the default refuses as a path, not a joint query', T,
     '            input: PATH_JOINT_DISTANCES_UNSUPPORTED,',
     '            input: CURVE_PATH_UNSUPPORTED,', [PATH]),
    ('the conformance suite accepts an interval missing the truth', F,
     '    if interval.lower > truth + rounding || interval.upper < truth - rounding {',
     '    if false {', [PATH]),
    # The graph.
    ('a negative window is accepted', V,
     '    if station.seam_window < 0.0 {', '    if false {', [STATION, WIRE]),
    ('a zero window is written', M,
     '    window.to_bits() == 0', '    false', [WIRE]),
    ('a window in a 1.0 payload is read', G,
     '    (version < SEAM_WINDOW_SINCE && station.seam_window.to_bits() != 0).then_some("seam_window")',
     '    (false && station.seam_window.to_bits() != 0).then_some("seam_window")', [WIRE]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
