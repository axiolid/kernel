"""Mutation probe for trims of trimmed, composite, offset and placed curves
compiled to polylines (#424).

Each mutant ignores a trim's sense, reads a window in the wrong parameter or
not modulo the period, accepts a trim outside its basis, reverses or places
a composite's piece wrongly, drops the certified cut error from the bound,
accepts a cut under an offset, a trim across a closed composite's joint or
a trim of a trim read as a distance, or breaks the evaluator's native span
or certified distance, and must turn a test red.
"""
import probe_lib

P = "crates/execution/compile/src/polyline.rs"
T = "crates/execution/compile/src/station/trim.rs"
D = "crates/execution/compile/src/directrix.rs"
B = "crates/execution/compile/src/station/basis.rs"
C = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
Z = "crates/algorithms/parametric/evaluate/src/station/certify.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "polyline_relation_trims"],
    ["-p", "axiolid-mesh-compile", "--test", "station_trim_of_trim"],
    ["-p", "axiolid-evaluate", "--test", "station_cut"],
]

MUTANTS = [
    # Trims of trims, offsets and placed curves.
    ('outer trim sense ignored', P, '    if !sense {\n        points.reverse();\n    }', '    if false && !sense {\n        points.reverse();\n    }', TESTS),
    ('own span used instead of the narrowed window', T, '        Some(window) => narrowed(own, carrier.period(), window, slack),', '        Some(_) => Ok(own),', TESTS),
    ('selectors not read modulo the period of a trim of a conic', T, '            narrowed((lo, hi), period, (a, b), slack)', '            narrowed((lo, hi), None, (a, b), slack)', TESTS),
    ('a trim of a closed conic loses its period', T, '            let period = inner.period();', '            let period: Option<Scalar> = None;', TESTS),
    ('window outside the trimmed curve accepted', T, '    if x < lo - slack || y > hi + slack {', '    if false && (x < lo - slack || y > hi + slack) {', TESTS),
    ('a full turn read at its seam keeps the near arc', D, '    if hi - lo >= period * (1.0 - PERIOD_SLACK) && !(s_lo && e_lo) {', '    if false && hi - lo >= period * (1.0 - PERIOD_SLACK) && !(s_lo && e_lo) {', TESTS),
    ('window slack as a length', T, '            Some(tolerance / speed)', '            Some(tolerance)', TESTS),
    ('parameter slack as a length', D, '        (Some(tolerance), Some(speed)) if speed.is_finite() && speed > 0.0 => tolerance / speed,', '        (Some(tolerance), Some(speed)) if speed.is_finite() && speed > 0.0 => tolerance + 0.0 * speed,', TESTS),
    # Trims of composites.
    ('piece sense ignored', P, '    if piece.reversed {\n        points.reverse();\n    }', '    if false && piece.reversed {\n        points.reverse();\n    }', TESTS),
    ('cut error not in the bound', P, '        let bound = worst + gap + cut;', '        let bound = worst + gap;', TESTS),
    ('cut under an offset accepted', P, '    if offset.is_some() && !(start.is_point() && end.is_point()) {', '    if false && offset.is_some() && !(start.is_point() && end.is_point()) {', TESTS),
    ('trim across a closing joint accepted', T, '            if closed && (b > a) != sense {', '            if false && closed && (b > a) != sense {', TESTS),
    ('piece offset not flipped by its sense', P, '    let law = sensed(offset, !piece.reversed);', '    let law = sensed(offset, true);', TESTS),
    ('piece placement ignored', P, '    match piece.placement {', '    match piece.placement.filter(|_| false) {', TESTS),
    # The evaluator's cut.
    ('reversed piece measured from its start', C, '                part.start = piece.end - into;', '                part.start = piece.start + into;', TESTS),
    ('native circle parameter is the measure', Z, '            return Ok(m / rate);', '            return Ok(m);', TESTS),
    ('partial piece left out of the distance interval', C, '            let length = piece_length(&part, self.tables[k].as_ref(), shares[k])?;', '            let length = DistanceInterval::point(0.0);', TESTS),
    # Stations.
    ('a trim of a trim read as a distance', B, '            let nested = is_trim_through_placements(graph, *basis);', '            let nested = false;', TESTS),
]

probe_lib.main(MUTANTS, timeout=1500)
