"""Mutation probe for the exact revolve of a profile touching the axis (#379).

Each mutant misreads the exact on-axis classification, the perpendicularity
of an edge ending on the axis, the admitted contacts, the disc cap and the
cone's apex (its vertex, its frame, the nesting of a cavity inside it), the
sphere's pole (the centre gate, the latitudes, the sweep check) and the
quadrant snap of the arc lowering, the axis frame of an axis in the profile
plane that is not along y (#450: its placement, its orientation, the exact
on-axis and perpendicular decisions kept through the rotation), or the
mirroring of a section on the other side, and must turn a test red. The tests
are `tests/revolve_axis_exact.rs` and `tests/revolve_exact.rs`; the
compile-level case in `boolean_deviation.rs` covers a solid cylinder
reaching a boolean as a certified bound.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/construction/construct/src/revolve_contour.rs"
L = "crates/algorithms/construction/construct/src/contour_lower.rs"
E = "crates/algorithms/construction/construct/src/revolve_exact.rs"
TESTS = [
    ["-p", "axiolid-construct", "--test", "revolve_axis_exact", "--test", "revolve_exact"],
    ["-p", "axiolid-construct", "--lib", "contour_lower"],
]

# Not mutated: `side_of_axis` itself, probed with the mesh revolve in
# `probe_revolve_axis_mutants.py`.
MUTANTS = [
    ('on-axis judged within the tolerance', C, '                on_axis: side == Sign::Zero,', '                on_axis: (vertex.point.x - axis_origin.x).abs() <= epsilon,', TESTS),
    ('nothing on the axis', C, '                on_axis: side == Sign::Zero,', '                on_axis: false,', TESTS),
    ('near-axis vertex accepted', C, '        .any(|station| !station.on_axis && station.radius.abs() <= epsilon)', '        .any(|_| false)', TESTS),
    ('crossing check dropped', C, '    if positive && negative {', '    if false {', TESTS),
    ('positive side read as negative', C, '                Sign::Positive => positive = true,', '                Sign::Positive => negative = true,', TESTS),
    ('vertex on the axis read as a side', C, '                _ => {}\n            }\n            Ok(Station {', '                _ => positive = true,\n            }\n            Ok(Station {', TESTS),
    ('any edge to the axis a disc', C, '        _ if here.height == there.height => Some(Contact::Cap),', '        _ => Some(Contact::Cap),', TESTS),
    ('disc judged within the tolerance', C, '        _ if here.height == there.height => Some(Contact::Cap),', '        _ if (here.height - there.height).abs() <= 1e-9 => Some(Contact::Cap),', TESTS),
    ('arc to the axis an axis edge', C, '        _ if here.bulge != 0.0 => Some(Contact::Arc),\n', '', TESTS),
    ('axis edge only before an axis edge', C, '(Some(Contact::Axis), Some(_)) | (Some(_), Some(Contact::Axis)) => {}', '(Some(Contact::Axis), Some(Contact::Axis)) | (Some(_), Some(Contact::Axis)) => {}', TESTS),
    ('axis edge only after a disc', C, '(Some(Contact::Axis), Some(_)) | (Some(_), Some(Contact::Axis)) => {}', '(Some(Contact::Axis), Some(_)) | (Some(Contact::Cap), Some(Contact::Axis)) => {}', TESTS),
    ('isolated vertex accepted', C, '            _ => {\n                return Err(unsupported(\n                    "exact revolution of a section touching', '            _ => {\n                continue;\n                #[allow(unreachable_code)]\n                return Err(unsupported(\n                    "exact revolution of a section touching', TESTS),
    ('axis edge sweeps a face', C, '            if here.on_axis && there.on_axis {\n                continue;\n            }\n', '', TESTS),
    ('disc faces the wrong way off the axis', C, '        (index, Orientation::Forward)\n', '        (index, Orientation::Reversed)\n', TESTS),
    ('disc faces the wrong way onto the axis', C, '        (next, Orientation::Reversed)\n', '        (next, Orientation::Forward)\n', TESTS),
    ('disc bounded by the axis station', C, '        slot(circles, rim)?,\n        stations[rim].radius,', '        slot(circles, index)?,\n        stations[rim].radius,', TESTS),
    ('section on the far side not mirrored', C, '    let stations: Vec<Station> = if positive {', '    let stations: Vec<Station> = if false {', TESTS),
    ('mirrored walk not reversed', C, '                let source = count - 1 - index;\n                Station {', '                let source = index;\n                Station {', TESTS),
    ('pieces not split', C, '        if stations[index].on_axis || stations[next].on_axis {\n            continue;\n        }\n        let (a, b)', '        let (a, b)', TESTS),
    ('cavity not nested in the solid', C, '    let contains = |outer: &Piece, inner: &Piece| outer.low < inner.low && inner.high < outer.high;', '    let contains = |_: &Piece, _: &Piece| false;', TESTS),
    ('solid nested in a cavity accepted', C, '    if depth.iter().any(|&d| d > 1) {', '    if false {', TESTS),
    ('arc without an axis edge accepted', C, '            (Some(Contact::Arc), _) | (_, Some(Contact::Arc)) => {', '            (Some(Contact::Arc), _) | (_, Some(Contact::Arc)) => {\n                continue;\n                #[allow(unreachable_code)]', TESTS),
    ('near-flat cone edge accepted', C, '                && (station.height - stations[other].height).abs() <= epsilon\n', '                && false\n', TESTS),
    ('apex without a vertex', C, '.any(|edge| matches!(contact(stations, edge), Some(Contact::Cone | Contact::Arc)));', '.any(|edge| matches!(contact(stations, edge), Some(Contact::Cone | Contact::Arc)) && false);', TESTS),
    ('cone framed at its start, not its apex', C, 'let basis = if there.on_axis { there } else { here };', 'let basis = here;', TESTS),
    ("piece range without the faces' heights", C, '            for station in [index, next] {', '            for station in [index, next].into_iter().take(0) {', TESTS),
    ('two faces refused', C, '    if faces.len() < 2 {', '    if faces.len() < 3 {', TESTS),
    ('lower circle omitted', C, '    if let Some(lower) = lower {', '    if let Some(lower) = lower.filter(|_| false) {', TESTS),
    ('pole without a vertex', C, 'Some(Contact::Cone | Contact::Arc)));', 'Some(Contact::Cone)));', TESTS),
    ('arc with both ends on the axis accepted', C, '        if station.bulge != 0.0 && stations[(index + 1) % count].on_axis {', '        if false {', TESTS),
    ('centre of the arc not gated', C, '    if offset.is_nan() || offset > CENTRE_ROUNDING * scale {', '    if false {', TESTS),
    ('centre gate has no rounding', C, 'const CENTRE_ROUNDING: Scalar = 64.0 * Scalar::EPSILON;', 'const CENTRE_ROUNDING: Scalar = 0.0;', TESTS),
    ('arc across the axis accepted', C, '    if (arc.sweep - (end - start)).abs() >= core::f64::consts::PI {', '    if false {', TESTS),
    ('pole latitude sign flipped', C, '                core::f64::consts::FRAC_PI_2\n            } else {\n                -core::f64::consts::FRAC_PI_2', '                -core::f64::consts::FRAC_PI_2\n            } else {\n                core::f64::consts::FRAC_PI_2', TESTS),
    ('no quadrant snap', L, [('    if cos.abs() <= slack {', '    if false {'), ('    } else if sin.abs() <= slack {', '    } else if false {')], TESTS),
    ('quadrant snap at a loose band', L, 'pub(crate) const QUADRANT_ROUNDING: Scalar = 2.0 * Scalar::EPSILON;', 'pub(crate) const QUADRANT_ROUNDING: Scalar = 1e-6;', TESTS),
    ('quadrant snap not scaled by the parameter', L, '    let slack = QUADRANT_ROUNDING * t.abs().max(1.0);', '    let slack = QUADRANT_ROUNDING;', TESTS),
    ('general axis built in place', E, '    let placement = (axis_direction.x != 0.0).then(|| {', '    let placement = (axis_direction.x != 0.0 && false).then(|| {', TESTS),
    ('axis out of the plane accepted', E, '    if axis_direction.z != 0.0 || axis_origin.z != 0.0 {', '    if false {', TESTS),
    ('axis origin out of the plane accepted by the ring entry', C, '    if axis_origin.z != 0.0 {', '    if false {', TESTS),
    ('frame\'s normal flipped', E, '        let normal = Vec2::new(unit.y, -unit.x);', '        let normal = Vec2::new(-unit.y, unit.x);', TESTS),
    ('partial angle negated in the axis frame', E, '        (Point3::ZERO, angle)', '        (Point3::ZERO, angle.map(|angle| -angle))', TESTS),
    ('vertex on the axis not put on it', E, '            if side == Sign::Zero { 0.0 } else { across },', '            across,', TESTS),
    ('disc height not copied to the axis end', E, '            points[on].y = points[off].y;', '            let _ = (on, off);', TESTS),
    ('hair off perpendicular read as a disc', E, '        } else if points[on].y == points[off].y {', '        } else if false {', TESTS),
    ('side of a vertex not kept by the rotation', E, '        if !consistent {', '        if false {', TESTS),
    ('near-flat cone carried by its angle', C, '        if carried.is_nan() || carried.abs() > epsilon.max(8.0 * Scalar::EPSILON * scale) {', '        if false {', TESTS),
    ('slope check at the tolerance only', C, 'carried.abs() > epsilon.max(8.0 * Scalar::EPSILON * scale)', 'carried.abs() > epsilon.max(1e3 * scale)', TESTS),
    ('sphere gate blind to the coordinates\' magnitude', C, '        .max(here.magnitude)\n        .max(there.magnitude)\n', '', TESTS),
    ('magnitude of a mapped section not threaded', E, '        if !placed {\n            return 0.0;', '        if true {\n            return 0.0;', TESTS),
    ('semicircle between axis points not split', C, '        if vertex.bulge.abs() == 1.0 && on_axis', '        if false && on_axis', TESTS),
    ('joint scale of an arc blind to its parameter', L, '                |t: Scalar| circle.frame.origin.length() + circle.radius.abs() * t.abs().max(1.0);', '                |_t: Scalar| circle.frame.origin.length() + circle.radius.abs();', TESTS),
    ('rectangle edge on the axis refused', E, '    if left == 0.0 || right == 0.0 {', '    if false {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
