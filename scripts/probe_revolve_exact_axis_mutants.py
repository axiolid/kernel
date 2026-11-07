"""Mutation probe for the exact revolve of a profile touching the axis (#379).

Each mutant misreads the exact on-axis classification, the perpendicularity
of an edge ending on the axis, the admitted contacts, or the disc cap and the
mirroring of a section on the other side, and must turn a test red. The tests
are `tests/revolve_axis_exact.rs` and `tests/revolve_exact.rs`; the
compile-level case in `boolean_deviation.rs` covers a solid cylinder
reaching a boolean as a certified bound.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/construction/construct/src/revolve_contour.rs"
E = "crates/algorithms/construction/construct/src/revolve_exact.rs"
TESTS = ["-p", "axiolid-construct", "--test", "revolve_axis_exact", "--test", "revolve_exact"]

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
    ('axis edge not followed by a disc', C, '(Some(Contact::Axis), Some(Contact::Axis | Contact::Cap))', '(Some(Contact::Axis), Some(Contact::Axis))', TESTS),
    ('disc not followed by an axis edge', C, '            | (Some(Contact::Cap), Some(Contact::Axis)) => {}', '            | (Some(Contact::Cap), Some(Contact::Axis)) if false => {}', TESTS),
    ('oblique end beside an axis edge accepted', C, '            (Some(Contact::Axis), _) | (_, Some(Contact::Axis)) => {', '            (Some(Contact::Axis), _) | (_, Some(Contact::Axis)) => {\n                return Ok(());', TESTS),
    ('isolated vertex accepted', C, '            _ => {\n                return Err(unsupported(\n                    "exact revolution of a section touching', '            _ => {\n                continue;\n                #[allow(unreachable_code)]\n                return Err(unsupported(\n                    "exact revolution of a section touching', TESTS),
    ('axis edge sweeps a face', C, '            if here.on_axis && there.on_axis {\n                continue;\n            }\n', '', TESTS),
    ('disc faces the wrong way off the axis', C, '        (index, Orientation::Forward)\n', '        (index, Orientation::Reversed)\n', TESTS),
    ('disc faces the wrong way onto the axis', C, '        (next, Orientation::Reversed)\n', '        (next, Orientation::Forward)\n', TESTS),
    ('disc bounded by the axis station', C, '        slot(circles, rim)?,\n        stations[rim].radius,', '        slot(circles, index)?,\n        stations[rim].radius,', TESTS),
    ('section on the far side not mirrored', C, '    let stations: Vec<Station> = if positive {', '    let stations: Vec<Station> = if false {', TESTS),
    ('mirrored walk not reversed', C, '                let source = count - 1 - index;\n                Station {', '                let source = index;\n                Station {', TESTS),
    ('tilted axis accepted', C, ' && !(exact_direction && axis_origin.z == 0.0)', ' && !(true && axis_origin.z == 0.0)', TESTS),
    ('axis off the plane accepted', C, ' && !(exact_direction && axis_origin.z == 0.0)', ' && !(exact_direction && axis_origin.z.abs() <= epsilon)', TESTS),
    ('axis direction judged within the tolerance', E, '    let exact_direction = axis_direction.x == 0.0 && axis_direction.z == 0.0;', '    let exact_direction = axis_direction.x.abs() <= tolerance.linear() && axis_direction.z.abs() <= tolerance.linear();', TESTS),
    ('pieces not split', C, '        if stations[index].on_axis || stations[next].on_axis {\n            continue;\n        }\n        let (a, b)', '        let (a, b)', TESTS),
    ('cavity not nested in the solid', C, '    let contains = |outer: &Piece, inner: &Piece| outer.low < inner.low && inner.high < outer.high;', '    let contains = |_: &Piece, _: &Piece| false;', TESTS),
    ('solid nested in a cavity accepted', C, '    if depth.iter().any(|&d| d > 1) {', '    if false {', TESTS),
    ('arc end named an oblique edge', C, '            (Some(Contact::Arc), _) | (_, Some(Contact::Arc)) => {', '            (Some(Contact::Arc), _) | (_, Some(Contact::Arc)) if false => {', TESTS),
    ('rectangle edge on the axis refused', E, '    if left == 0.0 || right == 0.0 {', '    if false {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
