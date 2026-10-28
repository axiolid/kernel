"""Mutation probe for the mesh revolve's axis handling (#351, #363).

Each mutant drops, weakens or misreads the check that a revolved profile does
not cross the axis, or the meshing of a profile that touches it, and must turn
a test red. The unit tests in
`tests/revolution.rs` cover the crossing, touching and one-sided profiles and
the compile-level test in `surface_deviation.rs` covers the typed error
and the closed mesh reaching a caller.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/construction/construct/src/revolve.rs"
S = "crates/algorithms/construction/construct/src/sweep.rs"
L = "crates/algorithms/construction/construct/src/loft.rs"
TESTS = ["-p", "axiolid-construct", "--test", "revolution", "--test", "sweeps"]
TESTS_SWEEP = ["-p", "axiolid-construct", "--test", "sweeps"]
COMPILE = ["-p", "axiolid-mesh-compile", "--test", "surface_deviation", "crosses_the_axis"]

# Not mutated: the start cap indexing station 0 raw. Station 0 allocates its
# ids in order, so `id[f] == f` there and the mapping is the identity.
MUTANTS = [
    ('crossing check dropped', R, '    if positive && negative {', '    if false {', TESTS),
    ('touching flags dropped', R, '        on_axis.push(side == Sign::Zero && !offset);', '        on_axis.push(false);', TESTS),
    ('touching flagged under an offset axis', R, '        on_axis.push(side == Sign::Zero && !offset);', '        on_axis.push(side == Sign::Zero);', TESTS),
    ('axis vertex duplicated per station', L, '                    (true, Some(v)) => v,', '                    (true, Some(_)) => positions.len() as u32,', TESTS),
    ('axis vertex placed at the first station', L, '                        positions.push(if welded_here { *last_point } else { *point });', '                        positions.push(*point);', TESTS),
    ('collapsed quad emits a degenerate triangle', L, '                    (true, false) => indices.extend([a, b, d]),', '                    (true, false) => indices.extend([a, b, d, a, d, c]),', TESTS),
    ('collapsed quad keeps the wrong triangle', L, '                    (false, true) => indices.extend([a, d, c]),', '                    (false, true) => indices.extend([a, b, d]),', TESTS),
    ('axis edge emits a quad', L, '                    (true, true) => {}', '                    (true, true) => indices.extend([a, b, d, a, d, c]),', TESTS),
    ('end cap does not share the axis vertices', L, '            let end = t.map(|f| id[last + f as usize]);', '            let end = t.map(|f| (last + f as usize) as u32);', TESTS),
    ('revolve does not weld', R, 'crate::loft::loft_welded(rings, &stations, full, &on_axis)?', 'crate::loft::loft(rings, &stations, full)?', TESTS),
    ('tapered touching no longer refused', S, '        return Err(crate::revolve::touching_unsupported());', '        let _ = crate::revolve::touching_unsupported();', TESTS_SWEEP),
    ('negative side read as positive', R, '            Sign::Negative => negative = true,', '            Sign::Negative => positive = true,', TESTS),
    ('positive side read as negative', R, '            Sign::Positive => positive = true,', '            Sign::Positive => negative = true,', TESTS),
    ('vertex on the axis read as a side', R, '            _ => {}', '            _ => positive = true,', TESTS),
    ('crossing needs a vertex on the axis too', R, '    if positive && negative {', '    if positive && negative && on_axis {', TESTS),
    ('check skipped for an in-plane axis', R, '    if axis_direction.z.abs() > tolerance.linear() * axis_direction.length() {', '    if true {', TESTS),
    ('leaning axis not skipped', R, '    if axis_direction.z.abs() > tolerance.linear() * axis_direction.length() {', '    if false {', TESTS),
    ('z noise switches the offset on', R, '    let offset = axis_origin.z.abs() > tolerance.linear();', '    let offset = axis_origin.z != 0.0;', TESTS),
    ('offset axis never offset', R, '    let offset = axis_origin.z.abs() > tolerance.linear();', '    let offset = false;', TESTS),
    ('axis direction ignored', R, '        let side = side_of_axis(a, axis_direction, *p);', '        let side = side_of_axis(a, Vec3::Y, *p);', TESTS),
    ('side decided against a rounded a + d', R, '    expansion_sign(&sum)', '    {\n        let b = Point2::new(a.x + d.x, a.y + d.y);\n        let c = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);\n        if c > 0.0 { Sign::Positive } else if c < 0.0 { Sign::Negative } else { Sign::Zero }\n    }', TESTS),
    ('holes not checked', R, '    for p in rings.outer.iter().chain(rings.holes.iter().flatten()) {\n        if !(p.x', '    for p in rings.outer.iter() {\n        if !(p.x', TESTS),
    ('check not called', R, '    let on_axis = check_axis_side(rings, axis_origin, axis_direction, tolerance)?;', '    let _ = check_axis_side;\n    let on_axis = vec![false; rings.outer.len() + rings.holes.iter().map(Vec::len).sum::<usize>()];', TESTS),
    ('taper end rings not checked', S, '    if crate::revolve::check_axis_side(&ends, axis_origin, axis_direction, tolerance)?', '    if { let _ = &ends; vec![false] }', TESTS_SWEEP),
    ('taper end judged apart from the start', S, '        outer: start.outer.iter().chain(&end.outer).copied().collect(),', '        outer: start.outer.clone(),', TESTS_SWEEP),
    ('crossing refused with the wrong variant', R, '        return Err(GeomError::InvalidInput(\n            "revolved profile crosses the axis', '        return Err(GeomError::Degenerate(\n            "revolved profile crosses the axis', COMPILE),
]

probe_lib.main(MUTANTS, timeout=1200)
