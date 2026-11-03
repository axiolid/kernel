"""Mutation probe for the mesh revolve's judgement of an axis leaning out of the profile plane (#362).

Each mutant drops, weakens or misplaces a step of `check_leaning_axis`: the
exact side test against the axis's projection into the plane, the refusal
of a profile that meets the point where the axis crosses the plane, the
margin that certifies a profile clear of it, and the reading of the whole
blend of a tapered revolution. Each must turn a test red.
`tests/revolution.rs` covers the plain revolve, `tests/sweeps.rs` the
tapered one, and the compile-level test in `surface_deviation.rs` the typed
error reaching a caller.

Not mutated: the early `return Ok(())` for rings with no vertex on the
projection (the contact loop then finds nothing, the same result) and the
non-finite guard on `t` (a non-finite offset is refused by the loop's own
finiteness test).
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/construction/construct/src/revolve.rs"
TESTS = ["-p", "axiolid-construct", "--test", "revolution", "--test", "sweeps"]
COMPILE = ["-p", "axiolid-mesh-compile", "--test", "surface_deviation", "leaning_axis"]
CROSSING = '(low - margin <= 0.0 && 0.0 <= high + margin)'
BOTH = '    if both(Sign::Positive) && both(Sign::Negative) {'

MUTANTS = [
    ('leaning axis not judged', R, '        check_leaning_axis(rings, rings, axis_origin, axis_direction, tolerance)?;', '', TESTS),
    ('tapered leaning axis not judged', R, '        check_leaning_axis(start, end, axis_origin, axis_direction, tolerance)?;', '', TESTS),
    ('perpendicular axis accepted', R, '    if axis_direction.x == 0.0 && axis_direction.y == 0.0 {', '    if false {', TESTS),
    ('projection crossing not refused', R, BOTH, '    if false {', TESTS),
    ('projection crossing needs a vertex on it too', R, BOTH, '    if both(Sign::Positive) && both(Sign::Negative) && both(Sign::Zero) {', TESTS),
    ('crossing between the ends not seen', R, '    let both = |sign: Sign| sides_start.contains(&sign) || sides_end.contains(&sign);', '    let both = |sign: Sign| sides_start.contains(&sign);', TESTS),
    ('projection crossing refused with the wrong variant', R, '        return Err(GeomError::InvalidInput(\n            "revolved profile crosses the axis: its rings lie on both sides of', '        return Err(GeomError::Degenerate(\n            "revolved profile crosses the axis: its rings lie on both sides of', COMPILE),
    ('projection judged from the wrong direction', R, '    let sides_start = ring_sides(start, a, axis_direction)?;', '    let sides_start = ring_sides(start, a, Vec3::X)?;', TESTS),
    ('contact margin dropped', R, CROSSING, '(low <= 0.0 && 0.0 <= high)', TESTS),
    ('tolerance left out of the margin', R, '        tolerance.linear() + 64.0 * Scalar::EPSILON * scale', '        64.0 * Scalar::EPSILON * scale', TESTS),
    ('contact only below the point', R, CROSSING, '(low - margin <= 0.0)', TESTS),
    ('contact only above the point', R, CROSSING, '(0.0 <= high + margin)', TESTS),
    ('start vertex on the projection not judged', R, '                vec![s0],\n', '', TESTS),
    ('end vertex on the projection not judged', R, '                vec![e0],\n', '', TESTS),
    ('sliding vertex span not judged', R, '                vec![s0, e0],\n', '', TESTS),
    ('start edge span not judged', R, '                vec![s0, s1],\n', '', TESTS),
    ('end edge span not judged', R, '                vec![e0, e1],\n', '', TESTS),
    ('vertex off the projection taken as contact', R, '        (sides[i] == Sign::Zero).then(|| offset_of(points[i]))', '        Some(offset_of(points[i]))', TESTS),
    ('crossing point on the wrong side of the origin', R, '    let t = -axis_origin.z / axis_direction.z;', '    let t = axis_origin.z / axis_direction.z;', TESTS),
    ('crossing point not scaled by the projection', R, '    let reach = t * along;', '    let reach = t;', TESTS),
    ('projection direction not normalised', R, '    let unit = (axis_direction.x / along, axis_direction.y / along);', '    let unit = (axis_direction.x, axis_direction.y);', TESTS),
    ('tapered leaning axis judged as the pooled rings', R, '        check_leaning_axis(start, end, axis_origin, axis_direction, tolerance)?;', '        check_leaning_axis(start, start, axis_origin, axis_direction, tolerance)?;', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
