"""Mutation probe for the mesh revolve's judgement of an axis leaning out of the profile plane (#362).

Each mutant drops, weakens or misplaces a step of `check_leaning_axis`: the
exact side test against the axis's projection into the plane, the refusal
of a profile that meets the point where the axis crosses the plane, and the
margin that certifies a profile clear of it. Each must turn a test red.
`tests/revolution.rs` covers the construct level and the compile-level test
in `surface_deviation.rs` the typed error reaching a caller.

Not mutated: the early `return Ok(vec![false; count])` for a profile with no
vertex on the projection (the contact loop then finds nothing, the same
result) and the non-finite guard on `t` (a non-finite offset is refused by
the loop's own finiteness test).
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/construction/construct/src/revolve.rs"
TESTS = ["-p", "axiolid-construct", "--test", "revolution"]
COMPILE = ["-p", "axiolid-mesh-compile", "--test", "surface_deviation", "leaning_axis"]
CROSSING = '(low - margin <= 0.0 && 0.0 <= high + margin)'

MUTANTS = [
    ('leaning axis not judged', R, '        return check_leaning_axis(rings, axis_origin, axis_direction, tolerance);', '        return Ok(vec![false; count]);', TESTS),
    ('perpendicular axis accepted', R, '    if axis_direction.x == 0.0 && axis_direction.y == 0.0 {', '    if false {', TESTS),
    ('projection crossing not refused', R, '    if sides.contains(&Sign::Positive) && sides.contains(&Sign::Negative) {', '    if false {', TESTS),
    ('projection crossing needs a vertex on it too', R, '    if sides.contains(&Sign::Positive) && sides.contains(&Sign::Negative) {', '    if sides.contains(&Sign::Positive) && sides.contains(&Sign::Negative) && sides.contains(&Sign::Zero) {', TESTS),
    ('projection crossing refused with the wrong variant', R, '        return Err(GeomError::InvalidInput(\n            "revolved profile crosses the axis: its rings lie on both sides of', '        return Err(GeomError::Degenerate(\n            "revolved profile crosses the axis: its rings lie on both sides of', COMPILE),
    ('projection judged from the wrong direction', R, '    let sides = ring_sides(rings, a, axis_direction)?;', '    let sides = ring_sides(rings, a, Vec3::X)?;', TESTS),
    ('contact margin dropped', R, CROSSING, '(low <= 0.0 && 0.0 <= high)', TESTS),
    ('tolerance left out of the margin', R, '        tolerance.linear() + 64.0 * Scalar::EPSILON * scale', '        64.0 * Scalar::EPSILON * scale', TESTS),
    ('contact only below the point', R, CROSSING, '(low - margin <= 0.0)', TESTS),
    ('contact only above the point', R, CROSSING, '(0.0 <= high + margin)', TESTS),
    ('edge span ignored', R, '(here.min(other), here.max(other))', '(here, here)', TESTS),
    ('vertex off the projection taken as contact', R, '            if sides[start + k] != Sign::Zero {\n                continue;\n            }', '', TESTS),
    ('crossing point on the wrong side of the origin', R, '    let t = -axis_origin.z / axis_direction.z;', '    let t = axis_origin.z / axis_direction.z;', TESTS),
    ('crossing point not scaled by the projection', R, '    let reach = t * along;', '    let reach = t;', TESTS),
    ('projection direction not normalised', R, '    let unit = (axis_direction.x / along, axis_direction.y / along);', '    let unit = (axis_direction.x, axis_direction.y);', TESTS),
    ('leaning axis welded', R, '    Ok(vec![false; count])\n}\n\n/// Refuse rings', '    Ok(sides.iter().map(|s| *s == Sign::Zero).collect())\n}\n\n/// Refuse rings', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
