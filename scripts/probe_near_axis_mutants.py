"""Mutation probe for the welding of near-axis vertices in mesh revolutions (#367).

Each mutant weakens the rule that a profile vertex within the linear
tolerance of an in-plane axis is projected onto it, welded, and charged to the
chord budget, and must turn a test red. The tests are in
`tests/revolution.rs` (the crescent, the tilted axis, the 2-tol vertex, the
boundary, the deviation) and `tests/sweeps.rs` (the tapered blend).
"""
import probe_lib

R = "crates/algorithms/construction/construct/src/revolve.rs"
S = "crates/algorithms/construction/construct/src/sweep.rs"
TESTS = ["-p", "axiolid-construct", "--test", "revolution", "--test", "sweeps"]

# Not mutated: the tapered end rings left unprojected. A welded vertex is
# placed from the start station and its along-axis position is unchanged by
# the projection, so the mesh is the same.
MUTANTS = [
    ('near vertex not projected', R, '(Sign::Zero, project_onto_axis(a, axis_direction, *p), true)', '(Sign::Zero, *p, true)', TESTS),
    ('projection distance not added to the deviation', R, '                    projection = projection.max(distance);', '                    projection = projection.max(0.0);', TESTS),
    ('budget ignores the projection', R, '    let budget = budget_after_projection(tolerance, projection)?;', '    let budget = tolerance.linear();', TESTS),
    ('taper budget ignores the projection', S, 'budget_after_projection(tolerance, contact.projection + slide)?', 'budget_after_projection(tolerance, 0.0)?', TESTS),
    ('tolerance multiplied by ten', R, '    let tol = tolerance.linear();\n    if tol <= 0.0 {', '    let tol = 10.0 * tolerance.linear();\n    if tol <= 0.0 {', TESTS),
    ('tilted axis ignored', R, '        let (side, p, near) = if side != Sign::Zero && exact && !offset {', '        let (side, p, near) = if side != Sign::Zero && exact && !offset && (axis_direction.x == 0.0 || axis_direction.y == 0.0) {', TESTS),
    ('boundary distance welded', R, '    if distance + error <= tol {', '    if distance <= tol {', TESTS),
    ('boundary distance kept', R, '    } else {\n        Err(touching_unsupported())\n    }\n}\n\n/// How far apart', '    } else {\n        Ok(None)\n    }\n}\n\n/// How far apart', TESTS),
    ('distance measured along the axis', R, '    let cross = d.x * (p.y - a.y) - d.y * (p.x - a.x);', '    let cross = d.x * (p.x - a.x) + d.y * (p.y - a.y);', TESTS),
    ('distance not divided by the direction length', R, '    let distance = cross.abs() / norm;', '    let distance = cross.abs();', TESTS),
    ('taper start not projected', R, '        start: Rings {\n            outer: contact.rings.outer,', '        start: Rings {\n            outer: start.outer.clone(),', TESTS),
    ('taper projection dropped', R, '        projection: contact.projection,\n    })\n}\n\n/// What', '        projection: 0.0,\n    })\n}\n\n/// What', TESTS),
    ('flattened ring accepted', R, '            if on_axis[first..end].iter().all(|&on| on) || area == 0.0 {', '            if false {', TESTS),
    ('flattened holes not judged', R, '    for ring in core::iter::once(&rings.outer).chain(&rings.holes) {\n        let end = first', '    for ring in core::iter::once(&rings.outer) {\n        let end = first', TESTS),
    ('far side vertex not welded', R, '        let (side, p, near) = if side != Sign::Zero && exact && !offset {', '        let (side, p, near) = if side == Sign::Positive && exact && !offset {', TESTS),
    ('taper slide compared bitwise', R, '    (shift <= tolerance.linear()).then_some(shift)', '    (shift == 0.0).then_some(shift)', TESTS),
    ('taper slide unbounded', R, '    (shift <= tolerance.linear()).then_some(shift)', '    Some(shift)', TESTS),
    ('taper slide not charged', S, 'budget_after_projection(tolerance, contact.projection + slide)?', 'budget_after_projection(tolerance, contact.projection)?', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
