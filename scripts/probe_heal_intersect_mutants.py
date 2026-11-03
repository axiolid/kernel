"""Mutation probe for the exact triangle-triangle test in heal (#399).

`crates/algorithms/repair/heal/src/intersect.rs` decides every pair with
certified `orient3d` / `orient2d` signs only: the plane-side rejection, the
edge-against-triangle sign test and the coplanar separating-edge test. Each
mutant breaks one orientation decision, one rejection or the coplanar branch
and must turn the tests red: `tests/exact_intersection.rs` (the issue's
rotated repro, the slotted face, the intersecting and touching cases and the
seeded comparison with a rational oracle) or `tests/self_intersection.rs`.

Equivalent mutants, deliberately not listed: dropping one side of the
whole-triangle plane rejection (`separated`). The rejection is a shortcut:
a triangle wholly off the other's plane has every edge skipped as lying on
one side, and the other triangle's edges cannot reach it either, so the edge
tests give the same answer. The test that contact with the plane is not a
rejection is listed.

Runs on `scripts/probe_lib.py`: each mutant in its own process group under a
timeout that kills the group, the source always restored, and a missing or
ambiguous anchor aborts the probe. A mutant that does not compile counts as
BROKEN, not killed.
"""
import probe_lib

SRC = "crates/algorithms/repair/heal/src/intersect.rs"
EXACT = ["-p", "axiolid-heal", "--test", "exact_intersection"]
CLASSIC = ["-p", "axiolid-heal", "--test", "self_intersection"]
BOTH = [EXACT, CLASSIC]
TIMEOUT = 1800

SIDE = ("    match orient3d(a, b, c, d).sign() {\n"
        "        Some(Sign::Positive) => 1,\n"
        "        Some(Sign::Negative) => -1,\n"
        "        _ => 0,")
TURN = ("    match orient2d(a, b, c).sign() {\n"
        "        Some(Sign::Positive) => 1,\n"
        "        Some(Sign::Negative) => -1,\n"
        "        _ => 0,")

MUTANTS = [
    # orient3d side decision
    ('a vertex in the plane read as above it', SRC, SIDE,
     SIDE.replace("_ => 0,", "_ => 1,"), BOTH),
    ('the plane sides negated', SRC, SIDE,
     SIDE.replace("=> 1,", "=> 0,").replace("=> -1,", "=> 1,").replace("_ => 0,", "_ => -1,"),
     BOTH),
    # plane rejection
    ('contact with the plane rejected', SRC,
     'sides.iter().all(|&s| s > 0) || sides.iter().all(|&s| s < 0)',
     'sides.iter().all(|&s| s >= 0) || sides.iter().all(|&s| s <= 0)', BOTH),
    # coplanar dispatch
    ('coplanar only when no vertex is off the plane', SRC,
     'sides.iter().any(|&s| s != 0).then_some(sides)',
     'sides.iter().all(|&s| s != 0).then_some(sides)', BOTH),
    ('the second plane test dropped', SRC,
     '    let Some(q_sides) = plane_sides(p, q) else {\n'
     '        return coplanar_triangles_overlap(p, q);\n'
     '    };\n'
     '    if separated(q_sides) {\n'
     '        return false;\n'
     '    }',
     '    let q_sides = plane_sides(p, q).unwrap_or([0; 3]);', BOTH),
    # edge against triangle
    ('an edge with an end in the plane rejected', SRC,
     'ends[0] * ends[1] > 0 ||', 'ends[0] * ends[1] >= 0 ||', BOTH),
    ('an edge on one side accepted', SRC,
     'ends[0] * ends[1] > 0 ||', 'false ||', BOTH),
    ('an edge in the plane accepted', SRC,
     '(ends[0] == 0 && ends[1] == 0)', 'false', BOTH),
    ('the closing edge test reversed', SRC,
     'side(a, b, t[2], t[0]),', 'side(a, b, t[0], t[2]),', BOTH),
    ('the first edge test reversed', SRC,
     'side(a, b, t[0], t[1]),', 'side(a, b, t[1], t[0]),', BOTH),
    ('a point on an edge line rejected', SRC,
     '!(around.iter().any(|&s| s > 0) && around.iter().any(|&s| s < 0))',
     '!(around.iter().any(|&s| s >= 0) && around.iter().any(|&s| s <= 0))', BOTH),
    ('any one-signed test accepted', SRC,
     '!(around.iter().any(|&s| s > 0) && around.iter().any(|&s| s < 0))',
     'around.iter().all(|&s| s >= 0)', BOTH),
    ('the edges of the first triangle only', SRC,
     '    (0..3).any(|i| edge_meets_triangle(p[i], p[(i + 1) % 3], [p_sides[i], p_sides[(i + 1) % 3]], q))\n'
     '        || (0..3)',
     '    (0..3).any(|i| edge_meets_triangle(p[i], p[(i + 1) % 3], [p_sides[i], p_sides[(i + 1) % 3]], q))\n'
     '        && (0..3)', BOTH),
    ('the edges of the first triangle skipped', SRC,
     '    (0..3).any(|i| edge_meets_triangle(p[i], p[(i + 1) % 3], [p_sides[i], p_sides[(i + 1) % 3]], q))\n'
     '        || (0..3)',
     '    false\n'
     '        || (0..3)', BOTH),
    # coplanar branch
    ('the first coplanar projection always taken', SRC,
     '        if p_turn == 0 {\n            continue;\n        }',
     '        let _ = p_turn == 0;', BOTH),
    ('a flat second triangle overlapping', SRC,
     '        if q_turn == 0 {\n            return false;\n        }',
     '        if q_turn == 0 {\n            return true;\n        }', BOTH),
    ('a flat first triangle overlapping', SRC,
     '        return !has_separating_edge(pp, p_turn, qq) && !has_separating_edge(qq, q_turn, pp);\n'
     '    }\n'
     '    false',
     '        return !has_separating_edge(pp, p_turn, qq) && !has_separating_edge(qq, q_turn, pp);\n'
     '    }\n'
     '    true', BOTH),
    ('separation by the first triangle ignored', SRC,
     '!has_separating_edge(pp, p_turn, qq) && !has_separating_edge(qq, q_turn, pp)',
     '!has_separating_edge(qq, q_turn, pp)', BOTH),
    ('separation by the second triangle ignored', SRC,
     '!has_separating_edge(pp, p_turn, qq) && !has_separating_edge(qq, q_turn, pp)',
     '!has_separating_edge(pp, p_turn, qq)', BOTH),
    ('coplanar contact reported', SRC,
     '* turn <= 0)', '* turn < 0)', BOTH),
    ('every edge line a separator needs all vertices', SRC,
     '    (0..3).any(|i| {\n        other', '    (0..3).all(|i| {\n        other', BOTH),
    ('the winding of the triangle ignored', SRC,
     '* turn <= 0)', '<= 0)', BOTH),
    ('a collinear triple read as turning left', SRC, TURN,
     TURN.replace("_ => 0,", "_ => 1,"), BOTH),
    ('the coplanar turn negated', SRC, TURN,
     TURN.replace("=> 1,", "=> 0,").replace("=> -1,", "=> 1,").replace("_ => 0,", "_ => -1,"),
     BOTH),
    ('the dropped coordinate kept', SRC,
     '        _ => Point2::new(p.x, p.y),', '        _ => Point2::new(p.x, p.z),', BOTH),
]


probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
