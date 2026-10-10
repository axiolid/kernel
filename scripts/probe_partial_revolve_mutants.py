"""Mutation probe for exact partial-turn revolution (#172, C4).

Each mutant winds, places or trims the partial revolution wrongly, accepts
an angle or section that must be refused, or loses the axis handling, and
must turn a test red. `revolve_partial.rs` holds the Pappus volume and area
oracles, the audits and the cross-check against the mesh path.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/construction/construct/src/revolve_partial.rs"
X = "crates/algorithms/construction/construct/src/revolve_exact.rs"
TESTS = ["-p", "axiolid-construct", "--test", "revolve_partial"]

MUTANTS = [
    ('end cap used forward', P, '        rings.iter().map(|ring| ring.end.as_slice()),\n        Orientation::Reversed,', '        rings.iter().map(|ring| ring.end.as_slice()),\n        Orientation::Forward,', TESTS),
    ('start cap used reversed', P, '        rings.iter().map(|ring| ring.start.as_slice()),\n        Orientation::Forward,', '        rings.iter().map(|ring| ring.start.as_slice()),\n        Orientation::Reversed,', TESTS),
    ('sweep range reversed', P, '    let (u0, u1) = (base + angle.min(0.0), base + angle.max(0.0));', '    let (u0, u1) = (base + angle.max(0.0), base + angle.min(0.0));', TESTS),
    ('far side swept from the near side', P, '    let base = if side > 0.0 { 0.0 } else { PI };', '    let base = 0.0;', TESTS),
    ('far-side arcs not mirrored', P, '        let bulge = side * vertex.bulge;', '        let bulge = vertex.bulge;', TESTS),
    ('downward axis read as upward', X, '    let angle = angle.map(|angle| if axis.y < 0.0 { -angle } else { angle });', '    let angle = angle.map(|angle| angle);', TESTS),
    ('radial turns left-handed', P, '        Vec3::new(cos, 0.0, -sin)', '        Vec3::new(cos, 0.0, sin)', TESTS),
    ('vertex near the axis kept off it', P, '        let radius = if radius <= epsilon { 0.0 } else { radius };', '        let radius = radius.max(0.0);', TESTS),
    ('segment on the axis given a wall', P, '            if here.bulge == 0.0 && here.radius == 0.0 && there.radius == 0.0 {', '            if false {', TESTS),
    ('hole walked like the outer ring', P, '    let keep = (area > 0.0) == outer;', '    let keep = area > 0.0;', TESTS),
    ('cap arc trimmed backwards', P, '                (curve3, curve2, Interval::new(start, end))', '                (curve3, curve2, Interval::new(end, start))', TESTS),
    ('wall arc trimmed to a full turn', P, '            builder.set_edge_interval(arc, Interval::new(u0, u1));', '            builder.set_edge_interval(arc, Interval::new(u0, u0 + TAU));', TESTS),
    ('torus tube centred on the axis side', P, '                    major_radius: centre.x,', '                    major_radius: centre.x - radius,', TESTS),
    ('cone slope inverted', P, '                    semi_angle: ((there.radius - here.radius) / (there.height - here.height))', '                    semi_angle: ((here.radius - there.radius) / (there.height - here.height))', TESTS),
    ('sub-tolerance sweep accepted', P, '    if reach * angle.abs() <= epsilon {', '    if false {', TESTS),
    ('crossing section accepted', P, '    if min < -epsilon && max > epsilon {', '    if false {', TESTS),
    ('tube reaching the axis accepted', P, '                if radius >= centre.x.abs() - Scalar::EPSILON {', '                if false {', TESTS),
    ('turn beyond a full turn not named', X, '    if angle.abs() > TAU + tolerance.linear() {', '    if false {', TESTS),
    ('partial turn built as a full turn', X, '        return revolve_partial(profile, axis_origin, axis_direction, angle, tolerance);', '        let _ = revolve_partial;', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
