"""Mutation probe for trimmed, composite, offset and intrinsic curves
compiled to polylines (#317).

Each mutant traverses a trim or a segment the wrong way, offsets to the
wrong side, loses or duplicates a joint, under-reports a joint's gap,
accepts a gap, an out-of-domain trim, a corner, an offset of an offset or a
twisted space curve, chords across a break of an offset's basis or a seam
of its curvature law, or under-reports an offset's or an intrinsic curve's
chord bound, and must turn a test red.
"""
import probe_lib

P = "crates/execution/compile/src/polyline.rs"
O = "crates/algorithms/parametric/evaluate/src/bound/offset.rs"
I = "crates/algorithms/parametric/evaluate/src/bound/intrinsic.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "polyline_relations"],
    ["-p", "axiolid-mesh-compile", "--test", "polyline"],
]

MUTANTS = [
    # Trims.
    ('trim run against its sense', P, '    if !sense {\n        points.reverse();\n    }', '    if false && !sense {\n        points.reverse();\n    }', TESTS),
    ('trim offset not flipped by its sense', P, '    let offset = sensed(offset, sense);', '    let offset = sensed(offset, true);', TESTS),
    ('periodic trim read against the domain', P, '    if periodic || unbounded {', '    if unbounded {', TESTS),
    ('trim past the domain accepted', P, [('    if (start.is_finite() && lo < start - slack(start))\n', '    if false && ((start.is_finite() && lo < start - slack(start))\n'), ('        || (end.is_finite() && hi > end + slack(end))\n', '        || (end.is_finite() && hi > end + slack(end)))\n')], TESTS),
    # Composites.
    ('joint vertex duplicated', P, '                child.remove(0);', '                child.remove(0);\n                child.insert(0, last);', TESTS),
    ('segment sense ignored', P, '            if !segment.same_sense {\n                child.reverse();\n            }', '            if false && !segment.same_sense {\n                child.reverse();\n            }', TESTS),
    ('segment offset not flipped by its sense', P, '            let law = sensed(offset, segment.same_sense);', '            let law = sensed(offset, true);', TESTS),
    ('joint gap not in the bound', P, '        let bound = if gap > 0.0 { worst + gap } else { worst };', '        let bound = worst;', TESTS),
    ('gap beyond the tolerance accepted', P, '                if joint.is_nan() || joint > tolerance {', '                if false && (joint.is_nan() || joint > tolerance) {', TESTS),
    ('closing joint not shared', P, '            if closing > 0.0 && closing <= tolerance {', '            if false && closing > 0.0 && closing <= tolerance {', TESTS),
    ('no room left for the gap', P, '        if bound <= chord || budget < chord || gap >= chord {', '        if true {', TESTS),
    # Offsets.
    ('2D offset distance dropped', P, '    let distance = offset.distance;\n    bisected(\n        |t| bound::offset_point2(', '    let distance = 0.0 * offset.distance;\n    bisected(\n        |t| bound::offset_point2(', TESTS),
    ('3D offset distance dropped', P, '    let distance = offset.distance;\n    bisected(\n        |t| bound::offset_point3(', '    let distance = 0.0 * offset.distance;\n    bisected(\n        |t| bound::offset_point3(', TESTS),
    ('2D offset bounded by its basis', P, '        |a, b| Some(bound::offset_chord_bound2(curve, distance, a, b).unwrap_or(Scalar::INFINITY)),', '        |a, b| bound::chord_bound2(curve, a, b),', TESTS),
    ('offset across a corner accepted', P, '    if inside(&bound::continuity_breaks2(curve, 1), span) {', '    if false && inside(&bound::continuity_breaks2(curve, 1), span) {', TESTS),
    ('C2 breaks of an offset basis not cut', P, '        &with_seams(bound::continuity_breaks2(curve, 2), seams),', '        &with_seams(Vec::new(), seams),', TESTS),
    ('seams under an offset not cut', P, '            Some(intrinsic.curvature.seams_within(intrinsic.length))', '            None', TESTS),
    ('offset of an offset accepted', P, '            if offset.is_some() {', '            if false && offset.is_some() {', TESTS),
    ('left normal flipped', O, '    let left = unit(Vec3::new(-d.y, d.x, 0.0))?;', '    let left = unit(Vec3::new(d.y, -d.x, 0.0))?;', TESTS),
    ('3D offset direction reversed', O, '    Ok(p + distance * unit(v.cross(d))?)', '    Ok(p + distance * unit(d.cross(v))?)', TESTS),
    ('unit normal bend not in the bound', O, '    let second = bounds.second + distance.abs() * (y + 3.0 * x * x);', '    let second = bounds.second;', TESTS),
    # Intrinsic curves.
    ('2D intrinsic bound halved', I, '    let out = inflate(taylor(b - a, bounds.second));', '    let out = inflate(taylor(b - a, 0.5 * bounds.second));', TESTS),
    ('3D intrinsic bound halved', I, '    let out = inflate(taylor(hi - lo, frame_stretch3(&curve.start) * k));', '    let out = inflate(taylor(hi - lo, 0.5 * frame_stretch3(&curve.start) * k));', TESTS),
    ('intrinsic frame stretch ignored', I, '        frame_stretch2(&curve.start),', '        1.0,', TESTS),
    ('torsion accepted', P, '            if !intrinsic.is_planar()', '            if false && !intrinsic.is_planar()', TESTS),
]

probe_lib.main(MUTANTS, timeout=1500)
