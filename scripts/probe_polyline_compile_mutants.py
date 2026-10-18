"""Mutation probe for curves compiled to polylines (#311).

Each mutant under-reports the chord bound, accepts a span over budget,
chords across a corner, flattens a curve it cannot bound, misplaces or
mis-scales a placed curve, or loses a curve's closure, and must turn a
test red.
"""
import probe_lib

P = "crates/execution/compile/src/polyline.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "polyline"],
    ["-p", "axiolid-mesh-compile", "--test", "polyline_relations"],
]

MUTANTS = [
    ('bound reported as zero', P, '                worst = worst.max(span);', '                worst = worst.max(0.0 * span);', TESTS),
    ('span accepted over the budget', P, '            if span <= chord {', '            if span <= 4.0 * chord {', TESTS),
    ('corners not cut', P, '    ends.extend(cuts);', '    ends.extend(cuts.into_iter().filter(|_| false));', TESTS),
    ('uncertified 3D family flattened', P, '        if !bound::certifies_flattening3(curve) {\n            return Err(refused(uncertified3(curve)));\n        }\n        return bisected(', '        if false && !bound::certifies_flattening3(curve) {\n            return Err(refused(uncertified3(curve)));\n        }\n        return bisected(', TESTS),
    ('placed bound not grown by the stretch', P, '        let grown = crate::deviation::grown(bound, stretch, rounding);', '        let grown = bound + rounding;', TESTS),
    ('placed budget not shrunk by the stretch', P, '        let local = crate::deviation::shrunk_budget_within(chord, stretch, room);', '        let local = chord - room;', TESTS),
    ('placement not applied', P, '            .map(|point| transform.transform_point3(point))', '            .map(|point| point)', TESTS),
    ('unbounded 3D line flattened over its unit window', P, '        (Curve3::Line(_), None, _) => {', '        (Curve3::Line(_), None, _) if false => {', TESTS),
    ('closed polyline left open', P, '    if closed && points.first() != points.last() {', '    if false && closed && points.first() != points.last() {', TESTS),
    ('full turn not closed exactly', P, '    if full_turn {', '    if false && full_turn {', TESTS),
    ('numerically resolved station accepted', P, '            if !placement.exact {', '            if false && !placement.exact {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
