"""Mutation probe for the mesh revolve's axis-side refusal (#351).

Each mutant drops, weakens or misreads the check that a revolved profile lies
strictly on one side of the axis, and must turn a test red. The unit tests in
`tests/revolution.rs` cover the crossing, touching and one-sided profiles and
the compile-level test in `surface_deviation.rs` covers the typed error
reaching a caller.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/construction/construct/src/revolve.rs"
S = "crates/algorithms/construction/construct/src/sweep.rs"
TESTS = ["-p", "axiolid-construct", "--test", "revolution", "--test", "sweeps"]
TESTS_SWEEP = ["-p", "axiolid-construct", "--test", "sweeps"]
COMPILE = ["-p", "axiolid-mesh-compile", "--test", "surface_deviation", "crosses_or_touches"]

MUTANTS = [
    ('crossing check dropped', R, '    if positive && negative {', '    if false {', TESTS),
    ('touching check dropped', R, '    if on_axis && !offset {', '    if false {', TESTS),
    ('negative side read as positive', R, '            Sign::Negative => negative = true,', '            Sign::Negative => positive = true,', TESTS),
    ('positive side read as negative', R, '            Sign::Positive => positive = true,', '            Sign::Positive => negative = true,', TESTS),
    ('vertex on the axis read as a side', R, '            _ => on_axis = true,', '            _ => positive = true,', TESTS),
    ('crossing needs a vertex on the axis too', R, '    if positive && negative {', '    if positive && negative && on_axis {', TESTS),
    ('touching refused only without a crossing side', R, '    if on_axis && !offset {', '    if on_axis && !offset && !positive && !negative {', TESTS),
    ('check skipped for an in-plane axis', R, '    if axis_direction.z.abs() > tolerance.linear() * axis_direction.length() {', '    if true {', TESTS),
    ('leaning axis not skipped', R, '    if axis_direction.z.abs() > tolerance.linear() * axis_direction.length() {', '    if false {', TESTS),
    ('z noise switches the offset on', R, '    let offset = axis_origin.z.abs() > tolerance.linear();', '    let offset = axis_origin.z != 0.0;', TESTS),
    ('offset axis never offset', R, '    let offset = axis_origin.z.abs() > tolerance.linear();', '    let offset = false;', TESTS),
    ('axis direction ignored', R, '        match side_of_axis(a, axis_direction, *p) {', '        match side_of_axis(a, Vec3::Y, *p) {', TESTS),
    ('side decided against a rounded a + d', R, '    expansion_sign(&sum)', '    {\n        let b = Point2::new(a.x + d.x, a.y + d.y);\n        let c = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);\n        if c > 0.0 { Sign::Positive } else if c < 0.0 { Sign::Negative } else { Sign::Zero }\n    }', TESTS),
    ('holes not checked', R, '    for p in rings.outer.iter().chain(rings.holes.iter().flatten()) {\n        if !(p.x', '    for p in rings.outer.iter() {\n        if !(p.x', TESTS),
    ('check not called', R, '    check_axis_side(rings, axis_origin, axis_direction, tolerance)?;', '    let _ = check_axis_side;', TESTS),
    ('taper end rings not checked', S, '    crate::revolve::check_axis_side(&ends, axis_origin, axis_direction, tolerance)?;', '    let _ = &ends;', TESTS_SWEEP),
    ('taper end judged apart from the start', S, '        outer: start.outer.iter().chain(&end.outer).copied().collect(),', '        outer: start.outer.clone(),', TESTS_SWEEP),
    ('crossing refused with the wrong variant', R, '        return Err(GeomError::InvalidInput(\n            "revolved profile crosses the axis', '        return Err(GeomError::Degenerate(\n            "revolved profile crosses the axis', COMPILE),
]

probe_lib.main(MUTANTS, timeout=1200)
