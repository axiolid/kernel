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
TESTS = ["-p", "axiolid-construct", "--test", "revolution"]
COMPILE = ["-p", "axiolid-mesh-compile", "--test", "surface_deviation", "crosses_or_touches"]

MUTANTS = [
    ('crossing check dropped', R, '    if positive && negative {', '    if false {', TESTS),
    ('touching check dropped', R, '    if on_axis {', '    if false {', TESTS),
    ('negative side read as positive', R, '            } => negative = true,', '            } => positive = true,', TESTS),
    ('positive side read as negative', R, '            } => positive = true,', '            } => negative = true,', TESTS),
    ('vertex on the axis read as a side', R, '            } => on_axis = true,', '            } => positive = true,', TESTS),
    ('crossing needs a vertex on the axis too', R, '    if positive && negative {', '    if positive && negative && on_axis {', TESTS),
    ('touching refused only without a crossing side', R, '    if on_axis {', '    if on_axis && !positive && !negative {', TESTS),
    ('check skipped for an in-plane axis', R, '    if axis_direction.z != 0.0 || axis_origin.z != 0.0 {', '    if true {', TESTS),
    ('axis direction ignored', R, '    let b = Point2::new(a.x + axis_direction.x, a.y + axis_direction.y);', '    let b = Point2::new(a.x, a.y + 1.0);', TESTS),
    ('holes not checked', R, '    for p in rings.outer.iter().chain(rings.holes.iter().flatten()) {\n        match orient2d', '    for p in rings.outer.iter() {\n        match orient2d', TESTS),
    ('check not called', R, '    check_axis_side(rings, axis_origin, axis_direction)?;', '    let _ = check_axis_side;', TESTS),
    ('crossing refused with the wrong variant', R, '        return Err(GeomError::InvalidInput(\n            "revolved profile crosses the axis', '        return Err(GeomError::Degenerate(\n            "revolved profile crosses the axis', COMPILE),
]

probe_lib.main(MUTANTS, timeout=1200)
