"""Mutation probe for the planarity of B-rep faces without a surface (#343).

Each mutant brings back the tolerance rule (a warp within the linear
tolerance read `0`), decides coplanarity in floating point or against the
wrong plane, or loses the slab, and must turn a test red: a corner lifted
`9e-4`, `1e-6` and `4e-10` m reports half its lift, and an exactly
coplanar tilted top reports `0`.
"""
import probe_lib

P = "crates/execution/compile/src/planar.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "brep_warped_faces"],
    ["-p", "axiolid-mesh-compile", "--lib", "planar::"],
]

MUTANTS = [
    ('tolerance rule again: a warp within it reads planar', P,
     '    Ok((!exactly_coplanar(rings)).then(|| spread.slab_width()))',
     '    Ok((spread.largest > 1e-3).then(|| spread.slab_width()))', TESTS),
    ('coplanarity ignored: every face reports a slab', P,
     '    Ok((!exactly_coplanar(rings)).then(|| spread.slab_width()))',
     '    Ok(Some(spread.slab_width()))', TESTS),
    ('coplanarity read from the computed slab', P,
     '    Ok((!exactly_coplanar(rings)).then(|| spread.slab_width()))',
     '    Ok((spread.above - spread.below > 0.0).then(|| spread.slab_width()))', TESTS),
    ('only the first corners are tested against the plane', P,
     '    corners\n        .iter()\n        .all(|&p| orient3d(a, b, c, p).sign() == Some(Sign::Zero))',
     '    corners\n        .iter()\n        .take(3)\n        .all(|&p| orient3d(a, b, c, p).sign() == Some(Sign::Zero))', TESTS),
    ('collinearity read from one projection', P,
     '    [(0, 1), (1, 2), (2, 0)].into_iter().all(|(i, j)| {',
     '    [(0, 1)].into_iter().all(|(i, j)| {', TESTS),
]

# The slab's rounding allowance is #261's and probed by
# probe_polygon_faces_mutants.py.
probe_lib.main(MUTANTS, timeout=2400)
