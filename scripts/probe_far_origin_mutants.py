"""Mutation probe for the exact recentring of far operands (#418, #419).

`axiolid-brep-boolean` moves operands lying far from the origin, relative
to their size, by `-c` before it cuts them, and moves the result back by
`+c` (`recentre.rs`). The move is taken only where `ExactBRep::
translated_exactly` proves every coordinate moves without rounding. Each
mutant disables the recentring, accepts a rounded translation, drops the
rounding of `c` to the ulp, the skip near the origin (per axis), the move back, or the cut of the
operands as given after a refused recentred cut,
and must turn a test red: the far column-and-corbel sweep, the recentring
unit tests, or the translation's exactness test in `axiolid-brep`.

Each mutant runs in its own process group under a timeout that kills the
whole group; the mutated file is always restored. A missing or ambiguous
anchor aborts the probe.
"""
import probe_lib

L = "crates/algorithms/construction/brep-boolean/src/lib.rs"
R = "crates/algorithms/construction/brep-boolean/src/recentre.rs"
T = "crates/representations/brep/src/transform.rs"
BOOLEAN = [
    ["-p", "axiolid-brep-boolean", "--test", "far_origin"],
    ["-p", "axiolid-brep-boolean", "--lib", "recentre::"],
]
BREP = [["-p", "axiolid-brep", "--test", "transform"]]

# Dropping `exact_sum`'s check that the sum is `a + b` is equivalent: the
# two-sum error of a rounded sum is nonzero, and of a non-finite one NaN.
MUTANTS = [
    ('recentring disabled', L,
     '.then(|| recentre::operands(a, b)).flatten();',
     '.then(|| recentre::operands(a, b).filter(|_| false)).flatten();', BOOLEAN),
    ('the two-sum error not checked', T,
     '    error == 0.0\n',
     '    error.is_finite()\n', BREP + BOOLEAN),
    ('a translation that rounds is not recorded', T,
     '                lossless.set(false);',
     '                let _ = lossless;', BREP + BOOLEAN),
    ('the centre not rounded to the ulp', R,
     '            (mid / quantum).round() * quantum',
     '            mid', BOOLEAN),
    ('an axis near zero recentred too', R,
     '        if mid.abs() <= near {',
     '        if false {', BOOLEAN),
    ('operands near the origin on every axis recentred', R,
     'if !rounded.is_finite() || rounded == Vec3::ZERO {',
     'if !rounded.is_finite() {', BOOLEAN),
    ('no cut of the operands as given after a refused recentred cut', L,
     '        session.restart(tolerance);\n',
     '        return Err(BooleanError::Undecided);\n', BOOLEAN),
    ('the session not reset for the cut as given',
     'crates/algorithms/construction/brep-boolean/src/report.rs',
     '                s.finer = false;\n                s.extends_ends = false;\n',
     '                s.finer = false;\n',
     [["-p", "axiolid-brep-boolean", "--lib", "report::"]]),
    ('the result not moved back', R,
     '        .transformed(&Transform3::from_translation(centre))\n',
     '        .transformed(&Transform3::IDENTITY)\n', BOOLEAN),
    ('the recentred centre taken from one operand only', R,
     '    let centre = centre(&[a, b])?;',
     '    let centre = centre(&[a])?;', BOOLEAN),
]

probe_lib.main(MUTANTS, timeout=1200)
