"""Mutation probe for least-squares curve and surface fitting (#150, B11).

Each mutant weakens a refusal (coincident points, invalid degree, an out of
range control count, an unmet tolerance budget) or the certified deviation
report itself, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
F = "crates/algorithms/parametric/nurbs/src/approximate.rs"
TEST = ["-p", "axiolid-nurbs", "--test", "least_squares_fit"]

MUTANTS = [
    (
        'centripetal coincidence unchecked',
        F,
        '        if step <= 0.0 {\n            return Err(GeomError::Degenerate(\n'
        '                "consecutive fitting points coincide, so centripetal spacing is undefined"\n'
        '                    .to_owned(),\n            ));\n        }',
        '',
        TEST,
    ),
    (
        'degree zero accepted',
        F,
        'fn validate_degree(degree: u16) -> FitResult<usize> {\n    if degree == 0 {',
        'fn validate_degree(degree: u16) -> FitResult<usize> {\n    if false {',
        TEST,
    ),
    (
        'control count upper bound dropped',
        F,
        '    if control_count < degree + 1 || control_count > points.len() {',
        '    if control_count < degree + 1 {',
        TEST,
    ),
    (
        'smoothing weight sign unchecked',
        F,
        '    if !options.smoothing.is_finite() || options.smoothing < 0.0 {',
        '    if !options.smoothing.is_finite() {',
        TEST,
    ),
    (
        'rank deficiency relabelled as success message only',
        F,
        '        } => GeomError::Degenerate(format!(\n'
        '            "{name} has numerical rank {rank}, needs {required} for a unique least-squares fit"\n'
        '        )),',
        '        } => GeomError::Degenerate(String::new()),',
        TEST,
    ),
    (
        'endpoint constraint dropped',
        F,
        '    if interpolate_endpoints {',
        '    if false {',
        TEST,
    ),
    (
        'tolerance budget ceiling ignored',
        F,
        '        if control_count > tolerance.max_control_points || control_count > points.len() {',
        '        if control_count > points.len() {',
        TEST,
    ),
    (
        'tolerance budget never checked against the data itself',
        F,
        '        if control_count > tolerance.max_control_points || control_count > points.len() {',
        '        if control_count > tolerance.max_control_points {',
        TEST,
    ),
    (
        'curve deviation reports zero regardless of fit',
        F,
        '        let fitted = evaluate(curve, knots, degree, tk);\n        let d = (fitted - *point).length();',
        '        let fitted = evaluate(curve, knots, degree, tk);\n        let d = 0.0 * (fitted - *point).length();',
        TEST,
    ),
    (
        'ragged surface rows unchecked',
        F,
        '        if row.len() != cols {',
        '        if false {',
        TEST,
    ),
]

probe_lib.main(MUTANTS, timeout=1200)
