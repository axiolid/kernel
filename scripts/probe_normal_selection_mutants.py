"""Mutation probe for faces selected by normal (#417).

Each mutant breaks one rule of the cone test, the area sum, or the extent,
and must turn a test of `face_selection` red. All 10 are killed.
"""
import probe_lib

F = "crates/algorithms/query/measure/src/mesh_faces.rs"
M = "crates/algorithms/query/measure/src/mesh.rs"
TESTS = ["-p", "axiolid-measure", "--test", "face_selection"]

MUTANTS = [
    ('angular slack dropped', F, 'let reach = cone.half_angle + tolerance.angular();', 'let reach = cone.half_angle;', TESTS),
    ('inclusive bound made exclusive', F, '        if angle <= reach {', '        if angle < reach {', TESTS),
    ('area factor one half dropped', M, '    (b - a).cross(c - a).length() * 0.5\n}', '    (b - a).cross(c - a).length()\n}', TESTS),
    ('axis sign ignored', F, 'normal.dot(cone.axis));', 'normal.dot(cone.axis).abs());', TESTS),
    ('degenerate gate removed', F, 'if !health.is_surface_usable() || health.degenerate_triangles != 0 {', 'if !health.is_surface_usable() {', TESTS),
    ('direction not normalised', F, 'let direction = unit(along).ok_or(FaceSelectionError::InvalidDirection)?;', 'let direction = {\n        unit(along).ok_or(FaceSelectionError::InvalidDirection)?;\n        along\n    };', TESTS),
    ('extent reads only the first corner', F, 'for corner in mesh.triangle(triangle) {', 'for corner in mesh.triangle(triangle).into_iter().take(1) {', TESTS),
    ('overflowed weight not refused', F, '        if !normal.is_finite() || !weight.is_finite() {', '        if !normal.is_finite() {', TESTS),
    ('unit vector not pre-scaled', F, '    let scaled = direction / largest;', '    let scaled = direction;', TESTS),
    ('half-angle above pi accepted', F, '(0.0..=core::f64::consts::PI).contains(&half_angle)', '(0.0..=4.0).contains(&half_angle)', TESTS),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
