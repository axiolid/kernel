"""Mutation probe for barycentric and mean-value coordinates (#143).

Each mutant weights the wrong corner, drops a refusal, or breaks the
mean-value formula or its boundary cases, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/algorithms/query/spatial/src/barycentric.rs"
TESTS = ["-p", "axiolid-spatial", "--test", "barycentric"]
# Wachspress and discrete-harmonic coordinates share the mean-value
# preamble and its vertex/edge cases, so these anchors name the mean-value
# function (#282).
MV = ('pub fn mean_value_coordinates2(\n    polygon: &Polygon2,\n    point: Point2,\n'
      '    tolerance: Tolerance,\n) -> Result<Vec<Scalar>, BarycentricError> {\n'
      '    let v = &polygon.vertices;\n    let n = v.len();\n')
MV_GUARD = '    if !point.is_finite() || !v.iter().all(|p| p.is_finite()) {\n'
MV_N = '        return Err(BarycentricError::NonFinite);\n    }\n    if n < 3 {'
MV_R = '    let r: Vec<Scalar> = s.iter().map(|q| q.length()).collect();\n    if let Some(i) = (0..n)\n'
MV_EDGE = '        if (0.0..=1.0).contains(&t) && s[i].perp_dot(s[j]).abs() / length <= linear {\n'

MUTANTS = [
    ('2D weights of b and c swapped', B, '        area2(a, point, c) / total,\n        area2(a, b, point) / total,', '        area2(a, b, point) / total,\n        area2(a, point, c) / total,', TESTS),
    ('thin shapes not refused', B, '    if thickness > tolerance.linear() {', '    if thickness >= 0.0 {', TESTS),
    ('3D weights not normalised by the squared normal', B, '    let total = normal.dot(normal);', '    let total = normal.length();', TESTS),
    ('3D weight of b taken the wrong way round', B, '        normal.dot((c - point).cross(a - point)) / total,', '        normal.dot((a - point).cross(c - point)) / total,', TESTS),
    ('tetrahedron weight of c with the wrong sign', B, '        volume6(a, b, point, d) / total,', '        volume6(a, point, b, d) / total,', TESTS),
    ('mean-value weight uses one half-angle only', B, '        let w = (half_tangent[(i + n - 1) % n] + half_tangent[i]) / r[i];', '        let w = (half_tangent[i] + half_tangent[i]) / r[i];', TESTS),
    ('mean-value weight not divided by the distance', B, '        let w = (half_tangent[(i + n - 1) % n] + half_tangent[i]) / r[i];', '        let w = half_tangent[(i + n - 1) % n] + half_tangent[i];', TESTS),
    ('half-angle tangent from the unsigned area', B, '            s[i].perp_dot(s[j]) / (r[i] * r[j] + s[i].dot(s[j]))', '            s[i].perp_dot(s[j]).abs() / (r[i] * r[j] + s[i].dot(s[j]))', TESTS),
    ('points on a vertex not given to it', B, MV_R + '        .filter(|&i| r[i] <= linear)', MV_R + '        .filter(|&i| r[i] < 0.0)', TESTS),
    ('points on an edge not interpolated along it', B, '        if (0.0..=1.0).contains(&t) && s[i].perp_dot(s[j]).abs() / length <= linear {', '        if false && (0.0..=1.0).contains(&t) {', TESTS),
    ('edge weights swapped', B, MV_EDGE + '            weights[i] = 1.0 - t;\n            weights[j] = t;', MV_EDGE + '            weights[i] = t;\n            weights[j] = 1.0 - t;', TESTS),
    ('edges two apart never tested', B, '        for j in i + 2..last {', '        for j in i + 3..last {', TESTS),
    ('crossing polygons accepted', B, '            if segment_distance(p, q, r, s) <= linear {', '            if segment_distance(p, q, r, s) < 0.0 {', TESTS),
    ('proper crossings missed', B, '    if o1 * o2 < 0.0 && o3 * o4 < 0.0 {\n        return 0.0;\n    }', '', TESTS),
    ('short edges accepted', B, '        if (q - p).length() <= linear {', '        if (q - p).length() < 0.0 {', TESTS),
    ('too few vertices accepted', B, MV + MV_GUARD + MV_N, MV + MV_GUARD + MV_N.replace('n < 3', 'n < 2'), TESTS),
    ('non-finite points accepted', B, MV + MV_GUARD, MV + '    if !point.is_finite() {\n', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
