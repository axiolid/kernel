"""Mutation probe for the torus and wedge primitives (#142).

Each mutant accepts what must be refused, winds or places the mesh wrongly,
or joins its faces badly, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/reference/src/primitive.rs"
TESTS = ["-p", "axiolid-reference", "--test", "torus_wedge"]

MUTANTS = [
    ('horn torus accepted', P, '    if r >= big {', '    if r > big {', TESTS),
    ('spindle torus accepted', P, '    if r >= big {', '    if r >= big && r == big {', TESTS),
    ('torus cells wound inward', P, '            idx.extend([a, b, c]);\n            idx.extend([a, c, d]);', '            idx.extend([a, c, b]);\n            idx.extend([a, d, c]);', TESTS),
    ('one torus cell half flipped', P, '            idx.extend([a, c, d]);', '            if i + j == 0 {\n                idx.extend([a, d, c]);\n            } else {\n                idx.extend([a, c, d]);\n            }', TESTS),
    ('torus seam not closed round the axis', P, '    let at = |i: usize, j: usize| ((i % n) * m + (j % m)) as u32;', '    let at = |i: usize, j: usize| (i.min(n - 1) * m + (j % m)) as u32;', TESTS),
    ('tube circle sized by the major radius', P, '            let rho = big + r * phi.cos();', '            let rho = big + big * phi.cos();', TESTS),
    ('torus grid round the axis sized by the tube', P, '    let n = segments(big + r, tol);', '    let n = segments(r, tol);', TESTS),
    ('wedge top x range reversed accepted', P, '    if x0 > x1 || y0 > y1 {', '    if y0 > y1 {', TESTS),
    ('wedge top y range reversed accepted', P, '    if x0 > x1 || y0 > y1 {', '    if x0 > x1 {', TESTS),
    ('non-finite wedge top accepted', P, '        if !value.is_finite() {', '        if false && !value.is_finite() {', TESTS),
    ('collapsed wedge top vertices kept apart', P, '        .map(|&c| match p.iter().position(|&q| q == c) {', '        .map(|&c| match p.iter().position(|&q| q == c && false) {', TESTS),
    ('wedge top face dropped', P, '        [4, 5, 6, 7],\n', '', TESTS),
    ('wedge top taken as the base extent', P, '        Point3::new(x1, y1, h),', '        Point3::new(x, y1, h),', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
