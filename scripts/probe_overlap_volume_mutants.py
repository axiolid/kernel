"""Mutation probe for certified overlap volumes (#183).

Each mutant narrows or skews the interval -- drops the outward rounding,
misjudges which plane is lower or which side of a line a vertex is on,
double-counts coplanar faces, loses the sign of a face -- and must fail an
exact containment test.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
O = "crates/algorithms/query/inspect/src/overlap.rs"
TESTS = ["-p", "axiolid-inspect", "--test", "overlap"]

MUTANTS = [
    ('interval not rounded outward', O, '        Self::new(lo.next_down(), hi.next_up())', '        Self::new(lo, hi)', TESTS),
    ('interval narrowed by a sliver', O, '        lower: v.lo.max(low).min(high),', '        lower: (v.lo + 1e-12).max(low).min(high),', TESTS),
    ('clip keeps the right side', O, '        if si != Sign::Negative {\n            // From', '        if si != Sign::Positive {\n            // From', TESTS),
    ('crossing sign ignores the denominator', O, '        let s = dx.mul(&ry).sub(&dy.mul(&rx)).sign()?;\n        Some(times(s, w.sign()?))', '        let s = dx.mul(&ry).sub(&dy.mul(&rx)).sign()?;\n        Some(s)', TESTS),
    ('lower plane misjudged', O, '        Some(times(difference.sign()?, w.sign()?))', '        Some(times(difference.sign()?, w.sign()?).flip())', TESTS),
    ('coplanar overlap counted twice', O, '    if signs.iter().all(|s| *s == Sign::Zero) {', '    if false && signs.iter().all(|s| *s == Sign::Zero) {', TESTS),
    ('split drops the crossing points', O, '        if opposite {\n            // Where', '        if false && opposite {\n            // Where', TESTS),
    ('face signs ignored', O, '            total = if f.up == g.up {', '            total = if true {', TESTS),
    ('inward meshes not turned', O, '        if volume.hi < 0.0 {\n            volume', '        if false && volume.hi < 0.0 {\n            volume', TESTS),
    ('self-intersection not refused', O, '        if let Some(pair) = self_intersections(mesh).first() {', '        if let Some(pair) = self_intersections(mesh).first().filter(|_| false) {', TESTS),
    ('result not clamped at zero', O, '    clamp(total, 0.0, a.volume.hi.min(b.volume.hi))', '    clamp(total, f64::NEG_INFINITY, a.volume.hi.min(b.volume.hi))', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
