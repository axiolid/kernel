"""Mutation probe for certified overlap volumes (#183, #313).

Each mutant narrows or skews the interval -- drops the outward rounding,
misjudges which plane is lower or which side of a line a vertex is on,
double-counts coplanar faces, loses the sign of a face -- and must fail an
exact containment test. The #313 mutants drop or weaken the a priori bound
on each term, which near-vertical slivers and the exported cabinet in its
room need to stay tight, or make it unsound, which containment catches.
Taking the higher of the two tops instead of the lower is left out: it is
still sound, and on a sliver it loosens the bound by under 1e-15, below
what any width test can tell apart.
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
    ('term bound dropped', O, '            let term = overlap(f, g, a.base).within(0.0, cap(f, g, a.base));', '            let term = overlap(f, g, a.base);', TESTS),
    ('term bound from the larger shadow', O, '    let area = f.area().hi.min(g.area().hi);', '    let area = f.area().hi.max(g.area().hi);', TESTS),
    ('term bound from the lowest corner', O, '        self.p[0].z.max(self.p[1].z).max(self.p[2].z)', '        self.p[0].z.min(self.p[1].z).min(self.p[2].z)', TESTS),
    ('term bound ignores the base', O, '    Iv::point(area).mul(Iv::point(top).sub(Iv::point(base))).hi', '    Iv::point(area).mul(Iv::point(top)).hi', TESTS),
    ('term not clamped at zero', O, '.within(0.0, cap(f, g, a.base));', '.within(f64::NEG_INFINITY, cap(f, g, a.base));', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
