"""Mutation probe for the unit-winding certificate (#420, ADR 0099 D6).

Each mutant breaks one rule of `certify_unit_winding` and must turn a test
of `axiolid-inspect` red. All 9 are killed.
"""
import probe_lib

W = "crates/algorithms/query/inspect/src/winding.rs"
O = "crates/algorithms/query/inspect/src/overlap.rs"

TEST = ["-p", "axiolid-inspect", "--test", "winding", "--test", "overlap"]

MUTANTS = [
    ('pair sign inverted', O, 'if f.up == g.up {\n                total.add(term)', 'if f.up != g.up {\n                total.add(term)', TEST),
    ('the half dropped', W, 'second.sub(first).mul(Iv::point(0.5))', 'second.sub(first)', TEST),
    ('faces normalised as in Solid::new', W, '    let faces = faces_as_given(mesh);\n', '    let mut faces = faces_as_given(mesh);\n    if signed_volume(&faces, base).hi < 0.0 {\n        for f in &mut faces {\n            f.up = !f.up;\n        }\n    }\n', TEST),
    ('accept when the lower end fits', W, 'if excess.upper <= max_excess {', 'if excess.lower <= max_excess {', TEST),
    ('clamp lower end to -infinity', W, 'clamp(second.sub(first).mul(Iv::point(0.5)), 0.0, f64::INFINITY)', 'clamp(second.sub(first).mul(Iv::point(0.5)), f64::NEG_INFINITY, f64::INFINITY)', TEST),
    ('budget guard removed', W, '    if !max_excess.is_finite() || max_excess <= 0.0 {', '    if false && (!max_excess.is_finite() || max_excess <= 0.0) {', TEST),
    ('finite guard removed', W, '.any(|i| !mesh.position(i).is_finite())', '.any(|i| false && !mesh.position(i).is_finite())', TEST),
    ('signed volume sign flipped', O, 'volume = if f.up {\n            volume.add(term)\n        } else {\n            volume.sub(term)', 'volume = if f.up {\n            volume.sub(term)\n        } else {\n            volume.add(term)', TEST),
    ('audit skipped', W, '    if !audit.is_closed_oriented_cycle() {', '    if false && !audit.is_closed_oriented_cycle() {', TEST),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
