"""Mutation probe for exact ruled sections (#119, ADR 0076).

Each mutant weakens one step -- the coefficients of the substituted
quadratic, the exact span decisions, the branch evaluation, or a refusal --
and must turn a test red.

Equivalent mutant, deliberately not listed: testing the discriminant with
`!= Negative` instead of `== Positive` at a span's sample point. Samples lie
strictly between isolated roots of the same square-free polynomial, so the
discriminant is never zero there; touching surfaces are pinned instead by
`apart_and_touching_pairs_are_named` through the empty-span path.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
X = "crates/algorithms/parametric/nurbs/src/exact_surface_intersection.rs"
R = "crates/algorithms/parametric/nurbs/src/ruled_section.rs"
Q = "crates/representations/analytic/curve/src/quadric_section.rs"
T = "crates/algorithms/parametric/nurbs/src/torus_section.rs"
A = "crates/representations/analytic/curve/src/torus_section.rs"
RULED = ["-p", "axiolid-nurbs", "--test", "ruled_section"]
PAIRS = ["-p", "axiolid-nurbs", "--test", "exact_intersection"]
EVAL = ["-p", "axiolid-evaluate", "--test", "quadratic_graph"]

MUTANTS = [
    ('trig derivative sign', Q, '        -self.cos * s + self.sin * c - 2.0 * self.cos2 * s2 + 2.0 * self.sin2 * c2', '        self.cos * s + self.sin * c - 2.0 * self.cos2 * s2 + 2.0 * self.sin2 * c2', EVAL),
    ('branches swapped', Q, '        let denominator = -b - self.branch.sign() * root;', '        let denominator = -b + self.branch.sign() * root;', EVAL),
    ('span-end rounding refused', Q, '        if d < -1e-12 * scale {', '        if d < 0.0 {', RULED),
    ('squares not halved', R, '        f(&x.e0, &y.e0).add(&cc.add(&ss).mul(&half())),', '        f(&x.e0, &y.e0).add(&cc.add(&ss)),', RULED),
    ('discriminant 2ac', R, '        let discriminant = psub(&pmul(&tb, &tb), &pmul(&trig_scale_poly(&ta, 4), &tc));', '        let discriminant = psub(&pmul(&tb, &tb), &pmul(&trig_scale_poly(&ta, 2), &tc));', RULED),
    ('nappe side inverted', R, '            Some(sn == Sign::Zero || (sn == sb))', '            Some(sn == Sign::Zero || (sn != sb))', RULED),
    ('root at pi ignored', R, '        for (start, end) in spans(&breaks, &good, at_pi_ok, vanishes) {', '        for (start, end) in spans(&breaks, &good, at_pi_ok, false) {', RULED),
    ('wrap span dropped', R, '        out.push((last.2, first.2 + core::f64::consts::TAU));', '', RULED),
    ('apex rulings on the other nappe', X, '    let z = axis * slope.signum();', '    let z = -axis * slope.signum();', RULED),
    ('other cone keeps both nappes', R, '        out.push((h0, linear(&sz, a1)));', '', RULED),
    ('carrier cone keeps both nappes', R, "        out.push((constant(d(c.radius)?), constant(d(c.semi_angle.tan())?)));", '', RULED),
    ('sign of a ignored on nappes', R, "                    (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)", "                    (_, Sign::Negative)", RULED),
    ('torus angle branch swapped', A, '        let u = (b * c + s * a * root).atan2(a * c - s * b * root);', '        let u = (b * c - s * a * root).atan2(a * c - s * b * root);', EVAL),
    ('torus sphere cross term dropped', T, '                    two.mul(&small).mul(&dz).neg(),', '                    zero(),', RULED),
    ('torus pieces cross the cut', T, '        &[e.clone(), cut],', '        &[e.clone()],', RULED),
    ('degenerate frame accepted', R, '    if dot(&cross, &z).sign() == Some(Sign::Zero) {\n        return Err(refuse());\n    }', '', PAIRS),
]

probe_lib.main(MUTANTS, timeout=1200)
