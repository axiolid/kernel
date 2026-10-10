"""Mutation probe for sections of two B-spline surfaces (#119, ADR 0077).

Each mutant weakens one step -- the loop-free test, the window seeding,
the Krawczyk certificates and their enclosures, the check for curves
already traced, the rates, the edge exits, the pcurve read from the solve,
the domain breaks, or the refusal of touching surfaces -- and must turn a
test red.

Accepting every crossing that cannot be proven, tangent or not, is not
listed on its own: the tests hold no regular crossing that fails its proof
unless the certificate is broken, and `partials not scaled to the box`
breaks it -- and is killed only because such crossings are refused.

Equivalent mutants, deliberately not listed: accepting a chord whose box
was not proven. Every chord the tests trace is provable, so only a march
that strayed would tell, and the tests' curves are far apart; the
certificate is pinned instead through its enclosures and exclusions,
which the tests do reach.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/parametric/nurbs/src/pair_trace.rs"
C = "crates/algorithms/parametric/nurbs/src/pair_certify.rs"
S = "crates/representations/analytic/curve/src/pair_section.rs"
E = "crates/algorithms/parametric/evaluate/src/curve.rs"
M = "crates/algorithms/query/measure/src/exact_domain.rs"
PAIR = ["-p", "axiolid-nurbs", "--test", "pair_section"]
BOOL = ["-p", "axiolid-brep-boolean", "--test", "splines"]

MUTANTS = [
    ('normal cones never compared', P, 'if angle - ha - hb > 1e-9 && (core::f64::consts::PI - angle) - ha - hb > 1e-9 {', 'if true {', PAIR),
    ('touching pair dropped, not refused', P, '        if depth > 24 {\n            return Err(PairRefusal::Unresolved);', '        if depth > 24 {\n            continue;', PAIR),
    ('windows not seeded', P, '        clipped(patches(b1).ok_or(PairRefusal::Unsupported)?, w1),', '        patches(b1).ok_or(PairRefusal::Unsupported)?,', PAIR),
    ('seeds on traced curves traced again', P, '        if boxes.iter().any(|(x, a, b)| in_chord(x, a, b, &seed)) {', '        if false {', PAIR),
    ('every box excluded', C, '    if (0..3).any(|i| !k[i].meets(x[i])) {', '    if true {', PAIR),
    ('partials not scaled to the box', C, '.map(|r| r.scale(1.0 / width))', '.map(|r| r)', PAIR),
    ('rates not per chord', S, '        let x = solve4(matrix(&j1, &j2, d), [0.0, 0.0, 0.0, d.dot(d)])?;', '        let x = solve4(matrix(&j1, &j2, d), [0.0, 0.0, 0.0, 1.0])?;', BOOL),
    ('window edges never crossed', P, '                let outside = if value == w.0[axis] {\n                    x1 < value', '                let outside = if false {\n                    x1 < value', BOOL),
    ('a seed on an edge heading out fails', P, '                    && (step.length() <= 1e-12 * scale', '                    && (false', PAIR),
    ('a tangent edge crossing refused', C, '                _ if tangent(&root, curve, enclosure) => {', '                _ if false => {', PAIR),
    ('pcurve reads the other surface', E, '                return finite2(if first { a } else { b }, "curve point");', '                return finite2(if first { b } else { a }, "curve point");', BOOL),
    ('no domain breaks at nodes', M, 'axiolid_curve::Curve3::ImplicitSection(_) | axiolid_curve::Curve3::PairSection(_)', 'axiolid_curve::Curve3::ImplicitSection(_)', BOOL),
]

probe_lib.main(MUTANTS, timeout=1200)
