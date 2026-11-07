"""Mutation probe for the chord bounds of lifted pcurves, pair sections and
intrinsic pcurves (#414).

Each mutant breaks one step the certificates rest on: the parameter
separation (Weyl's spread, the Gram matrix's off-diagonal, the scaling, the
crease refusal), the lifted continuation (the reach check against the box,
its growth, the end error, the terms and weights of the bend), the
derivative bounds of the section curves (the carrier's partials, the graph
axes, the implicit rate and the bridge's cubic), the pair system (the
Neumann check, the matrix perturbation, the chord length, the second
surface's curvature, the 4D residual), the intrinsic ends (the position
table's error and the quadrature's distance from it) and the wiring into
`chord_bound2` and the continuity breaks. Each must turn a test of
`axiolid-evaluate` red.

Not probed, because no test can tell them apart: the jet pad (`256 eps`) of
the separation and the matrix pad of the pair system, a few parts in `1e14`
of a bound; the end error's division by the separation where the carrier's
columns are orthogonal (the true separation is then `1`, and the certified
one only a Weyl spread under it, which no sampled deviation exposes); the
Neumann check `|I - R M'| < 1` of the pair system, which `mu > 0` subsumes
(a failed check makes `(1 - eps) / |R|` negative); the quadrature's own
distance of an intrinsic end from the exact point, about `1e-14` where the
position table's error, an eighth of the law's bound, is orders larger; and
the rounding direction of the interval layer, which `bound::graph`'s probe
(`probe_graph_chord_bound_mutants.py`) owns.
"""
import probe_lib

B = "crates/algorithms/parametric/evaluate/src/bound.rs"
D = "crates/algorithms/parametric/evaluate/src/bound/"
SEP = D + "separation.rs"
LIF = D + "lifted.rs"
CON = D + "continuation.rs"
SEC = D + "section.rs"
GRA = D + "graph.rs"
IMP = D + "implicit.rs"
PAI = D + "pair.rs"
INT = D + "intrinsic.rs"

E = ["-p", "axiolid-evaluate"]
SEPT = E + ["--test", "separation"]
LIFT = E + ["--test", "lifted_chord_bound"]
LIFU = E + ["--lib", "bound::lifted"]
SECU = E + ["--lib", "bound::section"]
PAIT = E + ["--test", "pair_chord_bound"]
PAIU = E + ["--lib", "bound::pair"]
INTT = E + ["--test", "intrinsic_chord_bound"]
BREAKS = E + ["--lib", "implicit_and_pair_sections"]

MUTANTS = [
    # The separation.
    ('the Weyl spread is dropped', SEP,
     '        let sigma = Iv::point(lambda.sqrt().lo).add(Iv::point(-jet.spread)).lo;',
     '        let sigma = lambda.sqrt().lo;', [SEPT, LIFT]),
    ('the off-diagonal of the Gram matrix is added', SEP,
     '        let lambda = Iv::point(guu.lo.min(gvv.lo)).add(Iv::point(-guv.mag()));',
     '        let lambda = Iv::point(guu.lo.min(gvv.lo)).add(Iv::point(guv.mag()));', [SEPT]),
    ('the scales are swapped', SEP,
     '        let scale = [su.length(), sv.length()];',
     '        let scale = [sv.length(), su.length()];', [SEPT]),
    ('a crease is not refused', SEP,
     '        if !bounds.smooth {\n            return None;\n        }',
     '        if false && !bounds.smooth {\n            return None;\n        }', [SEPT]),
    # The lifted continuation.
    ('rho is not checked against the reach', LIF,
     '        if reach < rho {', '        if true {', [LIFT, LIFU]),
    ('rho is not grown', LIF,
     '        rho = (1.1 * reach).max(1.1 * rho).max(Scalar::MIN_POSITIVE);',
     '        rho = Scalar::MIN_POSITIVE;', [LIFT]),
    ('the ends carry no error', LIF,
     '                start: Iv::point(start.residual).div(floor).hi,\n'
     '                end: Iv::point(finish.residual).div(floor).hi,',
     '                start: 0.0,\n                end: 0.0,', [LIFT, LIFU]),
    ('the speed is not divided by the separation', LIF,
     '        let slope = Iv::point(space.first).div(floor).hi;',
     '        let slope = Iv::point(space.first).div(Iv::point(m)).hi;', [LIFU]),
    ('the curvature along the speed is dropped', LIF,
     '        .mul(coefficient)\n        .add(Iv::point(second))',
     '        .mul(Iv::point(0.0))\n        .add(Iv::point(second))', [LIFU]),
    ('the space curve\'s own bend is dropped', LIF,
     '        .mul(coefficient)\n        .add(Iv::point(second))',
     '        .mul(coefficient)\n        .add(Iv::point(0.0))', [LIFU]),
    ('the parameters are scaled by the larger factor', LIF,
     '        let m = sep.scale[0].min(sep.scale[1]);',
     '        let m = sep.scale[0].max(sep.scale[1]);', [LIFT, LIFU]),
    ('a guide is not held to half a turn', LIF,
     '    .all(|(periodic, g_lo, g_hi, (x_lo, x_hi))| {\n        !periodic',
     '    .all(|(periodic, g_lo, g_hi, (x_lo, x_hi))| {\n        true || !periodic', [LIFT]),
    ('a pair section is read as a plain lift', LIF,
     '        return super::pair::chord_bound(section, first, lo, hi);',
     '        let _ = (section, first);', [PAIT]),
    # The continuation's combination of pieces.
    ('the first-order chord is span / 4', CON,
     '    let first = Iv::point(span).mul(Iv::point(slope)).scale(0.5).hi;',
     '    let first = Iv::point(span).mul(Iv::point(slope)).scale(0.25).hi;', [LIFU, LIFT]),
    ('the second-order chord is span^2 / 16', CON,
     '            .scale(0.125)\n            .mul(Iv::point(bend))',
     '            .scale(0.0625)\n            .mul(Iv::point(bend))', [LIFU, LIFT]),
    ('the end errors are not added', CON,
     '    Some(Iv::point(within).add(Iv::point(start.max(end))).hi)',
     '    Some(within)', [LIFU, LIFT]),
    ('the loosest split is taken', CON,
     '            if bound.is_finite() && best.is_none_or(|b| bound < b) {',
     '            if bound.is_finite() && best.is_none_or(|b| bound > b) {', [LIFT]),
    # The section curves.
    ('a ruled section has the axes of a torus section', GRA,
     '        first: [1.0, e.slope],\n        second: [0.0, e.bend],',
     '        first: [e.slope, 1.0],\n        second: [e.bend, 0.0],', [SECU]),
    ('a torus section has the axes of a ruled section', GRA,
     '        first: [e.slope, 1.0],\n        second: [e.bend, 0.0],',
     '        first: [1.0, e.slope],\n        second: [0.0, e.bend],', [SECU]),
    ('the second partial along v drops its speed term', SEC,
     '    let first = up(&[(b.du, fu), (b.dv, fv)]);',
     '    let first = up(&[(b.du, fu)]);', [SECU]),
    ('the mixed partial loses its factor 2', SEC,
     '            (2.0 * b.duv, up(&[(fu, fv)])),',
     '            (b.duv, up(&[(fu, fv)])),', [SECU]),
    ('the curve\'s own bend along v is dropped', SEC,
     '            (b.dv, sv),\n', '            (0.0, sv),\n', [SECU]),
    ('an implicit cell\'s rate is not squared', IMP,
     '        (rate * lipschitz, rate * rate * bend, range)',
     '        (rate * lipschitz, rate * bend, range)', [SECU]),
    ('a bridge\'s speed loses its factor 3', IMP,
     '        (3.0 * d1 / dt, 6.0 * d2 / (dt * dt), part.solved_range())',
     '        (d1 / dt, 6.0 * d2 / (dt * dt), part.solved_range())', [SECU]),
    ('a bridge\'s bend loses its factor 6', IMP,
     '        (3.0 * d1 / dt, 6.0 * d2 / (dt * dt), part.solved_range())',
     '        (3.0 * d1 / dt, d2 / (dt * dt), part.solved_range())', [SECU]),
    # The pair system.
    ('the matrix perturbation is dropped', PAI,
     '        let mu = Iv::point(centre).add(Iv::point(-delta)).lo;',
     '        let mu = centre;', [PAIU, PAIT]),
    ('the chord length is dropped from the speed', PAI,
     '            let speed = self.length.div(Iv::point(mu)).hi;',
     '            let speed = Iv::point(1.0).div(Iv::point(mu)).hi;', [PAIU, PAIT]),
    ('the second surface\'s curvature is dropped', PAI,
     '        .mul(q[0])\n        .add(q[1])',
     '        .mul(q[0])\n        .add(Iv::point(0.0))', [PAIU]),
    ('the residual of an end is dropped', PAI,
     '        let residual = squares.add(plane.sq()).sqrt().hi;',
     '        let residual = 0.0;', [PAIU]),
    # The intrinsic ends.
    ('the table\'s error is dropped', INT,
     '        Some(Iv::point(d.hi).add(Iv::point(error)).hi)',
     '        Some(d.hi)', [INTT]),
    # The wiring.
    ('a lifted pcurve has no bound', B,
     '            lifted::chord_bound(l, lo, hi).map(inflate)',
     '            Some(0.0).filter(|_| l.start.is_nan())', [LIFT]),
    ('a lifted pcurve is bounded by zero', B,
     '            lifted::chord_bound(l, lo, hi).map(inflate)',
     '            Some(0.0)', [LIFT]),
    ('a lifted pcurve has no breaks', B,
     '        Curve2::Lifted(l) => continuity_breaks3(&l.curve, k),',
     '        Curve2::Lifted(_) => Vec::new(),', [BREAKS]),
    ('an implicit section has no breaks', B,
     '        Curve3::ImplicitSection(s) if k >= 1 => implicit::joins(&s.curve),',
     '        Curve3::ImplicitSection(_) if k >= 1 => Vec::new(),', [BREAKS]),
    ('a pair section has no breaks', B,
     '        Curve3::PairSection(p) if k >= 1 => (1..p.nodes.len().saturating_sub(1))',
     '        Curve3::PairSection(p) if k >= 1 => (1..p.nodes.len().saturating_sub(1)).filter(|_| false)', [BREAKS]),
    ('an intrinsic pcurve has no bound', B,
     '            intrinsic::evaluated_chord_bound2(i, lo, hi)',
     '            intrinsic::evaluated_chord_bound2(i, lo, hi).filter(|_| i.length.is_nan())', [INTT]),
]

probe_lib.main(MUTANTS, timeout=1800)
