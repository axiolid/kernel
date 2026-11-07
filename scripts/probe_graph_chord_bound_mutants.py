"""Mutation probe for the chord bounds of quadratic-graph and angle-graph
pcurves (#414).

Each mutant breaks one step the certificate rests on: the outward rounding
of the interval layer and of the sine and cosine ranges, the height range
`rho` about the ends and its growth, the refusal where `F_w` may vanish,
the terms and weights of the implicit second derivative `w''`, the partials
of the two fields, the choice among splits, and the wiring of both families
into `chord_bound2`. Each must turn the unit tests of `bound::graph` or
`evaluate/tests/graph_chord_bound.rs` red.

The second-derivative terms enter as magnitudes (`|F_tt| + 2 |F_tw| L +
|F_ww| L^2`), so a flipped sign of the formula's terms cannot undershoot and
is not a mutant; dropping a term or its weight is. Not probed, because no
test can tell them apart: the early refusal where the discriminant (or `a^2
+ b^2 - c^2`) is certainly negative, which the `F_w` refusal subsumes, and
the added rounding of the computed ends, a few ulps of the bound; and the
factor 1/2 of the reach `L (hi - lo) / 2`, a safety margin that only a
graph climbing at its full slope bound to the middle of a piece and falling
again would use, which no smooth quadratic or angle graph does; and the
end distance `d` in the acceptance `rho > reach + d`, which only widens the
box by `d`, too little for the bound to change visibly.
"""
import probe_lib

G = "crates/algorithms/parametric/evaluate/src/bound/graph.rs"
I = "crates/algorithms/parametric/evaluate/src/bound/interval.rs"
B = "crates/algorithms/parametric/evaluate/src/bound.rs"
UNIT = ["-p", "axiolid-evaluate", "--lib", "bound::graph"]
PUBLIC = ["-p", "axiolid-evaluate", "--test", "graph_chord_bound"]

FW_CHECK = (
    "if ![p.t, p.w, p.tt, p.tw, p.ww].iter().all(|i| i.finite()) || p.w.straddles_zero() {"
)

MUTANTS = [
    # Outward rounding.
    ('widen rounds inward', I,
     '            lo: self.lo.next_down(),\n            hi: self.hi.next_up(),',
     '            lo: self.lo.next_up(),\n            hi: self.hi.next_down(),', [UNIT]),
    ('a sum is not widened', I,
     '            hi: self.hi + o.hi,\n        }\n        .widen()',
     '            hi: self.hi + o.hi,\n        }', [UNIT]),
    ('a product is not widened', I,
     '            hi: p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),\n        }\n        .widen()',
     '            hi: p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),\n        }', [UNIT]),
    ('a square is not widened above', I,
     '            hi: a.max(b).next_up(),',
     '            hi: a.max(b),', [UNIT]),
    ('a cosine range is padded inward', I,
     '    let (min, max) = cos_range(lo, hi);\n    Iv {\n        lo: min - TRIG_PAD,\n        hi: max + TRIG_PAD,',
     '    let (min, max) = cos_range(lo, hi);\n    Iv {\n        lo: min + TRIG_PAD,\n        hi: max - TRIG_PAD,', [UNIT]),
    ('a sine range is padded inward', I,
     '        lo.sin().min(hi.sin())\n    };\n    Iv {\n        lo: min - TRIG_PAD,\n        hi: max + TRIG_PAD,',
     '        lo.sin().min(hi.sin())\n    };\n    Iv {\n        lo: min + TRIG_PAD,\n        hi: max - TRIG_PAD,', [UNIT]),
    ('a sine maximum inside is not searched', I,
     '    let max = if contains_lattice(lo, hi, FRAC_PI_2, TAU) {',
     '    let max = if false {', [UNIT]),
    ('a sine minimum inside is not searched', I,
     '    let min = if contains_lattice(lo, hi, -FRAC_PI_2, TAU) {',
     '    let min = if false {', [UNIT]),
    # The height range.
    ('rho is not checked against the reach', G,
     '        if Iv::point(reach).add(Iv::point(d)).hi < rho {', '        if true {', [UNIT, PUBLIC]),
    ('rho is not grown', G,
     '        rho = (1.1 * reach).max(1.1 * rho);', '        rho = 1.1 * rho;', [UNIT, PUBLIC]),
    # The ends and the cut.
    ('the angle cut is ignored', G,
     '        if field.cut(heights) {', '        if false && field.cut(heights) {', [UNIT, PUBLIC]),
    ('the cut margin is dropped', G,
     '        heights.hi >= core::f64::consts::PI - margin\n            || heights.lo <= -core::f64::consts::PI + margin',
     '        heights.hi > core::f64::consts::PI || heights.lo < -core::f64::consts::PI', [UNIT]),
    ('an end is not certified', G,
     '        if (below.hi < 0.0 && above.lo > 0.0) || (below.lo > 0.0 && above.hi < 0.0) {',
     '        if true {', [UNIT, PUBLIC]),
    ('an end needs one sign change either way', G,
     '(below.hi < 0.0 && above.lo > 0.0) || (below.lo > 0.0 && above.hi < 0.0)',
     '(below.hi < 0.0 && above.lo > 0.0) || (below.lo > 0.0 && above.hi < 1e300)', [UNIT]),
    ('the end distance is not carried into the bound', G,
     '    Some((bend, start.max(end)))', '    Some((bend, 0.0))', [UNIT, PUBLIC]),
    # The refusal near a vertical tangent.
    ('F_w may vanish', G,
     FW_CHECK,
     'if ![p.t, p.w, p.tt, p.tw, p.ww].iter().all(|i| i.finite()) {', [UNIT, PUBLIC]),
    ('the floor of F_w is its ceiling', G,
     '        let floor = p.w.mig();', '        let floor = p.w.mag();', [UNIT, PUBLIC]),
    # The second derivative.
    ('the slope bound is not carried', G,
     '        let lipschitz = (p.t.mag() / floor).next_up();',
     '        let lipschitz = (p.t.mag() / floor).next_up() * 0.5;', [UNIT, PUBLIC]),
    ('F_tt is dropped', G,
     '            let numerator = Iv::point(p.tt.mag())',
     '            let numerator = Iv::point(0.0)', [UNIT, PUBLIC]),
    ('the cross term loses its factor 2', G,
     '                .add(Iv::point(p.tw.mag()).mul(l).scale(2.0))',
     '                .add(Iv::point(p.tw.mag()).mul(l))', [UNIT, PUBLIC]),
    ('F_ww is dropped', G,
     '                .add(Iv::point(p.ww.mag()).mul(l).mul(l));',
     '                .add(Iv::point(0.0));', [UNIT, PUBLIC]),
    ('the bend is not divided by F_w', G,
     '            let bend = (numerator.hi / floor).next_up();',
     '            let bend = numerator.hi.next_up();', [UNIT, PUBLIC]),
    ('the chord bound is span^2 / 4', G,
     '                .scale(0.125)\n                .mul(Iv::point(bend))',
     '                .scale(0.25)\n                .mul(Iv::point(bend))', [PUBLIC]),
    ('the loosest split is taken', G,
     'best.is_none_or(|b| bound < b)', 'best.is_none_or(|b| bound > b)', [PUBLIC]),
    ('a piece is not certified from its own box', G,
     '        let piece = piece_enclosure(field, a, b)?;',
     '        let piece = piece_enclosure(field, lo, hi)?;', [PUBLIC]),
    # The partials of the two fields.
    ('the quadratic F_t loses its linear term', G,
     '            t: poly(k.a.first, k.b.first, k.c.first),',
     '            t: poly(k.a.first, Iv::point(0.0), k.c.first),', [UNIT, PUBLIC]),
    ('the quadratic F_w has no constant coefficient', G,
     '            w: k.a.value.mul(w).scale(2.0).add(k.b.value),',
     '            w: k.a.value.mul(w).scale(2.0),', [UNIT, PUBLIC]),
    ('the quadratic F_tw has no b', G,
     '            tw: k.a.first.mul(w).scale(2.0).add(k.b.first),',
     '            tw: k.a.first.mul(w).scale(2.0),', [UNIT, PUBLIC]),
    ('the angle F_w adds its sine term', G,
     '        let turn = |a: Iv, b: Iv| b.mul(cw).add(a.mul(sw).neg());',
     '        let turn = |a: Iv, b: Iv| b.mul(cw).add(a.mul(sw));', [UNIT, PUBLIC]),
    ('the angle F_t keeps the right-hand side', G,
     '            t: mix(k.a.first, k.b.first).add(k.c.first.neg()),',
     '            t: mix(k.a.first, k.b.first),', [UNIT, PUBLIC]),
    ('the angle F_tt keeps the right-hand side', G,
     '            tt: mix(k.a.second, k.b.second).add(k.c.second.neg()),',
     '            tt: mix(k.a.second, k.b.second),', [UNIT, PUBLIC]),
    # The enclosure of a coefficient.
    ('a first derivative loses its sine term', G,
     '            s1.scale(-t.cos),\n            c1.scale(t.sin),',
     '            s1.scale(0.0),\n            c1.scale(t.sin),', [UNIT, PUBLIC]),
    ('a second derivative loses its double angle', G,
     '            c2.scale(-4.0 * t.cos2),\n            s2.scale(-4.0 * t.sin2),',
     '            c2.scale(-t.cos2),\n            s2.scale(-t.sin2),', [UNIT, PUBLIC]),
    # The wiring.
    ('a quadratic graph has no bound', B,
     '            graph::quadratic_chord_bound(g, lo, hi).map(inflate)',
     '            None::<Scalar>.map(inflate)', [PUBLIC]),
    ('an angle graph has no bound', B,
     '            graph::angle_chord_bound(g, lo, hi).map(inflate)',
     '            None::<Scalar>.map(inflate)', [PUBLIC]),
]

probe_lib.main(MUTANTS, timeout=1800)
