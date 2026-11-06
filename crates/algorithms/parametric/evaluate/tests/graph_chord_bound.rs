//! Certified chord bounds of quadratic-graph and angle-graph pcurves (#414).
//!
//! Each bound is checked against a dense sample of the curve: no sampled
//! point may lie further from the chord of the span's ends than the bound,
//! and the bound must stay within a small factor of the sampled maximum, so
//! it is a usable certificate and not merely a safe one. The quadratic
//! graph's points are computed here from the textbook root formula, which
//! shares nothing with the curve's rationalised evaluation; the angle
//! graph's from the curve, each point checked against its own equation.

use axiolid_core::{Point2, Scalar};
use axiolid_curve::{AngleGraph2, Branch, Curve2, QuadraticGraph2, Trig2};
use axiolid_evaluate::bound::chord_bound2;
use axiolid_evaluate::curve::evaluate2;

const SAMPLES: usize = 4000;

/// How far above the sampled maximum a bound may sit and stay useful.
const USEFUL: Scalar = 10.0;

fn trig(constant: Scalar, cos: Scalar, sin: Scalar, cos2: Scalar, sin2: Scalar) -> Trig2 {
    Trig2 {
        constant,
        cos,
        sin,
        cos2,
        sin2,
    }
}

fn quadratic(branch: Branch) -> QuadraticGraph2 {
    QuadraticGraph2 {
        a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
        b: trig(0.5, 0.0, -1.5, 0.25, 0.0),
        c: trig(-2.0, 0.75, 0.0, 0.0, 0.5),
        branch,
    }
}

fn angle(branch: Branch) -> AngleGraph2 {
    AngleGraph2 {
        a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
        b: trig(0.0, 0.0, 0.5, 0.1, 0.0),
        c: trig(0.3, 0.0, 0.0, 0.0, 0.2),
        branch,
    }
}

/// The textbook root `(-b +- sqrt(b^2 - 4ac)) / 2a` of the quadratic at `t`.
fn root(g: &QuadraticGraph2, t: Scalar) -> Scalar {
    let (a, b, c) = (g.a.value(t), g.b.value(t), g.c.value(t));
    let d = (b * b - 4.0 * a * c).sqrt();
    (-b + g.branch.sign() * d) / (2.0 * a)
}

/// The largest distance of the sampled points from the chord of the ends.
fn sampled(point: impl Fn(Scalar) -> Point2, lo: Scalar, hi: Scalar) -> Scalar {
    let (p, q) = (point(lo), point(hi));
    let chord = q - p;
    (1..SAMPLES)
        .map(|i| {
            let t = lo + (hi - lo) * i as Scalar / SAMPLES as Scalar;
            let r = point(t) - p;
            let s = (r.dot(chord) / chord.length_squared()).clamp(0.0, 1.0);
            (r - chord * s).length()
        })
        .fold(0.0, Scalar::max)
}

fn quadratic_point(g: &QuadraticGraph2) -> impl Fn(Scalar) -> Point2 + '_ {
    move |t| Point2::new(t, root(g, t))
}

fn angle_point(g: AngleGraph2) -> impl Fn(Scalar) -> Point2 {
    move |t| {
        let p = evaluate2(&Curve2::AngleGraph(g), t).expect("defined");
        let residual = g.a.value(t) * p.x.cos() + g.b.value(t) * p.x.sin() - g.c.value(t);
        assert!(
            residual.abs() < 1e-12,
            "off its equation at {t}: {residual}"
        );
        assert_eq!(p.y, t);
        p
    }
}

/// The sample must lie within the bound, and the bound within `USEFUL`
/// times the sample or, across an inflection where the curve is nearly
/// straight and any second-order bound is far above it, within `USEFUL`
/// times the chord deviation of a bend of 1.
fn assert_tight(bound: Option<Scalar>, sample: Scalar, (a, b): (Scalar, Scalar), what: &str) {
    let bound = bound.unwrap_or_else(|| panic!("{what}: no bound"));
    assert!(
        sample <= bound,
        "{what}: sampled {sample} above the bound {bound}"
    );
    let reference = sample.max((a - b) * (a - b) / 8.0);
    assert!(
        bound <= USEFUL * reference,
        "{what}: bound {bound} is over {USEFUL}x the reference {reference}"
    );
}

const SPANS: [(Scalar, Scalar); 6] = [
    (0.0, 0.05),
    (0.2, 0.4),
    (-0.5, 0.5),
    (-1.0, 1.0),
    (0.9, 0.2),
    (-0.9, -0.8),
];

#[test]
fn a_quadratic_graph_bound_holds_the_sampled_curve_and_stays_close() {
    for branch in [Branch::Plus, Branch::Minus] {
        let g = quadratic(branch);
        for (a, b) in SPANS {
            let (lo, hi) = (a.min(b), a.max(b));
            let sample = sampled(quadratic_point(&g), lo, hi);
            let bound = chord_bound2(&Curve2::QuadraticGraph(g), a, b);
            assert_tight(bound, sample, (a, b), &format!("{branch:?} [{a}, {b}]"));
        }
    }
}

#[test]
fn an_angle_graph_bound_holds_the_sampled_curve_and_stays_close() {
    for branch in [Branch::Plus, Branch::Minus] {
        let g = angle(branch);
        for (a, b) in SPANS {
            let (lo, hi) = (a.min(b), a.max(b));
            let sample = sampled(angle_point(g), lo, hi);
            let bound = chord_bound2(&Curve2::AngleGraph(g), a, b);
            assert_tight(bound, sample, (a, b), &format!("{branch:?} [{a}, {b}]"));
        }
    }
}

#[test]
fn halving_a_span_quarters_the_bound() {
    for curve in [
        Curve2::QuadraticGraph(quadratic(Branch::Plus)),
        Curve2::AngleGraph(angle(Branch::Minus)),
    ] {
        let whole = chord_bound2(&curve, 0.0, 0.4).unwrap();
        let half = chord_bound2(&curve, 0.0, 0.2).unwrap();
        assert!(half < 0.3 * whole, "{half} against {whole}");
    }
}

#[test]
fn a_bound_over_an_empty_span_is_zero() {
    for curve in [
        Curve2::QuadraticGraph(quadratic(Branch::Plus)),
        Curve2::AngleGraph(angle(Branch::Plus)),
    ] {
        assert_eq!(chord_bound2(&curve, 0.3, 0.3), Some(0.0));
        assert_eq!(chord_bound2(&curve, 0.3, Scalar::NAN), None);
    }
}

#[test]
fn a_height_that_dips_below_its_ends_is_still_bounded_or_refused() {
    // `w^2 = 0.5 - 0.45 cos t`: both ends sit at 0.51 and the middle dips to
    // 0.22, so the ends alone say nothing about `F_w = 2 w` between them.
    let g = QuadraticGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(-0.5, 0.45, 0.0, 0.0, 0.0),
        branch: Branch::Plus,
    };
    let sample = sampled(quadratic_point(&g), -1.0, 1.0);
    assert!(sample > 0.2, "{sample}");
    if let Some(bound) = chord_bound2(&Curve2::QuadraticGraph(g), -1.0, 1.0) {
        assert!(sample <= bound, "sampled {sample} above the bound {bound}");
    }
    // A narrower span with the same shape is bounded, and holds.
    let sample = sampled(quadratic_point(&g), -0.3, 0.3);
    let bound = chord_bound2(&Curve2::QuadraticGraph(g), -0.3, 0.3).expect("bounded");
    assert!(
        sample <= bound && bound <= USEFUL * sample,
        "{sample} {bound}"
    );
}

#[test]
fn a_span_where_the_graph_is_undefined_has_no_bound() {
    // `w^2 = 0.3 - 0.5 cos t` is undefined near `t = 0`, though both ends
    // of `[-1.5, 1.5]` are defined.
    let hole = QuadraticGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(-0.3, 0.5, 0.0, 0.0, 0.0),
        branch: Branch::Plus,
    };
    assert!(hole.height(-1.5).is_some() && hole.height(1.5).is_some());
    assert_eq!(chord_bound2(&Curve2::QuadraticGraph(hole), -1.5, 1.5), None);
    // `|c| > sqrt(a^2 + b^2)` in the middle only.
    let reach = AngleGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(0.2, 1.5, 0.0, 0.0, 0.0),
        branch: Branch::Plus,
    };
    assert!(reach.angle(-1.5).is_some() && reach.angle(1.5).is_some());
    assert_eq!(chord_bound2(&Curve2::AngleGraph(reach), -1.5, 1.5), None);
    // Undefined throughout.
    let never = QuadraticGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        branch: Branch::Minus,
    };
    assert_eq!(chord_bound2(&Curve2::QuadraticGraph(never), 0.0, 1.0), None);
}

#[test]
fn a_span_near_a_vertical_tangent_has_no_bound() {
    // `w^2 = 0.5 (1 - cos t)` has a branch end at `t = 0`, where `F_w = 2 w`
    // vanishes: the graph's slope is unbounded there.
    let g = QuadraticGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(-0.5, 0.5, 0.0, 0.0, 0.0),
        branch: Branch::Plus,
    };
    let curve = Curve2::QuadraticGraph(g);
    assert_eq!(chord_bound2(&curve, -0.5, 0.5), None);
    assert_eq!(chord_bound2(&curve, 0.0, 0.5), None);
    // `u = +-acos(cos t)` meets `u = 0` at `t = 0`, where `F_u` vanishes.
    let h = AngleGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(0.0, 1.0, 0.0, 0.0, 0.0),
        branch: Branch::Plus,
    };
    let curve = Curve2::AngleGraph(h);
    assert_eq!(chord_bound2(&curve, -0.3, 0.3), None);
    assert_eq!(chord_bound2(&curve, 0.0, 0.3), None);
}
