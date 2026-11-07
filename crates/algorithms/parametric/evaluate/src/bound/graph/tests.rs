//! The interval layer and the piece certificate of [`super`], below the
//! public chord bound, which the split into pieces would otherwise hide.

use super::*;
use axiolid_curve::Branch;

fn trig(constant: Scalar, cos: Scalar, sin: Scalar, cos2: Scalar, sin2: Scalar) -> Trig2 {
    Trig2 {
        constant,
        cos,
        sin,
        cos2,
        sin2,
    }
}

/// Whether `iv` holds `s + e`, the exact result of an operation whose
/// rounded result is `s` and whose rounding error is `e`.
fn holds(iv: Iv, s: Scalar, e: Scalar) -> bool {
    (e >= 0.0 || iv.lo < s) && (e <= 0.0 || iv.hi > s) && iv.lo <= s && s <= iv.hi
}

#[test]
fn interval_operations_round_outward() {
    let values = [0.1, 0.2, 0.3, 1.0 / 3.0, -0.7, 1e-3, 12.5, -3.1e-5];
    for &a in &values {
        for &b in &values {
            let s = a + b;
            let bb = s - a;
            let e = (a - (s - bb)) + (b - bb);
            assert!(holds(Iv::point(a).add(Iv::point(b)), s, e), "{a} + {b}");
            let p = a * b;
            let e = a.mul_add(b, -p);
            assert!(holds(Iv::point(a).mul(Iv::point(b)), p, e), "{a} * {b}");
        }
        let q = a * a;
        let e = a.mul_add(a, -q);
        assert!(holds(Iv::point(a).sq(), q, e), "{a}^2");
    }
}

#[test]
fn a_square_over_zero_starts_at_zero() {
    let iv = Iv { lo: -2.0, hi: 1.0 }.sq();
    assert_eq!(iv.lo, 0.0);
    assert!(iv.hi >= 4.0);
    assert!(Iv { lo: 1.0, hi: 3.0 }.sq().lo > 0.0);
}

#[test]
fn trigonometric_ranges_hold_their_extremes_and_every_sample() {
    for (lo, hi) in [(0.0, 0.5), (-0.3, 2.0), (1.0, 5.0), (-4.0, 4.0), (6.0, 6.1)] {
        let (c, s) = (cos_iv(lo, hi), sin_iv(lo, hi));
        for i in 0..=2000 {
            let t = lo + (hi - lo) * i as Scalar / 2000.0;
            assert!(c.lo <= t.cos() && t.cos() <= c.hi, "cos {t}");
            assert!(s.lo <= t.sin() && t.sin() <= s.hi, "sin {t}");
        }
    }
    // A turning point inside is the extreme, rounded outward.
    assert!(cos_iv(-0.1, 0.1).hi >= 1.0);
    assert!(cos_iv(3.0, 3.3).lo <= -1.0);
    assert!(sin_iv(1.4, 1.7).hi >= 1.0);
    assert!(sin_iv(4.6, 4.8).lo <= -1.0);
}

#[test]
fn enclosures_hold_the_value_and_both_derivatives() {
    let t = trig(0.4, -1.2, 0.8, 0.3, -0.6);
    for (lo, hi) in [(0.0, 0.3), (-1.0, 1.0), (2.0, 3.5)] {
        let e = trig_enclose(&t, lo, hi);
        for i in 0..=500 {
            let x = lo + (hi - lo) * i as Scalar / 500.0;
            for (iv, v) in [
                (e.value, t.value(x)),
                (e.first, t.derivative(x)),
                (e.second, t.second(x)),
            ] {
                assert!(
                    iv.lo <= v && v <= iv.hi,
                    "[{lo}, {hi}] at {x}: {v} in {iv:?}"
                );
            }
        }
    }
}

/// `w^2 = 0.5 - 0.45 cos t`, whose height dips from 0.51 at `t = +-1` to
/// 0.22 at `t = 0`: the ends alone say nothing between them.
fn dip() -> QuadraticGraph2 {
    QuadraticGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: trig(-0.5, 0.45, 0.0, 0.0, 0.0),
        branch: Branch::Plus,
    }
}

#[test]
fn a_piece_holds_a_height_that_leaves_the_hull_of_its_ends() {
    let g = dip();
    // `w'' = q''/(2w) - q'^2/(4 w^3)` at the bottom of the dip.
    let at_bottom = g.bend(0.0).unwrap().abs();
    assert!(at_bottom > 0.9, "{at_bottom}");
    if let Some(bend) = piece_enclosure(&g, -1.0, 1.0).map(|p| p.bend) {
        assert!(bend >= at_bottom, "{bend} under {at_bottom}");
    }
}

#[test]
fn a_piece_refuses_a_zero_of_the_height_derivative() {
    // `F_w = 2 w` vanishes at the branch end `t = 0` of `w^2 = 0.5 (1 - cos t)`.
    let g = QuadraticGraph2 {
        c: trig(-0.5, 0.5, 0.0, 0.0, 0.0),
        ..dip()
    };
    assert!(piece_enclosure(&g, -0.5, 0.5).is_none());
    assert!(piece_enclosure(&g, 0.0, 0.2).is_none());
}

#[test]
fn a_tiny_span_is_bounded_as_tightly_as_the_curve_bends() {
    // On a span this short the enclosures are nearly points, so the bound is
    // the sampled deviation of a parabola, `h^2/8 |w''|`, to a few percent:
    // dropping a term of the second derivative, or its weight, undershoots
    // it where that term is not cancelled by another.
    let h = 1e-3;
    let g = QuadraticGraph2 {
        a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
        b: trig(0.5, 0.0, -1.5, 0.25, 0.0),
        c: trig(-2.0, 0.75, 0.0, 0.0, 0.5),
        branch: Branch::Plus,
    };
    let mut tight = 0;
    for i in 0..36 {
        let t = -0.9 + 0.05 * i as Scalar;
        let exact = h * h / 8.0 * g.bend(t + 0.5 * h).unwrap().abs();
        let bound = certify(&g, t, t + h).unwrap();
        assert!(exact <= bound, "{t}: {exact} above {bound}");
        if bound <= 1.05 * exact + 1e-14 {
            tight += 1;
        }
    }
    assert!(tight >= 12, "only {tight} spans bounded within 5%");
}

#[test]
fn a_tiny_angle_span_is_bounded_as_tightly_as_the_curve_bends() {
    let h = 1e-3;
    for branch in [Branch::Plus, Branch::Minus] {
        let g = AngleGraph2 {
            a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
            b: trig(0.0, 0.0, 0.5, 0.1, 0.0),
            c: trig(0.3, 0.0, 0.0, 0.0, 0.2),
            branch,
        };
        let mut tight = 0;
        for i in 0..36 {
            let t = -0.9 + 0.05 * i as Scalar;
            let exact = h * h / 8.0 * g.bend(t + 0.5 * h).unwrap().abs();
            let bound = certify(&g, t, t + h).unwrap();
            assert!(exact <= bound, "{branch:?} {t}: {exact} above {bound}");
            if bound <= 1.05 * exact + 1e-14 {
                tight += 1;
            }
        }
        assert!(
            tight >= 6,
            "{branch:?}: only {tight} spans bounded within 5%"
        );
    }
}

#[test]
fn an_end_is_certified_within_the_distance_it_is_off_its_root() {
    let g = QuadraticGraph2 {
        a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
        b: trig(0.5, 0.0, -1.5, 0.25, 0.0),
        c: trig(-2.0, 0.75, 0.0, 0.0, 0.5),
        branch: Branch::Plus,
    };
    let (t, w) = (0.3, g.height(0.3).unwrap());
    let exact = end_error(&g, t, w).unwrap();
    assert!(exact < 1e-12, "{exact}");
    // A height 1e-9 off is certified no closer than that, and not wildly
    // further.
    let off = end_error(&g, t, w + 1e-9).unwrap();
    assert!((1e-9..=1e-8).contains(&off), "{off}");
}

#[test]
fn a_box_reaching_the_cut_of_the_angle_is_refused() {
    let g = AngleGraph2 {
        a: trig(1.0, 0.0, 0.0, 0.0, 0.0),
        b: Trig2::default(),
        c: Trig2::default(),
        branch: Branch::Plus,
    };
    let pi = core::f64::consts::PI;
    assert!(g.cut(Iv { lo: 3.0, hi: pi }));
    assert!(g.cut(Iv { lo: -pi, hi: -3.0 }));
    assert!(g.cut(Iv { lo: 3.0, hi: 3.3 }));
    assert!(!g.cut(Iv { lo: -3.0, hi: 3.0 }));
    assert!(!QuadraticGraph2::cut(&dip(), Iv { lo: 3.0, hi: 4.0 }));
}

/// A quadratic graph whose evaluator is off by `shift`, as one with a loose
/// closed form would be.
struct Shifted {
    graph: QuadraticGraph2,
    shift: Scalar,
}

impl Field for Shifted {
    type Piece = Coefficients;
    fn is_finite(&self) -> bool {
        self.graph.is_finite()
    }
    fn piece(&self, lo: Scalar, hi: Scalar) -> Option<Coefficients> {
        self.graph.piece(lo, hi)
    }
    fn height(&self, t: Scalar) -> Option<Scalar> {
        self.graph.height(t).map(|w| w + self.shift)
    }
    fn cut(&self, heights: Iv) -> bool {
        self.graph.cut(heights)
    }
    fn value(&self, k: &Coefficients, heights: Iv) -> Iv {
        self.graph.value(k, heights)
    }
    fn partials(&self, k: &Coefficients, heights: Iv) -> Partials {
        self.graph.partials(k, heights)
    }
}

#[test]
fn the_distance_of_the_computed_ends_from_their_roots_is_in_the_bound() {
    let graph = QuadraticGraph2 {
        a: trig(1.0, 0.2, 0.0, 0.0, 0.0),
        b: trig(0.5, 0.0, -1.5, 0.25, 0.0),
        c: trig(-2.0, 0.75, 0.0, 0.0, 0.5),
        branch: Branch::Plus,
    };
    let shift = 1e-6;
    // On a span this short the curve's own bend is below the end error, so
    // the bound is the end error to within a little.
    let bound = certify(&Shifted { graph, shift }, 0.3, 0.3001).unwrap();
    assert!((shift..2.0 * shift).contains(&bound), "{bound}");
}
