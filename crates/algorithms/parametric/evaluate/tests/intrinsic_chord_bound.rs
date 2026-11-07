//! Certified chord bounds of intrinsic curves between their evaluated ends
//! (#414).
//!
//! The evaluator reads an intrinsic curve's points by quadrature, so the
//! chord of its ends is not the chord of the exact ends. The bound must
//! cover every exact point's distance from the chord of the *evaluated*
//! ends, and stay close to it. The exact points come from closed forms: a
//! clothoid from the Fresnel integrals' power series, an arc, and a
//! clothoid running into an arc across a seam of the law.

use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{CurvatureLaw, Curve2, Intrinsic2};
use axiolid_evaluate::bound::{chord_bound2, continuity_breaks2, intrinsic_chord_bound2};
use axiolid_evaluate::curve::evaluate2;

const SAMPLES: usize = 2000;

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    }
}

fn place(local: Vec2) -> Point2 {
    let f = frame();
    f.origin + f.x * local.x + f.y * local.y
}

/// `int_0^s (cos(a u^2), sin(a u^2)) du` by its power series.
fn fresnel(a: Scalar, s: Scalar) -> Vec2 {
    let (mut x, mut y) = (0.0, 0.0);
    let mut factorial = 1.0;
    for n in 0..40 {
        if n > 0 {
            factorial *= n as Scalar;
        }
        let value = a.powi(n) * s.powi(2 * n + 1) / (factorial * (2 * n + 1) as Scalar);
        match n % 4 {
            0 => x += value,
            1 => y += value,
            2 => x -= value,
            _ => y -= value,
        }
    }
    Vec2::new(x, y)
}

fn rotate(angle: Scalar, v: Vec2) -> Vec2 {
    let (s, c) = angle.sin_cos();
    Vec2::new(c * v.x - s * v.y, s * v.x + c * v.y)
}

/// The distance of the exact points from the chord of the evaluated ends.
fn deviation(curve: &Curve2, exact: &impl Fn(Scalar) -> Point2, lo: Scalar, hi: Scalar) -> Scalar {
    let (p, q) = (evaluate2(curve, lo).unwrap(), evaluate2(curve, hi).unwrap());
    let chord = q - p;
    (0..=SAMPLES)
        .map(|i| {
            let s = lo + (hi - lo) * i as Scalar / SAMPLES as Scalar;
            let r = exact(s) - p;
            let t = (r.dot(chord) / chord.length_squared()).clamp(0.0, 1.0);
            (r - chord * t).length()
        })
        .fold(0.0, Scalar::max)
}

/// Sound and close: the bound holds the deviation, and is within ten times
/// it (or a nanometre).
fn check(curve: &Curve2, exact: impl Fn(Scalar) -> Point2, spans: &[(Scalar, Scalar)]) {
    for &(lo, hi) in spans {
        let sample = deviation(curve, &exact, lo, hi);
        let bound = chord_bound2(curve, lo, hi).unwrap_or_else(|| panic!("[{lo}, {hi}]: none"));
        assert!(sample <= bound, "[{lo}, {hi}]: {sample} above {bound}");
        assert!(
            bound <= 10.0 * sample + 1e-9,
            "[{lo}, {hi}]: {bound} over ten times {sample}"
        );
        // Reversed ends bound the same chord.
        assert_eq!(chord_bound2(curve, hi, lo), Some(bound));
        // The law's own bound for the exact ends is no larger.
        let Curve2::Intrinsic(i) = curve else {
            unreachable!()
        };
        assert!(intrinsic_chord_bound2(i, lo, hi).unwrap() <= bound);
    }
}

#[test]
fn a_clothoid_against_the_fresnel_integrals() {
    let (l, r) = (60.0, 300.0);
    let curve = Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::clothoid(0.0, 1.0 / r, l),
        l,
    ));
    let exact = |s: Scalar| place(fresnel(0.5 / (r * l), s));
    check(
        &curve,
        exact,
        &[(0.0, 60.0), (10.0, 25.0), (40.0, 41.0), (0.0, 0.5)],
    );
}

#[test]
fn a_constant_arc_in_closed_form() {
    let r = 12.0;
    let curve = Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::Constant { curvature: 1.0 / r },
        40.0,
    ));
    let exact = |s: Scalar| {
        let a = s / r;
        place(Vec2::new(r * a.sin(), r * (1.0 - a.cos())))
    };
    check(
        &curve,
        exact,
        &[(0.0, 40.0), (5.0, 12.0), (30.0, 30.5), (0.0, 0.001)],
    );
}

#[test]
fn a_clothoid_into_an_arc_across_the_seam() {
    let (l, r, arc) = (60.0, 300.0, 80.0);
    let law = CurvatureLaw::piecewise(
        vec![l],
        vec![
            CurvatureLaw::clothoid(0.0, 1.0 / r, l),
            CurvatureLaw::Constant { curvature: 1.0 / r },
        ],
    );
    let curve = Curve2::Intrinsic(Intrinsic2::new(frame(), law, l + arc));
    let turn = 0.5 * l / r;
    let exact = move |s: Scalar| {
        if s <= l {
            place(fresnel(0.5 / (r * l), s))
        } else {
            let a = (s - l) / r;
            let end = fresnel(0.5 / (r * l), l);
            place(end + rotate(turn, Vec2::new(r * a.sin(), r * (1.0 - a.cos()))))
        }
    };
    // The seam is not a break: the heading is continuous across it.
    assert!(continuity_breaks2(&curve, 1).is_empty());
    check(
        &curve,
        exact,
        &[(50.0, 70.0), (0.0, 140.0), (59.0, 61.0), (61.0, 100.0)],
    );
}

#[test]
fn a_span_outside_the_length_or_without_a_law_bound_has_none() {
    let curve = Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::Constant { curvature: 0.1 },
        10.0,
    ));
    assert!(chord_bound2(&curve, -1.0, 5.0).is_none());
    assert!(chord_bound2(&curve, 5.0, 10.5).is_none());
    assert!(chord_bound2(&curve, 5.0, Scalar::NAN).is_none());
    assert_eq!(chord_bound2(&curve, 3.0, 3.0), Some(0.0));
    assert!(chord_bound2(&curve, 2.0, 3.0).is_some());
}
