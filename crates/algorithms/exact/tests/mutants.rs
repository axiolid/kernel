//! Targeted tests for behaviour the broad tests left unpinned.
//!
//! Each case was added because mutation testing showed a mutant surviving
//! the rest of the suite. The comment names the mutant it kills.

use axiolid_core::Point2;
use axiolid_exact::{line_circle_hits, Arith, Circle, Dyadic, HitCount, Interval, Line, Root2};
use axiolid_guarantees::Sign;

fn exact(value: i64) -> Dyadic {
    Dyadic::from_f64(value as f64)
}

#[test]
fn the_filter_proves_an_exact_zero_input() {
    // Kills: "[0,0] not recognised as zero". The exact tier would still
    // answer, so only the fast tier's own answer shows the difference.
    assert_eq!(Interval::point(0.0).sign(), Some(Sign::Zero));
    assert_eq!(Interval::point(-0.0).neg().sign(), Some(Sign::Zero));
    // A widened result is never [0, 0], even when the true value is zero:
    // that is the case the exact tier exists for.
    let widened = Interval::point(0.0).mul(&Interval::point(3.0));
    assert_eq!(widened.sign(), None);
}

#[test]
fn a_square_is_never_negative_even_across_zero() {
    // Kills: "square lower bound 0 dropped". An interval straddling zero,
    // multiplied by itself as a general product, gets a negative lower
    // bound (lo * hi < 0). As a square the true lower bound is 0. The
    // fixture is 1 - 1: exactly zero, but widened outward to straddle it.
    let span = Interval::point(1.0).sub(&Interval::point(1.0));
    assert!(
        span.lo() < 0.0 && span.hi() > 0.0,
        "a widened zero straddles zero"
    );
    let squared = span.square();
    assert_eq!(
        squared.lo(),
        0.0,
        "a square's lower bound is zero, not negative"
    );
    assert!(squared.hi() > 0.0);
    let general = span.mul(&span);
    assert!(general.lo() < 0.0, "the general product really is looser");
}

#[test]
fn comparison_respects_negative_denominators() {
    // Kills: "cmp_sign denominator signs ignored". x = (1 + sqrt 2) / -2 is
    // negative (-1.207..); y = 0 / 1. So x < y, even though the numerator
    // difference (1 + sqrt 2) is positive.
    let x = Root2 {
        a: exact(1),
        b: exact(1),
        c: exact(2),
        d: exact(-2),
    };
    let y = Root2 {
        a: exact(0),
        b: exact(0),
        c: exact(0),
        d: exact(1),
    };
    assert_eq!(x.sign(), Some(Sign::Negative));
    assert_eq!(x.cmp_sign(&y), Some(Sign::Negative));
    assert_eq!(y.cmp_sign(&x), Some(Sign::Positive));
    // Both negative: (3 - sqrt 3) / -1 = -1.267.. < (1 + sqrt 2) / -2.
    let z = Root2 {
        a: exact(3),
        b: exact(-1),
        c: exact(3),
        d: exact(-1),
    };
    assert_eq!(z.cmp_sign(&x), Some(Sign::Negative));
}

#[test]
fn the_two_hits_lie_on_opposite_sides_of_a_line_through_the_centre() {
    // Kills: "hit orientation: branch ignored". The x axis crosses the
    // circle of radius 2 about the origin at x = -2 (first) and x = 2
    // (second). The y axis (upwards) has x < 0 on its left.
    let axis = Line::new(Point2::new(-5.0, 0.0), Point2::new(5.0, 0.0)).unwrap();
    let circle = Circle::new(Point2::new(0.0, 0.0), 2.0).unwrap();
    let Ok(HitCount::Secant(first, second)) = line_circle_hits(axis, circle) else {
        panic!("expected two hits");
    };
    let (up_from, up_to) = (Point2::new(0.0, -1.0), Point2::new(0.0, 1.0));
    assert_eq!(first.orientation(up_from, up_to), Ok(Sign::Positive));
    assert_eq!(second.orientation(up_from, up_to), Ok(Sign::Negative));
    // And a line through x = 2 exactly: the second hit lies on it.
    let (on_from, on_to) = (Point2::new(2.0, -1.0), Point2::new(2.0, 1.0));
    assert_eq!(second.orientation(on_from, on_to), Ok(Sign::Zero));
    assert_eq!(first.orientation(on_from, on_to), Ok(Sign::Positive));
}
