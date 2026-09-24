//! Constructions: exact answers on configurations built to be degenerate.
//!
//! Every expected answer comes from geometry chosen by hand (a crossing
//! placed exactly on a line, a line exactly tangent to a circle), never from
//! another call into the code under test.

use axiolid_core::Point2;
use axiolid_exact::{
    compare_along, crossing_orientation, line_circle_hits, Circle, ExactError, HitCount, Line,
};
use axiolid_guarantees::Sign;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn line(ax: f64, ay: f64, bx: f64, by: f64) -> Line {
    Line::new(p(ax, ay), p(bx, by)).expect("valid line")
}

fn circle(x: f64, y: f64, r: f64) -> Circle {
    Circle::new(p(x, y), r).expect("valid circle")
}

#[test]
fn a_crossing_exactly_on_the_query_line_is_zero() {
    // The diagonals of the square [0,3]^2 cross at (1.5, 1.5), which lies
    // exactly on y = x. Shift everything by a large offset so doubles can
    // no longer represent the arithmetic along the way.
    let o = 1e9 + 0.1;
    let first = line(o, o, o + 3.0, o + 3.0);
    let second = line(o, o + 3.0, o + 3.0, o);
    let on = crossing_orientation(first, second, p(o - 7.0, o - 7.0), p(o + 11.0, o + 11.0));
    assert_eq!(on, Ok(Some(Sign::Zero)));
    let left = crossing_orientation(first, second, p(o, o - 1.0), p(o + 1.0, o));
    assert_eq!(left, Ok(Some(Sign::Positive)), "above y = x - 1");
}

#[test]
fn parallel_lines_have_no_crossing() {
    let first = line(0.0, 0.0, 1.0, 1.0);
    let second = line(0.0, 1.0, 1.0, 2.0);
    assert_eq!(
        crossing_orientation(first, second, p(0.0, 0.0), p(1.0, 0.0)),
        Ok(None)
    );
}

#[test]
fn tangency_is_decided_exactly() {
    // y = 5 touches the circle of radius 5 at the origin exactly once.
    let unit = circle(0.0, 0.0, 5.0);
    match line_circle_hits(line(-10.0, 5.0, 10.0, 5.0), unit) {
        Ok(HitCount::Tangent(hit)) => {
            assert_eq!(
                hit.cmp_param(0.5),
                Ok(Sign::Zero),
                "touches at the midpoint"
            );
            let touch = hit.approx_point();
            assert!(touch.x.abs() < 1e-12 && (touch.y - 5.0).abs() < 1e-12);
        }
        other => panic!("expected tangent, got {other:?}"),
    }
    // One ulp above y = 5 misses; one below crosses twice. No tolerance.
    let above = 5f64.next_up();
    let below = 5f64.next_down();
    assert!(matches!(
        line_circle_hits(line(-10.0, above, 10.0, above), unit),
        Ok(HitCount::Missed)
    ));
    assert!(matches!(
        line_circle_hits(line(-10.0, below, 10.0, below), unit),
        Ok(HitCount::Secant(..))
    ));
}

#[test]
fn secant_hits_come_in_order_and_sit_on_the_circle() {
    // The x axis crosses the circle about (1, 0) with radius 2 at x = -1
    // and x = 3. On the line from (-5,0) to (3,0), t = (x + 5) / 8, so the
    // hits are at t = 0.5 and 1.0: both exact doubles, so "equal" is exact.
    let axis = line(-5.0, 0.0, 3.0, 0.0);
    let Ok(HitCount::Secant(first, second)) = line_circle_hits(axis, circle(1.0, 0.0, 2.0)) else {
        panic!("expected two hits");
    };
    assert_eq!(first.cmp_param(0.5), Ok(Sign::Zero));
    assert_eq!(second.cmp_param(1.0), Ok(Sign::Zero));
    // 0.4 as a double is not 0.4: the exact comparison must not say equal.
    assert_eq!(first.cmp_param(0.4), Ok(Sign::Positive));
    assert_eq!(compare_along(first, second), Ok(Sign::Negative));
    assert_eq!(compare_along(second, first), Ok(Sign::Positive));
    assert_eq!(compare_along(first, first), Ok(Sign::Zero));
}

#[test]
fn hits_on_different_circles_order_exactly() {
    // Circles about (0,0) with radii 3 and 5 meet the x axis at x = +-3
    // and +-5: rational hits, so ordering has exact answers.
    let axis = line(-10.0, 0.0, 10.0, 0.0);
    let Ok(HitCount::Secant(_, three)) = line_circle_hits(axis, circle(0.0, 0.0, 3.0)) else {
        panic!();
    };
    let Ok(HitCount::Secant(minus_five, _)) = line_circle_hits(axis, circle(0.0, 0.0, 5.0)) else {
        panic!();
    };
    assert_eq!(compare_along(minus_five, three), Ok(Sign::Negative));
    // Circles of radius 1 about (1,0) and (-1,0) both pass through the
    // origin exactly, so their hits there are the same point.
    let right = circle(1.0, 0.0, 1.0);
    let left = circle(-1.0, 0.0, 1.0);
    let (Ok(HitCount::Secant(right_first, _)), Ok(HitCount::Secant(_, left_second))) =
        (line_circle_hits(axis, right), line_circle_hits(axis, left))
    else {
        panic!();
    };
    assert_eq!(compare_along(right_first, left_second), Ok(Sign::Zero));
}

#[test]
fn irrational_hits_are_compared_across_radicands() {
    // The diagonal y = x meets the unit circle at +-1/sqrt(2) and the circle
    // of radius 2 about the origin at +-sqrt(2). t on the line (0,0)->(1,1)
    // is x itself, so the hits are at -0.707.., 0.707.., -1.414.., 1.414...
    let diagonal = line(0.0, 0.0, 1.0, 1.0);
    let Ok(HitCount::Secant(u_minus, u_plus)) = line_circle_hits(diagonal, circle(0.0, 0.0, 1.0))
    else {
        panic!();
    };
    let Ok(HitCount::Secant(w_minus, w_plus)) = line_circle_hits(diagonal, circle(0.0, 0.0, 2.0))
    else {
        panic!();
    };
    assert_eq!(compare_along(w_minus, u_minus), Ok(Sign::Negative));
    assert_eq!(compare_along(u_plus, w_plus), Ok(Sign::Negative));
    // 1/sqrt(2) against 0.7071067811865476 (the double nearest it): the
    // double is above, and only exact arithmetic can say so.
    assert_eq!(
        u_plus.cmp_param(std::f64::consts::FRAC_1_SQRT_2),
        Ok(Sign::Negative)
    );
    // Which side of the vertical line x = 1/sqrt(2)-ish each hit lies on.
    let x = std::f64::consts::FRAC_1_SQRT_2;
    assert_eq!(
        u_plus.orientation(p(x, -1.0), p(x, 1.0)),
        Ok(Sign::Positive)
    );
}

#[test]
fn invalid_input_is_refused_by_reason() {
    assert_eq!(
        Line::new(p(1.0, 1.0), p(1.0, 1.0)),
        Err(ExactError::DegenerateLine)
    );
    assert_eq!(
        Line::new(p(f64::NAN, 0.0), p(1.0, 1.0)),
        Err(ExactError::NonFinite)
    );
    assert_eq!(
        Circle::new(p(0.0, 0.0), -1.0),
        Err(ExactError::NegativeRadius)
    );
    let axis = line(-5.0, 0.0, 5.0, 0.0);
    let other = line(-5.0, 1.0, 5.0, 1.0);
    let unit = circle(0.0, 0.0, 2.0);
    let (Ok(HitCount::Secant(a, _)), Ok(HitCount::Secant(b, _))) =
        (line_circle_hits(axis, unit), line_circle_hits(other, unit))
    else {
        panic!();
    };
    assert_eq!(compare_along(a, b), Err(ExactError::DifferentLines));
    assert_eq!(a.cmp_param(f64::INFINITY), Err(ExactError::NonFinite));
}
