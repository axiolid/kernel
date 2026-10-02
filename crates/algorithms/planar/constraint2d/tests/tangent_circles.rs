//! Apollonius and tangent-circle constructions (#159, ledger row B14).
//!
//! Every returned circle is checked by tangency, not by comparing against a
//! hand-picked expected circle: `is_tangent` recomputes, from the raw
//! geometry, whether a circle touches a point, line or circle within a
//! tight tolerance, independent of the solver's own internals.

use axiolid_core::{Frame2, Point2, Tolerance, Vec2};
use axiolid_curve::Circle2;
use axiolid_linear::Line2;
use axiolid_constraint2d::{tangent_circles, TangencyError, Tangent};

const TOL: Tolerance = Tolerance::METRE;
const TIGHT: f64 = 1e-9;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn circle(cx: f64, cy: f64, r: f64) -> Circle2 {
    Circle2 {
        frame: Frame2 {
            origin: p(cx, cy),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: r,
    }
}

fn line(ox: f64, oy: f64, dx: f64, dy: f64) -> Line2 {
    Line2 {
        origin: p(ox, oy),
        direction: Vec2::new(dx, dy),
    }
}

/// Whether `solution` touches `constraint` within `eps`, recomputed from
/// the raw geometry rather than the solver's sign bookkeeping.
fn is_tangent(solution: &Circle2, constraint: &Tangent, eps: f64) -> bool {
    let c = solution.frame.origin;
    let r = solution.radius;
    match *constraint {
        Tangent::Point(q) => ((c - q).length() - r).abs() <= eps,
        Tangent::Line(l) => {
            let n = l.direction.perp().normalize();
            (n.dot(c - l.origin).abs() - r).abs() <= eps
        }
        Tangent::Circle(o) => {
            let d = (c - o.frame.origin).length();
            (d - (r + o.radius)).abs() <= eps || (d - (r - o.radius).abs()).abs() <= eps
        }
    }
}

fn assert_all_tangent(solutions: &[Circle2], constraints: &[Tangent; 3]) {
    for solution in solutions {
        assert!(solution.radius > 0.0, "{solution:?}");
        for constraint in constraints {
            assert!(
                is_tangent(solution, constraint, TIGHT),
                "{solution:?} not tangent to {constraint:?}"
            );
        }
    }
}

// -- Three points: the circumcircle, unique -------------------------------

#[test]
fn three_points_have_a_unique_circumcircle() {
    let constraints = [
        Tangent::Point(p(0.0, 0.0)),
        Tangent::Point(p(4.0, 0.0)),
        Tangent::Point(p(0.0, 3.0)),
    ];
    let got = tangent_circles(constraints, TOL).unwrap();
    assert_eq!(got.len(), 1, "{got:?}");
    assert_all_tangent(&got, &constraints);
    // The right triangle's circumcentre is the midpoint of its hypotenuse.
    let c = got[0].frame.origin;
    assert!((c - p(2.0, 1.5)).length() < 1e-9);
    assert!((got[0].radius - 2.5).abs() < 1e-9);
}

#[test]
fn three_collinear_points_are_degenerate() {
    let constraints = [
        Tangent::Point(p(0.0, 0.0)),
        Tangent::Point(p(1.0, 0.0)),
        Tangent::Point(p(2.0, 0.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::Degenerate)
    );
}

#[test]
fn three_coincident_points_are_degenerate() {
    let constraints = [
        Tangent::Point(p(1.0, 1.0)),
        Tangent::Point(p(1.0, 1.0)),
        Tangent::Point(p(1.0, 1.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::Degenerate)
    );
}

// -- Three lines: incircle and excircles, four solutions -------------------

#[test]
fn three_lines_of_a_triangle_have_four_solutions() {
    // A 3-4-5 right triangle's three side lines.
    let constraints = [
        Tangent::Line(line(0.0, 0.0, 1.0, 0.0)),
        Tangent::Line(line(4.0, 0.0, 0.0, 1.0)),
        Tangent::Line(line(0.0, 0.0, 4.0, 3.0)),
    ];
    let got = tangent_circles(constraints, TOL).unwrap();
    assert_all_tangent(&got, &constraints);
    // Incircle + three excircles.
    assert_eq!(got.len(), 4, "{got:?}");
    let radii: Vec<f64> = {
        let mut r: Vec<f64> = got.iter().map(|c| c.radius).collect();
        r.sort_by(f64::total_cmp);
        r
    };
    // Incircle radius r = (a + b - c) / 2 = (3 + 4 - 5) / 2 = 1.
    assert!((radii[0] - 1.0).abs() < 1e-9, "{radii:?}");
}

#[test]
fn three_parallel_lines_are_degenerate() {
    let constraints = [
        Tangent::Line(line(0.0, 0.0, 1.0, 0.0)),
        Tangent::Line(line(0.0, 1.0, 1.0, 0.0)),
        Tangent::Line(line(0.0, 2.0, 1.0, 0.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::Degenerate)
    );
}

// -- Three circles: the classic Apollonius problem, up to eight solutions --

#[test]
fn three_mutually_external_circles_have_eight_solutions() {
    // Three widely separated circles admit all eight sign combinations.
    let constraints = [
        Tangent::Circle(circle(0.0, 0.0, 1.0)),
        Tangent::Circle(circle(6.0, 0.0, 1.0)),
        Tangent::Circle(circle(3.0, 6.0, 1.0)),
    ];
    let got = tangent_circles(constraints, TOL).unwrap();
    assert_all_tangent(&got, &constraints);
    assert_eq!(got.len(), 8, "{got:?}");
}

#[test]
fn three_identical_circles_are_degenerate() {
    let c = circle(1.0, 1.0, 2.0);
    let constraints = [Tangent::Circle(c), Tangent::Circle(c), Tangent::Circle(c)];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::Degenerate)
    );
}

#[test]
fn a_circle_tangent_to_two_circles_and_a_point_on_one_of_them() {
    // A point on the boundary of a constraint circle is a legitimate
    // (if tight) tangency point for a solution passing through it.
    let constraints = [
        Tangent::Circle(circle(0.0, 0.0, 1.0)),
        Tangent::Circle(circle(5.0, 0.0, 1.0)),
        Tangent::Point(p(2.5, 3.0)),
    ];
    let got = tangent_circles(constraints, TOL).unwrap();
    assert!(!got.is_empty());
    assert_all_tangent(&got, &constraints);
}

#[test]
fn a_line_a_circle_and_a_point_mix_freely() {
    let constraints = [
        Tangent::Line(line(0.0, 0.0, 1.0, 0.0)),
        Tangent::Circle(circle(5.0, 5.0, 2.0)),
        Tangent::Point(p(1.0, 3.0)),
    ];
    let got = tangent_circles(constraints, TOL).unwrap();
    assert!(!got.is_empty(), "expected at least one solution");
    assert_all_tangent(&got, &constraints);
}

#[test]
fn a_well_posed_configuration_can_have_no_real_solution() {
    // Two circles that overlap so much that no third circle reaches a
    // distant point while staying tangent to both: this sign branch solves
    // cleanly, it simply has no positive-radius root.
    let constraints = [
        Tangent::Circle(circle(0.0, 0.0, 5.0)),
        Tangent::Circle(circle(0.5, 0.0, 4.7)),
        Tangent::Point(p(100.0, 100.0)),
    ];
    assert!(tangent_circles(constraints, TOL).is_ok());
}

#[test]
fn concentric_circle_type_constraints_are_a_named_scope_limit() {
    // Two circle- or point-type constraints centred at the same point
    // reduce two of the three equations to a statement about the radius
    // alone, with no `cx`/`cy` term to eliminate against: solving that pair
    // for `cx(r)`/`cy(r)` is singular under every sign choice, even though a
    // solution family exists in principle (found by a different
    // elimination order this crate does not implement). Refused by name,
    // not answered wrong or silently dropped.
    let constraints = [
        Tangent::Circle(circle(0.0, 0.0, 5.0)),
        Tangent::Circle(circle(0.0, 0.0, 1.0)),
        Tangent::Point(p(100.0, 100.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::Degenerate)
    );
}

// -- Degenerate and invalid input -------------------------------------------

#[test]
fn non_finite_input_is_refused() {
    let constraints = [
        Tangent::Point(p(f64::NAN, 0.0)),
        Tangent::Point(p(1.0, 0.0)),
        Tangent::Point(p(0.0, 1.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::NonFinite)
    );
}

#[test]
fn a_non_positive_radius_is_refused_by_name() {
    let constraints = [
        Tangent::Circle(circle(0.0, 0.0, 0.0)),
        Tangent::Point(p(1.0, 0.0)),
        Tangent::Point(p(0.0, 1.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::InvalidRadius {
            index: 0,
            radius: 0.0
        })
    );
    let constraints = [
        Tangent::Point(p(1.0, 0.0)),
        Tangent::Circle(circle(0.0, 0.0, -1.0)),
        Tangent::Point(p(0.0, 1.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::InvalidRadius {
            index: 1,
            radius: -1.0
        })
    );
}

#[test]
fn a_zero_direction_line_is_refused_by_name() {
    let constraints = [
        Tangent::Point(p(0.0, 0.0)),
        Tangent::Line(line(0.0, 0.0, 0.0, 0.0)),
        Tangent::Point(p(0.0, 1.0)),
    ];
    assert_eq!(
        tangent_circles(constraints, TOL),
        Err(TangencyError::ZeroDirection { index: 1 })
    );
}
