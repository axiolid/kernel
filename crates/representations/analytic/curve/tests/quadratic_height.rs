//! The height of a quadratic graph is exact to rounding for either sign of
//! `b`, where either textbook form of the root loses digits to cancellation.
//!
//! The coefficients are built from roots that are powers of two, so `b` (the
//! negated sum of the roots) and `c` (their product) are exact and the roots
//! are known without any arithmetic of the height's own.

use axiolid_curve::{Branch, QuadraticGraph2, Trig2};

fn constant(x: f64) -> Trig2 {
    Trig2 {
        constant: x,
        ..Trig2::default()
    }
}

#[test]
fn both_roots_are_exact_to_rounding_for_either_sign_of_b() {
    let big = 2.0_f64.powi(24);
    let small = 2.0_f64.powi(-24);
    for sign in [1.0, -1.0] {
        let (r1, r2) = (sign * small, sign * big);
        let g = |branch| QuadraticGraph2 {
            a: constant(1.0),
            b: constant(-(r1 + r2)),
            c: constant(r1 * r2),
            branch,
        };
        for branch in [Branch::Plus, Branch::Minus] {
            let v = g(branch).height(0.3).expect("two real roots");
            // `(-b +- sqrt(D)) / 2` is the root of larger magnitude for the
            // sign that follows `-b`, and the smaller for the other.
            let exact = if (branch == Branch::Plus) == (sign > 0.0) {
                r2
            } else {
                r1
            };
            assert!(
                (v - exact).abs() <= 4.0 * f64::EPSILON * exact.abs(),
                "{branch:?}, b {}: {v} against {exact}",
                -(r1 + r2)
            );
        }
        // The two branches are the two roots.
        let (p, m) = (
            g(Branch::Plus).height(0.3).unwrap(),
            g(Branch::Minus).height(0.3).unwrap(),
        );
        assert!((p - m).abs() > 1.0 && p.signum() == m.signum());
    }
}

#[test]
fn a_branch_that_diverges_where_a_vanishes_has_no_height() {
    // `b v + c = 0` with `b > 0`: the plus root is `-c / b`, the minus root
    // runs off to infinity.
    let g = |branch| QuadraticGraph2 {
        a: Trig2::default(),
        b: constant(2.0),
        c: constant(1.0),
        branch,
    };
    assert_eq!(g(Branch::Plus).height(0.0), Some(-0.5));
    assert_eq!(g(Branch::Minus).height(0.0), None);
}
