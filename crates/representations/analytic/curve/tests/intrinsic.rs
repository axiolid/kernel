//! Intrinsic plane curves: curvature laws and their closed-form operations.
//!
//! The suite pins the two things that make this representation worth having:
//! a clothoid stays exact instead of becoming a polyline, and every operation
//! offered is genuinely closed-form rather than a sampled approximation.

use axiolid_core::{Frame2, Point2, Vec2};
use axiolid_curve::{CurvatureLaw, Curve2, Intrinsic2};

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

#[test]
fn a_clothoid_derives_its_rate_from_endpoint_curvatures() {
    // Straight into a 50 m radius over 30 m of transition.
    let law = CurvatureLaw::clothoid(0.0, 1.0 / 50.0, 30.0);
    let CurvatureLaw::Polynomial { coefficients } = &law else {
        panic!("a clothoid is a degree-one polynomial law");
    };
    assert_eq!(coefficients[0], 0.0);
    // Rate is exactly (end - start) / length, stored, not sampled.
    assert_eq!(coefficients[1], (1.0 / 50.0) / 30.0);
}

#[test]
fn a_clothoid_of_zero_length_degrades_to_its_start_curvature() {
    // No division by zero smuggled into a coefficient.
    let law = CurvatureLaw::clothoid(0.25, 1.0, 0.0);
    assert_eq!(law, CurvatureLaw::Constant { curvature: 0.25 });
}

#[test]
fn the_derivative_of_a_clothoid_is_its_constant_sharpness() {
    let law = CurvatureLaw::clothoid(0.0, 0.02, 30.0);
    let rate = law.derivative();
    assert!(rate.is_constant());
    let CurvatureLaw::Polynomial { coefficients } = &rate else {
        panic!("differentiating a polynomial law stays polynomial");
    };
    assert_eq!(coefficients, &[0.02 / 30.0]);
}

#[test]
fn a_circular_arc_turns_its_length_over_its_radius() {
    // A quarter circle of radius 2 turns exactly pi/2, in closed form.
    let radius = 2.0;
    let quarter = core::f64::consts::FRAC_PI_2 * radius;
    let curve = Intrinsic2::new(frame(), CurvatureLaw::circular(1.0 / radius), quarter);
    let turning = curve.total_turning().expect("finite length");
    assert!((turning - core::f64::consts::FRAC_PI_2).abs() < 1e-12);
}

#[test]
fn a_clothoid_turns_the_mean_of_its_endpoint_curvatures() {
    // For a linear law the integral is the average curvature times length.
    // Checking against that independent formula, not against itself.
    let (start, end, length) = (0.0, 0.04, 25.0);
    let curve = Intrinsic2::new(frame(), CurvatureLaw::clothoid(start, end, length), length);
    let turning = curve.total_turning().expect("finite length");
    assert!((turning - (start + end) / 2.0 * length).abs() < 1e-12);
}

#[test]
fn a_straight_law_accumulates_no_turning() {
    let curve = Intrinsic2::new(frame(), CurvatureLaw::straight(), 100.0);
    assert!(curve.is_straight());
    assert_eq!(curve.total_turning(), Some(0.0));
}

#[test]
fn an_infinite_length_has_no_defined_turning() {
    // Undefined, not merely large: refuse rather than return a number.
    let curve = Intrinsic2::new(frame(), CurvatureLaw::circular(1.0), f64::INFINITY);
    assert_eq!(curve.total_turning(), None);
}

#[test]
fn a_sinusoid_over_a_whole_period_turns_only_by_its_mean() {
    // The oscillation cancels exactly over one period; only the mean survives.
    let w = 0.5;
    let period = 2.0 * core::f64::consts::PI / w;
    let law = CurvatureLaw::Sinusoid {
        mean: 0.01,
        amplitude: 0.3,
        angular_frequency: w,
        phase: 0.7,
    };
    let curve = Intrinsic2::new(frame(), law, period);
    let turning = curve.total_turning().expect("finite length");
    assert!((turning - 0.01 * period).abs() < 1e-12);
}

#[test]
fn a_zero_frequency_sinusoid_still_integrates() {
    // The closed form divides by the frequency; the degenerate case is
    // handled by integrating the frozen constant instead of dividing by zero.
    let law = CurvatureLaw::Sinusoid {
        mean: 0.2,
        amplitude: 1.0,
        angular_frequency: 0.0,
        phase: core::f64::consts::FRAC_PI_2,
    };
    let curve = Intrinsic2::new(frame(), law, 10.0);
    let turning = curve.total_turning().expect("finite length");
    // sin(pi/2) = 1, so the integrand is the constant 0.2 + 1.0.
    assert!((turning - 12.0).abs() < 1e-12);
    assert!(turning.is_finite());
}

#[test]
fn differentiating_a_sinusoid_stays_in_the_family() {
    // d/ds of a sine is a cosine, folded back to a sine by a phase shift,
    // so the law family is closed under differentiation.
    let law = CurvatureLaw::Sinusoid {
        mean: 5.0,
        amplitude: 2.0,
        angular_frequency: 3.0,
        phase: 0.0,
    };
    let CurvatureLaw::Sinusoid {
        mean,
        amplitude,
        angular_frequency,
        phase,
    } = law.derivative()
    else {
        panic!("the derivative of a sinusoid is a sinusoid");
    };
    assert_eq!(mean, 0.0);
    assert_eq!(amplitude, 6.0);
    assert_eq!(angular_frequency, 3.0);
    assert_eq!(phase, core::f64::consts::FRAC_PI_2);
}

#[test]
fn mirroring_negates_the_turning_exactly() {
    let law = CurvatureLaw::clothoid(0.01, 0.05, 20.0);
    let forward = Intrinsic2::new(frame(), law.clone(), 20.0);
    let mirrored = Intrinsic2::new(frame(), law.reversed_orientation(), 20.0);
    let a = forward.total_turning().expect("finite");
    let b = mirrored.total_turning().expect("finite");
    assert_eq!(a, -b);
}

#[test]
fn an_empty_polynomial_is_the_straight_law() {
    let law = CurvatureLaw::Polynomial {
        coefficients: Vec::new(),
    };
    assert!(law.is_straight());
    assert!(law.is_constant());
}

#[test]
fn a_curvature_law_is_reachable_as_a_plane_curve_variant() {
    // The point of the variant: a clothoid is a Curve2 like any other.
    let curve = Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::clothoid(0.0, 0.02, 30.0),
        30.0,
    ));
    let Curve2::Intrinsic(intrinsic) = &curve else {
        panic!("intrinsic curve round-trips through the atomic enum");
    };
    assert_eq!(intrinsic.length, 30.0);
    assert!(!intrinsic.is_straight());
}

#[test]
fn differentiating_a_higher_degree_law_applies_the_power_rule() {
    // Degree 1 hides the power rule: the factor there is 1, so a derivative
    // that forgot to multiply by the power would still look correct. A cubic
    // law (the Bloss/cubic-parabola family) is what actually pins it.
    // d/ds [1 + 2s + 3s^2 + 4s^3] = 2 + 6s + 12s^2
    let law = CurvatureLaw::Polynomial {
        coefficients: vec![1.0, 2.0, 3.0, 4.0],
    };
    let CurvatureLaw::Polynomial { coefficients } = law.derivative() else {
        panic!("differentiating a polynomial law stays polynomial");
    };
    assert_eq!(coefficients, vec![2.0, 6.0, 12.0]);
}

#[test]
fn a_quadratic_law_integrates_with_descending_weights() {
    // Integral of 6s^2 over [0, 2] is 2s^3 = 16, not 6 * 2^3 = 48.
    // Pins the 1/(i+1) weight that a constant-weight integral would drop.
    let law = CurvatureLaw::Polynomial {
        coefficients: vec![0.0, 0.0, 6.0],
    };
    let curve = Intrinsic2::new(frame(), law, 2.0);
    let turning = curve.total_turning().expect("finite length");
    assert!((turning - 16.0).abs() < 1e-12);
}
