//! Intrinsic plane curves: curvature laws and their closed-form operations.
//!
//! The suite pins the two things that make this representation worth having:
//! a clothoid stays exact instead of becoming a polyline, and every operation
//! offered is genuinely closed-form rather than a sampled approximation.

use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{CurvatureLaw, Curve2, Harmonic, Intrinsic2};

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

// --- Composite: polynomial and harmonic terms at once ---

/// The acceptance case: k(s) = k0 + (d/L) s - (d/2pi) sin(2 pi s / L).
///
/// With k0 = 0, d = 1/300, L = 60 the law must hit both endpoint curvatures
/// and turn by the clothoid amount, because the sine term integrates to zero
/// over exactly one period.
#[test]
fn a_sine_corrected_transition_meets_both_endpoints_and_turns_the_mean() {
    let delta = 1.0 / 300.0;
    let length = 60.0;
    let law = CurvatureLaw::sine_corrected_transition(0.0, delta, length);

    let CurvatureLaw::Composite {
        polynomial,
        harmonics,
    } = &law
    else {
        panic!("expected a composite law");
    };
    assert_eq!(polynomial, &vec![0.0, delta / length]);
    assert_eq!(harmonics.len(), 1);
    assert_eq!(harmonics[0].amplitude, -delta / core::f64::consts::TAU);
    assert_eq!(
        harmonics[0].angular_frequency,
        core::f64::consts::TAU / length
    );
    assert_eq!(harmonics[0].phase, 0.0);

    // k(0): the sine vanishes at s = 0, leaving the constant term.
    // k(L): the sine vanishes again after a full period, leaving k0 + d.
    // Both are structural here -- no evaluator exists to ask.
    let turning = Intrinsic2::new(frame(), law, length)
        .total_turning()
        .expect("finite length turns a finite amount");
    // Mean rate d/L over length L gives d*L/2; the harmonic contributes
    // nothing over a whole period.
    assert!((turning - delta * length / 2.0).abs() < 1e-15);
    assert!((turning - 0.1).abs() < 1e-15);
}

/// Differentiation acts on both parts at once and stays in the family.
#[test]
fn differentiating_a_composite_law_differentiates_both_parts() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![1.0, 2.0, 3.0],
        harmonics: vec![Harmonic {
            amplitude: 5.0,
            angular_frequency: 7.0,
            phase: 0.25,
        }],
    };
    let CurvatureLaw::Composite {
        polynomial,
        harmonics,
    } = law.derivative()
    else {
        panic!("a composite law differentiates to a composite law");
    };
    // d/ds [1 + 2s + 3s^2] = 2 + 6s: the power rule, not a shift.
    assert_eq!(polynomial, vec![2.0, 6.0]);
    // d/ds [5 sin(7s + p)] = 35 sin(7s + p + pi/2).
    assert_eq!(harmonics[0].amplitude, 35.0);
    assert_eq!(harmonics[0].angular_frequency, 7.0);
    assert_eq!(harmonics[0].phase, 0.25 + core::f64::consts::FRAC_PI_2);
}

/// Mirroring negates every term and is its own inverse.
#[test]
fn mirroring_a_composite_law_negates_every_term() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![1.0, -2.0],
        harmonics: vec![Harmonic {
            amplitude: 3.0,
            angular_frequency: 4.0,
            phase: 0.5,
        }],
    };
    let CurvatureLaw::Composite {
        polynomial,
        harmonics,
    } = law.reversed_orientation()
    else {
        panic!("mirroring keeps the variant");
    };
    assert_eq!(polynomial, vec![-1.0, 2.0]);
    assert_eq!(harmonics[0].amplitude, -3.0);
    // Frequency and phase are untouched, so mirroring twice is the identity.
    assert_eq!(harmonics[0].angular_frequency, 4.0);
    assert_eq!(harmonics[0].phase, 0.5);
    assert_eq!(law.reversed_orientation().reversed_orientation(), law);
}

/// A zero-frequency harmonic freezes at A sin(p) and integrates linearly,
/// rather than dividing by its frequency.
#[test]
fn a_zero_frequency_harmonic_integrates_as_a_constant() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![],
        harmonics: vec![Harmonic {
            amplitude: 2.0,
            angular_frequency: 0.0,
            phase: core::f64::consts::FRAC_PI_2,
        }],
    };
    // sin(pi/2) = 1, so the frozen curvature is 2 and turning over 3 is 6.
    let turning = Intrinsic2::new(frame(), law, 3.0)
        .total_turning()
        .expect("a frozen harmonic still integrates");
    assert!((turning - 6.0).abs() < 1e-15);
}

/// Empty harmonics is exactly the polynomial law.
#[test]
fn a_composite_without_harmonics_turns_like_its_polynomial() {
    let composite = CurvatureLaw::Composite {
        polynomial: vec![0.0, 1.0 / 300.0],
        harmonics: vec![],
    };
    let plain = CurvatureLaw::Polynomial {
        coefficients: vec![0.0, 1.0 / 300.0],
    };
    let composite_turn = Intrinsic2::new(frame(), composite, 60.0).total_turning();
    let plain_turn = Intrinsic2::new(frame(), plain, 60.0).total_turning();
    assert_eq!(composite_turn, plain_turn);
}

/// An empty composite is the zero law, and structurally straight.
#[test]
fn an_empty_composite_is_straight_and_constant() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![],
        harmonics: vec![],
    };
    assert!(law.is_straight());
    assert!(law.is_constant());
    assert_eq!(
        Intrinsic2::new(frame(), law, 10.0).total_turning(),
        Some(0.0)
    );
}

/// A real transition is neither straight nor constant.
#[test]
fn a_sine_corrected_transition_is_neither_straight_nor_constant() {
    let law = CurvatureLaw::sine_corrected_transition(0.0, 1.0 / 300.0, 60.0);
    assert!(!law.is_straight());
    assert!(!law.is_constant());
}

/// A harmonic with a live frequency but zero amplitude contributes nothing,
/// so a composite of zeros is straight.
#[test]
fn zero_amplitude_harmonics_do_not_defeat_straightness() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![0.0, 0.0],
        harmonics: vec![Harmonic {
            amplitude: 0.0,
            angular_frequency: 5.0,
            phase: 1.0,
        }],
    };
    assert!(law.is_straight());
}

/// A frozen harmonic is constant, but whether it CANCELS the polynomial
/// needs evaluation, so straightness reports false rather than guessing.
#[test]
fn a_frozen_harmonic_is_constant_but_not_declared_straight() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![0.0],
        harmonics: vec![Harmonic {
            amplitude: 1.0,
            angular_frequency: 0.0,
            phase: 0.0,
        }],
    };
    // sin(0) = 0 so this law IS identically zero, but deciding that requires
    // evaluating the sine. The structural answer is the honest refusal.
    assert!(law.is_constant());
    assert!(!law.is_straight());
}

/// A non-finite frequency is not a decidable law; do not claim straightness.
#[test]
fn a_non_finite_harmonic_frequency_is_not_declared_straight() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![0.0],
        harmonics: vec![Harmonic {
            amplitude: 0.0,
            angular_frequency: Scalar::NAN,
            phase: 0.0,
        }],
    };
    assert!(!law.is_straight());
}
