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

/// Curvature of a composite law at one arc length, computed in the TEST
/// only. This is direct closed-form arithmetic on stored coefficients --
/// not position evaluation, quadrature, or sampling -- and it lives here
/// rather than in the crate so the representation keeps refusing to
/// evaluate. It exists to pin endpoint values the spec names explicitly.
fn composite_curvature_at(law: &CurvatureLaw, s: Scalar) -> Scalar {
    let CurvatureLaw::Composite {
        polynomial,
        harmonics,
    } = law
    else {
        panic!("expected a composite law");
    };
    let poly: Scalar = polynomial
        .iter()
        .enumerate()
        .map(|(power, c)| c * s.powi(power as i32))
        .sum();
    let harm: Scalar = harmonics
        .iter()
        .map(|h| h.amplitude * (h.angular_frequency * s + h.phase).sin())
        .sum();
    poly + harm
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

    // The spec names both endpoint curvatures, so assert them rather than
    // asserting the coefficients that imply them. The sine vanishes at s = 0
    // and again after a whole period, so the endpoints are exactly k0 and k0+d.
    assert!(composite_curvature_at(&law, 0.0).abs() < 1e-18);
    assert!((composite_curvature_at(&law, length) - delta).abs() < 1e-18);

    // A quarter along, the sine is at its extreme rather than a zero, so this
    // is the point that actually discriminates the harmonic term: a law that
    // dropped the sine would still pass the endpoints and the total turning.
    let quarter = composite_curvature_at(&law, length / 4.0);
    let expected_quarter = delta / 4.0 - delta / core::f64::consts::TAU;
    assert!((quarter - expected_quarter).abs() < 1e-15);
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

// --- Piecewise: several laws over one arc-length domain ---

/// The alignment shape the variant exists for: straight, then a clothoid
/// transition, then a circular arc -- all under ONE start frame, because
/// no interior frame can be computed without a Fresnel integral.
#[test]
fn a_three_piece_alignment_turns_the_sum_of_its_pieces() {
    let radius = 300.0;
    let curvature = 1.0 / radius;
    let (tangent, transition, arc) = (50.0, 60.0, 100.0);

    let law = CurvatureLaw::piecewise(
        vec![tangent, tangent + transition],
        vec![
            CurvatureLaw::straight(),
            CurvatureLaw::clothoid(0.0, curvature, transition),
            CurvatureLaw::circular(curvature),
        ],
    );
    assert!(law.is_well_formed());

    let total = tangent + transition + arc;
    let turning = Intrinsic2::new(frame(), law, total)
        .total_turning()
        .expect("a well-formed piecewise law over a finite length turns");

    // Straight contributes nothing; the clothoid turns the mean of its
    // endpoint curvatures; the arc turns length/radius.
    let expected = 0.0 + curvature * transition / 2.0 + curvature * arc;
    assert!((turning - expected).abs() < 1e-15);
}

/// Each piece integrates over its OWN subinterval, not the whole domain.
/// Two arcs of equal length and opposite curvature must cancel exactly;
/// integrating either over the full span would not cancel.
#[test]
fn piecewise_pieces_integrate_over_their_own_subintervals() {
    let law = CurvatureLaw::piecewise(
        vec![10.0],
        vec![CurvatureLaw::circular(0.2), CurvatureLaw::circular(-0.2)],
    );
    let turning = Intrinsic2::new(frame(), law, 20.0)
        .total_turning()
        .expect("well formed");
    assert!(turning.abs() < 1e-15);
}

/// Seams are positions in arc length, so an off-centre break must shift
/// the balance. This is what pins the subinterval arithmetic.
#[test]
fn an_off_centre_seam_shifts_the_turning() {
    let law = CurvatureLaw::piecewise(
        vec![5.0],
        vec![CurvatureLaw::circular(0.2), CurvatureLaw::circular(-0.2)],
    );
    let turning = Intrinsic2::new(frame(), law, 20.0)
        .total_turning()
        .expect("well formed");
    // 0.2*5 - 0.2*15 = -2.0
    assert!((turning - -2.0).abs() < 1e-15);
}

/// Straight exactly when every piece is straight -- seams are irrelevant.
#[test]
fn a_piecewise_law_is_straight_only_when_every_piece_is() {
    let all_straight = CurvatureLaw::piecewise(
        vec![5.0],
        vec![CurvatureLaw::straight(), CurvatureLaw::straight()],
    );
    assert!(all_straight.is_straight());
    assert!(all_straight.is_constant());

    let one_bent = CurvatureLaw::piecewise(
        vec![5.0],
        vec![CurvatureLaw::straight(), CurvatureLaw::circular(0.1)],
    );
    assert!(!one_bent.is_straight());
}

/// Constant needs every piece constant AND mutually equal: two arcs of
/// different radius are each constant, but the law that joins them is not.
#[test]
fn differing_constant_pieces_are_not_a_constant_law() {
    let same = CurvatureLaw::piecewise(
        vec![5.0],
        vec![CurvatureLaw::circular(0.1), CurvatureLaw::circular(0.1)],
    );
    assert!(same.is_constant());
    assert!(!same.is_straight());

    let differing = CurvatureLaw::piecewise(
        vec![5.0],
        vec![CurvatureLaw::circular(0.1), CurvatureLaw::circular(0.2)],
    );
    assert!(!differing.is_constant());
}

/// Differentiation is per piece and keeps the seams, so the derivative of
/// a straight-then-clothoid law is 0 then the constant sharpness -- a real
/// jump at the seam, which is mathematically correct for this family.
#[test]
fn differentiating_a_piecewise_law_differentiates_each_piece() {
    let sharpness = 0.5;
    let law = CurvatureLaw::piecewise(
        vec![10.0],
        vec![
            CurvatureLaw::straight(),
            CurvatureLaw::Polynomial {
                coefficients: vec![0.0, sharpness],
            },
        ],
    );

    let CurvatureLaw::Piecewise { breaks, laws } = law.derivative() else {
        panic!("the derivative of a piecewise law stays piecewise");
    };
    assert_eq!(breaks, vec![10.0]);
    assert_eq!(laws[0], CurvatureLaw::Constant { curvature: 0.0 });
    assert_eq!(
        laws[1],
        CurvatureLaw::Polynomial {
            coefficients: vec![sharpness],
        }
    );
}

/// Mirroring negates every piece and leaves the seams alone, so the
/// mirrored law turns exactly the opposite amount.
#[test]
fn mirroring_a_piecewise_law_negates_every_piece() {
    let law = CurvatureLaw::piecewise(
        vec![10.0],
        vec![
            CurvatureLaw::clothoid(0.0, 0.01, 10.0),
            CurvatureLaw::circular(0.01),
        ],
    );
    let forward = Intrinsic2::new(frame(), law.clone(), 25.0)
        .total_turning()
        .expect("well formed");
    let mirrored = Intrinsic2::new(frame(), law.reversed_orientation(), 25.0)
        .total_turning()
        .expect("well formed");
    assert!((forward + mirrored).abs() < 1e-15);
}

/// A malformed law refuses to report turning rather than guessing which
/// piece the missing seam belonged to.
#[test]
fn a_malformed_piecewise_law_refuses_to_turn() {
    // Three pieces need two seams; one seam leaves the tiling ambiguous.
    let law = CurvatureLaw::Piecewise {
        breaks: vec![5.0],
        laws: vec![
            CurvatureLaw::straight(),
            CurvatureLaw::circular(0.1),
            CurvatureLaw::circular(0.2),
        ],
    };
    assert!(!law.is_well_formed());
    assert_eq!(Intrinsic2::new(frame(), law, 20.0).total_turning(), None);
}

/// Unordered seams are malformed: a descending break would make a piece
/// span a negative length.
#[test]
fn descending_seams_are_malformed() {
    let law = CurvatureLaw::Piecewise {
        breaks: vec![10.0, 5.0],
        laws: vec![
            CurvatureLaw::straight(),
            CurvatureLaw::circular(0.1),
            CurvatureLaw::circular(0.2),
        ],
    };
    assert!(!law.is_well_formed());
    assert_eq!(Intrinsic2::new(frame(), law, 20.0).total_turning(), None);
}

/// A seam beyond the curve length does not tile the domain, so the stored
/// integral is not the one being asked for. Refuse instead of clamping.
#[test]
fn a_seam_outside_the_curve_length_refuses() {
    let law = CurvatureLaw::piecewise(
        vec![50.0],
        vec![CurvatureLaw::straight(), CurvatureLaw::circular(0.1)],
    );
    assert!(law.is_well_formed());
    assert_eq!(Intrinsic2::new(frame(), law, 20.0).total_turning(), None);
}

/// A single piece with no seams is the law itself; the wrapper must not
/// change the answer.
#[test]
fn a_single_piece_law_matches_the_bare_law() {
    let bare = CurvatureLaw::circular(0.05);
    let wrapped = CurvatureLaw::piecewise(vec![], vec![bare.clone()]);
    assert!(wrapped.is_well_formed());

    let a = Intrinsic2::new(frame(), bare, 12.0)
        .total_turning()
        .expect("finite");
    let b = Intrinsic2::new(frame(), wrapped, 12.0)
        .total_turning()
        .expect("finite");
    assert_eq!(a, b);
}

/// An empty piecewise law has nothing to integrate. It is well formed
/// (zero pieces need zero seams) and turns nothing.
#[test]
fn an_empty_piecewise_law_turns_nothing() {
    let law = CurvatureLaw::piecewise(vec![], vec![]);
    assert!(law.is_well_formed());
    assert!(law.is_straight());
    assert!(law.is_constant());
    let turning = Intrinsic2::new(frame(), law, 10.0)
        .total_turning()
        .expect("an empty law still integrates, to nothing");
    assert_eq!(turning, 0.0);
}

/// Pieces may themselves be composite laws: a straight, then the
/// sine-corrected transition, then an arc. Nesting must stay closed.
#[test]
fn a_piece_may_itself_be_a_composite_law() {
    let delta = 1.0 / 300.0;
    let (transition, arc) = (60.0, 40.0);
    let law = CurvatureLaw::piecewise(
        vec![transition],
        vec![
            CurvatureLaw::sine_corrected_transition(0.0, delta, transition),
            CurvatureLaw::circular(delta),
        ],
    );
    let turning = Intrinsic2::new(frame(), law, transition + arc)
        .total_turning()
        .expect("well formed");
    let expected = delta * transition / 2.0 + delta * arc;
    assert!((turning - expected).abs() < 1e-15);
}
