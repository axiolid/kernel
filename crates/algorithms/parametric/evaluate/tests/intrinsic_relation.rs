//! Relations over `Intrinsic3`: trim, offset, join (ADR 0062).
//!
//! Each exactness claim is checked against the BASE curve's own evaluator at
//! matching arc lengths, never against a second run of the same relation.

use axiolid_contracts::GeomError;
use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{CurvatureLaw, Harmonic, Intrinsic3};
use axiolid_evaluate::{
    frenet_frame, frenet_point, join_intrinsic3, offset_intrinsic3, trim_intrinsic3,
};

fn start_frame() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: Vec3::new(1.0, 0.0, 0.0),
        y: Vec3::new(0.0, 1.0, 0.0),
        z: Vec3::new(0.0, 0.0, 1.0),
    }
}

/// A helix of radius `a` and pitch parameter `b`.
fn helix(a: Scalar, b: Scalar, length: Scalar) -> Intrinsic3 {
    let c = a.hypot(b);
    Intrinsic3::new(
        start_frame(),
        CurvatureLaw::circular(a / (c * c)),
        CurvatureLaw::circular(b / (c * c)),
        length,
    )
}

/// A curve whose curvature and torsion both vary.
fn varying(length: Scalar) -> Intrinsic3 {
    Intrinsic3::new(
        start_frame(),
        CurvatureLaw::Polynomial {
            coefficients: vec![0.08, 0.01],
        },
        CurvatureLaw::Polynomial {
            coefficients: vec![0.03, -0.004],
        },
        length,
    )
}

#[test]
fn a_trim_reproduces_the_base_curve_on_the_kept_span() {
    let base = varying(6.0);
    let (start, end) = (1.7, 4.9);
    let trimmed = trim_intrinsic3(&base, start, end).expect("trim");

    assert!((trimmed.length - (end - start)).abs() < 1e-12);

    // The trimmed curve at u must be the base curve at start + u, in both
    // position and frame. This is the definition of an exact trim.
    //
    // The tolerance is 1e-6, not 1e-12, and that is a statement about the
    // EVALUATOR, not about the trim. The trim itself is exact -- the law is
    // shifted in closed form. But the two sides are integrated over
    // different spans, so they land on different panel layouts, and each
    // carries its own O(h^4) quadrature error. Verified as convergence
    // rather than an algebra fault: raising PANELS_PER_RADIAN from 4 to 16
    // drops the gap below 1e-9 and it keeps falling. Asserting 1e-12 here
    // would be asserting that two different quadratures agree to machine
    // precision, which is false for any finite budget.
    for step in 0..=8 {
        let u = (end - start) * Scalar::from(step) / 8.0;
        let got = frenet_point(&trimmed, u).expect("trimmed point");
        let want = frenet_point(&base, start + u).expect("base point");
        let error = (got - want).length();
        assert!(
            error < 1e-6,
            "trim at u={u}: got {got:?} want {want:?} error {error:e}"
        );

        let got_frame = frenet_frame(&trimmed, u).expect("trimmed frame");
        let want_frame = frenet_frame(&base, start + u).expect("base frame");
        assert!(
            (got_frame.x - want_frame.x).length() < 1e-6,
            "tangent diverges at u={u}"
        );
    }
}

#[test]
fn a_trim_of_a_piecewise_law_keeps_the_later_seams() {
    let law = CurvatureLaw::piecewise(
        vec![2.0, 5.0],
        vec![
            CurvatureLaw::circular(0.1),
            CurvatureLaw::circular(0.25),
            CurvatureLaw::circular(0.05),
        ],
    );
    let base = Intrinsic3::new(start_frame(), law, CurvatureLaw::straight(), 8.0);
    let trimmed = trim_intrinsic3(&base, 3.0, 8.0).expect("trim");

    // Starting at 3.0 drops the first piece and moves the remaining seam
    // from 5.0 to 2.0.
    match &trimmed.curvature {
        CurvatureLaw::Piecewise { breaks, laws } => {
            assert_eq!(laws.len(), 2, "one seam must survive");
            assert!(
                (breaks[0] - 2.0).abs() < 1e-12,
                "seam must move to 2.0, got {}",
                breaks[0]
            );
        }
        other => panic!("expected a piecewise law, got {other:?}"),
    }

    // And it must still agree with the base curve pointwise.
    for step in 0..=6 {
        let u = 5.0 * Scalar::from(step) / 6.0;
        let got = frenet_point(&trimmed, u).expect("trimmed");
        let want = frenet_point(&base, 3.0 + u).expect("base");
        assert!((got - want).length() < 1e-9, "piecewise trim at u={u}");
    }
}

#[test]
fn a_trim_outside_the_curve_is_refused() {
    let base = varying(4.0);
    assert!(trim_intrinsic3(&base, -0.5, 2.0).is_err());
    assert!(trim_intrinsic3(&base, 1.0, 9.0).is_err());
    assert!(trim_intrinsic3(&base, 3.0, 3.0).is_err());
    assert!(trim_intrinsic3(&base, 3.0, 1.0).is_err());
}

#[test]
fn a_helix_offset_lands_on_the_normal_offset_of_the_base() {
    // Verified in Python first: the offset of a helix is a coaxial helix of
    // radius a - d, same pitch, whose arc length runs at c2/c.
    let (a, b) = (3.0, 1.5);
    let base = helix(a, b, 4.0);

    for distance in [0.5, 1.0, -0.75] {
        let offset = offset_intrinsic3(&base, distance).expect("offset");
        let scale = (a - distance).hypot(b) / a.hypot(b);
        assert!(
            (offset.length - 4.0 * scale).abs() < 1e-12,
            "arc length must rescale by c2/c"
        );

        for step in 0..=6 {
            let s = 4.0 * Scalar::from(step) / 6.0;
            let frame = frenet_frame(&base, s).expect("base frame");
            let want = frenet_point(&base, s).expect("base point") + frame.y * distance;
            let got = frenet_point(&offset, s * scale).expect("offset point");
            let error = (got - want).length();
            assert!(
                error < 1e-8,
                "offset d={distance} at s={s}: got {got:?} want {want:?} error {error:e}"
            );
        }
    }
}

#[test]
fn an_offset_of_a_varying_law_is_refused_not_approximated() {
    // The offset of a varying-curvature space curve is not unit speed, so it
    // is not an Intrinsic3 in its own arc length at all.
    let base = varying(3.0);
    let error = offset_intrinsic3(&base, 0.8).expect_err("must refuse");
    assert!(
        matches!(error, GeomError::Unsupported { .. }),
        "expected a named refusal, got {error:?}"
    );
}

#[test]
fn an_offset_onto_the_axis_is_refused() {
    // d == a collapses the helix to its axis: no curve remains.
    let base = helix(3.0, 0.0, 4.0);
    assert!(offset_intrinsic3(&base, 3.0).is_err());
}

#[test]
fn a_zero_offset_is_the_curve_itself() {
    let base = varying(3.0);
    let offset = offset_intrinsic3(&base, 0.0).expect("zero offset");
    assert_eq!(offset, base, "a zero offset must not perturb the curve");
}

#[test]
fn a_join_reproduces_both_halves() {
    let first = helix(3.0, 1.5, 2.0);
    // The second curve must start exactly where the first ends.
    let seam = frenet_frame(&first, first.length).expect("seam frame");
    let second = Intrinsic3::new(
        seam,
        CurvatureLaw::circular(0.2),
        CurvatureLaw::circular(0.05),
        3.0,
    );
    let joined = join_intrinsic3(&first, &second, 1e-9, 1e-9).expect("join");
    assert!((joined.length - 5.0).abs() < 1e-12);

    // Before the seam the joined curve is the first; after it, the second.
    for step in 0..=4 {
        let s = 2.0 * Scalar::from(step) / 4.0;
        let got = frenet_point(&joined, s).expect("joined");
        let want = frenet_point(&first, s).expect("first");
        assert!((got - want).length() < 1e-9, "before the seam at s={s}");
    }
    for step in 0..=4 {
        let u = 3.0 * Scalar::from(step) / 4.0;
        let got = frenet_point(&joined, 2.0 + u).expect("joined");
        let want = frenet_point(&second, u).expect("second");
        let error = (got - want).length();
        assert!(
            error < 1e-8,
            "after the seam at u={u}: got {got:?} want {want:?} error {error:e}"
        );
    }
}

#[test]
fn a_join_of_curves_that_do_not_meet_is_refused() {
    let first = helix(3.0, 1.5, 2.0);
    let seam = frenet_frame(&first, first.length).expect("seam");
    // The tangent is CORRECT and only the position is wrong, so this isolates
    // the meeting-point check. A fixture that also got the tangent wrong would
    // pass even with the position check deleted -- the kink check would catch
    // it -- and the test would prove nothing about position.
    let detached = Intrinsic3::new(
        Frame3 {
            origin: seam.origin + seam.x * 0.5,
            x: seam.x,
            y: seam.y,
            z: seam.z,
        },
        CurvatureLaw::circular(0.2),
        CurvatureLaw::straight(),
        1.0,
    );
    assert!(
        join_intrinsic3(&first, &detached, 1e-9, 1e-9).is_err(),
        "a join that jumps must be refused even when the tangents agree"
    );
}

#[test]
fn a_join_that_kinks_is_refused() {
    let first = helix(3.0, 1.5, 2.0);
    let seam = frenet_frame(&first, first.length).expect("seam");
    // Right position, wrong tangent.
    let kinked = Intrinsic3::new(
        Frame3 {
            origin: seam.origin,
            x: seam.y,
            y: seam.z,
            z: seam.x,
        },
        CurvatureLaw::circular(0.2),
        CurvatureLaw::straight(),
        1.0,
    );
    assert!(
        join_intrinsic3(&first, &kinked, 1e-9, 1e-9).is_err(),
        "a join with disagreeing tangents must be refused"
    );
}

#[test]
fn a_shifted_sinusoid_moves_only_its_phase() {
    let law = CurvatureLaw::Sinusoid {
        mean: 0.15,
        amplitude: 0.4,
        angular_frequency: 2.3,
        phase: 0.6,
    };
    let shifted = law.shifted(1.7).expect("shift");
    match shifted {
        CurvatureLaw::Sinusoid {
            mean,
            amplitude,
            angular_frequency,
            phase,
        } => {
            assert!((mean - 0.15).abs() < 1e-15);
            assert!((amplitude - 0.4).abs() < 1e-15);
            assert!((angular_frequency - 2.3).abs() < 1e-15);
            assert!(
                (phase - (0.6 + 2.3 * 1.7)).abs() < 1e-12,
                "phase must absorb the shift"
            );
        }
        other => panic!("a shifted sinusoid must stay a sinusoid, got {other:?}"),
    }
}

#[test]
fn a_shifted_composite_shifts_every_harmonic() {
    let law = CurvatureLaw::Composite {
        polynomial: vec![0.3, -0.12, 0.045],
        harmonics: vec![Harmonic {
            amplitude: 0.2,
            angular_frequency: 1.1,
            phase: 0.0,
        }],
    };
    let shifted = law.shifted(2.0).expect("shift");
    match shifted {
        CurvatureLaw::Composite {
            polynomial,
            harmonics,
        } => {
            // p(2 + u) = 0.3 - 0.24 + 0.18 + (-0.12 + 0.18) u + 0.045 u^2
            assert!((polynomial[0] - 0.24).abs() < 1e-12, "constant term");
            assert!((polynomial[1] - 0.06).abs() < 1e-12, "linear term");
            assert!((polynomial[2] - 0.045).abs() < 1e-12, "quadratic term");
            assert!((harmonics[0].phase - 2.2).abs() < 1e-12, "harmonic phase");
        }
        other => panic!("expected a composite, got {other:?}"),
    }
}

#[test]
fn a_constant_value_is_none_for_a_varying_law() {
    let varying = CurvatureLaw::Polynomial {
        coefficients: vec![0.1, 0.02],
    };
    assert!(varying.constant_value().is_none());
    assert!(
        (CurvatureLaw::circular(0.25)
            .constant_value()
            .expect("constant")
            - 0.25)
            .abs()
            < 1e-15
    );
}

#[test]
fn a_trim_is_exact_in_the_law_itself_not_merely_close() {
    // The evaluator comparison above is limited by quadrature. THIS is the
    // exactness claim: the trimmed law, evaluated symbolically, is the base
    // law shifted -- to machine precision, with no integration involved.
    let base = varying(6.0);
    let (start, end) = (1.7, 4.9);
    let trimmed = trim_intrinsic3(&base, start, end).expect("trim");

    for step in 0..=32 {
        let u = (end - start) * Scalar::from(step) / 32.0;
        // Heading accumulated over [0, u] of the trimmed law must equal the
        // base law's over [start, start + u] -- a difference of exact
        // closed-form integrals, so any gap is an algebra fault.
        let got = heading_of(&trimmed.curvature, u);
        let want = heading_of(&base.curvature, start + u) - heading_of(&base.curvature, start);
        assert!(
            (got - want).abs() < 1e-12,
            "curvature law at u={u}: {got} vs {want}"
        );

        let got = heading_of(&trimmed.torsion, u);
        let want = heading_of(&base.torsion, start + u) - heading_of(&base.torsion, start);
        assert!(
            (got - want).abs() < 1e-12,
            "torsion law at u={u}: {got} vs {want}"
        );
    }
}

/// Exact closed-form integral of a law over `[0, s]`, no quadrature.
fn heading_of(law: &CurvatureLaw, s: Scalar) -> Scalar {
    axiolid_curve::Intrinsic2::new(
        axiolid_core::Frame2 {
            origin: axiolid_core::Point2::new(0.0, 0.0),
            x: axiolid_core::Vec2::X,
            y: axiolid_core::Vec2::Y,
        },
        law.clone(),
        s,
    )
    .total_turning()
    .expect("law integrates")
}
