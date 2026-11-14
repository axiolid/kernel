//! The rates of a section frame law (#406, ADR 0082 amendment):
//! `FrameLaw::rates` bounds the forward rate, the rise, the two bends and
//! the twist of `E = [T, L, U]` over a span.
//!
//! The references are independent of the jets the rates are built from:
//! central differences (`h = 1e-5`) of the kernel's closed-form frames,
//! `tau = U(s) . (L(s + h) - L(s - h)) / 2h` and `rise = (O(s + h) - O(s -
//! h)) / 2h . U(s)`. Each sampled value must lie inside the reported range,
//! and the middle value within its own error of the reference there.

use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BSplineCurve2, BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2,
    Circle2, Circle3, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Ellipse2, Intrinsic2,
    KnotSpec, Line2, OffsetFrame, OffsetLaw, PathOffsets, SeamSide,
};
use axiolid_evaluate::station::{
    offset_pieces, station_section2, station_section3, CompositeBasis, FrameLaw, FrameRates,
    FrameReading, RateBounds, SectionFrame, StationCurve, StationPiece,
};

const TOLERANCE: Scalar = 1e-7;
const H: Scalar = 1e-5;
/// The step of the differences of frames read by quadrature, whose own
/// quadrature errors (about `1e-9`) a smaller one would amplify.
const H_READ: Scalar = 1e-3;

fn contains(rate: &RateBounds, value: Scalar) -> bool {
    rate.range[0] - 1e-7 <= value && value <= rate.range[1] + 1e-7
}

/// `U . L'` by a central difference of the frame's lateral axis.
fn twist_at(frame: &dyn Fn(Scalar) -> SectionFrame, s: Scalar, h: Scalar) -> Scalar {
    let rate = (frame(s + h).lateral - frame(s - h).lateral) / (2.0 * h);
    frame(s).up.dot(rate)
}

/// The forward rate, the bends, the twist and the rise sampled at `samples + 1` points over
/// `[lo, hi]` lie within the rates' ranges, and their middle values within
/// their errors (and `slack`) of the references.
fn check_sampled(
    rates: &FrameRates,
    frame: &dyn Fn(Scalar) -> SectionFrame,
    (lo, hi): (Scalar, Scalar),
    (h, slack): (Scalar, Scalar),
    samples: usize,
) {
    let at = |s: Scalar| {
        let (before, here, after) = (frame(s - h), frame(s), frame(s + h));
        let bend = (after.tangent - before.tangent) / (2.0 * h);
        let velocity = (after.point - before.point) / (2.0 * h);
        [
            (&rates.forward, velocity.dot(here.tangent), "forward"),
            (&rates.lateral, here.lateral.dot(bend), "lateral"),
            (&rates.up, here.up.dot(bend), "up"),
            (
                &rates.twist,
                here.up.dot((after.lateral - before.lateral) / (2.0 * h)),
                "twist",
            ),
            (&rates.rise, velocity.dot(here.up), "rise"),
        ]
    };
    for i in 0..=samples {
        let s = lo + (hi - lo) * i as Scalar / samples as Scalar;
        for (rate, reference, name) in at(s) {
            assert!(
                rate.range[0] - slack <= reference && reference <= rate.range[1] + slack,
                "{name} {reference} at {s}: {:?}",
                rate.range
            );
        }
    }
    for (rate, reference, name) in at(0.5 * (lo + hi)) {
        assert!(
            (rate.middle - reference).abs() <= rate.error + slack,
            "{name} {} against {reference}, error {}",
            rate.middle,
            rate.error
        );
    }
}

fn check_with(
    rates: &FrameRates,
    frame: &dyn Fn(Scalar) -> SectionFrame,
    span: (Scalar, Scalar),
    steps: (Scalar, Scalar),
) {
    check_sampled(rates, frame, span, steps, 40);
}

/// [`check_with`] for the frames read in closed form.
fn check(rates: &FrameRates, frame: &dyn Fn(Scalar) -> SectionFrame, span: (Scalar, Scalar)) {
    check_with(rates, frame, span, (H, 1e-7));
}

fn plan_of(frame: SectionFrame) -> SectionFrame {
    frame.plan().unwrap()
}

fn circle(radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(1.0, 2.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

#[test]
fn a_level_circle_neither_twists_nor_rises() {
    let curve = circle(10.0);
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let span = (2.0, 6.0);
        let rates = law.rates(span.0, span.1, reading).unwrap();
        for rate in [&rates.twist, &rates.rise] {
            assert!(contains(rate, 0.0));
            assert!(rate.middle.abs() < 1e-9);
        }
        // The bend is the curvature.
        assert!(contains(&rates.lateral, 0.1) && rates.lateral.range[0] > 0.0);
        check(&rates, &|s| station_section2(&curve, s).unwrap(), span);
    }
}

/// A circle of radius 5 whose plane is tilted 30 degrees about `x`.
fn tilted_circle() -> Curve3 {
    let (sin, cos) = (30.0_f64).to_radians().sin_cos();
    let (x, y) = (Vec3::X, Vec3::new(0.0, cos, sin));
    Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 1.0),
            x,
            y,
            z: x.cross(y),
        },
        radius: 5.0,
    })
}

#[test]
fn a_tilted_circle_twists_in_the_section_reading() {
    let curve = tilted_circle();
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let span = (1.0, 3.0);
    let rates = law.rates(span.0, span.1, FrameReading::Section).unwrap();
    // Narrow enough to tell the sign of the twist.
    let narrow = law.rates(1.9, 2.1, FrameReading::Section).unwrap();
    assert!(narrow.twist.range[0] > 0.0 || narrow.twist.range[1] < 0.0);
    assert!(contains(&rates.rise, 0.0) && rates.rise.range == [0.0, 0.0]);
    check(&rates, &|s| station_section3(&curve, s).unwrap(), span);
    // The plan reading has the plan frame: no twist.
    let plan = law.rates(span.0, span.1, FrameReading::Plan).unwrap();
    assert_eq!(plan.twist.range, [0.0, 0.0]);
    assert_eq!(plan.up.range, [0.0, 0.0]);
    check(
        &plan,
        &|s| plan_of(station_section3(&curve, s).unwrap()),
        span,
    );
}

const CANT_B: Scalar = 1.5;

/// A banked curve over `plan` that rises on the profile `elevation` and
/// cants linearly from level to `0.1` over its first 60.
fn canted(plan: Curve2, elevation: ElevationLaw, convention: BankConvention) -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(plan, elevation),
        CantLaw::new(vec![CantPiece::linear(60.0, 0.0, 0.1)]),
        CantLaw::zero(60.0),
        CANT_B,
        convention,
    ))
}

fn steep() -> ElevationLaw {
    ElevationLaw::Polynomial {
        coefficients: vec![0.0, 0.3],
    }
}

/// The twist of a canted level circle is the roll's rate in closed form.
fn roll_rate_on_level(convention: BankConvention, grade: Scalar, s: Scalar) -> Scalar {
    let rate = 0.1 / 60.0 / CANT_B;
    let scale = match convention {
        BankConvention::TangentRotation => 1.0,
        _ => grade.hypot(1.0),
    };
    let g = 0.1 * s / 60.0 / CANT_B * scale;
    rate * scale / (1.0 - g * g).sqrt()
}

#[test]
fn a_canted_transition_twists_by_its_roll_rate() {
    let plan = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 100.0,
    });
    for (convention, elevation, grade) in [
        (
            BankConvention::TangentRotation,
            ElevationLaw::level(0.0),
            0.0,
        ),
        (BankConvention::TangentRotation, steep(), 0.3),
        (BankConvention::VerticalRise, steep(), 0.3),
    ] {
        let curve = canted(plan.clone(), elevation, convention);
        let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
        let span = (14.0, 16.0);
        let rates = law.rates(span.0, span.1, FrameReading::Section).unwrap();
        // Narrow enough to tell the sign of the twist.
        let narrow = law.rates(14.9, 15.1, FrameReading::Section).unwrap();
        assert!(
            narrow.twist.range[0] > 0.0,
            "{convention:?} grade {grade}: {:?}",
            narrow.twist
        );
        check(&rates, &|s| station_section3(&curve, s).unwrap(), span);
        if grade == 0.0 {
            // Level: the twist is the roll's rate alone.
            let closed = roll_rate_on_level(convention, grade, 15.0);
            assert!((rates.twist.middle - closed).abs() <= rates.twist.error + 1e-9);
        }
        let plan = law.rates(span.0, span.1, FrameReading::Plan).unwrap();
        assert_eq!(plan.twist.range, [0.0, 0.0]);
    }
}

fn clothoid() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::new(3.0, -1.0),
            x: Vec2::new(0.6, 0.8),
            y: Vec2::new(-0.8, 0.6),
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
        60.0,
    ))
}

#[test]
fn a_canted_clothoid_twists_by_its_roll_rate() {
    for (convention, elevation) in [
        (BankConvention::TangentRotation, ElevationLaw::level(0.0)),
        (BankConvention::VerticalRise, steep()),
    ] {
        let curve = canted(clothoid(), elevation, convention);
        let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
        let span = (14.0, 16.0);
        let rates = law.rates(span.0, span.1, FrameReading::Section).unwrap();
        let narrow = law.rates(14.9, 15.1, FrameReading::Section).unwrap();
        assert!(narrow.twist.range[0] > 0.0, "{convention:?}");
        check(&rates, &|s| station_section3(&curve, s).unwrap(), span);
    }
}

/// `z = s^2 / 100` over a line along `x`.
fn sag() -> Curve3 {
    Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::Polynomial {
            coefficients: vec![0.0, 0.0, 0.01],
        },
    ))
}

fn graded() -> Curve3 {
    Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::Polynomial {
            coefficients: vec![0.0, 0.5],
        },
    ))
}

#[test]
fn a_graded_curve_rises_in_the_plan_reading_only() {
    // The sag's rise is `z' = s / 50` and the graded line's is its grade.
    for (curve, rise) in [
        (sag(), (&|s: Scalar| 0.02 * s) as &dyn Fn(Scalar) -> Scalar),
        (graded(), &|_| 0.5),
    ] {
        let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
        let span = (5.0, 15.0);
        let plan = law.rates(span.0, span.1, FrameReading::Plan).unwrap();
        for i in 0..=20 {
            let s = span.0 + (span.1 - span.0) * Scalar::from(i) / 20.0;
            assert!(contains(&plan.rise, rise(s)), "rise at {s}");
        }
        assert!((plan.rise.middle - rise(10.0)).abs() <= plan.rise.error + 1e-12);
        check(
            &plan,
            &|s| plan_of(station_section3(&curve, s).unwrap()),
            span,
        );
        // The section's frame follows the tangent: nothing rises.
        let section = law.rates(span.0, span.1, FrameReading::Section).unwrap();
        assert_eq!(section.rise.range, [0.0, 0.0]);
        assert_eq!(section.rise.middle, 0.0);
        assert_eq!(section.rise.error, 0.0);
        check(&section, &|s| station_section3(&curve, s).unwrap(), span);
    }
    // The sag rises, narrowly enough to tell the sign.
    let curve = sag();
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let narrow = law.rates(9.0, 11.0, FrameReading::Plan).unwrap();
    assert!(narrow.rise.range[0] > 0.0);
}

#[test]
fn a_reversed_and_a_placed_piece_keep_their_twist() {
    let curve = tilted_circle();
    let canted_curve = canted(
        Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 100.0,
        }),
        ElevationLaw::level(0.0),
        BankConvention::TangentRotation,
    );
    // Reversed: the tilted circle's, and the canted transition's.
    for curve in [&curve, &canted_curve] {
        let piece = StationPiece::between(StationCurve::Three(curve), 0.0, 20.0)
            .unwrap()
            .reversed();
        let basis = CompositeBasis::new(vec![piece]).unwrap();
        let law = FrameLaw::composite(&basis, TOLERANCE).unwrap();
        let span = (7.0, 9.0);
        let rates = law.rates(span.0, span.1, FrameReading::Section).unwrap();
        let narrow = law.rates(7.9, 8.1, FrameReading::Section).unwrap();
        assert!(
            narrow.twist.range[0] > 0.0 || narrow.twist.range[1] < 0.0,
            "{:?}",
            narrow.twist
        );
        check(
            &rates,
            &|s| piece.section_on(s, SeamSide::Outgoing).unwrap(),
            span,
        );
    }
    // Placed by a turn about Z (and moved): `+Z` is kept.
    let (sin, cos) = (35.0_f64).to_radians().sin_cos();
    let turn = Transform3::from_cols(
        Vec3::new(cos, sin, 0.0),
        Vec3::new(-sin, cos, 0.0),
        Vec3::Z,
        Point3::new(4.0, -2.0, 1.0),
    );
    for (curve, at) in [(&curve, 2.0), (&canted_curve, 8.0)] {
        let piece = StationPiece::between(StationCurve::Three(curve), 0.0, 20.0)
            .unwrap()
            .placed(turn, true);
        let law = FrameLaw::atomic(StationCurve::Three(curve), Some(turn), TOLERANCE).unwrap();
        let span = (at - 1.0, at + 1.0);
        let rates = law.rates(span.0, span.1, FrameReading::Section).unwrap();
        let narrow = law
            .rates(at - 0.1, at + 0.1, FrameReading::Section)
            .unwrap();
        assert!(narrow.twist.range[0] > 0.0 || narrow.twist.range[1] < 0.0);
        check(
            &rates,
            &|s| piece.section_on(s, SeamSide::Outgoing).unwrap(),
            span,
        );
    }
}

#[test]
fn bend_is_the_projection_of_the_rates() {
    let curve = tilted_circle();
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let rates = law.rates(1.0, 3.0, reading).unwrap();
        let bend = law.bend(1.0, 3.0, reading).unwrap();
        assert_eq!(bend.forward, rates.forward.range);
        assert_eq!(bend.lateral, rates.lateral.range);
        assert_eq!(bend.up, rates.up.range);
        assert_eq!(
            bend.middle,
            [rates.forward.middle, rates.lateral.middle, rates.up.middle]
        );
        assert_eq!(
            bend.middle_error,
            [rates.forward.error, rates.lateral.error, rates.up.error]
        );
    }
}

// --- frames read by quadrature ----------------------------------------------------

/// The ellipse of semi-axes 20 and 10 about `(30, 0)`.
fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(30.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 20.0,
        semi_axis_y: 10.0,
    })
}

/// The ellipse's parameter at arc length `s` from `t = 0`, by RK4 on
/// `dt/ds = 1 / |P'(t)|` with step `1e-3`.
fn ellipse_parameter(s: Scalar) -> Scalar {
    let speed = |t: Scalar| (400.0 * t.sin().powi(2) + 100.0 * t.cos().powi(2)).sqrt();
    let steps = (s / 1e-3).ceil().max(1.0) as usize;
    let h = s / steps as Scalar;
    let mut t = 0.0;
    for _ in 0..steps {
        let k1 = 1.0 / speed(t);
        let k2 = 1.0 / speed(t + 0.5 * h * k1);
        let k3 = 1.0 / speed(t + 0.5 * h * k2);
        let k4 = 1.0 / speed(t + h * k3);
        t += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    t
}

/// The ellipse's signed curvature `a b / (a^2 sin^2 t + b^2 cos^2 t)^(3/2)`.
fn ellipse_curvature(t: Scalar) -> Scalar {
    200.0 / (400.0 * t.sin().powi(2) + 100.0 * t.cos().powi(2)).powf(1.5)
}

#[test]
fn an_ellipse_bends_at_its_curvature() {
    let curve = ellipse();
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        // Narrow enough to tell the sign of the bend, at the sharpest end.
        let sharp = law.rates(0.0, 1.0, reading).unwrap();
        assert!(sharp.lateral.range[0] > 0.0, "{:?}", sharp.lateral);
        assert!(contains(&sharp.lateral, 0.2) && contains(&sharp.forward, 1.0));
        for span in [(0.05, 1.05), (10.0, 14.0), (30.0, 33.0)] {
            let rates = law.rates(span.0, span.1, reading).unwrap();
            for i in 0..=20 {
                let s = span.0 + (span.1 - span.0) * Scalar::from(i) / 20.0;
                let kappa = ellipse_curvature(ellipse_parameter(s));
                assert!(contains(&rates.lateral, kappa), "{kappa} at {s}");
            }
            let mid = 0.5 * (span.0 + span.1);
            let kappa = ellipse_curvature(ellipse_parameter(mid));
            assert!(
                (rates.lateral.middle - kappa).abs() <= rates.lateral.error + 1e-9,
                "{} against {kappa}, error {}",
                rates.lateral.middle,
                rates.lateral.error
            );
            // A plane curve neither twists, rises nor bends upward, and
            // its error is the certified one, not a loose estimate.
            for rate in [&rates.twist, &rates.rise, &rates.up] {
                assert!(contains(rate, 0.0));
            }
            assert!(rates.lateral.error < 1e-6, "{}", rates.lateral.error);
            check_with(
                &rates,
                &|s| station_section2(&curve, s).unwrap(),
                span,
                (H_READ, 1e-6),
            );
        }
    }
}

#[test]
fn a_reversed_ellipse_piece_bends_the_other_way() {
    let curve = ellipse();
    let piece = StationPiece::between(StationCurve::Two(&curve), 0.0, 20.0)
        .unwrap()
        .reversed();
    let basis = CompositeBasis::new(vec![piece]).unwrap();
    let law = FrameLaw::composite(&basis, TOLERANCE).unwrap();
    // The reversed piece's distance `s` is the forward arc length `20 - s`.
    let span = (18.9, 19.9);
    let rates = law.rates(span.0, span.1, FrameReading::Section).unwrap();
    assert!(rates.lateral.range[1] < 0.0, "{:?}", rates.lateral);
    for i in 0..=20 {
        let s = span.0 + (span.1 - span.0) * Scalar::from(i) / 20.0;
        let kappa = -ellipse_curvature(ellipse_parameter(20.0 - s));
        assert!(contains(&rates.lateral, kappa), "{kappa} at {s}");
    }
    check_with(
        &rates,
        &|s| piece.section_on(s, SeamSide::Outgoing).unwrap(),
        span,
        (H_READ, 1e-6),
    );
}

/// A 3D cubic B-spline climbing gently, about 30 long.
fn spline3() -> Curve3 {
    Curve3::BSpline(BSplineCurve {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(8.0, 2.0, 0.5),
            Point3::new(16.0, -3.0, 1.5),
            Point3::new(24.0, 2.0, 2.0),
            Point3::new(30.0, 6.0, 2.2),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

#[test]
fn a_3d_b_spline_bends_and_twists_as_its_frames_turn() {
    let curve = spline3();
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    for (span, reading) in [
        ((3.0, 5.0), FrameReading::Section),
        ((12.0, 14.0), FrameReading::Section),
        ((20.0, 22.0), FrameReading::Section),
        ((12.0, 14.0), FrameReading::Plan),
    ] {
        let rates = law.rates(span.0, span.1, reading).unwrap();
        let plan = reading == FrameReading::Plan;
        let frame = |s: Scalar| {
            let frame = station_section3(&curve, s).unwrap();
            if plan {
                plan_of(frame)
            } else {
                frame
            }
        };
        check_with(&rates, &frame, span, (H_READ, 1e-6));
        if !plan {
            assert!(contains(&rates.rise, 0.0) && rates.rise.range == [0.0, 0.0]);
        }
    }
    // A space curve twists its reference-up frame: narrowly enough to see,
    // where it twists most.
    let frame = |s: Scalar| station_section3(&curve, s).unwrap();
    let (at, _) = (1..290)
        .map(|i| {
            (
                0.1 * Scalar::from(i),
                twist_at(&frame, 0.1 * Scalar::from(i), H_READ),
            )
        })
        .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        .unwrap();
    let twisting = law
        .rates(at - 0.004, at + 0.004, FrameReading::Section)
        .unwrap();
    assert!(
        twisting.twist.range[0] > 0.0 || twisting.twist.range[1] < 0.0,
        "{:?} at {at}",
        twisting.twist
    );
}

/// A quadratic plan B-spline whose curvature jumps at its interior knot.
fn quadratic() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(12.0, 0.0),
            Point2::new(12.0, 12.0),
            Point2::new(24.0, 14.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 1, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

#[test]
fn a_span_across_a_curvature_break_has_no_rates() {
    let curve = quadratic();
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    // Spans of 0.05 along the curve: the ones across the knot have bounds
    // but no rates, every other has both.
    let mut across = Vec::new();
    let mut s = 1.0;
    while law.bounds(s, s + 0.05, FrameReading::Section).is_some() {
        for reading in [FrameReading::Section, FrameReading::Plan] {
            let rates = law.rates(s, s + 0.05, reading);
            if rates.is_none() {
                across.push(s);
            }
            assert_eq!(rates.is_some(), law.bend(s, s + 0.05, reading).is_some());
        }
        s += 0.05;
    }
    assert!(s > 20.0, "scanned to {s}");
    // The knot: one span, found by both readings.
    assert!(
        (2..=4).contains(&across.len()) && across[across.len() - 1] - across[0] < 0.06,
        "{across:?}"
    );
    assert!(law.has_bend(across[0], across[0] + 0.05));
}

#[test]
fn chains_have_no_rates() {
    let chain = Curve2::Chain(Chain2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        vec![ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
            length: 60.0,
        }],
    ));
    let law = FrameLaw::atomic(StationCurve::Two(&chain), None, TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        assert!(law.rates(5.0, 8.0, reading).is_none());
    }
}

// --- offsets beside curves read by quadrature -----------------------------------

/// The ellipse's parameter at arc length `s` of its offset by `distance`
/// (positive to the left), by RK4 on `dt/ds = 1 / |P'(t)|` with
/// `|P'| = |B'| - distance a b / |B'|^2`.
fn offset_parameter(distance: Scalar, s: Scalar) -> Scalar {
    let speed = |t: Scalar| {
        let base = (400.0 * t.sin().powi(2) + 100.0 * t.cos().powi(2)).sqrt();
        base - distance * 200.0 / (base * base)
    };
    let steps = (s / 1e-3).ceil().max(1.0) as usize;
    let h = s / steps as Scalar;
    let mut t = 0.0;
    for _ in 0..steps {
        let k1 = 1.0 / speed(t);
        let k2 = 1.0 / speed(t + 0.5 * h * k1);
        let k3 = 1.0 / speed(t + 0.5 * h * k2);
        let k4 = 1.0 / speed(t + h * k3);
        t += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    t
}

#[test]
fn an_offset_beside_an_ellipse_bends_at_its_own_radius() {
    let plan = ellipse();
    let distance = -2.0;
    let base = StationPiece::between(StationCurve::Two(&plan), 0.0, 40.0).unwrap();
    let basis =
        CompositeBasis::new(offset_pieces(&[base], OffsetLaw::Planar { distance }).unwrap())
            .unwrap();
    let law = basis.frame_law(TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        // The curvature of the offset is `kappa / (1 - d kappa)`.
        let offset_curvature = |s: Scalar| {
            let kappa = ellipse_curvature(offset_parameter(distance, s));
            kappa / (1.0 - distance * kappa)
        };
        // Near the sharpest end the radius is `5 + 2 = 7`.
        let sharp = law.rates(0.0, 0.05, reading).unwrap();
        assert!(sharp.lateral.range[0] > 0.0, "{:?}", sharp.lateral);
        assert!((offset_curvature(0.0) - 1.0 / 7.0).abs() < 1e-12);
        for span in [(0.0, 1.0), (10.0, 13.0), (30.0, 33.0)] {
            let rates = law.rates(span.0, span.1, reading).unwrap();
            for i in 0..=20 {
                let s = span.0 + (span.1 - span.0) * Scalar::from(i) / 20.0;
                assert!(contains(&rates.lateral, offset_curvature(s)), "at {s}");
            }
            let mid = 0.5 * (span.0 + span.1);
            assert!(
                (rates.lateral.middle - offset_curvature(mid)).abs() <= rates.lateral.error + 1e-9,
                "{} against {}, error {}",
                rates.lateral.middle,
                offset_curvature(mid),
                rates.lateral.error
            );
            assert!(rates.lateral.error < 1e-6, "{}", rates.lateral.error);
            for rate in [&rates.twist, &rates.rise, &rates.up] {
                assert!(contains(rate, 0.0));
            }
        }
    }
}

fn widening() -> OffsetLaw {
    OffsetLaw::Linear {
        start: PathOffsets::new(0.5, 0.0, 0.0),
        end: PathOffsets::new(0.75, 0.2, 0.0),
        frame: OffsetFrame::Section,
    }
}

#[test]
fn offsets_by_distances_beside_an_ellipse_and_a_b_spline_bend() {
    let (plan, solid) = (ellipse(), spline3());
    let cases: [(StationCurve<'_>, OffsetLaw, Scalar, Scalar); 4] = [
        (StationCurve::Two(&plan), widening(), 2.0, 20.0),
        (
            StationCurve::Two(&plan),
            OffsetLaw::Linear {
                start: PathOffsets::new(0.5, 0.0, 0.0),
                end: PathOffsets::new(0.75, 0.2, 0.0),
                frame: OffsetFrame::Plan,
            },
            2.0,
            20.0,
        ),
        (StationCurve::Three(&solid), widening(), 1.0, 28.0),
        (
            StationCurve::Three(&solid),
            OffsetLaw::Directed {
                distance: 1.0,
                reference_direction: Vec3::Z,
            },
            0.0,
            28.0,
        ),
    ];
    for (k, (curve, offset_law, start, end)) in cases.into_iter().enumerate() {
        let base = StationPiece::between(curve, start, end).unwrap();
        // The first offset by distances runs along its base backwards.
        let base = if k == 0 { base.reversed() } else { base };
        let basis = CompositeBasis::new(offset_pieces(&[base], offset_law).unwrap()).unwrap();
        let law = basis.frame_law(TOLERANCE).unwrap();
        let section = |s: Scalar| basis.section_on(s, SeamSide::Outgoing).unwrap();
        // The plan reading differs from the section's where the curve
        // climbs or the law is read in the plan frame.
        let readings: &[FrameReading] = if k == 1 || k == 2 {
            &[FrameReading::Section, FrameReading::Plan]
        } else {
            &[FrameReading::Section]
        };
        for &reading in readings {
            let span = (3.0, 4.5);
            let rates = law
                .rates(span.0, span.1, reading)
                .unwrap_or_else(|| panic!("no rates for case {k} {reading:?}"));
            if reading == FrameReading::Plan {
                check_sampled(&rates, &|s| plan_of(section(s)), span, (H_READ, 1e-6), 4);
            } else {
                check_sampled(&rates, &section, span, (H_READ, 1e-6), 4);
            }
        }
    }
}

#[test]
fn a_constant_offset_beside_a_cubic_has_rates_only_between_its_knots() {
    let solid = spline3();
    let base = StationPiece::between(StationCurve::Three(&solid), 0.0, 28.0).unwrap();
    let basis = CompositeBasis::new(
        offset_pieces(
            &[base],
            OffsetLaw::Directed {
                distance: 1.0,
                reference_direction: Vec3::Z,
            },
        )
        .unwrap(),
    )
    .unwrap();
    let law = basis.frame_law(TOLERANCE).unwrap();
    // The offset's `P''` holds the base's `B'''`, which jumps at the
    // cubic's interior knot: a span across it has bounds but no rates,
    // however narrow.
    let none = |lo: Scalar, hi: Scalar| law.rates(lo, hi, FrameReading::Section).is_none();
    // A span too wide for the bounds has none either; halving tells it from
    // the knot, which stays without as narrow as the halving goes.
    let mut knots = Vec::new();
    for start in (3..13).map(|k| 2.0 * Scalar::from(k)) {
        let (mut lo, mut hi) = (start, start + 2.0);
        while none(lo, hi) && hi - lo > 0.01 {
            let mid = 0.5 * (lo + hi);
            if none(lo, mid) {
                hi = mid;
            } else if none(mid, hi) {
                lo = mid;
            } else {
                break;
            }
        }
        if hi - lo <= 0.01 && none(lo, hi) {
            knots.push((lo, hi));
        }
    }
    assert_eq!(knots.len(), 1, "{knots:?}");
    let (lo, hi) = knots[0];
    // The bounds stay, and the sides have rates again.
    assert!(law.bounds(lo, hi, FrameReading::Section).is_some());
    assert!(law.rates(lo - 0.2, lo, FrameReading::Section).is_some());
    assert!(law.rates(hi, hi + 0.2, FrameReading::Section).is_some());

    // The padded range of a span is held to the break too: spans 0.2 long
    // lose their rates before the span itself reaches the knot, and keep
    // none after it, by the padding on each side, `2 measure` of arc length
    // and more of distance at the speed floor.
    let knot = 0.5 * (lo + hi);
    // The smallest end `x` of a span `[x - 0.2, x]` without rates, and the
    // largest start `x` of a span `[x, x + 0.2]` without them.
    let (mut some, mut absent) = (knot - 0.1, knot + 0.1);
    for _ in 0..40 {
        let mid = 0.5 * (some + absent);
        if none(mid - 0.2, mid) {
            absent = mid;
        } else {
            some = mid;
        }
    }
    let first = absent;
    let (mut absent, mut some) = (knot - 0.1, knot + 0.1);
    for _ in 0..40 {
        let mid = 0.5 * (some + absent);
        if none(mid, mid + 0.2) {
            absent = mid;
        } else {
            some = mid;
        }
    }
    let last = absent;
    let measure = law
        .bounds(first - 0.2, first, FrameReading::Section)
        .unwrap()
        .measure;
    assert!(
        last - first >= 4.0 * measure,
        "no rates between {first} and {last}, measure {measure}"
    );
}

#[test]
fn rates_over_a_span_do_not_depend_on_what_was_asked_before() {
    // The reads at a span's ends are kept for its halves: asking the halves
    // first, or the whole twice, changes nothing.
    let plan = ellipse();
    let base = StationPiece::between(StationCurve::Two(&plan), 0.0, 40.0).unwrap();
    let basis =
        CompositeBasis::new(offset_pieces(&[base], OffsetLaw::Planar { distance: -2.0 }).unwrap())
            .unwrap();
    let ask = |law: &FrameLaw<'_, '_>, a: Scalar, b: Scalar| {
        law.rates(a, b, FrameReading::Section).unwrap()
    };
    let cold = basis.frame_law(TOLERANCE).unwrap();
    let first = ask(&cold, 4.0, 6.0);
    let warm = basis.frame_law(TOLERANCE).unwrap();
    let _ = (ask(&warm, 4.0, 5.0), ask(&warm, 5.0, 6.0));
    assert_eq!(first, ask(&warm, 4.0, 6.0));
    assert_eq!(first, ask(&warm, 4.0, 6.0));
}
