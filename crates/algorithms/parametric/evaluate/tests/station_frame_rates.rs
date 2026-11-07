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
    BankConvention, Banked3, CantLaw, CantPiece, Circle2, Circle3, CurvatureLaw, Curve2, Curve3,
    Elevated3, ElevationLaw, Intrinsic2, Line2, SeamSide,
};
use axiolid_evaluate::station::{
    station_section2, station_section3, CompositeBasis, FrameLaw, FrameRates, FrameReading,
    RateBounds, SectionFrame, StationCurve, StationPiece,
};

const TOLERANCE: Scalar = 1e-7;
const H: Scalar = 1e-5;

fn contains(rate: &RateBounds, value: Scalar) -> bool {
    rate.range[0] - 1e-7 <= value && value <= rate.range[1] + 1e-7
}

/// `U . L'` by a central difference of the frame's lateral axis.
fn twist_at(frame: &dyn Fn(Scalar) -> SectionFrame, s: Scalar) -> Scalar {
    let rate = (frame(s + H).lateral - frame(s - H).lateral) / (2.0 * H);
    frame(s).up.dot(rate)
}

/// `O' . U` by a central difference of the frame's point.
fn rise_at(frame: &dyn Fn(Scalar) -> SectionFrame, s: Scalar) -> Scalar {
    let rate = (frame(s + H).point - frame(s - H).point) / (2.0 * H);
    rate.dot(frame(s).up)
}

/// The twist and the rise sampled over `[lo, hi]` lie within the rates'
/// ranges, and their middle values within their errors of the references.
fn check(rates: &FrameRates, frame: &dyn Fn(Scalar) -> SectionFrame, (lo, hi): (Scalar, Scalar)) {
    for i in 0..=40 {
        let s = lo + (hi - lo) * Scalar::from(i) / 40.0;
        let (twist, rise) = (twist_at(frame, s), rise_at(frame, s));
        assert!(contains(&rates.twist, twist), "twist {twist} at {s}");
        assert!(contains(&rates.rise, rise), "rise {rise} at {s}");
    }
    let mid = 0.5 * (lo + hi);
    for (rate, reference, name) in [
        (&rates.twist, twist_at(frame, mid), "twist"),
        (&rates.rise, rise_at(frame, mid), "rise"),
    ] {
        assert!(
            (rate.middle - reference).abs() <= rate.error + 1e-7,
            "{name} {} against {reference}, error {}",
            rate.middle,
            rate.error
        );
    }
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
