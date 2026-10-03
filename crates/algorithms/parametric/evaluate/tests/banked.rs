//! Banked curves: points, tangents and section frames at chosen stations
//! against hand-computed values, under both bank conventions (ADR 0081).

use axiolid_core::{Frame2, Interval, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Circle2, Curve2, Curve3, Elevated3, ElevationLaw,
    Line2,
};
use axiolid_curve_evaluate_contract::{CurveEvaluator, CurveMeasure, DistanceConvention};
use axiolid_evaluate::banked::{banked_point, banked_section, banked_tangent, BankedSection};
use axiolid_evaluate::bound::certifies_flattening3;
use axiolid_evaluate::curve::{derivative3, domain3, evaluate3, flatten3};
use axiolid_evaluate::ReferenceCurveEvaluator;

const EPS: Scalar = 1e-12;

fn close(actual: Scalar, expected: Scalar, what: &str) {
    assert!(
        (actual - expected).abs() <= EPS,
        "{what}: {actual} != {expected} (off by {:e})",
        actual - expected
    );
}

fn close3(actual: Vec3, expected: Vec3, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= EPS,
        "{what}: {actual:?} != {expected:?}"
    );
}

/// Along `+x` from the origin, climbing 2% from 100 m.
fn straight_two_percent() -> Elevated3 {
    Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::new(0.0, 0.0),
            direction: Vec2::new(1.0, 0.0),
        }),
        ElevationLaw::constant_grade(100.0, 0.02),
    )
}

/// 0 -> 150 mm linear over 100 m, then 150 mm for 50 m.
fn ramp_then_hold() -> CantLaw {
    CantLaw::new(vec![
        CantPiece::linear(100.0, 0.0, 0.15),
        CantPiece::constant(50.0, 0.15),
    ])
}

fn banked(base: Elevated3, cant: CantLaw, pivot: CantLaw, convention: BankConvention) -> Banked3 {
    Banked3::new(base, cant, pivot, 1.5, convention)
}

fn assert_orthonormal(section: &BankedSection, what: &str) {
    let (t, l, v) = (section.tangent, section.lateral, section.up);
    for (axis, name) in [(t, "tangent"), (l, "lateral"), (v, "up")] {
        close(axis.length(), 1.0, &format!("{what}: |{name}|"));
    }
    close(t.dot(l), 0.0, &format!("{what}: t.l"));
    close(t.dot(v), 0.0, &format!("{what}: t.v"));
    close(l.dot(v), 0.0, &format!("{what}: l.v"));
    close3(t.cross(l), v, &format!("{what}: t x l = v"));
    let frame = section.frame();
    close3(
        frame.x.cross(frame.y),
        frame.z,
        &format!("{what}: right-handed"),
    );
}

#[test]
fn a_tangent_rotation_frame_at_a_two_percent_grade() {
    let curve = banked(
        straight_two_percent(),
        ramp_then_hold(),
        CantLaw::zero(150.0),
        BankConvention::TangentRotation,
    );
    let section = banked_section(&curve, 100.0).unwrap();
    let s = 1.0004_f64.sqrt();
    let cos_psi = 0.99_f64.sqrt();
    // t = (1, 0, g) / s, n = +y, u = t x n = (-g, 0, 1) / s, psi = asin(0.1).
    close3(section.point, Point3::new(100.0, 0.0, 102.0), "point");
    close3(
        section.tangent,
        Vec3::new(1.0 / s, 0.0, 0.02 / s),
        "tangent",
    );
    close3(
        section.lateral,
        Vec3::new(-0.1 * 0.02 / s, cos_psi, 0.1 / s),
        "lateral = cos psi n + sin psi u",
    );
    close3(
        section.up,
        Vec3::new(-cos_psi * 0.02 / s, -0.1, cos_psi / s),
        "up = -sin psi n + cos psi u",
    );
    close(section.cant, 0.15, "cant");
    close(section.bank_angle, 0.100_167_421_161_559_8, "psi");
    close(section.roll, section.bank_angle, "rotation rolls by psi");
    assert_orthonormal(&section, "tangent rotation");
    // The rail heads rise D cos(theta) across, not D.
    let (left, right) = section.rail_heads();
    close(left.z - right.z, 0.15 / s, "rise is D cos theta");
}

#[test]
fn a_vertical_rise_frame_at_a_two_percent_grade_and_the_difference() {
    let rotation = banked(
        straight_two_percent(),
        ramp_then_hold(),
        CantLaw::zero(150.0),
        BankConvention::TangentRotation,
    );
    let rise = Banked3 {
        convention: BankConvention::VerticalRise,
        ..rotation.clone()
    };
    let s = 1.0004_f64.sqrt();
    let section = banked_section(&rise, 100.0).unwrap();
    // rho = asin(D / (b cos theta)) = asin(0.1 s).
    close(section.roll, (0.1 * s).asin(), "vertical rise roll");
    close(
        section.bank_angle,
        0.1_f64.asin(),
        "nominal psi is unchanged",
    );
    let sin_rho = 0.1 * s;
    let cos_rho = (1.0 - sin_rho * sin_rho).sqrt();
    close3(
        section.lateral,
        Vec3::new(-sin_rho * 0.02 / s, cos_rho, sin_rho / s),
        "lateral",
    );
    assert_orthonormal(&section, "vertical rise");
    let (left, right) = section.rail_heads();
    close(left.z - right.z, 0.15, "rise is exactly D");
    // The two conventions differ by D (1 - cos theta) = 1.9994e-4 D.
    let (rotated_left, rotated_right) = banked_section(&rotation, 100.0).unwrap().rail_heads();
    let difference = (left.z - right.z) - (rotated_left.z - rotated_right.z);
    close(
        difference,
        0.15 * (1.0 - 1.0 / s),
        "difference D (1 - cos theta)",
    );
    close(
        difference / 0.15,
        1.999_400_199_929_057_8e-4,
        "about 2e-4 D at 2%",
    );
    // On level track the two agree.
    let level = |convention| {
        let mut curve = banked(
            straight_two_percent(),
            ramp_then_hold(),
            CantLaw::zero(150.0),
            convention,
        );
        curve.base.elevation = ElevationLaw::level(5.0);
        banked_section(&curve, 120.0).unwrap()
    };
    let (a, b) = (
        level(BankConvention::TangentRotation),
        level(BankConvention::VerticalRise),
    );
    close3(a.lateral, b.lateral, "level: same lateral");
    close3(a.up, b.up, "level: same up");
}

#[test]
fn a_centreline_pivot_splits_the_cant_between_the_rails() {
    let curve = banked(
        straight_two_percent(),
        ramp_then_hold(),
        CantLaw::zero(150.0),
        BankConvention::VerticalRise,
    );
    // d = 30: D = 0.045, centreline at 100.6.
    let section = banked_section(&curve, 30.0).unwrap();
    close3(
        section.point,
        Point3::new(30.0, 0.0, 100.6),
        "rotation point",
    );
    close(section.pivot, 0.0, "pivot");
    let (left, right) = section.rail_heads();
    close(left.z, 100.6 + 0.0225, "left head up D/2");
    close(right.z, 100.6 - 0.0225, "right head down D/2");
}

#[test]
fn a_low_rail_pivot_keeps_the_low_rail_on_the_profile() {
    // e = D / 2 throughout: the left rail is high, the right rail the pivot.
    let curve = banked(
        straight_two_percent(),
        ramp_then_hold(),
        CantLaw::new(vec![
            CantPiece::linear(100.0, 0.0, 0.075),
            CantPiece::constant(50.0, 0.075),
        ]),
        BankConvention::VerticalRise,
    );
    // In the transition: e = 0.0375 rising at 7.5e-4 per metre, so the
    // curve's own grade is 0.02075 and its tangent tilts with it.
    let section = banked_section(&curve, 50.0).unwrap();
    close3(
        section.point,
        Point3::new(50.0, 0.0, 101.0375),
        "raised point",
    );
    close(section.grade, 0.020_75, "grade plus pivot rate");
    let g: Scalar = 0.020_75;
    close3(
        section.tangent,
        Vec3::new(1.0, 0.0, g) / (1.0 + g * g).sqrt(),
        "tangent of the raised curve",
    );
    let (left, right) = section.rail_heads();
    close(right.z, 101.0, "low rail on the profile");
    close(left.z, 101.075, "high rail D above it");
    assert_orthonormal(&section, "low-rail pivot");
    // In the hold: e constant, tangent back to the profile's.
    let held = banked_section(&curve, 120.0).unwrap();
    close(held.grade, 0.02, "constant pivot leaves the grade");
    close(
        held.rail_heads().1.z,
        102.4,
        "low rail on the profile in the hold",
    );
}

#[test]
fn a_viennese_bend_rolls_by_its_own_angle() {
    let mut base = straight_two_percent();
    base.elevation = ElevationLaw::level(0.0);
    let cant = CantLaw::new(vec![CantPiece::viennese_bend(100.0, 0.0, 0.1)]);
    for convention in [
        BankConvention::TangentRotation,
        BankConvention::VerticalRise,
    ] {
        let curve = banked(base.clone(), cant.clone(), CantLaw::zero(100.0), convention);
        let section = banked_section(&curve, 50.0).unwrap();
        close(section.bank_angle, 0.05, "psi at xi = 1/2");
        close(section.roll, 0.05, "level: roll is psi");
        close(section.cant, 1.5 * 0.05_f64.sin(), "D = b sin psi");
        close3(
            section.lateral,
            Vec3::new(0.0, 0.05_f64.cos(), 0.05_f64.sin()),
            "lateral",
        );
        let quarter = banked_section(&curve, 25.0).unwrap();
        close(quarter.bank_angle, 0.007_055_664_062_5, "psi at xi = 1/4");
    }
    // On a grade, a vertical rise rolls further than psi.
    let curve = banked(
        straight_two_percent(),
        cant,
        CantLaw::zero(100.0),
        BankConvention::VerticalRise,
    );
    let section = banked_section(&curve, 50.0).unwrap();
    close(
        section.roll,
        (0.05_f64.sin() * 1.0004_f64.sqrt()).asin(),
        "rho = asin(sin(psi) / cos(theta))",
    );
    let (left, right) = section.rail_heads();
    close(
        left.z - right.z,
        1.5 * 0.05_f64.sin(),
        "rise b sin psi exactly",
    );
}

#[test]
fn a_negative_cant_on_a_left_curve_raises_the_outer_rail() {
    let plan = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 300.0,
    });
    let curve = banked(
        Elevated3::new(plan, ElevationLaw::level(10.0)),
        CantLaw::new(vec![CantPiece::constant(500.0, -0.15)]),
        CantLaw::zero(500.0),
        BankConvention::TangentRotation,
    );
    let quarter = 300.0 * core::f64::consts::FRAC_PI_2;
    let section = banked_section(&curve, quarter).unwrap();
    close3(
        section.point,
        Point3::new(0.0, 300.0, 10.0),
        "quarter point",
    );
    close3(section.tangent, Vec3::new(-1.0, 0.0, 0.0), "tangent");
    // n points at the centre; psi = -asin(0.1).
    close3(
        section.lateral,
        Vec3::new(0.0, -(0.99_f64.sqrt()), -0.1),
        "lateral",
    );
    let (inner, outer) = section.rail_heads();
    close(inner.z, 10.0 - 0.075, "inner (left) rail low");
    close(outer.z, 10.0 + 0.075, "outer (right) rail high");
    assert!(outer.y > 300.0, "the right rail is outside the curve");
}

#[test]
fn frames_are_orthonormal_everywhere_and_continuous_at_seams() {
    let plan = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 250.0,
    });
    let cant = CantLaw::new(vec![
        CantPiece::cosine(60.0, 0.0, 0.16),
        CantPiece::constant(40.0, 0.16),
        CantPiece::sine(60.0, 0.16, 0.02),
    ]);
    let pivot = CantLaw::new(vec![
        CantPiece::cosine(60.0, 0.0, 0.08),
        CantPiece::constant(40.0, 0.08),
        CantPiece::sine(60.0, 0.08, 0.01),
    ]);
    for convention in [
        BankConvention::TangentRotation,
        BankConvention::VerticalRise,
    ] {
        let curve = banked(
            Elevated3::new(
                plan.clone(),
                ElevationLaw::parabolic(20.0, 0.04, -0.03, 160.0),
            ),
            cant.clone(),
            pivot.clone(),
            convention,
        );
        for step in 0..=32 {
            let d = 160.0 * Scalar::from(step) / 32.0;
            let section = banked_section(&curve, d).unwrap();
            assert_orthonormal(&section, &format!("{convention:?} at {d}"));
            close3(
                section.tangent,
                banked_tangent(&curve, d).unwrap(),
                "section tangent is the curve's",
            );
        }
        for seam in cant.seams() {
            let before = banked_section(&curve, seam - 1e-9).unwrap();
            let at = banked_section(&curve, seam).unwrap();
            assert!(
                (before.lateral - at.lateral).length() < 1e-9,
                "lateral jumps at {seam}"
            );
            assert!((before.up - at.up).length() < 1e-9, "up jumps at {seam}");
            // The point moves 1e-9 along the curve, at speed just over one.
            assert!(
                (before.point - at.point).length() < 2e-9,
                "point jumps at {seam}"
            );
        }
    }
}

#[test]
fn a_cant_beyond_the_rail_heads_or_their_vertical_span_is_refused_by_name() {
    let too_much = banked(
        straight_two_percent(),
        CantLaw::new(vec![CantPiece::constant(10.0, 1.6)]),
        CantLaw::zero(10.0),
        BankConvention::TangentRotation,
    );
    let message = banked_section(&too_much, 5.0).unwrap_err().to_string();
    assert!(
        message.contains("exceeds the rail-head distance"),
        "{message}"
    );
    // 1.4999 fits b = 1.5 but not b cos(theta) = 1.4925 at 10%.
    let mut steep = banked(
        straight_two_percent(),
        CantLaw::new(vec![CantPiece::constant(10.0, 1.4999)]),
        CantLaw::zero(10.0),
        BankConvention::VerticalRise,
    );
    steep.base.elevation = ElevationLaw::constant_grade(0.0, 0.1);
    let message = banked_section(&steep, 5.0).unwrap_err().to_string();
    assert!(message.contains("vertical span"), "{message}");
    steep.convention = BankConvention::TangentRotation;
    assert!(banked_section(&steep, 5.0).is_ok());
    // Off the law, and an angle piece as a pivot.
    let message = banked_point(&steep, 10.5).unwrap_err().to_string();
    assert!(message.contains("cant law"), "{message}");
    steep.pivot = CantLaw::new(vec![CantPiece::viennese_bend(10.0, 0.0, 0.1)]);
    let message = banked_point(&steep, 5.0).unwrap_err().to_string();
    assert!(
        message.contains("pivot law has an angle-form piece"),
        "{message}"
    );
}

#[test]
fn the_reference_evaluator_frames_a_banked_curve_by_its_section() {
    let curve = Curve3::Banked(banked(
        straight_two_percent(),
        ramp_then_hold(),
        CantLaw::zero(150.0),
        BankConvention::VerticalRise,
    ));
    let Curve3::Banked(inner) = &curve else {
        unreachable!()
    };
    let evaluator = ReferenceCurveEvaluator::new();
    assert_eq!(
        evaluator.distance_convention(&curve),
        DistanceConvention::PlanDistance
    );
    let section = banked_section(inner, 75.0).unwrap();
    for at in [CurveMeasure::Distance(75.0), CurveMeasure::Parameter(75.0)] {
        assert_eq!(evaluator.frame_at(&curve, at).unwrap(), section.frame());
        assert_eq!(evaluator.point_at(&curve, at).unwrap(), section.point);
        close3(
            evaluator.tangent_at(&curve, at).unwrap(),
            section.tangent,
            "tangent",
        );
    }
    // The cant is against +Z; another reference up is refused.
    let tilted = ReferenceCurveEvaluator::with_up(Vec3::Y).unwrap();
    assert!(tilted
        .frame_at(&curve, CurveMeasure::Distance(75.0))
        .is_err());
    // At zero cant the section frame is the reference-up frame.
    let flat = Curve3::Banked(banked(
        straight_two_percent(),
        CantLaw::zero(150.0),
        CantLaw::zero(150.0),
        BankConvention::TangentRotation,
    ));
    let elevated = Curve3::Elevated(straight_two_percent());
    let at = CurveMeasure::Distance(40.0);
    let (a, b) = (
        evaluator.frame_at(&flat, at).unwrap(),
        evaluator.frame_at(&elevated, at).unwrap(),
    );
    close3(a.x, b.x, "x");
    close3(a.y, b.y, "y");
    close3(a.z, b.z, "z");
}

#[test]
fn a_banked_curve_evaluates_and_flattens_by_plan_distance_uncertified() {
    let curve = Curve3::Banked(banked(
        straight_two_percent(),
        ramp_then_hold(),
        CantLaw::new(vec![
            CantPiece::linear(100.0, 0.0, 0.075),
            CantPiece::constant(50.0, 0.075),
        ]),
        BankConvention::TangentRotation,
    ));
    assert_eq!(domain3(&curve), Interval::new(0.0, 150.0));
    close3(
        evaluate3(&curve, 50.0).unwrap(),
        Point3::new(50.0, 0.0, 101.0375),
        "evaluate3",
    );
    close3(
        derivative3(&curve, 50.0).unwrap(),
        Vec3::new(1.0, 0.0, 0.020_75),
        "derivative3 in plan distance",
    );
    assert!(evaluate3(&curve, 150.5).is_err());
    let points = flatten3(&curve, domain3(&curve), 1e-3, 20).unwrap();
    assert!(points.len() >= 2);
    close3(points[0], Point3::new(0.0, 0.0, 100.0), "first");
    close3(
        *points.last().unwrap(),
        Point3::new(150.0, 0.0, 103.075),
        "last",
    );
    // Midpoint-sagitta flattening, not a certified bound (#232).
    assert!(!certifies_flattening3(&curve));
}
