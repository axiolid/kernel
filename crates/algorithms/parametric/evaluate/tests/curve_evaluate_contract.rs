//! The reference evaluator against the curve-evaluation contract (#106).

use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{
    Circle3, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Ellipse3, Intrinsic2,
    Intrinsic3, Line3,
};
use axiolid_curve_evaluate_contract::{conformance, CurveEvaluator, DistanceConvention};
use axiolid_evaluate::ReferenceCurveEvaluator;

fn frame() -> Frame3 {
    Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

/// The provider must satisfy every published conformance expectation.
#[test]
fn the_reference_evaluator_is_conformant() {
    let failures = conformance::check(&ReferenceCurveEvaluator::new());
    assert!(failures.is_empty(), "conformance failures: {failures:#?}");
}

/// A non-unit line direction must not leak into the distance.
///
/// Import adapters preserve a non-unit direction, so a provider that
/// passed the distance through as a parameter would place an object at
/// twice the requested station here while looking correct on a unit line.
#[test]
fn a_non_unit_line_direction_does_not_scale_the_distance() {
    let e = ReferenceCurveEvaluator::new();
    let origin = Point3::new(1.0, 2.0, 3.0);
    let line = Curve3::Line(Line3 {
        origin,
        direction: Vec3::new(0.0, 3.0, 4.0),
    });
    for distance in [0.0, 1.0, 7.5, 100.0] {
        let p = e.point_at(&line, distance).expect("line point");
        let moved = (p - origin).length();
        assert!(
            (moved - distance).abs() < 1e-9,
            "distance {distance} moved {moved}"
        );
    }
}

/// A circle must advance arc length, not angle.
#[test]
fn a_circle_advances_arc_length_not_angle() {
    let e = ReferenceCurveEvaluator::new();
    let radius = 7.0;
    let circle = Curve3::Circle(Circle3 {
        frame: frame(),
        radius,
    });
    // Half the circumference must land diametrically opposite.
    let half = core::f64::consts::PI * radius;
    let p = e.point_at(&circle, half).expect("circle point");
    let want = Point3::new(-radius, 0.0, 0.0);
    assert!((p - want).length() < 1e-9, "half circumference gave {p:?}");
}

/// Straight plan, crest then sag: an alignment with an inflection.
fn crest_then_sag() -> Curve3 {
    // Plan is a straight line along +x, 200 m long.
    let plan = Curve2::Intrinsic(Intrinsic2::new(
        axiolid_core::Frame2 {
            origin: axiolid_core::Point2::new(0.0, 0.0),
            x: axiolid_core::Vec2::X,
            y: axiolid_core::Vec2::Y,
        },
        CurvatureLaw::straight(),
        200.0,
    ));
    // Crest: +5% into -3% over 100 m. Sag: -3% into +4% over 100 m.
    // Each piece is written in its own distance, restarting at its seam.
    let elevation = ElevationLaw::Piecewise {
        breaks: vec![100.0],
        laws: vec![
            ElevationLaw::parabolic(0.0, 0.05, -0.03, 100.0),
            ElevationLaw::parabolic(1.0, -0.03, 0.04, 100.0),
        ],
    };
    Curve3::Elevated(Elevated3 {
        plan: Box::new(plan),
        elevation,
    })
}

/// The placement frame must stay upright through an inflection.
///
/// This is the test that rules out the Frenet frame. On this alignment
/// the Frenet normal points DOWN on the crest and UP in the sag, so an
/// object placed by it would be upside down on one half. The
/// reference-up frame keeps a positive up component throughout.
#[test]
fn the_placement_frame_stays_upright_across_an_inflection() {
    let e = ReferenceCurveEvaluator::new();
    let curve = crest_then_sag();
    for step in 0..=20 {
        let d = 200.0 * Scalar::from(step) / 20.0;
        let f = e.frame_at(&curve, d).expect("frame");
        assert!(
            f.y.z > 0.9,
            "frame tipped at distance {d}: up component {}",
            f.y.z
        );
    }
}

/// Plan distance is reported as such, never as 3D arc length.
///
/// On a graded alignment the two differ by the grade factor. Reporting
/// PlanDistance is what lets a caller know it must not treat a station
/// as a travelled distance.
#[test]
fn an_elevated_curve_reports_plan_distance_not_arc_length() {
    let e = ReferenceCurveEvaluator::new();
    assert_eq!(
        e.distance_convention(&crest_then_sag()),
        DistanceConvention::PlanDistance
    );
}

/// Families with no closed-form arc length are refused, not approximated.
///
/// An ellipse's arc length is an elliptic integral. Returning the angle
/// as though it were a distance would be a lie the caller cannot see.
#[test]
fn a_family_without_closed_form_arc_length_is_refused() {
    let e = ReferenceCurveEvaluator::new();
    let ellipse = Curve3::Ellipse(Ellipse3 {
        frame: frame(),
        semi_axis_x: 5.0,
        semi_axis_y: 2.0,
    });
    assert_eq!(
        e.distance_convention(&ellipse),
        DistanceConvention::Unsupported
    );
    assert!(e.point_at(&ellipse, 1.0).is_err());
    assert!(e.tangent_at(&ellipse, 1.0).is_err());
    assert!(e.frame_at(&ellipse, 1.0).is_err());
}

/// A vertical tangent leaves roll undetermined, so the frame is refused.
///
/// Returning an arbitrary roll would rotate whatever is placed by it,
/// silently. The point and tangent are still available.
#[test]
fn a_tangent_parallel_to_up_refuses_a_frame() {
    let e = ReferenceCurveEvaluator::new();
    let vertical = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::Z,
    });
    assert!(e.point_at(&vertical, 3.0).is_ok(), "point is still defined");
    assert!(e.tangent_at(&vertical, 3.0).is_ok(), "tangent is defined");
    assert!(
        e.frame_at(&vertical, 3.0).is_err(),
        "roll is undetermined, so the frame must be refused"
    );
}

/// A helix is a torsion curve; distance along it is its native parameter.
#[test]
fn a_torsion_curve_uses_arc_length_directly() {
    let e = ReferenceCurveEvaluator::new();
    let helix = Curve3::Intrinsic(Intrinsic3::new(
        frame(),
        CurvatureLaw::circular(0.12),
        CurvatureLaw::circular(0.05),
        40.0,
    ));
    assert_eq!(
        e.distance_convention(&helix),
        DistanceConvention::ArcLength3d
    );
    // Arc length between two stations must equal their separation,
    // measured along the curve by dense sampling.
    let mut travelled = 0.0;
    let steps = 4000;
    let mut previous = e.point_at(&helix, 0.0).expect("start");
    for step in 1..=steps {
        let d = 40.0 * Scalar::from(step) / Scalar::from(steps);
        let current = e.point_at(&helix, d).expect("point");
        travelled += (current - previous).length();
        previous = current;
    }
    assert!(
        (travelled - 40.0).abs() < 1e-3,
        "40 m of arc length measured {travelled} m"
    );
}

/// The frame must sit on the curve with its x axis on the tangent.
#[test]
fn the_frame_agrees_with_the_point_and_tangent() {
    let e = ReferenceCurveEvaluator::new();
    let curve = crest_then_sag();
    for step in 0..=10 {
        let d = 200.0 * Scalar::from(step) / 10.0;
        let f = e.frame_at(&curve, d).expect("frame");
        let p = e.point_at(&curve, d).expect("point");
        let t = e.tangent_at(&curve, d).expect("tangent");
        assert!((f.origin - p).length() < 1e-12, "origin off curve at {d}");
        assert!((f.x - t).length() < 1e-12, "x is not the tangent at {d}");
        // Orthonormal and right-handed.
        assert!((f.x.length() - 1.0).abs() < 1e-12);
        assert!(f.x.dot(f.y).abs() < 1e-12);
        assert!((f.x.cross(f.y).dot(f.z) - 1.0).abs() < 1e-12);
    }
}

/// A non-finite distance is refused by name, not passed downstream.
///
/// The parameter conversion divides by a speed or radius, so a NaN
/// would survive as a NaN parameter and be caught further down with a
/// message about a `curve parameter`. A caller asked about a DISTANCE,
/// so the refusal must name the distance -- otherwise the error points
/// at an internal concept the caller never supplied.
#[test]
fn a_non_finite_distance_is_refused_as_a_distance() {
    let e = ReferenceCurveEvaluator::new();
    let line = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::new(2.0, 0.0, 0.0),
    });
    for bad in [Scalar::NAN, Scalar::INFINITY, Scalar::NEG_INFINITY] {
        let err = e.point_at(&line, bad).expect_err("must refuse");
        let text = format!("{err:?}");
        assert!(
            text.contains("distance"),
            "refusal for {bad} must name the distance, got: {text}"
        );
    }
}
