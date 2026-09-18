//! The reference evaluator against the curve-evaluation contract (#106).

use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{
    Circle3, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Ellipse3, Intrinsic2,
    Intrinsic3, Line3, Polyline3,
};
use axiolid_curve_evaluate_contract::{
    conformance, CurveEvaluator, CurveMeasure, DistanceConvention,
};
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
        let p = e
            .point_at(&line, CurveMeasure::Distance(distance))
            .expect("line point");
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
    let p = e
        .point_at(&circle, CurveMeasure::Distance(half))
        .expect("circle point");
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
        let f = e
            .frame_at(&curve, CurveMeasure::Distance(d))
            .expect("frame");
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
    assert!(e.point_at(&ellipse, CurveMeasure::Distance(1.0)).is_err());
    assert!(e.tangent_at(&ellipse, CurveMeasure::Distance(1.0)).is_err());
    assert!(e.frame_at(&ellipse, CurveMeasure::Distance(1.0)).is_err());
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
    assert!(
        e.point_at(&vertical, CurveMeasure::Distance(3.0)).is_ok(),
        "point is still defined"
    );
    assert!(
        e.tangent_at(&vertical, CurveMeasure::Distance(3.0)).is_ok(),
        "tangent is defined"
    );
    assert!(
        e.frame_at(&vertical, CurveMeasure::Distance(3.0)).is_err(),
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
    let mut previous = e
        .point_at(&helix, CurveMeasure::Distance(0.0))
        .expect("start");
    for step in 1..=steps {
        let d = 40.0 * Scalar::from(step) / Scalar::from(steps);
        let current = e
            .point_at(&helix, CurveMeasure::Distance(d))
            .expect("point");
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
        let f = e
            .frame_at(&curve, CurveMeasure::Distance(d))
            .expect("frame");
        let p = e
            .point_at(&curve, CurveMeasure::Distance(d))
            .expect("point");
        let t = e
            .tangent_at(&curve, CurveMeasure::Distance(d))
            .expect("tangent");
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
        let err = e
            .point_at(&line, CurveMeasure::Distance(bad))
            .expect_err("must refuse");
        let text = format!("{err:?}");
        assert!(
            text.contains("distance"),
            "refusal for {bad} must name the distance, got: {text}"
        );
    }
}

/// A refused DISTANCE must not also close the parameter route.
///
/// This is the case that matters for a consumer holding an authored
/// `IfcParameterValue` on an ellipse: arc length needs elliptic integrals
/// and is refused, but the native parameter is perfectly meaningful and
/// must stay reachable through the contract -- otherwise the consumer is
/// pushed back onto the engine crate their architecture gate forbids.
#[test]
fn a_parameter_works_where_a_distance_is_refused() {
    let e = ReferenceCurveEvaluator::new();
    let ellipse = Curve3::Ellipse(Ellipse3 {
        frame: frame(),
        semi_axis_x: 5.0,
        semi_axis_y: 2.0,
    });
    assert!(
        e.point_at(&ellipse, CurveMeasure::Distance(1.0)).is_err(),
        "arc length on an ellipse has no closed form"
    );
    // Quarter turn: the parameter IS the angle, so this is the +y apex.
    let quarter = core::f64::consts::FRAC_PI_2;
    let p = e
        .point_at(&ellipse, CurveMeasure::Parameter(quarter))
        .expect("the parameter route must stay open");
    let want = Point3::new(0.0, 2.0, 0.0);
    assert!((p - want).length() < 1e-12, "got {p:?}, want {want:?}");
    // And the frame follows, which is what a placement actually needs.
    assert!(e
        .frame_at(&ellipse, CurveMeasure::Parameter(quarter))
        .is_ok());
}

/// The reporter's exact hazard: the same number, two meanings.
///
/// On a circle of radius 4, `1.5` as a parameter is 1.5 rad round; as a
/// distance it is 1.5 m along, i.e. 0.375 rad. Roughly 86 degrees apart.
/// Both are finite, plausible points -- nothing downstream could detect
/// the confusion, which is why the method of measurement is carried in
/// the value rather than left to the caller to remember.
#[test]
fn the_same_number_means_different_places() {
    let e = ReferenceCurveEvaluator::new();
    let radius = 4.0;
    let circle = Curve3::Circle(Circle3 {
        frame: frame(),
        radius,
    });
    let value = 1.5;
    let as_parameter = e
        .point_at(&circle, CurveMeasure::Parameter(value))
        .expect("parameter");
    let as_distance = e
        .point_at(&circle, CurveMeasure::Distance(value))
        .expect("distance");
    // The distance lands at angle d/r; the parameter lands at angle d.
    let separation = (as_parameter - as_distance).length();
    assert!(
        separation > 1.0,
        "parameter and distance must not collapse: {separation} apart"
    );
    // Pin each against its own closed form, so the test fails loudly if
    // either route silently changes meaning.
    let want_parameter = Point3::new(radius * value.cos(), radius * value.sin(), 0.0);
    let angle = value / radius;
    let want_distance = Point3::new(radius * angle.cos(), radius * angle.sin(), 0.0);
    assert!((as_parameter - want_parameter).length() < 1e-12);
    assert!((as_distance - want_distance).length() < 1e-12);
}

/// An intrinsic curve is arc-length parameterised, so the two routes
/// legitimately AGREE. Pinned so the distinction is not over-enforced.
#[test]
fn the_routes_agree_where_the_parameter_is_arc_length() {
    let e = ReferenceCurveEvaluator::new();
    let helix = Curve3::Intrinsic(Intrinsic3::new(
        frame(),
        CurvatureLaw::circular(0.12),
        CurvatureLaw::circular(0.05),
        40.0,
    ));
    for s in [0.0, 7.5, 21.0] {
        let by_p = e.point_at(&helix, CurveMeasure::Parameter(s)).expect("p");
        let by_d = e.point_at(&helix, CurveMeasure::Distance(s)).expect("d");
        assert!((by_p - by_d).length() < 1e-12, "disagreed at {s}");
    }
}

/// A non-finite PARAMETER is refused as a curve measure.
///
/// The parameter route skips the distance conversion entirely, so it
/// needs its own guard: without one a NaN reaches the evaluators and is
/// refused as a `curve parameter`, naming an internal concept instead of
/// the value the caller actually handed over.
#[test]
fn a_non_finite_parameter_is_refused_as_a_measure() {
    let e = ReferenceCurveEvaluator::new();
    let ellipse = Curve3::Ellipse(Ellipse3 {
        frame: frame(),
        semi_axis_x: 5.0,
        semi_axis_y: 2.0,
    });
    for bad in [Scalar::NAN, Scalar::INFINITY, Scalar::NEG_INFINITY] {
        let err = e
            .point_at(&ellipse, CurveMeasure::Parameter(bad))
            .expect_err("must refuse");
        let text = format!("{err:?}");
        assert!(
            text.contains("measure"),
            "refusal for {bad} must name the measure, got: {text}"
        );
        assert!(e
            .tangent_at(&ellipse, CurveMeasure::Parameter(bad))
            .is_err());
        assert!(e.frame_at(&ellipse, CurveMeasure::Parameter(bad)).is_err());
    }
}

/// A polyline's arc length is an exact finite sum, so distance works.
///
/// The reporter's example (kernel#107): segments of 5 and 12, so a
/// distance of 5 lands exactly on the interior vertex. No integral
/// and no iteration are involved.
#[test]
fn a_polyline_is_evaluable_by_distance() {
    let e = ReferenceCurveEvaluator::new();
    let curve = Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 0.0),
            Point3::new(3.0, 4.0, 12.0),
        ],
        closed: false,
    });
    assert_eq!(
        e.distance_convention(&curve),
        DistanceConvention::ArcLength3d
    );
    let at_vertex = e
        .point_at(&curve, CurveMeasure::Distance(5.0))
        .expect("5 m along a 5 + 12 polyline is exactly the vertex");
    assert_eq!(at_vertex, Point3::new(3.0, 4.0, 0.0));
}

/// A distance landing on a seam reads the OUTGOING heading.
///
/// At an interior vertex the tangent is two-valued. Leaving that to
/// chance would silently pick one of two headings, so the convention
/// is pinned: the curve is about to travel +z, not +x/+y.
#[test]
fn a_seam_reports_the_outgoing_tangent() {
    let e = ReferenceCurveEvaluator::new();
    let curve = Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 0.0),
            Point3::new(3.0, 4.0, 12.0),
        ],
        closed: false,
    });
    let tangent = e
        .tangent_at(&curve, CurveMeasure::Distance(5.0))
        .expect("the seam resolves to the outgoing segment");
    assert_eq!(tangent, Vec3::new(0.0, 0.0, 1.0));
}

/// A repeated point has no direction, so it is refused, not skipped.
///
/// Skipping a zero-length segment would silently change the
/// parameterisation; normalising its zero tangent would invent a
/// heading. The whole curve is ill-defined, so even a distance that
/// stops short of the repeat is refused.
#[test]
fn a_zero_length_segment_is_refused() {
    let e = ReferenceCurveEvaluator::new();
    let curve = Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        closed: false,
    });
    assert!(e.point_at(&curve, CurveMeasure::Distance(0.5)).is_err());
    assert!(e.point_at(&curve, CurveMeasure::Distance(1.5)).is_err());
}

/// A closed polyline's wrap segment is length, but does not wrap.
///
/// The closing segment counts toward the total, so a unit square is
/// 4 long rather than 3. Distance past the end is refused rather
/// than wrapped around, so a caller cannot silently lap the curve.
#[test]
fn a_closed_polyline_counts_the_wrap_but_does_not_lap() {
    let e = ReferenceCurveEvaluator::new();
    let square = Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        closed: true,
    });
    // 3.5 is on the closing segment, half way back to the start.
    let on_wrap = e
        .point_at(&square, CurveMeasure::Distance(3.5))
        .expect("the wrap segment is evaluable");
    assert_eq!(on_wrap, Point3::new(0.0, 0.5, 0.0));
    assert!(e.point_at(&square, CurveMeasure::Distance(4.5)).is_err());
}
