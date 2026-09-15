//! Arc-length evaluation of spirals and of the plan-plus-elevation
//! composition (ADR 0060, issue #105).
//!
//! The clothoid's position is checked against the FRESNEL integral computed
//! independently in the test, not against another run of the kernel: a
//! quadrature checked against itself proves nothing.

use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2};
use axiolid_curve::{CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2, Line2};
use axiolid_evaluate::{elevated_point, elevated_tangent, intrinsic_point, intrinsic_tangent};

fn origin_frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

/// Fresnel integrals by series, independent of the kernel's quadrature.
///
/// C(x) = int_0^x cos(pi/2 t^2) dt, S(x) = int_0^x sin(pi/2 t^2) dt.
/// The series converges fast for the moderate arguments used here.
fn fresnel(x: Scalar) -> (Scalar, Scalar) {
    let mut c = 0.0;
    let mut s = 0.0;
    for n in 0..60 {
        let n_f = n as Scalar;
        // C term: (-1)^n (pi/2)^{2n} x^{4n+1} / ((2n)! (4n+1))
        let mut term_c = 1.0;
        for k in 1..=(2 * n) {
            term_c *= x * x * core::f64::consts::FRAC_PI_2 / k as Scalar;
        }
        term_c *= x / (4.0 * n_f + 1.0);
        c += if n % 2 == 0 { term_c } else { -term_c };

        // S term: (-1)^n (pi/2)^{2n+1} x^{4n+3} / ((2n+1)! (4n+3))
        let mut term_s = 1.0;
        for k in 1..=(2 * n + 1) {
            term_s *= x * x * core::f64::consts::FRAC_PI_2 / k as Scalar;
        }
        term_s *= x / (4.0 * n_f + 3.0);
        s += if n % 2 == 0 { term_s } else { -term_s };
    }
    (c, s)
}

#[test]
fn a_clothoid_position_matches_the_fresnel_closed_form() {
    // Transition from straight into a 300 m radius over 120 m: k(s) = s / (R L).
    let (radius, length) = (300.0, 120.0);
    let rate = 1.0 / (radius * length);
    let curve = Intrinsic2::new(
        origin_frame(),
        CurvatureLaw::clothoid(0.0, 1.0 / radius, length),
        length,
    );

    // Fresnel reference: with k = rate * s, heading = rate * s^2 / 2, so
    // x = a C(L/a), y = a S(L/a) for a = sqrt(pi / rate).
    let a = (core::f64::consts::PI / rate).sqrt();
    let (c, s) = fresnel(length / a);
    let (expected_x, expected_y) = (a * c, a * s);

    let got = intrinsic_point(&curve, length).expect("a clothoid evaluates");
    assert!(
        (got.x - expected_x).abs() < 1e-9 && (got.y - expected_y).abs() < 1e-9,
        "expected ({expected_x}, {expected_y}), got ({}, {})",
        got.x,
        got.y
    );
}

#[test]
fn a_circular_law_reproduces_the_exact_arc() {
    // A constant curvature law is a circle, which DOES have a closed form:
    // it pins the quadrature against elementary geometry.
    let radius = 50.0;
    let quarter = core::f64::consts::FRAC_PI_2 * radius;
    let curve = Intrinsic2::new(
        origin_frame(),
        CurvatureLaw::circular(1.0 / radius),
        quarter,
    );

    let got = intrinsic_point(&curve, quarter).expect("a circular law evaluates");
    // Starting at the origin heading +x, curving left: quarter turn ends at
    // (r, r).
    assert!(
        (got.x - radius).abs() < 1e-9 && (got.y - radius).abs() < 1e-9,
        "expected ({radius}, {radius}), got ({}, {})",
        got.x,
        got.y
    );

    let tangent = intrinsic_tangent(&curve, quarter).expect("tangent");
    assert!(
        tangent.x.abs() < 1e-12 && (tangent.y - 1.0).abs() < 1e-12,
        "a quarter turn must end heading +y, got ({}, {})",
        tangent.x,
        tangent.y
    );
}

#[test]
fn a_straight_law_is_exactly_the_line() {
    let curve = Intrinsic2::new(origin_frame(), CurvatureLaw::straight(), 100.0);
    let got = intrinsic_point(&curve, 42.0).expect("a straight law evaluates");
    assert!(
        (got.x - 42.0).abs() < 1e-12 && got.y.abs() < 1e-12,
        "expected (42, 0), got ({}, {})",
        got.x,
        got.y
    );
}

#[test]
fn an_elevated_spiral_keeps_both_halves_exact() {
    // The composition the issue asks for: an exact clothoid plan paired with
    // an exact parabolic vertical profile, neither approximated.
    let (radius, length) = (300.0, 120.0);
    let plan = Curve2::Intrinsic(Intrinsic2::new(
        origin_frame(),
        CurvatureLaw::clothoid(0.0, 1.0 / radius, length),
        length,
    ));
    // z(d) = 100 + 0.02 d - 0.0002 d^2 : entry grade 2%, falling.
    let elevation = ElevationLaw::parabolic(100.0, 0.02, -0.02, length);
    let curve = Elevated3::new(plan.clone(), elevation);

    // The plan half must equal the bare plan curve, unchanged.
    let at = 75.0;
    let spatial = elevated_point(&curve, at).expect("elevated evaluates");
    let planar = intrinsic_point(
        match &plan {
            Curve2::Intrinsic(i) => i,
            _ => unreachable!(),
        },
        at,
    )
    .expect("plan evaluates");
    assert!(
        (spatial.x - planar.x).abs() < 1e-12 && (spatial.y - planar.y).abs() < 1e-12,
        "pairing must not disturb the plan"
    );

    // The vertical half must equal the law, unchanged.
    let expected_z = 100.0 + 0.02 * at - (0.04 / (2.0 * length)) * at * at;
    assert!(
        (spatial.z - expected_z).abs() < 1e-12,
        "expected z {expected_z}, got {}",
        spatial.z
    );
}

#[test]
fn the_elevated_tangent_carries_the_grade() {
    // Flat straight plan, constant 3% grade: the tangent must be
    // (1, 0, 0.03) normalised. This is the sqrt(1 + g^2) factor by which 3D
    // arc length runs ahead of plan distance.
    let plan = Curve2::Line(Line2 {
        origin: Point2::new(0.0, 0.0),
        direction: Vec2::new(1.0, 0.0),
    });
    let curve = Elevated3::new(plan, ElevationLaw::constant_grade(10.0, 0.03));

    let tangent = elevated_tangent(&curve, 500.0).expect("tangent");
    let scale = (1.0 + 0.03 * 0.03_f64).sqrt();
    assert!(
        (tangent.x - 1.0 / scale).abs() < 1e-12
            && tangent.y.abs() < 1e-12
            && (tangent.z - 0.03 / scale).abs() < 1e-12,
        "got ({}, {}, {})",
        tangent.x,
        tangent.y,
        tangent.z
    );
    let norm = (tangent.x * tangent.x + tangent.y * tangent.y + tangent.z * tangent.z).sqrt();
    assert!(
        (norm - 1.0).abs() < 1e-12,
        "tangent must be unit, got {norm}"
    );
}

#[test]
fn a_piecewise_profile_restarts_each_piece_at_its_own_zero() {
    // Two 100 m pieces: level at 50, then climbing at 5% from 50.
    // At d = 150 the second piece is 50 m in, so z = 50 + 0.05*50 = 52.5.
    let elevation = ElevationLaw::Piecewise {
        breaks: vec![100.0],
        laws: vec![
            ElevationLaw::level(50.0),
            ElevationLaw::constant_grade(50.0, 0.05),
        ],
    };
    assert_eq!(elevation.height_at(50.0), Some(50.0));
    assert_eq!(elevation.height_at(150.0), Some(52.5));
    // The seam belongs to the piece that starts there.
    assert_eq!(elevation.height_at(100.0), Some(50.0));
    assert!(elevation.is_well_formed());
}

#[test]
fn a_mismatched_piecewise_profile_refuses_rather_than_guessing() {
    let broken = ElevationLaw::Piecewise {
        breaks: vec![10.0, 20.0],
        laws: vec![ElevationLaw::level(0.0)],
    };
    assert!(!broken.is_well_formed());
    assert_eq!(broken.height_at(5.0), None);
}

#[test]
fn a_bspline_plan_is_refused_because_its_parameter_is_not_arc_length() {
    // An elevation law is written against DISTANCE. A B-spline's parameter is
    // not arc length, so pairing one would silently mean something else.
    let plan = Curve2::BSpline(axiolid_curve::BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: axiolid_curve::KnotSpec::Unspecified,
    });
    let curve = Elevated3::new(plan, ElevationLaw::level(0.0));
    assert!(
        elevated_point(&curve, 1.0).is_err(),
        "a non-arc-length parameterisation must refuse, not be reinterpreted"
    );
}

#[test]
fn an_elevated_curve_is_a_curve3_value() {
    // The composition must be storable as an ordinary Curve3, which is what
    // lets it flow through the graph like any other curve.
    let plan = Curve2::Line(Line2 {
        origin: Point2::new(0.0, 0.0),
        direction: Vec2::new(1.0, 0.0),
    });
    let curve = Curve3::Elevated(Elevated3::new(plan, ElevationLaw::level(7.0)));
    let Curve3::Elevated(elevated) = &curve else {
        panic!("must round trip as the same variant");
    };
    let point = elevated_point(elevated, 3.0).expect("evaluates");
    assert!((point.z - 7.0).abs() < 1e-12);
    let _ = Tolerance::METRE;
}
