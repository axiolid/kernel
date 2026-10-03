//! Stations (#241): points and section frames at a distance along lines,
//! arcs, polylines, clothoids, chains, elevated and banked curves, against
//! values computed by hand -- closed forms, and the clothoid by its Fresnel
//! power series, never by the quadrature under test.

use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, Circle2, CurvatureLaw,
    Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2, Line2, Line3, Polyline2,
};
use axiolid_curve_evaluate_contract::{CurveEvaluator, CurveMeasure};
use axiolid_evaluate::station::{
    station_length2, station_length3, station_section2, station_section3, SectionFrame,
};
use axiolid_evaluate::ReferenceCurveEvaluator;

const EPS: Scalar = 1e-12;

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?} (off by {:e})",
        (actual - expected).abs().max_element()
    );
}

fn assert_frame(
    section: &SectionFrame,
    point: Point3,
    tangent: Vec3,
    lateral: Vec3,
    eps: Scalar,
    what: &str,
) {
    close3(section.point, point, eps, &format!("{what}: point"));
    close3(section.tangent, tangent, eps, &format!("{what}: tangent"));
    close3(section.lateral, lateral, eps, &format!("{what}: lateral"));
    close3(
        section.up,
        tangent.cross(lateral),
        eps,
        &format!("{what}: up = t x l"),
    );
}

fn refused(error: GeomError, needle: &str) {
    assert!(
        error.to_string().contains(needle),
        "expected a refusal naming {needle:?}, got {error}"
    );
}

fn frame2(origin: Point2, heading: Scalar) -> Frame2 {
    let (sin, cos) = heading.sin_cos();
    Frame2 {
        origin,
        x: Vec2::new(cos, sin),
        y: Vec2::new(-sin, cos),
    }
}

#[test]
fn a_line_is_measured_by_length_not_by_parameter() {
    // Speed 5: distance 10 is parameter 2.
    let line = Curve2::Line(Line2 {
        origin: Point2::new(1.0, 2.0),
        direction: Vec2::new(3.0, 4.0),
    });
    let section = station_section2(&line, 10.0).unwrap();
    assert_frame(
        &section,
        Point3::new(7.0, 10.0, 0.0),
        Vec3::new(0.6, 0.8, 0.0),
        Vec3::new(-0.8, 0.6, 0.0),
        EPS,
        "line",
    );
    close3(section.up, Vec3::Z, EPS, "a 2D curve's up is +Z");
    assert_eq!(station_length2(&line).unwrap(), None, "a line is unbounded");
    // Offsets: 2 left, 1.5 up, 0.5 ahead.
    close3(
        section.place(2.0, 1.5, 0.5),
        Point3::new(7.0 - 1.6 + 0.3, 10.0 + 1.2 + 0.4, 1.5),
        EPS,
        "place",
    );
}

#[test]
fn an_arc_station_turns_by_distance_over_radius() {
    let radius = 20.0;
    let circle = Curve2::Circle(Circle2 {
        frame: frame2(Point2::new(5.0, -3.0), 0.0),
        radius,
    });
    for s in [0.0, 7.5, 31.0, 100.0] {
        let angle: Scalar = s / radius;
        let (sin, cos) = angle.sin_cos();
        let section = station_section2(&circle, s).unwrap();
        assert_frame(
            &section,
            Point3::new(5.0 + radius * cos, -3.0 + radius * sin, 0.0),
            Vec3::new(-sin, cos, 0.0),
            // Left of a counter-clockwise circle is its centre.
            Vec3::new(-cos, -sin, 0.0),
            1e-10,
            &format!("arc at {s}"),
        );
    }
    let length = station_length2(&circle).unwrap().unwrap();
    assert!((length - core::f64::consts::TAU * radius).abs() <= 1e-10);
    refused(
        station_section2(&circle, length + 1e-3).unwrap_err(),
        "beyond the curve's length",
    );
}

#[test]
fn a_polyline_station_walks_its_segments() {
    let polyline = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 4.0),
        ],
        closed: false,
    });
    assert_eq!(station_length2(&polyline).unwrap(), Some(7.0));
    let section = station_section2(&polyline, 5.0).unwrap();
    assert_frame(
        &section,
        Point3::new(3.0, 2.0, 0.0),
        Vec3::Y,
        -Vec3::X,
        EPS,
        "second leg",
    );
    // The end itself, and a hair past it within tolerance, read the end.
    let end = station_section2(&polyline, 7.0 + 1e-14).unwrap();
    close3(end.point, Point3::new(3.0, 4.0, 0.0), EPS, "end");
    refused(
        station_section2(&polyline, 7.001).unwrap_err(),
        "beyond the curve's length",
    );
}

/// The clothoid from rest with `kappa = s / A^2` by its Fresnel power
/// series: `x = sum (-1)^n s^(4n+1) / ((4n+1) (2n)! (2A^2)^(2n))`, `y` the
/// odd terms.
fn clothoid_reference(a2: Scalar, s: Scalar) -> (Point2, Scalar) {
    let q = s * s / (2.0 * a2);
    let (mut x, mut y) = (0.0, 0.0);
    let mut term = s; // s q^k / k!
    for k in 0..60 {
        let sign = if (k / 2) % 2 == 0 { 1.0 } else { -1.0 };
        let piece = sign * term / (2 * k + 1) as Scalar;
        if k % 2 == 0 {
            x += piece;
        } else {
            y += piece;
        }
        term *= q / (k + 1) as Scalar;
    }
    (Point2::new(x, y), q)
}

#[test]
fn a_clothoid_station_matches_the_fresnel_series() {
    let (radius, length) = (300.0, 120.0);
    let a2 = radius * length;
    let clothoid = Curve2::Intrinsic(Intrinsic2::new(
        frame2(Point2::ZERO, 0.0),
        CurvatureLaw::clothoid(0.0, 1.0 / radius, length),
        length,
    ));
    for s in [0.0, 40.0, 90.0, 120.0] {
        let (point, heading) = clothoid_reference(a2, s);
        let (sin, cos) = heading.sin_cos();
        let section = station_section2(&clothoid, s).unwrap();
        assert_frame(
            &section,
            Point3::new(point.x, point.y, 0.0),
            Vec3::new(cos, sin, 0.0),
            Vec3::new(-sin, cos, 0.0),
            1e-10,
            &format!("clothoid at {s}"),
        );
    }
    refused(
        station_section2(&clothoid, 121.0).unwrap_err(),
        "beyond the curve's length",
    );
}

#[test]
fn a_chain_station_crosses_into_its_clothoid() {
    // 10 m straight heading 0.4 rad from (1, 2), then the clothoid above.
    let (radius, length) = (300.0, 120.0);
    let start = frame2(Point2::new(1.0, 2.0), 0.4);
    let chain = Curve2::Chain(Chain2::new(
        start,
        vec![
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::circular(0.0),
                length: 10.0,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::clothoid(0.0, 1.0 / radius, length),
                length,
            },
        ],
    ));
    assert_eq!(station_length2(&chain).unwrap(), Some(130.0));
    let (local, heading) = clothoid_reference(radius * length, 70.0);
    let world = |v: Vec2| start.x * v.x + start.y * v.y;
    let origin = start.origin + world(Vec2::new(10.0 + local.x, local.y));
    let tangent = world(Vec2::new(heading.cos(), heading.sin()));
    let section = station_section2(&chain, 80.0).unwrap();
    assert_frame(
        &section,
        Point3::new(origin.x, origin.y, 0.0),
        Vec3::new(tangent.x, tangent.y, 0.0),
        Vec3::new(-tangent.y, tangent.x, 0.0),
        1e-10,
        "chain",
    );
}

fn straight_two_percent() -> Elevated3 {
    Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::new(0.0, 0.0),
            direction: Vec2::new(1.0, 0.0),
        }),
        ElevationLaw::constant_grade(100.0, 0.02),
    )
}

#[test]
fn an_elevated_station_is_a_plan_distance_with_an_unbanked_frame() {
    let curve = Curve3::Elevated(straight_two_percent());
    let section = station_section3(&curve, 50.0).unwrap();
    let k = 1.0004_f64.sqrt();
    // Plan distance 50 is x = 50, not 50 m along the slope.
    assert_frame(
        &section,
        Point3::new(50.0, 0.0, 101.0),
        Vec3::new(1.0 / k, 0.0, 0.02 / k),
        Vec3::Y,
        EPS,
        "elevated",
    );
    close3(
        section.up,
        Vec3::new(-0.02 / k, 0.0, 1.0 / k),
        EPS,
        "up leans back with the grade",
    );
    // The plan frame stands upright.
    let plan = section.plan().unwrap();
    assert_frame(&plan, section.point, Vec3::X, Vec3::Y, EPS, "plan");
    close3(plan.up, Vec3::Z, EPS, "plan up");
    // The provider's layout: x tangent, y up, z right.
    let frame = section.frame();
    let provider = ReferenceCurveEvaluator::new()
        .frame_at(&curve, CurveMeasure::Distance(50.0))
        .unwrap();
    close3(frame.x, provider.x, EPS, "x matches frame_at");
    close3(frame.y, provider.y, EPS, "y matches frame_at");
    close3(frame.z, provider.z, EPS, "z matches frame_at");
    close3(
        frame.origin,
        provider.origin,
        EPS,
        "origin matches frame_at",
    );
}

#[test]
fn an_elevated_arc_station_reads_plan_and_profile_at_one_distance() {
    let radius = 200.0;
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Circle(Circle2 {
            frame: frame2(Point2::ZERO, 0.0),
            radius,
        }),
        ElevationLaw::constant_grade(10.0, -0.03),
    ));
    let s = 60.0;
    let (sin, cos) = (s / radius).sin_cos();
    let k = (1.0_f64 + 0.03 * 0.03).sqrt();
    let section = station_section3(&curve, s).unwrap();
    assert_frame(
        &section,
        Point3::new(radius * cos, radius * sin, 10.0 - 0.03 * s),
        Vec3::new(-sin / k, cos / k, -0.03 / k),
        Vec3::new(-cos, -sin, 0.0),
        1e-10,
        "elevated arc",
    );
}

#[test]
fn a_banked_station_is_framed_by_its_rolled_section() {
    // As the banked tests: 0 -> 150 mm over 100 m, rotation about the
    // tangent, rail heads 1.5 m apart, read at the end of the ramp.
    let curve = Curve3::Banked(Banked3::new(
        straight_two_percent(),
        CantLaw::new(vec![
            CantPiece::linear(100.0, 0.0, 0.15),
            CantPiece::constant(50.0, 0.15),
        ]),
        CantLaw::zero(150.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    assert_eq!(station_length3(&curve).unwrap(), Some(150.0));
    let section = station_section3(&curve, 100.0).unwrap();
    let s = 1.0004_f64.sqrt();
    let cos_psi = 0.99_f64.sqrt();
    let tangent = Vec3::new(1.0 / s, 0.0, 0.02 / s);
    let lateral = Vec3::new(-0.1 * 0.02 / s, cos_psi, 0.1 / s);
    assert_frame(
        &section,
        Point3::new(100.0, 0.0, 102.0),
        tangent,
        lateral,
        EPS,
        "banked",
    );
    close3(
        section.up,
        Vec3::new(-cos_psi * 0.02 / s, -0.1, cos_psi / s),
        EPS,
        "banked up",
    );
    // A lateral offset of b/2 is the left rail head, raised by D cos(theta) / 2.
    let left = section.place(0.75, 0.0, 0.0);
    assert!((left.z - (102.0 + 0.075 / s)).abs() <= EPS);
    refused(
        station_section3(&curve, 150.5).unwrap_err(),
        "beyond the curve's length",
    );
}

#[test]
fn malformed_stations_and_vertical_tangents_are_refused_by_name() {
    let line = Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    });
    refused(
        station_section2(&line, -1.0).unwrap_err(),
        "before the curve's start",
    );
    refused(
        station_section2(&line, Scalar::NAN).unwrap_err(),
        "not finite",
    );
    refused(
        station_section2(&line, Scalar::INFINITY).unwrap_err(),
        "not finite",
    );
    let vertical = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::Z,
    });
    refused(
        station_section3(&vertical, 1.0).unwrap_err(),
        "the tangent is vertical",
    );
    let sloped = station_section3(
        &Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::new(0.0, 3.0, 4.0),
        }),
        5.0,
    )
    .unwrap();
    close3(sloped.point, Point3::new(0.0, 3.0, 4.0), EPS, "3D line");
    let frame = Frame3 {
        origin: sloped.point,
        x: sloped.tangent,
        y: sloped.up,
        z: -sloped.lateral,
    };
    close3(frame.x.cross(frame.y), frame.z, EPS, "right-handed");
}

/// Arc length of `y = a x^3` from 0 to `x` by its binomial series (as
/// `arc_length_chain.rs` pins it), independent of the quadrature.
fn cubic_length(a: Scalar, x: Scalar) -> Scalar {
    let z = 9.0 * a * a * x.powi(4);
    assert!(z < 0.5);
    let (mut binomial, mut power, mut sum) = (1.0, 1.0, 0.0);
    for n in 0..200 {
        let term = binomial * power / (4.0 * n as Scalar + 1.0);
        sum += term;
        if term.abs() < 1e-30 {
            break;
        }
        binomial *= (0.5 - n as Scalar) / (n as Scalar + 1.0);
        power *= z;
    }
    x * sum
}

#[test]
fn a_cubic_parabola_station_inverts_its_arc_length() {
    use axiolid_curve::{BSplineCurve2, KnotSpec};
    let (a, reach) = (1.0 / (6.0 * 300.0 * 60.0), 60.0);
    let cubic = Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(reach / 3.0, 0.0),
            Point2::new(2.0 * reach / 3.0, 0.0),
            Point2::new(reach, a * reach.powi(3)),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    });
    let total = station_length2(&cubic).unwrap().unwrap();
    assert!((total - cubic_length(a, reach)).abs() <= 1e-10);
    let s = 45.0;
    // Invert the series by bisection.
    let (mut lo, mut hi) = (0.0, reach);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if cubic_length(a, mid) < s {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let x = 0.5 * (lo + hi);
    let slope = 3.0 * a * x * x;
    let norm = slope.hypot(1.0);
    let section = station_section2(&cubic, s).unwrap();
    assert_frame(
        &section,
        Point3::new(x, a * x.powi(3), 0.0),
        Vec3::new(1.0 / norm, slope / norm, 0.0),
        Vec3::new(-slope / norm, 1.0 / norm, 0.0),
        1e-10,
        "cubic parabola",
    );
}

#[test]
fn an_explicit_orientation_turns_the_section_in_its_own_axes() {
    // Line heading (0.6, 0.8): t = (0.6, 0.8, 0), l = (-0.8, 0.6, 0), u = Z.
    let line = Curve2::Line(Line2 {
        origin: Point2::new(1.0, 2.0),
        direction: Vec2::new(3.0, 4.0),
    });
    let section = station_section2(&line, 10.0).unwrap();
    let (t, l, u) = (section.tangent, section.lateral, section.up);
    // Defaults: the section itself.
    let same = section.oriented(None, None).unwrap();
    assert_frame(&same, section.point, t, l, EPS, "unturned");
    // Axis (0, 1, 1) and reference (1, 0, 5): up' = (l + u) / sqrt 2, the
    // reference loses its component 5 / sqrt 2 along up' ... in components
    // (1, 0, 5) - 2.5 (0, 1, 1) = (1, -2.5, 2.5), normalised.
    let r = core::f64::consts::FRAC_1_SQRT_2;
    let turned = section
        .oriented(
            Some(Vec3::new(0.0, 1.0, 1.0)),
            Some(Vec3::new(1.0, 0.0, 5.0)),
        )
        .unwrap();
    let n = (1.0_f64 + 12.5).sqrt();
    let tangent = (t - 2.5 * l + 2.5 * u) / n;
    let up = r * (l + u);
    assert_frame(
        &turned,
        section.point,
        tangent,
        up.cross(tangent),
        EPS,
        "turned",
    );
    close3(turned.up, up, EPS, "the axis is exact");
    // Parallel, anti-parallel and degenerate pairs are refused by name.
    refused(
        section.oriented(Some(Vec3::Y), Some(-Vec3::Y)).unwrap_err(),
        "parallel",
    );
    refused(
        section.oriented(None, Some(Vec3::Z)).unwrap_err(),
        "parallel",
    );
    refused(
        section.oriented(Some(Vec3::ZERO), None).unwrap_err(),
        "axis",
    );
    refused(
        section
            .oriented(None, Some(Vec3::splat(Scalar::INFINITY)))
            .unwrap_err(),
        "reference direction",
    );
}
