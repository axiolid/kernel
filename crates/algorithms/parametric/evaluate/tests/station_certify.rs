//! Certified joint distances and seam-snapping windows (#294).
//!
//! Every expected length here is independent of the code under test: a
//! closed form where one exists (the offset of a clothoid by `d` is `L - d
//! theta`), else the offset curve written out here in components and its
//! speed integrated by an adaptive Simpson rule to `1e-13`, far inside the
//! tolerances asked of the kernel.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, Chain2, ChainPiece2, Circle2, CurvatureLaw, Curve2, Curve3, CurvePath, Elevated3,
    ElevationLaw, Intrinsic2, KnotSpec, Line2, Line3, OffsetFrame, OffsetLaw, PathCurve,
    PathOffset, PathOffsets, PathPiece, Polyline2, SeamSide,
};
use axiolid_curve_evaluate_contract::{CurveEvaluator, CurveMeasure, DistanceInterval};
use axiolid_evaluate::station::{
    offset_pieces, snapped_distance2, CompositeBasis, StationCurve, StationPiece,
    SEAM_WINDOW_AMBIGUOUS, SEAM_WINDOW_STRADDLED, UNCERTIFIED_LENGTH, WINDOW_CERTIFICATION,
};
use axiolid_evaluate::ReferenceCurveEvaluator;

// --- an independent reference -------------------------------------------------

/// Adaptive Simpson, each panel refined until its Richardson estimate is
/// within its share of `tolerance`.
fn simpson(f: &dyn Fn(Scalar) -> Scalar, a: Scalar, b: Scalar, tolerance: Scalar) -> Scalar {
    fn rule(f: &dyn Fn(Scalar) -> Scalar, a: Scalar, b: Scalar) -> Scalar {
        (b - a) / 6.0 * (f(a) + 4.0 * f(0.5 * (a + b)) + f(b))
    }
    fn refine(
        f: &dyn Fn(Scalar) -> Scalar,
        a: Scalar,
        b: Scalar,
        whole: Scalar,
        tolerance: Scalar,
        depth: u32,
    ) -> Scalar {
        let mid = 0.5 * (a + b);
        let (left, right) = (rule(f, a, mid), rule(f, mid, b));
        let delta = left + right - whole;
        if depth > 40 || delta.abs() <= 15.0 * tolerance {
            return left + right + delta / 15.0;
        }
        refine(f, a, mid, left, 0.5 * tolerance, depth + 1)
            + refine(f, mid, b, right, 0.5 * tolerance, depth + 1)
    }
    // Start from many panels, so no feature hides between the first nodes.
    let panels = 64;
    let h = (b - a) / panels as Scalar;
    (0..panels)
        .map(|i| {
            let (lo, hi) = (a + h * i as Scalar, a + h * (i + 1) as Scalar);
            refine(f, lo, hi, rule(f, lo, hi), tolerance / panels as Scalar, 0)
        })
        .sum()
}

fn certifies(interval: DistanceInterval, truth: Scalar, tolerance: Scalar, what: &str) {
    assert!(
        interval.lower <= interval.upper && interval.width() <= tolerance,
        "{what}: {interval:?} is wider than {tolerance:e}"
    );
    assert!(
        interval.contains(truth),
        "{what}: {interval:?} does not contain {truth} (off by {:e})",
        if truth < interval.lower {
            interval.lower - truth
        } else {
            truth - interval.upper
        }
    );
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, input: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input: found, .. }) => assert_eq!(found, input),
        other => panic!("expected a refusal naming {input:?}, got {other:?}"),
    }
}

// --- bases --------------------------------------------------------------------

fn frame(origin: Point2, x: Vec2, y: Vec2) -> Frame2 {
    Frame2 { origin, x, y }
}

fn line2(origin: Point2, direction: Vec2) -> Curve2 {
    Curve2::Line(Line2 { origin, direction })
}

/// `z = 10 + 0.02 v - 0.0002 v^2`: a crest from +2% to -2% over 100 m.
fn crest() -> ElevationLaw {
    ElevationLaw::Polynomial {
        coefficients: vec![10.0, 0.02, -0.0002],
    }
}

fn crest_grade(v: Scalar) -> (Scalar, Scalar) {
    (0.02 - 0.0004 * v, -0.0004)
}

/// An offset by distances at stations `(distance, lateral, vertical)` along
/// `curve`, in its section frame: one law per interval, as the compiler
/// flattens it.
fn by_stations<'c>(
    curve: StationCurve<'c>,
    stations: &[(Scalar, Scalar, Scalar)],
) -> Vec<StationPiece<'c>> {
    let mut out = Vec::new();
    for pair in stations.windows(2) {
        let ((a, la, va), (b, lb, vb)) = (pair[0], pair[1]);
        let span = StationPiece::between(curve, a, b).unwrap();
        let law = OffsetLaw::Linear {
            start: PathOffsets::new(la, va, 0.0),
            end: PathOffsets::new(lb, vb, 0.0),
            frame: OffsetFrame::Section,
        };
        out.extend(offset_pieces(&[span], law).unwrap());
    }
    out
}

/// The plan of an offset by distances beside an elevated curve, in its
/// section frame, written out: `P(v) = B(v) + a N + b T`, `b = -vertical
/// sin(atan z')`, its plan speed from the components' derivatives.
fn plan_speed_beside(
    plan_velocity: &dyn Fn(Scalar, Scalar, Scalar, Scalar, Scalar) -> Vec2,
    station: [(Scalar, Scalar, Scalar); 2],
) -> impl Fn(Scalar) -> Scalar + '_ {
    let ((a0, l0, w0), (a1, l1, w1)) = (station[0], station[1]);
    move |v: Scalar| {
        let u = (v - a0) / (a1 - a0);
        let (lateral, dlateral) = (l0 + (l1 - l0) * u, (l1 - l0) / (a1 - a0));
        let (vertical, dvertical) = (w0 + (w1 - w0) * u, (w1 - w0) / (a1 - a0));
        let (g, dg) = crest_grade(v);
        let sin = g / (1.0 + g * g).sqrt();
        let dsin = dg / (1.0 + g * g).powf(1.5);
        let b = -vertical * sin;
        let db = -(dvertical * sin + vertical * dsin);
        plan_velocity(v, lateral, dlateral, b, db).length()
    }
}

#[test]
fn an_offset_by_distances_beside_an_elevated_line_is_certified() {
    let elevated = Curve3::Elevated(Elevated3::new(line2(Point2::ZERO, Vec2::X), crest()));
    let stations = [(0.0, 3.0, -0.2), (40.0, 3.5, 0.1), (100.0, 4.5, 0.3)];
    let basis =
        CompositeBasis::new(by_stations(StationCurve::Three(&elevated), &stations)).unwrap();
    // Along +x: B = (v, 0), N = +y, T = +x, no curvature.
    let velocity =
        |_: Scalar, _: Scalar, dl: Scalar, _: Scalar, db: Scalar| Vec2::new(1.0 + db, dl);
    let speed = plan_speed_beside(&velocity, [stations[0], stations[1]]);
    let truth = simpson(&speed, 0.0, 40.0, 1e-13);
    for tolerance in [1e-6, 1e-9] {
        let joints = basis.joint_intervals(tolerance).unwrap();
        assert_eq!(joints.len(), 1);
        certifies(joints[0], truth, tolerance, "elevated line");
        assert!(!joints[0].is_point());
    }
    // The kernel's own joint, an estimate, agrees with it.
    assert!((basis.pieces()[0].length() - truth).abs() <= 1e-9 * truth);
}

#[test]
fn an_offset_by_distances_beside_an_elevated_arc_is_certified() {
    let r = 200.0;
    // Starts at the origin heading +x, turning left about (0, 200).
    let plan = Curve2::Circle(Circle2 {
        frame: frame(Point2::new(0.0, r), Vec2::new(0.0, -1.0), Vec2::X),
        radius: r,
    });
    let elevated = Curve3::Elevated(Elevated3::new(plan, crest()));
    let stations = [(0.0, -2.0, 0.15), (40.0, -2.5, -0.1), (100.0, -1.5, 0.2)];
    let basis =
        CompositeBasis::new(by_stations(StationCurve::Three(&elevated), &stations)).unwrap();
    // B = (r sin t, r - r cos t), T = (cos t, sin t), N = (-sin t, cos t),
    // t = v / r: P = B + a N + b T, differentiated component by component.
    let velocity = move |v: Scalar, a: Scalar, da: Scalar, b: Scalar, db: Scalar| {
        let t = v / r;
        let (sin, cos) = t.sin_cos();
        let dt = 1.0 / r;
        Vec2::new(
            r * cos * dt - da * sin - a * cos * dt + db * cos - b * sin * dt,
            r * sin * dt + da * cos - a * sin * dt + db * sin + b * cos * dt,
        )
    };
    let speed = plan_speed_beside(&velocity, [stations[0], stations[1]]);
    let truth = simpson(&speed, 0.0, 40.0, 1e-13);
    for tolerance in [1e-6, 1e-9] {
        let joints = basis.joint_intervals(tolerance).unwrap();
        certifies(joints[0], truth, tolerance, "elevated arc");
    }
    // Outside the arc the offset is longer than its base.
    assert!(truth > 40.0);
}

/// A clothoid from straight to `1 / 300` over 100 m, then the arc it
/// reaches: the seam at 100 splits the offset into two pieces.
fn clothoid_then_arc() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        frame(Point2::new(5.0, -3.0), Vec2::X, Vec2::Y),
        CurvatureLaw::piecewise(
            vec![100.0],
            vec![
                CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 100.0),
                CurvatureLaw::Constant {
                    curvature: 1.0 / 300.0,
                },
            ],
        ),
        150.0,
    ))
}

#[test]
fn a_constant_offset_of_a_clothoid_is_certified() {
    let spiral = clothoid_then_arc();
    let whole = StationPiece::whole(StationCurve::Two(&spiral)).unwrap();
    // The spiral turns by 100 / 600 rad; beside it at `d` the offset runs
    // `1 - d k` times as fast, so its length is `100 - d / 6`.
    for distance in [5.0, -5.0] {
        let basis =
            CompositeBasis::new(offset_pieces(&[whole], OffsetLaw::Planar { distance }).unwrap())
                .unwrap();
        let truth = 100.0 - distance / 6.0;
        let kappa = |s: Scalar| s / 30_000.0;
        let quadrature = simpson(&|s| 1.0 - distance * kappa(s), 0.0, 100.0, 1e-13);
        assert!((quadrature - truth).abs() <= 1e-11);
        for tolerance in [1e-6, 1e-10] {
            let joints = basis.joint_intervals(tolerance).unwrap();
            certifies(joints[0], truth, tolerance, "clothoid");
        }
    }
}

/// A cubic Bezier as a B-spline, then the line continuing its end tangent.
fn bezier() -> ([Point2; 4], Curve2) {
    let points = [
        Point2::new(0.0, 0.0),
        Point2::new(30.0, 0.0),
        Point2::new(60.0, 20.0),
        Point2::new(90.0, 20.0),
    ];
    (
        points,
        Curve2::BSpline(BSplineCurve {
            degree: 3,
            control_points: points.to_vec(),
            knots: vec![0.0, 1.0],
            multiplicities: vec![4, 4],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }),
    )
}

#[test]
fn an_offset_of_a_b_spline_is_certified() {
    let (p, spline) = bezier();
    let tail = line2(Point2::new(90.0, 20.0), Vec2::X);
    let base = [
        StationPiece::whole(StationCurve::Two(&spline)).unwrap(),
        StationPiece::between(StationCurve::Two(&tail), 0.0, 10.0).unwrap(),
    ];
    let distance = 2.0;
    let basis =
        CompositeBasis::new(offset_pieces(&base, OffsetLaw::Planar { distance }).unwrap()).unwrap();
    // Bernstein derivatives, and the offset's velocity B' + d (rot B'/|B'|)'.
    let first = |t: Scalar| {
        let u = 1.0 - t;
        3.0 * (u * u * (p[1] - p[0]) + 2.0 * u * t * (p[2] - p[1]) + t * t * (p[3] - p[2]))
    };
    let second =
        |t: Scalar| 6.0 * ((1.0 - t) * (p[2] - 2.0 * p[1] + p[0]) + t * (p[3] - 2.0 * p[2] + p[1]));
    let speed = |t: Scalar| {
        let (d1, d2) = (first(t), second(t));
        let n = d1.length();
        let unit_rate = (d2 * (n * n) - d1 * d1.dot(d2)) / (n * n * n);
        (d1 + distance * unit_rate.perp()).length()
    };
    let truth = simpson(&speed, 0.0, 1.0, 1e-13);
    for tolerance in [1e-4, 1e-6, 1e-7] {
        let joints = basis.joint_intervals(tolerance).unwrap();
        certifies(joints[0], truth, tolerance, "B-spline");
    }
    // The plain B-spline's own length is certified the same way.
    let plain = CompositeBasis::new(base.to_vec()).unwrap();
    let length = simpson(&|t| first(t).length(), 0.0, 1.0, 1e-13);
    certifies(
        plain.joint_intervals(1e-7).unwrap()[0],
        length,
        1e-7,
        "plain B-spline",
    );
}

#[test]
fn an_offset_of_a_3d_b_spline_in_its_leaning_section_is_certified() {
    let p = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(30.0, 0.0, 2.0),
        Point3::new(60.0, 20.0, 5.0),
        Point3::new(90.0, 20.0, 4.0),
    ];
    let spline = Curve3::BSpline(BSplineCurve {
        degree: 3,
        control_points: p.to_vec(),
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let tail = Curve3::Line(Line3 {
        origin: p[3],
        direction: p[3] - p[2],
    });
    let base = [
        StationPiece::whole(StationCurve::Three(&spline)).unwrap(),
        StationPiece::between(StationCurve::Three(&tail), 0.0, 1.0).unwrap(),
    ];
    let (lateral, vertical) = (1.5, 0.4);
    let law = OffsetLaw::Linear {
        start: PathOffsets::new(lateral, vertical, 0.0),
        end: PathOffsets::new(lateral, vertical, 0.0),
        frame: OffsetFrame::Section,
    };
    let basis = CompositeBasis::new(offset_pieces(&base, law).unwrap()).unwrap();
    // The offset written out (Bernstein point and derivative, the
    // reference-up frame by hand), its speed by a five-point difference.
    let point = |t: Scalar| {
        let u = 1.0 - t;
        let b = p[0] * (u * u * u)
            + p[1] * (3.0 * u * u * t)
            + p[2] * (3.0 * u * t * t)
            + p[3] * (t * t * t);
        let d = 3.0 * (u * u * (p[1] - p[0]) + 2.0 * u * t * (p[2] - p[1]) + t * t * (p[3] - p[2]));
        let tangent = d.normalize();
        let left = Vec3::Z.cross(d).normalize();
        b + lateral * left + vertical * tangent.cross(left)
    };
    let h = 1e-3;
    let speed = |t: Scalar| {
        ((point(t - 2.0 * h) - point(t + 2.0 * h)) + 8.0 * (point(t + h) - point(t - h))).length()
            / (12.0 * h)
    };
    let truth = simpson(&speed, 0.0, 1.0, 1e-11);
    let joints = basis.joint_intervals(1e-6).unwrap();
    certifies(joints[0], truth, 1e-6, "3D B-spline");
}

#[test]
fn exact_lengths_are_point_intervals() {
    let line = line2(Point2::ZERO, Vec2::new(2.0, 0.0));
    let arc = Curve2::Circle(Circle2 {
        frame: frame(Point2::new(15.0, 5.0), Vec2::new(0.0, -1.0), Vec2::X),
        radius: 5.0,
    });
    let quarter = core::f64::consts::FRAC_PI_2 * 5.0;
    let up = line2(Point2::new(20.0, 5.0), Vec2::Y);
    let base = [
        StationPiece::between(StationCurve::Two(&line), 0.0, 15.0).unwrap(),
        StationPiece::between(StationCurve::Two(&arc), 0.0, quarter).unwrap(),
        StationPiece::between(StationCurve::Two(&up), 0.0, 1.0).unwrap(),
    ];
    let plain = CompositeBasis::new(base.to_vec()).unwrap();
    let joints = plain.joint_intervals(1e-9).unwrap();
    assert_eq!(
        joints,
        vec![
            DistanceInterval::point(15.0),
            DistanceInterval::point(15.0 + quarter),
        ]
    );
    // A constant offset of a line is a line, of a circle a circle of
    // radius r - d: both closed forms.
    let offset =
        CompositeBasis::new(offset_pieces(&base, OffsetLaw::Planar { distance: 1.0 }).unwrap())
            .unwrap();
    let joints = offset.joint_intervals(1e-9).unwrap();
    assert!(joints.iter().all(DistanceInterval::is_point), "{joints:?}");
    assert_eq!(joints[0].lower, 15.0);
    assert!((joints[1].lower - (15.0 + quarter * 4.0 / 5.0)).abs() <= 1e-12);
    // Through the contract, the same.
    let path = offset.path();
    assert_eq!(
        ReferenceCurveEvaluator::new()
            .path_joint_distances(&path, 1e-9)
            .unwrap(),
        joints
    );
}

#[test]
fn an_uncertifiable_length_and_a_bad_tolerance_are_refused_by_name() {
    let chain = parametric_chain(Vec2::ZERO);
    let whole = StationPiece::whole(StationCurve::Two(&chain)).unwrap();
    let tail = line2(Point2::new(90.0, 20.0), Vec2::X);
    let next = StationPiece::between(StationCurve::Two(&tail), 0.0, 10.0).unwrap();
    // An offset by distances along a chain's parametric piece.
    let widening = OffsetLaw::Linear {
        start: PathOffsets::new(1.0, 0.0, 0.0),
        end: PathOffsets::new(2.0, 0.0, 0.0),
        frame: OffsetFrame::Section,
    };
    let mut pieces = offset_pieces(&[whole], widening).unwrap();
    pieces.extend(offset_pieces(&[next], OffsetLaw::Planar { distance: 2.0 }).unwrap());
    let basis = CompositeBasis::new(pieces).unwrap();
    refused(basis.joint_intervals(1e-6), UNCERTIFIED_LENGTH);
    // A 3D offset whose direction is not vertical beside an elevated curve.
    let elevated = Curve3::Elevated(Elevated3::new(line2(Point2::ZERO, Vec2::X), crest()));
    let span = [
        StationPiece::between(StationCurve::Three(&elevated), 0.0, 40.0).unwrap(),
        StationPiece::between(StationCurve::Three(&elevated), 40.0, 100.0).unwrap(),
    ];
    let tilted = OffsetLaw::Directed {
        distance: 1.0,
        reference_direction: Vec3::new(0.3, 0.0, 1.0),
    };
    let basis = CompositeBasis::new(offset_pieces(&span, tilted).unwrap()).unwrap();
    refused(basis.joint_intervals(1e-6), UNCERTIFIED_LENGTH);
    // An upright one is certified: a plan offset by 1.
    let upright = OffsetLaw::Directed {
        distance: 1.0,
        reference_direction: Vec3::Z,
    };
    let basis = CompositeBasis::new(offset_pieces(&span, upright).unwrap()).unwrap();
    certifies(
        basis.joint_intervals(1e-9).unwrap()[0],
        40.0,
        1e-9,
        "upright",
    );
    for bad in [0.0, -1.0, Scalar::NAN, Scalar::INFINITY] {
        assert!(matches!(
            basis.joint_intervals(bad),
            Err(GeomError::InvalidInput(_))
        ));
    }
}

// --- seam-snapping windows ------------------------------------------------------

/// A line along +x over [0, 10], then a quarter arc turning left at a
/// right angle.
fn line_then_arc() -> (Curve3, Curve3) {
    (
        Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::new(2.0, 0.0, 0.0),
        }),
        // About (15, 0), from (10, 0) heading +y: a right-angle corner.
        Curve3::Circle(axiolid_curve::Circle3 {
            frame: axiolid_core::Frame3 {
                origin: Point3::new(15.0, 0.0, 0.0),
                x: -Vec3::X,
                y: Vec3::Y,
                z: -Vec3::Z,
            },
            radius: 5.0,
        }),
    )
}

#[test]
fn a_station_within_the_window_reads_the_joint_from_its_side() {
    let (line, arc) = line_then_arc();
    let quarter = core::f64::consts::FRAC_PI_2 * 5.0;
    let basis = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Three(&line), 0.0, 10.0).unwrap(),
        StationPiece::between(StationCurve::Three(&arc), 0.0, quarter).unwrap(),
    ])
    .unwrap();
    let window = 1e-3;
    for near in [10.0 - 6e-4, 10.0 + 6e-4, 10.0] {
        let snapped = basis.snapped_distance(near, window).unwrap();
        assert_eq!(snapped.to_bits(), 10.0_f64.to_bits(), "{near}");
        for side in [SeamSide::Incoming, SeamSide::Outgoing] {
            assert_eq!(
                basis.section_on(snapped, side).unwrap(),
                basis.section_on(10.0, side).unwrap()
            );
        }
    }
    // Outgoing on the joint is the arc's start, heading +y; unsnapped, 0.6
    // mm before it, still the line.
    let on = basis.section_on(
        basis.snapped_distance(10.0 - 6e-4, window).unwrap(),
        SeamSide::Outgoing,
    );
    assert!((on.unwrap().tangent - Vec3::Y).length() < 1e-12);
    let off = basis.section_on(10.0 - 6e-4, SeamSide::Outgoing).unwrap();
    assert!((off.tangent - Vec3::X).length() < 1e-12);
    // Outside the window, and with no window, nothing moves.
    assert_eq!(
        basis.snapped_distance(10.0 - 2e-3, window).unwrap(),
        10.0 - 2e-3
    );
    assert_eq!(
        basis.snapped_distance(10.0 - 6e-4, 0.0).unwrap(),
        10.0 - 6e-4
    );
    for bad in [-1e-3, Scalar::NAN, Scalar::INFINITY] {
        assert!(matches!(
            basis.snapped_distance(10.0, bad),
            Err(GeomError::InvalidInput(_))
        ));
    }
}

#[test]
fn two_seams_in_one_window_are_refused_by_name() {
    let (line, _) = line_then_arc();
    let short = 5e-4;
    let basis = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Three(&line), 0.0, 10.0).unwrap(),
        StationPiece::between(StationCurve::Three(&line), 10.0, 10.0 + short).unwrap(),
        StationPiece::between(StationCurve::Three(&line), 10.0 + short, 20.0).unwrap(),
    ])
    .unwrap();
    refused(
        basis.snapped_distance(10.0 + 0.5 * short, 1e-3),
        SEAM_WINDOW_AMBIGUOUS,
    );
    // A window holding only one of them snaps onto it.
    assert_eq!(basis.snapped_distance(10.0 - 2e-4, 4e-4).unwrap(), 10.0);
}

#[test]
fn a_numerical_joint_snaps_inside_and_is_refused_on_the_window_edge() {
    let elevated = Curve3::Elevated(Elevated3::new(line2(Point2::ZERO, Vec2::X), crest()));
    let stations = [(0.0, 3.0, -0.2), (40.0, 3.5, 0.1), (100.0, 4.5, 0.3)];
    let basis =
        CompositeBasis::new(by_stations(StationCurve::Three(&elevated), &stations)).unwrap();
    let joint = basis.pieces()[0].length();
    let window = 1e-3;
    let certified = basis
        .joint_intervals(window * WINDOW_CERTIFICATION)
        .unwrap()[0];
    assert!(certified.width() > 0.0 && certified.contains(joint));
    for near in [joint - 6e-4, joint + 6e-4] {
        assert_eq!(basis.snapped_distance(near, window).unwrap(), joint);
    }
    // The window's edge through the middle of the joint's interval: inside
    // or outside cannot be told.
    refused(
        basis.snapped_distance(certified.midpoint() - window, window),
        SEAM_WINDOW_STRADDLED,
    );
    refused(
        basis.snapped_distance(certified.midpoint() + window, window),
        SEAM_WINDOW_STRADDLED,
    );
}

#[test]
fn an_atomic_curve_snaps_onto_its_vertex_or_corner_knot() {
    let corner = Curve2::Polyline(Polyline2 {
        points: vec![Point2::ZERO, Point2::new(3.0, 4.0), Point2::new(3.0, 10.0)],
        closed: false,
    });
    assert_eq!(snapped_distance2(&corner, 5.0 + 6e-4, 1e-3).unwrap(), 5.0);
    // A smooth seam (the clothoid meeting its arc) has one frame: no
    // snapping onto it.
    let spiral = clothoid_then_arc();
    assert_eq!(
        snapped_distance2(&spiral, 100.0 - 6e-4, 1e-3).unwrap(),
        100.0 - 6e-4
    );
    assert_eq!(
        snapped_distance2(&corner, 5.0 + 2e-3, 1e-3).unwrap(),
        5.0 + 2e-3
    );
    // A B-spline's corner knot, located by quadrature, certified.
    let spline = Curve2::BSpline(BSplineCurve {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
        ],
        knots: vec![0.0, 1.0, 2.0],
        multiplicities: vec![3, 2, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let knot = axiolid_evaluate::station::station_seams2(&spline).unwrap()[0].distance;
    assert!((knot - 2.0).abs() < 1e-9);
    assert_eq!(snapped_distance2(&spline, knot - 6e-4, 1e-3).unwrap(), knot);
}

#[test]
fn a_curve_path_answers_its_joints_through_the_contract() {
    let spiral = clothoid_then_arc();
    let whole = StationPiece::whole(StationCurve::Two(&spiral)).unwrap();
    let basis =
        CompositeBasis::new(offset_pieces(&[whole], OffsetLaw::Planar { distance: 5.0 }).unwrap())
            .unwrap();
    let path: CurvePath = basis.path();
    let joints = ReferenceCurveEvaluator::new()
        .path_joint_distances(&path, 1e-8)
        .unwrap();
    certifies(joints[0], 100.0 - 5.0 / 6.0, 1e-8, "path");
    // A path that does not meet is refused, as every path query refuses it.
    let gapped = CurvePath::new(vec![
        PathPiece::new(PathCurve::Two(line2(Point2::ZERO, Vec2::X)), 0.0, 1.0),
        PathPiece::new(
            PathCurve::Two(line2(Point2::new(2.0, 0.0), Vec2::X)),
            0.0,
            1.0,
        ),
    ]);
    assert!(ReferenceCurveEvaluator::new()
        .path_joint_distances(&gapped, 1e-8)
        .is_err());
}

// --- open ends (#298) -----------------------------------------------------------

/// Every reading of `open` through the contract equals, bitwise, the
/// station reading of `stated` (the composite a compiler would build) at
/// each of `distances`, on both sides; and `open` resolves to `stated`'s
/// own pieces.
fn reads_as_stated(open: &CurvePath, stated: &CompositeBasis<'_>, distances: &[Scalar]) {
    let resolved = CompositeBasis::from_path(open).unwrap();
    assert_eq!(resolved.path(), stated.path());
    assert_eq!(resolved.length().to_bits(), stated.length().to_bits());
    let provider = ReferenceCurveEvaluator::new();
    assert_eq!(provider.path_distance_convention(open), stated.convention());
    for &d in distances {
        let at = CurveMeasure::Distance(d);
        for side in [SeamSide::Outgoing, SeamSide::Incoming] {
            let want = stated.section_on(d, side).unwrap();
            let point = provider.path_point_at_on(open, at, side).unwrap();
            let tangent = provider.path_tangent_at_on(open, at, side).unwrap();
            assert_eq!(point, want.point, "point at {d} {side:?}");
            assert_eq!(tangent, want.tangent, "tangent at {d} {side:?}");
        }
    }
    let past = stated.length() * (1.0 + 1e-9) + 1e-9;
    assert!(provider
        .path_point_at(open, CurveMeasure::Distance(past))
        .is_err());
}

#[test]
fn an_open_offset_by_distances_beside_an_elevated_arc_needs_no_length() {
    let r = 200.0;
    let plan = Curve2::Circle(Circle2 {
        frame: frame(Point2::new(0.0, r), Vec2::new(0.0, -1.0), Vec2::X),
        radius: r,
    });
    let elevated = Curve3::Elevated(Elevated3::new(plan, crest()));
    // Dyadic offsets, so the compiler's interpolated laws are these exactly.
    let stations = [(0.0, -2.0, 0.25), (40.0, -2.5, -0.125), (100.0, -1.5, 0.5)];
    let stated =
        CompositeBasis::new(by_stations(StationCurve::Three(&elevated), &stations)).unwrap();
    // The caller states the stations along the alignment, which it authored,
    // and leaves every offset's own length to the kernel.
    let open: CurvePath = stations
        .windows(2)
        .map(|pair| {
            let ((a, la, va), (b, lb, vb)) = (pair[0], pair[1]);
            let base = PathPiece::new(PathCurve::Three(elevated.clone()), a, b);
            let law = OffsetLaw::Linear {
                start: PathOffsets::new(la, va, 0.0),
                end: PathOffsets::new(lb, vb, 0.0),
                frame: OffsetFrame::Section,
            };
            PathPiece::whole(PathCurve::Offset(Box::new(PathOffset::new(base, law))))
        })
        .collect();
    assert!(open.pieces().iter().all(|piece| !piece.is_stated()));
    let joint = stated.pieces()[0].length();
    reads_as_stated(
        &open,
        &stated,
        &[
            0.0,
            17.0,
            joint,
            joint + 1e-13 * joint,
            70.0,
            stated.length(),
        ],
    );
    // Where the joint truly lies, against the offset written out.
    let velocity = move |v: Scalar, a: Scalar, da: Scalar, b: Scalar, db: Scalar| {
        let t = v / r;
        let (sin, cos) = t.sin_cos();
        let dt = 1.0 / r;
        Vec2::new(
            r * cos * dt - da * sin - a * cos * dt + db * cos - b * sin * dt,
            r * sin * dt + da * cos - a * sin * dt + db * sin + b * cos * dt,
        )
    };
    let speed = plan_speed_beside(&velocity, [stations[0], stations[1]]);
    let truth = simpson(&speed, 0.0, 40.0, 1e-13);
    let provider = ReferenceCurveEvaluator::new();
    for tolerance in [1e-6, 1e-9] {
        let joints = provider.path_joint_distances(&open, tolerance).unwrap();
        assert_eq!(joints.len(), 1);
        certifies(joints[0], truth, tolerance, "open elevated arc");
    }
    // The kernel's own end, an estimate, is within the offset tolerance.
    assert!((joint - truth).abs() <= 1e-9 * truth);
}

#[test]
fn an_open_constant_offset_of_a_clothoid_needs_no_length() {
    let spiral = clothoid_then_arc();
    let whole = StationPiece::whole(StationCurve::Two(&spiral)).unwrap();
    let distance = 5.0;
    let law = OffsetLaw::Planar { distance };
    let stated = CompositeBasis::new(offset_pieces(&[whole], law).unwrap()).unwrap();
    // The spiral's seam at 100 splits it, as the compiler splits it; the
    // arc after it runs to the curve's end, unstated.
    let offset =
        |base: PathPiece| PathPiece::whole(PathCurve::Offset(Box::new(PathOffset::new(base, law))));
    let open = CurvePath::new(vec![
        offset(PathPiece::from_start(PathCurve::Two(spiral.clone()), 100.0)),
        offset(PathPiece::to_end(PathCurve::Two(spiral.clone()), 100.0)),
    ]);
    let truth = 100.0 - distance / 6.0;
    reads_as_stated(&open, &stated, &[0.0, 50.0, truth, 120.0, stated.length()]);
    // Beside the arc of radius 300 the offset runs `(300 - 5) / 300` as fast.
    let total = truth + 50.0 * (300.0 - distance) / 300.0;
    assert!((stated.length() - total).abs() <= 1e-9 * total);
    for tolerance in [1e-6, 1e-10] {
        let joints = ReferenceCurveEvaluator::new()
            .path_joint_distances(&open, tolerance)
            .unwrap();
        certifies(joints[0], truth, tolerance, "open clothoid");
    }
    // An open end still checks the start it states.
    for start in [151.0, Scalar::NAN] {
        let past = CurvePath::from(PathPiece::to_end(PathCurve::Two(spiral.clone()), start));
        assert!(CompositeBasis::from_path(&past).is_err(), "start {start}");
    }
}

#[test]
fn an_open_offset_of_a_b_spline_needs_no_length() {
    let (p, spline) = bezier();
    let tail = line2(Point2::new(90.0, 20.0), Vec2::X);
    let distance = 2.0;
    let law = OffsetLaw::Planar { distance };
    let base = [
        StationPiece::whole(StationCurve::Two(&spline)).unwrap(),
        StationPiece::between(StationCurve::Two(&tail), 0.0, 10.0).unwrap(),
    ];
    let stated = CompositeBasis::new(offset_pieces(&base, law).unwrap()).unwrap();
    // Neither the spline's arc length nor its offset's is stated.
    let offset =
        |base: PathPiece| PathPiece::whole(PathCurve::Offset(Box::new(PathOffset::new(base, law))));
    let open = CurvePath::new(vec![
        offset(PathPiece::whole(PathCurve::Two(spline.clone()))),
        offset(PathPiece::from_start(PathCurve::Two(tail.clone()), 10.0)),
    ]);
    let joint = stated.pieces()[0].length();
    reads_as_stated(
        &open,
        &stated,
        &[0.0, 33.3, joint, joint - 1e-13 * joint, joint + 4.0],
    );
    let first = |t: Scalar| {
        let u = 1.0 - t;
        3.0 * (u * u * (p[1] - p[0]) + 2.0 * u * t * (p[2] - p[1]) + t * t * (p[3] - p[2]))
    };
    let second =
        |t: Scalar| 6.0 * ((1.0 - t) * (p[2] - 2.0 * p[1] + p[0]) + t * (p[3] - 2.0 * p[2] + p[1]));
    let speed = |t: Scalar| {
        let (d1, d2) = (first(t), second(t));
        let n = d1.length();
        let unit_rate = (d2 * (n * n) - d1 * d1.dot(d2)) / (n * n * n);
        (d1 + distance * unit_rate.perp()).length()
    };
    let truth = simpson(&speed, 0.0, 1.0, 1e-13);
    for tolerance in [1e-4, 1e-7] {
        let joints = ReferenceCurveEvaluator::new()
            .path_joint_distances(&open, tolerance)
            .unwrap();
        certifies(joints[0], truth, tolerance, "open B-spline offset");
    }
    assert!((joint - truth).abs() <= 1e-9 * truth);
    // The plain spline left open: its own arc length, certified.
    let plain = CurvePath::new(vec![
        PathPiece::whole(PathCurve::Two(spline.clone())),
        PathPiece::from_start(PathCurve::Two(tail.clone()), 10.0),
    ]);
    let length = simpson(&|t| first(t).length(), 0.0, 1.0, 1e-13);
    certifies(
        ReferenceCurveEvaluator::new()
            .path_joint_distances(&plain, 1e-7)
            .unwrap()[0],
        length,
        1e-7,
        "open B-spline",
    );
    assert_eq!(
        CompositeBasis::from_path(&plain).unwrap().pieces()[0].end,
        base[0].end
    );
}

// --- joints near a window (#305) ------------------------------------------------

/// The Bezier of [`bezier`] moved by `shift`.
fn bezier_at(shift: Vec2) -> ([Point2; 4], Curve2) {
    let (mut p, mut curve) = bezier();
    for point in &mut p {
        *point += shift;
    }
    if let Curve2::BSpline(spline) = &mut curve {
        spline.control_points = p.to_vec();
    }
    (p, curve)
}

/// The Bezier of [`bezier`] moved by `shift`, read as the one parametric
/// piece of an arc-length chain: an offset of it has no certified length,
/// so it stands in for a piece that cannot be certified.
fn parametric_chain(shift: Vec2) -> Curve2 {
    let (p, curve) = bezier();
    Curve2::Chain(Chain2::new(
        frame(Point2::new(shift.x, shift.y), Vec2::X, Vec2::Y),
        vec![ChainPiece2::Parametric {
            curve,
            start: 0.0,
            length: bezier_offset_length(p, 0.0),
        }],
    ))
}

/// The arc length of the offset by `distance` of the Bezier with control
/// points `p`, written out and integrated independently.
fn bezier_offset_length(p: [Point2; 4], distance: Scalar) -> Scalar {
    let first = |t: Scalar| {
        let u = 1.0 - t;
        3.0 * (u * u * (p[1] - p[0]) + 2.0 * u * t * (p[2] - p[1]) + t * t * (p[3] - p[2]))
    };
    let second =
        |t: Scalar| 6.0 * ((1.0 - t) * (p[2] - 2.0 * p[1] + p[0]) + t * (p[3] - 2.0 * p[2] + p[1]));
    let speed = |t: Scalar| {
        let (d1, d2) = (first(t), second(t));
        let n = d1.length();
        let unit_rate = (d2 * (n * n) - d1 * d1.dot(d2)) / (n * n * n);
        (d1 + distance * unit_rate.perp()).length()
    };
    simpson(&speed, 0.0, 1.0, 1e-13)
}

/// The whole of an offset of the whole of `base`, no length stated.
fn open_offset(base: PathCurve, law: OffsetLaw) -> PathPiece {
    PathPiece::whole(PathCurve::Offset(Box::new(PathOffset::new(
        PathPiece::whole(base),
        law,
    ))))
}

/// Four pieces, each the offset of a whole curve, no length stated: beside
/// a Bezier by 1 (certified by chords), beside a line by 1 over 10, beside
/// a second Bezier by `middle` (`(1, 1)`: by 1, certified; `(1, 2)`: by
/// distances widening from 1 to 2 along the Bezier read as a chain's
/// parametric piece, whose length is not certified), beside a line by
/// `middle.1` over 10. The Beziers' ends run along `+x`, so the offsets
/// meet.
fn four_offsets(middle: (Scalar, Scalar)) -> CurvePath {
    let (_, first) = bezier();
    let (_, mut second) = bezier_at(Vec2::new(100.0, 20.0));
    let line = |x: Scalar, y: Scalar| line2(Point2::new(x, y), Vec2::new(10.0, 0.0));
    let by = |distance| OffsetLaw::Planar { distance };
    let law = if middle.0 == middle.1 {
        by(middle.0)
    } else {
        second = parametric_chain(Vec2::new(100.0, 20.0));
        OffsetLaw::Linear {
            start: PathOffsets::new(middle.0, 0.0, 0.0),
            end: PathOffsets::new(middle.1, 0.0, 0.0),
            frame: OffsetFrame::Section,
        }
    };
    CurvePath::new(vec![
        open_offset(PathCurve::Two(first), by(1.0)),
        open_offset(PathCurve::Two(line(90.0, 20.0)), by(1.0)),
        open_offset(PathCurve::Two(second), law),
        open_offset(PathCurve::Two(line(190.0, 40.0)), by(middle.1)),
    ])
}

/// The kernel's own (estimated) joint distances of `path`.
fn nominal_joints(path: &CurvePath) -> Vec<Scalar> {
    let basis = CompositeBasis::from_path(path).unwrap();
    let mut run = 0.0;
    let mut out = Vec::new();
    for piece in &basis.pieces()[..basis.pieces().len() - 1] {
        run += piece.length();
        out.push(run);
    }
    out
}

fn joints_of(found: &[(usize, DistanceInterval)]) -> Vec<usize> {
    found.iter().map(|(joint, _)| *joint).collect()
}

#[test]
fn an_uncertifiable_piece_after_the_window_leaves_the_joints_before_it_certified() {
    let path = four_offsets((1.0, 2.0));
    let provider = ReferenceCurveEvaluator::new();
    let tolerance = 1e-6;
    // Every joint at once: refused, since the widening offset's length is
    // not certified.
    refused(
        provider.path_joint_distances(&path, tolerance),
        UNCERTIFIED_LENGTH,
    );
    let (p, _) = bezier();
    let first = bezier_offset_length(p, 1.0);
    let truths = [first, first + 10.0];
    let nominal = nominal_joints(&path);
    // A window on the second joint, the one the widening offset starts at:
    // its interval sums a numerical length and an exact one, and the
    // widening offset's chord, 92 long, carries the next joint past it.
    let found = provider
        .path_joint_distances_within(&path, nominal[1], 0.5, tolerance)
        .unwrap();
    assert_eq!(joints_of(&found), [1]);
    certifies(found[0].1, truths[1], tolerance, "second joint");
    assert!(!found[0].1.is_point(), "the prefix's interval widens");
    // On the first joint, and over both.
    let found = provider
        .path_joint_distances_within(&path, nominal[0], 0.5, tolerance)
        .unwrap();
    assert_eq!(joints_of(&found), [0]);
    certifies(found[0].1, truths[0], tolerance, "first joint");
    let found = provider
        .path_joint_distances_within(&path, nominal[0] + 5.0, 6.0, tolerance)
        .unwrap();
    assert_eq!(joints_of(&found), [0, 1]);
    for (joint, interval) in &found {
        certifies(*interval, truths[*joint], tolerance, "both joints");
    }
    // Inside the widening offset, far from its ends: no joint can fall
    // there, so none is returned, and nothing is refused.
    assert_eq!(
        provider
            .path_joint_distances_within(&path, nominal[1] + 40.0, 1.0, tolerance)
            .unwrap(),
        vec![]
    );
    // Between the first two joints: none.
    assert_eq!(
        provider
            .path_joint_distances_within(&path, nominal[0] + 5.0, 1.0, tolerance)
            .unwrap(),
        vec![]
    );
    // Through the composite, the same answer.
    assert_eq!(
        CompositeBasis::from_path(&path)
            .unwrap()
            .joint_intervals_within(nominal[1], 0.5, tolerance)
            .unwrap(),
        provider
            .path_joint_distances_within(&path, nominal[1], 0.5, tolerance)
            .unwrap()
    );
}

#[test]
fn uncertainty_reaching_the_window_is_refused_by_name() {
    let path = four_offsets((1.0, 2.0));
    let provider = ReferenceCurveEvaluator::new();
    let nominal = nominal_joints(&path);
    // On the joint the widening offset ends at.
    refused(
        provider.path_joint_distances_within(&path, nominal[2], 0.5, 1e-6),
        UNCERTIFIED_LENGTH,
    );
    // Past it, on the last piece: the window's position against the joint
    // is what the widening offset's length leaves unknown.
    refused(
        provider.path_joint_distances_within(&path, nominal[2] + 5.0, 0.5, 1e-6),
        UNCERTIFIED_LENGTH,
    );
    // Before it, but within reach of its chord (92.2, of a length above
    // it): the joint could fall within the window.
    let chord = 90.0_f64.hypot(20.0);
    assert!(nominal[2] - nominal[1] > chord + 0.5);
    refused(
        provider.path_joint_distances_within(&path, nominal[1] + chord - 1.0, 1.5, 1e-6),
        UNCERTIFIED_LENGTH,
    );
}

#[test]
fn an_uncertifiable_piece_before_the_window_refuses_it_by_name() {
    // Reversed: the line by 2, the widening offset backwards, the line by
    // 1, the first offset backwards.
    let path = four_offsets((1.0, 2.0)).reversed();
    let provider = ReferenceCurveEvaluator::new();
    let nominal = nominal_joints(&path);
    assert_eq!(nominal[0], 10.0);
    // Every window after the widening offset is refused: no coarser bound
    // of its length is substituted.
    for (around, window) in [
        (nominal[1], 0.5),
        (nominal[1] + 5.0, 1.0),
        (nominal[2], 0.5),
    ] {
        refused(
            provider.path_joint_distances_within(&path, around, window, 1e-6),
            UNCERTIFIED_LENGTH,
        );
    }
    // Before it the joint is certified, and inside it none can fall.
    assert_eq!(
        provider
            .path_joint_distances_within(&path, 10.0, 0.5, 1e-6)
            .unwrap(),
        vec![(0, DistanceInterval::point(10.0))]
    );
    assert_eq!(
        provider
            .path_joint_distances_within(&path, 50.0, 1.0, 1e-6)
            .unwrap(),
        vec![]
    );
}

#[test]
fn joints_near_a_window_agree_with_every_joint() {
    let spiral = clothoid_then_arc();
    let line = line2(Point2::ZERO, Vec2::new(2.0, 0.0));
    let arc = Curve2::Circle(Circle2 {
        frame: frame(Point2::new(15.0, 5.0), Vec2::new(0.0, -1.0), Vec2::X),
        radius: 5.0,
    });
    let quarter = core::f64::consts::FRAC_PI_2 * 5.0;
    let up = line2(Point2::new(20.0, 5.0), Vec2::Y);
    let exact = CompositeBasis::new(
        offset_pieces(
            &[
                StationPiece::between(StationCurve::Two(&line), 0.0, 15.0).unwrap(),
                StationPiece::between(StationCurve::Two(&arc), 0.0, quarter).unwrap(),
                StationPiece::between(StationCurve::Two(&up), 0.0, 1.0).unwrap(),
            ],
            OffsetLaw::Planar { distance: 1.0 },
        )
        .unwrap(),
    )
    .unwrap()
    .path();
    let clothoid = CompositeBasis::new(
        offset_pieces(
            &[StationPiece::whole(StationCurve::Two(&spiral)).unwrap()],
            OffsetLaw::Planar { distance: 5.0 },
        )
        .unwrap(),
    )
    .unwrap()
    .path();
    let paths = [
        ("line, arc and line", exact),
        ("four offsets", four_offsets((1.0, 1.0))),
        ("four offsets reversed", four_offsets((1.0, 1.0)).reversed()),
        ("clothoid", clothoid),
    ];
    let provider = ReferenceCurveEvaluator::new();
    for (what, path) in &paths {
        let nominal = nominal_joints(path);
        let length = CompositeBasis::from_path(path).unwrap().length();
        for tolerance in [1e-4, 1e-6] {
            let all = provider.path_joint_distances(path, tolerance).unwrap();
            let mut windows = vec![(0.5 * length, 2.0 * length), (-5.0, 1.0)];
            for (k, joint) in nominal.iter().enumerate() {
                let next = nominal.get(k + 1).copied().unwrap_or(length);
                windows.extend([
                    (*joint, 0.1),
                    (joint - 0.3, 0.1),
                    (joint + 0.3, 0.1),
                    (0.5 * (joint + next), 0.25 * (next - joint)),
                    (0.5 * (joint + next), 0.5 * (next - joint) + 0.1),
                ]);
            }
            for (around, window) in windows {
                let span = DistanceInterval::new(around - window, around + window);
                let want: Vec<_> = all
                    .iter()
                    .copied()
                    .enumerate()
                    .filter(|(_, interval)| interval.meets(&span))
                    .collect();
                let found = provider
                    .path_joint_distances_within(path, around, window, tolerance)
                    .unwrap();
                assert_eq!(
                    found, want,
                    "{what} at {tolerance:e}, window {window} around {around}"
                );
            }
        }
    }
}

#[test]
fn a_window_that_is_not_finite_or_negative_is_refused() {
    let path = four_offsets((1.0, 2.0));
    let basis = CompositeBasis::from_path(&path).unwrap();
    for (around, window, tolerance) in [
        (Scalar::NAN, 1.0, 1e-6),
        (Scalar::NEG_INFINITY, 1.0, 1e-6),
        (10.0, -1.0, 1e-6),
        (10.0, Scalar::INFINITY, 1e-6),
        (10.0, 1.0, 0.0),
        (10.0, 1.0, Scalar::NAN),
    ] {
        assert!(
            matches!(
                basis.joint_intervals_within(around, window, tolerance),
                Err(GeomError::InvalidInput(_))
            ),
            "{around} {window} {tolerance}"
        );
    }
}

#[test]
fn a_plan_measured_piece_is_floored_by_its_plan_chord() {
    // At a grade of 100% the 3D chord of a piece is sqrt(2) times its plan
    // length, which is what a plan-measured path runs in: only the plan
    // chord bounds it from below.
    let steep = Curve3::Elevated(Elevated3::new(
        line2(Point2::ZERO, Vec2::X),
        ElevationLaw::Polynomial {
            coefficients: vec![0.0, 1.0],
        },
    ));
    let path = CurvePath::new(vec![
        PathPiece::new(PathCurve::Three(steep.clone()), 0.0, 10.0),
        PathPiece::new(PathCurve::Three(steep), 10.0, 20.0),
    ]);
    let provider = ReferenceCurveEvaluator::new();
    for (around, window) in [(10.0, 0.5), (11.0, 1.0), (12.0, 3.0)] {
        assert_eq!(
            provider
                .path_joint_distances_within(&path, around, window, 1e-9)
                .unwrap(),
            vec![(0, DistanceInterval::point(10.0))],
            "window {window} around {around}"
        );
    }
}
