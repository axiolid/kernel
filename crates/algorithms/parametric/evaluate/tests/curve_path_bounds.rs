//! Curve-path pieces bounded at a native parameter or a point (#304).
//!
//! A bounded piece resolves to the piece between those parameters, so a
//! path of them reads bitwise as the composite built from
//! `StationPiece::between_parameters`; its certified joints contain the
//! lengths of the curves written out here (Bernstein polynomials, an
//! ellipse's speed) and integrated by an adaptive Simpson rule to `1e-13`,
//! independent of the code under test.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionTarget, GeomError, GeomResult, Operation,
};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, CurvatureLaw, Curve2, Curve3, CurvePath, Ellipse2, Intrinsic2, KnotSpec, Line2,
    OffsetLaw, PathBound, PathCurve, PathOffset, PathPiece, Polyline2, SeamSide,
};
use axiolid_curve_evaluate_contract::{
    conformance, CurveEvaluator, CurveMeasure, DistanceConvention, DistanceInterval,
    PATH_PIECE_BOUND_UNSUPPORTED,
};
use axiolid_evaluate::station::{
    CompositeBasis, StationCurve, StationPiece, PATH_POINT_AMBIGUOUS, PATH_POINT_OFF_CURVE,
    POINT_BOUND_UNSUPPORTED,
};
use axiolid_evaluate::ReferenceCurveEvaluator;

const TAU: Scalar = core::f64::consts::TAU;

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
        "{what}: {interval:?} does not contain {truth}"
    );
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, input: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input: found, .. }) => assert_eq!(found, input),
        other => panic!("expected a refusal naming {input:?}, got {other:?}"),
    }
}

// --- curves -------------------------------------------------------------------

const BEZIER: [Point2; 4] = [
    Point2::new(0.0, 0.0),
    Point2::new(30.0, 0.0),
    Point2::new(60.0, 20.0),
    Point2::new(90.0, 20.0),
];

fn spline(points: &[Point2]) -> Curve2 {
    Curve2::BSpline(BSplineCurve {
        degree: 3,
        control_points: points.to_vec(),
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// The cubic Bezier through [`BEZIER`] and its derivative, written out.
fn bernstein(t: Scalar) -> Point2 {
    let p = BEZIER;
    let u = 1.0 - t;
    p[0] * (u * u * u) + p[1] * (3.0 * u * u * t) + p[2] * (3.0 * u * t * t) + p[3] * (t * t * t)
}

fn bernstein_speed(t: Scalar) -> Scalar {
    let p = BEZIER;
    let u = 1.0 - t;
    (3.0 * (u * u * (p[1] - p[0]) + 2.0 * u * t * (p[2] - p[1]) + t * t * (p[3] - p[2]))).length()
}

const SEMI: (Scalar, Scalar) = (20.0, 10.0);

fn ellipse(centre: Point2) -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: centre,
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: SEMI.0,
        semi_axis_y: SEMI.1,
    })
}

fn ellipse_point(centre: Point2, t: Scalar) -> Point2 {
    centre + Vec2::new(SEMI.0 * t.cos(), SEMI.1 * t.sin())
}

fn ellipse_speed(t: Scalar) -> Scalar {
    Vec2::new(SEMI.0 * t.sin(), SEMI.1 * t.cos()).length()
}

fn at(p: Point2) -> PathBound {
    PathBound::Point(Point3::new(p.x, p.y, 0.0))
}

fn parameter_of(curve: &Curve2, p: Point2) -> Scalar {
    StationCurve::Two(curve).parameter_at_bound(at(p)).unwrap()
}

/// Every reading of `bounded` through the contract equals, bitwise, the
/// station reading of `stated` at each of `distances`, on both sides; and
/// `bounded` resolves to `stated`'s own pieces.
fn reads_as_stated(bounded: &CurvePath, stated: &CompositeBasis<'_>, distances: &[Scalar]) {
    let resolved = CompositeBasis::from_path(bounded).unwrap();
    assert_eq!(resolved.path(), stated.path());
    assert_eq!(resolved.length().to_bits(), stated.length().to_bits());
    let provider = ReferenceCurveEvaluator::new();
    assert_eq!(
        provider.path_distance_convention(bounded),
        stated.convention()
    );
    for &d in distances {
        let at = CurveMeasure::Distance(d);
        for side in [SeamSide::Outgoing, SeamSide::Incoming] {
            let want = stated.section_on(d, side).unwrap();
            let point = provider.path_point_at_on(bounded, at, side).unwrap();
            let tangent = provider.path_tangent_at_on(bounded, at, side).unwrap();
            assert_eq!(point, want.point, "point at {d} {side:?}");
            assert_eq!(tangent, want.tangent, "tangent at {d} {side:?}");
        }
    }
    let past = stated.length() * (1.0 + 1e-9) + 1e-9;
    assert!(provider
        .path_point_at(bounded, CurveMeasure::Distance(past))
        .is_err());
}

// --- single pieces ------------------------------------------------------------

#[test]
fn an_ellipse_trimmed_at_parameters_and_at_points_needs_no_distance() {
    let centre = Point2::new(3.0, -2.0);
    let curve = ellipse(centre);
    let station = StationCurve::Two(&curve);
    let (t0, t1) = (0.3, 1.9);
    let by_parameters = CurvePath::from(PathPiece::between_parameters(
        PathCurve::Two(curve.clone()),
        t0,
        t1,
    ));
    // A closed conic's span runs round from `t0` by `t1 - t0`, as a
    // graph compiler's trim of it reads it.
    let stated = CompositeBasis::new(vec![StationPiece::between_parameters(
        station,
        t0,
        round(t0, t1),
    )
    .unwrap()])
    .unwrap();
    let length = stated.length();
    reads_as_stated(&by_parameters, &stated, &[0.0, 7.5, length / 2.0, length]);
    let truth = simpson(&ellipse_speed, t0, t1, 1e-13);
    assert!((length - truth).abs() <= 1e-9 * truth);
    // At points: the closed-form angle of each, then the same piece.
    let (p0, p1) = (ellipse_point(centre, t0), ellipse_point(centre, t1));
    let (u0, u1) = (parameter_of(&curve, p0), parameter_of(&curve, p1));
    assert!((u0 - t0).abs() <= 1e-12 && (u1 - t1).abs() <= 1e-12);
    let by_points = CurvePath::from(PathPiece::bounded(
        PathCurve::Two(curve.clone()),
        at(p0),
        at(p1),
    ));
    let stated_points = CompositeBasis::new(vec![StationPiece::between_parameters(
        station,
        u0,
        round(u0, u1),
    )
    .unwrap()])
    .unwrap();
    reads_as_stated(&by_points, &stated_points, &[0.0, 11.0, length]);
    // A reversed piece is the same span read the other way.
    let reversed = CurvePath::from(
        PathPiece::between_parameters(PathCurve::Two(curve.clone()), t0, t1).reversed(),
    );
    let stated_reversed = CompositeBasis::new(vec![StationPiece::between_parameters(
        station,
        t0,
        round(t0, t1),
    )
    .unwrap()
    .reversed()])
    .unwrap();
    reads_as_stated(&reversed, &stated_reversed, &[0.0, 7.5, length]);
    // Points across the seam run forwards from the start, through it.
    let (late, early) = (5.5, 0.7);
    let across = CurvePath::from(PathPiece::bounded(
        PathCurve::Two(curve.clone()),
        at(ellipse_point(centre, late)),
        at(ellipse_point(centre, early)),
    ));
    let span = CompositeBasis::from_path(&across).unwrap().length();
    let truth = simpson(&ellipse_speed, late, TAU + early, 1e-13);
    assert!(
        (span - truth).abs() <= 1e-9 * truth,
        "{span} against {truth}"
    );
    // An end on the seam ends the turn; a start on it starts it.
    let seam = at(ellipse_point(centre, 0.0));
    let to_seam = PathPiece::bounded(PathCurve::Two(curve.clone()), at(p1), seam);
    let from_seam = PathPiece::bounded(PathCurve::Two(curve.clone()), seam, at(p1));
    let whole = StationPiece::whole(station).unwrap();
    let to_seam = StationPiece::from_path_piece(&to_seam).unwrap();
    let from_seam = StationPiece::from_path_piece(&from_seam).unwrap();
    assert!((to_seam.end - whole.end).abs() <= 1e-12 * whole.end);
    assert_eq!(from_seam.start, 0.0);
    assert!((to_seam.length() + from_seam.length() - whole.length()).abs() <= 1e-9 * length);
    // Parameters a turn apart, or the same point twice: the whole turn.
    for piece in [
        PathPiece::between_parameters(PathCurve::Two(curve.clone()), 1.0, 1.0 + TAU),
        PathPiece::bounded(PathCurve::Two(curve.clone()), at(p1), at(p1)),
    ] {
        let turn = StationPiece::from_path_piece(&piece).unwrap();
        assert!((turn.length() - whole.length()).abs() <= 1e-12 * whole.length());
    }
    // Parameters given across the seam, or more than a turn apart, wrap.
    let wrapped = PathPiece::between_parameters(PathCurve::Two(curve.clone()), late, early);
    let wrapped = StationPiece::from_path_piece(&wrapped).unwrap();
    assert!((wrapped.length() - truth).abs() <= 1e-9 * truth);
}

/// Where a graph compiler's trim of a closed conic in its sense ends:
/// `t0` plus `t1 - t0`, wrapped into one turn when not within it.
fn round(t0: Scalar, t1: Scalar) -> Scalar {
    let travelled = t1 - t0;
    if travelled > 0.0 && travelled <= TAU {
        t0 + travelled
    } else {
        t0 + travelled.rem_euclid(TAU)
    }
}

#[test]
fn a_b_spline_trimmed_at_parameters_and_at_points_needs_no_distance() {
    let curve = spline(&BEZIER);
    let station = StationCurve::Two(&curve);
    let (t0, t1) = (0.2, 0.7);
    let by_parameters = CurvePath::from(PathPiece::between_parameters(
        PathCurve::Two(curve.clone()),
        t0,
        t1,
    ));
    let stated = CompositeBasis::new(vec![
        StationPiece::between_parameters(station, t0, t1).unwrap()
    ])
    .unwrap();
    let length = stated.length();
    reads_as_stated(&by_parameters, &stated, &[0.0, 20.0, length]);
    let truth = simpson(&bernstein_speed, t0, t1, 1e-13);
    assert!((length - truth).abs() <= 1e-9 * truth);
    // At points: Newton's parameter of each, within rounding of the true.
    let (p0, p1) = (bernstein(t0), bernstein(t1));
    let (u0, u1) = (parameter_of(&curve, p0), parameter_of(&curve, p1));
    assert!(
        (u0 - t0).abs() <= 1e-10 && (u1 - t1).abs() <= 1e-10,
        "{u0} {u1}"
    );
    let by_points = CurvePath::from(PathPiece::bounded(
        PathCurve::Two(curve.clone()),
        at(p0),
        at(p1),
    ));
    let stated_points =
        CompositeBasis::new(vec![
            StationPiece::between_parameters(station, u0, u1).unwrap()
        ])
        .unwrap();
    reads_as_stated(&by_points, &stated_points, &[0.0, 20.0, length]);
    // The domain's ends are points too.
    let ends = PathPiece::bounded(PathCurve::Two(curve.clone()), at(BEZIER[0]), at(BEZIER[3]));
    let whole = StationPiece::whole(station).unwrap();
    let ends = StationPiece::from_path_piece(&ends).unwrap();
    assert_eq!((ends.start, ends.end.to_bits()), (0.0, whole.end.to_bits()));
}

#[test]
fn one_bound_beside_a_stated_or_open_end_measures_as_between_parameters() {
    let curve = spline(&BEZIER);
    let station = StationCurve::Two(&curve);
    let between = StationPiece::between_parameters(station, 0.1, 0.6).unwrap();
    // A stated start and an end at a parameter.
    let end = PathPiece::new(PathCurve::Two(curve.clone()), between.start, 0.0)
        .with_end(PathBound::Parameter(0.6));
    let end = StationPiece::from_path_piece(&end).unwrap();
    assert_eq!(
        (end.start.to_bits(), end.end.to_bits()),
        (between.start.to_bits(), between.end.to_bits())
    );
    // A start at a parameter, the end open.
    let open =
        PathPiece::whole(PathCurve::Two(curve.clone())).with_start(PathBound::Parameter(0.1));
    let open = StationPiece::from_path_piece(&open).unwrap();
    let whole = StationPiece::whole(station).unwrap();
    assert_eq!(
        (open.start.to_bits(), open.end.to_bits()),
        (between.start.to_bits(), whole.end.to_bits())
    );
    // On an ellipse an end at a parameter follows the start round the seam.
    let centre = Point2::ZERO;
    let oval = ellipse(centre);
    let station = StationCurve::Two(&oval);
    let late = StationPiece::between_parameters(station, 5.0, TAU + 1.0).unwrap();
    let round = PathPiece::new(PathCurve::Two(oval.clone()), late.start, 0.0)
        .with_end(PathBound::Parameter(1.0));
    let round = StationPiece::from_path_piece(&round).unwrap();
    assert_eq!(round.start.to_bits(), late.start.to_bits());
    assert!((round.end - late.end).abs() <= 1e-12 * late.end);
    // A start at a parameter folds into the first turn.
    let folded =
        PathPiece::whole(PathCurve::Two(oval.clone())).with_start(PathBound::Parameter(TAU + 1.0));
    let folded = StationPiece::from_path_piece(&folded).unwrap();
    let first = StationPiece::between_parameters(station, 1.0, 2.0).unwrap();
    assert!((folded.start - first.start).abs() <= 1e-12 * first.start);
}

// --- a composite ----------------------------------------------------------------

/// A B-spline trimmed at two points, a straight line from where it ends,
/// and an ellipse trimmed at parameters from where that ends: no distance
/// stated anywhere.
struct Composite {
    spline: Curve2,
    line: Curve2,
    oval: Curve2,
    centre: Point2,
}

const ARC: (Scalar, Scalar) = (0.3, 1.9);
const STRAIGHT: Vec2 = Vec2::new(10.0, 5.0);

fn composite() -> Composite {
    let end = bernstein(0.75);
    let corner = end + STRAIGHT;
    let centre = corner - Vec2::new(SEMI.0 * ARC.0.cos(), SEMI.1 * ARC.0.sin());
    Composite {
        spline: spline(&BEZIER),
        line: Curve2::Line(Line2 {
            origin: end,
            direction: STRAIGHT,
        }),
        oval: ellipse(centre),
        centre,
    }
}

impl Composite {
    fn path(&self) -> CurvePath {
        CurvePath::new(vec![
            PathPiece::bounded(
                PathCurve::Two(self.spline.clone()),
                at(bernstein(0.25)),
                at(bernstein(0.75)),
            ),
            PathPiece::whole(PathCurve::Two(self.line.clone())),
            PathPiece::between_parameters(PathCurve::Two(self.oval.clone()), ARC.0, ARC.1),
        ])
    }

    fn stated(&self) -> Vec<StationPiece<'_>> {
        let spline = StationCurve::Two(&self.spline);
        let (u0, u1) = (
            parameter_of(&self.spline, bernstein(0.25)),
            parameter_of(&self.spline, bernstein(0.75)),
        );
        vec![
            StationPiece::between_parameters(spline, u0, u1).unwrap(),
            StationPiece::whole(StationCurve::Two(&self.line)).unwrap(),
            StationPiece::between_parameters(
                StationCurve::Two(&self.oval),
                ARC.0,
                round(ARC.0, ARC.1),
            )
            .unwrap(),
        ]
    }
}

#[test]
fn a_composite_of_bounded_pieces_reads_as_its_stated_composite_and_certifies() {
    let c = composite();
    let path = c.path();
    assert!(path.pieces().iter().all(|piece| !piece.is_stated()));
    let stated = CompositeBasis::new(c.stated()).unwrap();
    let first = stated.pieces()[0].length();
    let second = first + STRAIGHT.length();
    reads_as_stated(
        &path,
        &stated,
        &[
            0.0,
            12.5,
            first,
            first * (1.0 + 1e-13),
            first + 4.0,
            second,
            second + 9.0,
            stated.length(),
        ],
    );
    let spline_truth = simpson(&bernstein_speed, 0.25, 0.75, 1e-13);
    let truths = [spline_truth, spline_truth + STRAIGHT.length()];
    let provider = ReferenceCurveEvaluator::new();
    for tolerance in [1e-5, 1e-7] {
        let joints = provider.path_joint_distances(&path, tolerance).unwrap();
        assert_eq!(joints.len(), 2);
        for (joint, truth) in joints.iter().zip(truths) {
            certifies(*joint, truth, tolerance, "bounded composite");
        }
        // The joints near a window (#305) are the same certificates.
        for (k, truth) in truths.iter().enumerate() {
            let near = provider
                .path_joint_distances_within(&path, *truth, 0.5, tolerance)
                .unwrap();
            assert_eq!(near, vec![(k, joints[k])], "joint {k} at {tolerance:e}");
        }
    }
    // The ellipse ends where its closed form puts it.
    let end = provider
        .path_point_at(&path, CurveMeasure::Distance(stated.length()))
        .unwrap();
    let want = ellipse_point(c.centre, ARC.1);
    assert!((Point2::new(end.x, end.y) - want).length() <= 1e-9);
}

#[test]
fn a_reversed_composite_of_bounded_pieces_reads_and_certifies() {
    let c = composite();
    let path = c.path().reversed();
    let stated = CompositeBasis::new(
        c.stated()
            .into_iter()
            .rev()
            .map(StationPiece::reversed)
            .collect(),
    )
    .unwrap();
    let first = stated.pieces()[0].length();
    reads_as_stated(
        &path,
        &stated,
        &[0.0, 5.0, first, first + 3.0, stated.length()],
    );
    let arc = simpson(&ellipse_speed, ARC.0, ARC.1, 1e-13);
    let truths = [arc, arc + STRAIGHT.length()];
    for tolerance in [1e-5, 1e-7] {
        let joints = ReferenceCurveEvaluator::new()
            .path_joint_distances(&path, tolerance)
            .unwrap();
        for (joint, truth) in joints.iter().zip(truths) {
            certifies(*joint, truth, tolerance, "reversed bounded composite");
        }
    }
    // The reversed path starts where the ellipse's span ends.
    let start = ReferenceCurveEvaluator::new()
        .path_point_at(&path, CurveMeasure::Distance(0.0))
        .unwrap();
    let want = ellipse_point(c.centre, ARC.1);
    assert!((Point2::new(start.x, start.y) - want).length() <= 1e-9);
}

// --- refusals -------------------------------------------------------------------

fn resolve(piece: PathPiece) -> GeomResult<StationPiece<'static>> {
    let piece: &'static PathPiece = Box::leak(Box::new(piece));
    StationPiece::from_path_piece(piece)
}

#[test]
fn a_point_off_the_curve_is_refused_by_name() {
    let centre = Point2::new(3.0, -2.0);
    let oval = PathCurve::Two(ellipse(centre));
    let on = ellipse_point(centre, 1.0);
    // Off along the normal by a micrometre, and lifted off the plane.
    let normal = Vec2::new(SEMI.1 * 1.0_f64.cos(), SEMI.0 * 1.0_f64.sin()).normalize();
    refused(
        resolve(PathPiece::bounded(
            oval.clone(),
            at(on + normal * 1e-6),
            at(ellipse_point(centre, 2.0)),
        )),
        PATH_POINT_OFF_CURVE,
    );
    let lifted = PathBound::Point(Point3::new(on.x, on.y, 1e-6));
    refused(
        resolve(PathPiece::whole(oval.clone()).with_start(lifted)),
        PATH_POINT_OFF_CURVE,
    );
    // Within the joint tolerance (about 2.3e-8 here) it is on the curve.
    assert!(resolve(PathPiece::whole(oval).with_start(at(on + normal * 1e-8))).is_ok());
    let curve = PathCurve::Two(spline(&BEZIER));
    let beside = bernstein(0.4) + Vec2::new(0.0, 1e-5);
    refused(
        resolve(PathPiece::whole(curve.clone()).with_start(at(beside))),
        PATH_POINT_OFF_CURVE,
    );
    refused(
        resolve(PathPiece::whole(curve).with_start(at(Point2::new(200.0, 0.0)))),
        PATH_POINT_OFF_CURVE,
    );
}

#[test]
fn a_point_several_parameters_name_is_refused_by_name() {
    // A closed polyline's seam: its first vertex is parameter 0 and 4.
    let square = PathCurve::Two(Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        closed: true,
    }));
    refused(
        resolve(PathPiece::bounded(
            square.clone(),
            at(Point2::ZERO),
            at(Point2::new(1.0, 1.0)),
        )),
        PATH_POINT_AMBIGUOUS,
    );
    // An interior vertex is one place, though two segments reach it.
    let side = resolve(PathPiece::bounded(
        square.clone(),
        at(Point2::new(1.0, 0.0)),
        at(Point2::new(1.0, 1.0)),
    ))
    .unwrap();
    assert_eq!((side.start, side.end), (1.0, 2.0));
    let corner = resolve(PathPiece::bounded(
        square,
        at(Point2::new(0.5, 0.0)),
        at(Point2::new(1.0, 1.0 - 1e-12)),
    ))
    .unwrap();
    assert!((corner.length() - 1.5).abs() <= 1e-11);
    // A polyline crossing itself, at the crossing.
    let crossing = PathCurve::Two(Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 2.0),
            Point2::new(2.0, 0.0),
            Point2::new(0.0, 2.0),
        ],
        closed: false,
    }));
    refused(
        resolve(PathPiece::whole(crossing).with_start(at(Point2::new(1.0, 1.0)))),
        PATH_POINT_AMBIGUOUS,
    );
    // A B-spline loop: its start and its end are one point.
    let loop_points = [
        Point2::new(0.0, 0.0),
        Point2::new(40.0, 30.0),
        Point2::new(-40.0, 30.0),
        Point2::new(0.0, 0.0),
    ];
    let looped = PathCurve::Two(spline(&loop_points));
    refused(
        resolve(PathPiece::whole(looped.clone()).with_start(at(Point2::ZERO))),
        PATH_POINT_AMBIGUOUS,
    );
    let top = resolve(PathPiece::whole(looped).with_start(at(Point2::new(0.0, 22.5)))).unwrap();
    assert!(top.start > 0.0);
}

#[test]
fn bounds_a_curve_cannot_take_are_refused_by_name() {
    let spiral = PathCurve::Two(Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 100.0),
        100.0,
    )));
    refused(
        resolve(PathPiece::whole(spiral.clone()).with_start(at(Point2::ZERO))),
        POINT_BOUND_UNSUPPORTED,
    );
    // A curve parameterised by its measure takes a parameter as a distance.
    let tail = resolve(PathPiece::whole(spiral).with_start(PathBound::Parameter(40.0))).unwrap();
    assert_eq!((tail.start, tail.end), (40.0, 100.0));
    // An offset takes no bound of its own; its base piece may.
    let base = PathPiece::between_parameters(PathCurve::Two(spline(&BEZIER)), 0.2, 0.7);
    let offset = PathCurve::Offset(Box::new(PathOffset::new(
        base,
        OffsetLaw::Planar { distance: 1.0 },
    )));
    assert!(resolve(PathPiece::whole(offset.clone())).is_ok());
    refused(
        resolve(PathPiece::whole(offset.clone()).with_start(at(Point2::ZERO))),
        POINT_BOUND_UNSUPPORTED,
    );
    assert!(resolve(PathPiece::whole(offset).with_start(PathBound::Parameter(0.5))).is_err());
}

#[test]
fn malformed_bounds_are_refused() {
    let curve = PathCurve::Two(spline(&BEZIER));
    let station = StationCurve::Two(match &curve {
        PathCurve::Two(curve) => curve,
        _ => unreachable!(),
    });
    for (t0, t1) in [
        (0.7, 0.2),
        (0.4, 0.4),
        (Scalar::NAN, 0.5),
        (0.1, Scalar::INFINITY),
        (-0.5, 0.5),
    ] {
        assert!(
            resolve(PathPiece::between_parameters(curve.clone(), t0, t1)).is_err(),
            "[{t0}, {t1}]"
        );
    }
    let point = PathBound::Point(Point3::new(Scalar::NAN, 0.0, 0.0));
    assert!(matches!(
        station.parameter_at_bound(point),
        Err(GeomError::InvalidInput(_))
    ));
    // An end both open and bounded.
    let mut both = PathPiece::whole(curve.clone()).with_end(PathBound::Parameter(0.5));
    both.open_end = true;
    assert!(matches!(resolve(both), Err(GeomError::InvalidInput(_))));
    // The same point twice spans nothing on a curve that does not run
    // round.
    let p = at(bernstein(0.5));
    assert!(resolve(PathPiece::bounded(curve, p, p)).is_err());
}

// --- the conformance suite ------------------------------------------------------

/// A provider that mishandles bounds (#304): `Bounds::Shifted` reads a
/// bounded path a millimetre further along (its joints too),
/// `Bounds::Invalid` refuses a bound as bad input, `Bounds::Projected`
/// reads a point it should refuse at the curve's start instead,
/// `Bounds::Unsupported` refuses bounds typed, by name (which is
/// conformant).
#[derive(Debug)]
struct Bounded(ReferenceCurveEvaluator, Bounds);

#[derive(Debug, Clone, Copy)]
enum Bounds {
    Shifted,
    Invalid,
    Projected,
    Unsupported,
}

fn has_bounds(path: &CurvePath) -> bool {
    path.pieces()
        .iter()
        .any(|piece| piece.start_bound.is_some() || piece.end_bound.is_some())
}

impl Bounded {
    /// The path this provider actually reads.
    fn restated(&self, path: &CurvePath) -> GeomResult<CurvePath> {
        if !has_bounds(path) {
            return Ok(path.clone());
        }
        match self.1 {
            Bounds::Invalid => Err(GeomError::InvalidInput("bound".into())),
            Bounds::Unsupported => Err(GeomError::UnsupportedInput {
                backend: self.descriptor().id,
                operation: Operation::CurveEvaluation,
                input: PATH_PIECE_BOUND_UNSUPPORTED,
            }),
            Bounds::Shifted => Ok(path.clone()),
            Bounds::Projected => Ok(path
                .pieces()
                .iter()
                .map(|piece| {
                    let mut piece = piece.clone();
                    if StationPiece::from_path_piece(&piece.clone()).is_err() {
                        piece.start_bound = Some(PathBound::Parameter(0.0));
                        piece.end_bound = Some(PathBound::Parameter(1.0));
                    }
                    piece
                })
                .collect()),
        }
    }
}

impl Backend for Bounded {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new("bounded"), ExecutionTarget::PortableCpu)
    }
}

impl CurveEvaluator for Bounded {
    fn distance_convention(&self, curve: &Curve3) -> DistanceConvention {
        self.0.distance_convention(curve)
    }
    fn point_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Point3> {
        self.0.point_at(curve, at)
    }
    fn tangent_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Vec3> {
        self.0.tangent_at(curve, at)
    }
    fn frame_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Frame3> {
        self.0.frame_at(curve, at)
    }
    fn path_distance_convention(&self, path: &CurvePath) -> DistanceConvention {
        self.restated(path)
            .map_or(DistanceConvention::Unsupported, |path| {
                self.0.path_distance_convention(&path)
            })
    }
    fn path_point_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        side: SeamSide,
    ) -> GeomResult<Point3> {
        let at = match (self.1, at) {
            (Bounds::Shifted, CurveMeasure::Distance(d)) if has_bounds(path) => {
                CurveMeasure::Distance(d + 1e-3)
            }
            _ => at,
        };
        self.0.path_point_at_on(&self.restated(path)?, at, side)
    }
    fn path_tangent_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        side: SeamSide,
    ) -> GeomResult<Vec3> {
        self.0.path_tangent_at_on(&self.restated(path)?, at, side)
    }
    fn path_frame_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        side: SeamSide,
    ) -> GeomResult<Frame3> {
        self.0.path_frame_at_on(&self.restated(path)?, at, side)
    }
    fn path_joint_distances(
        &self,
        path: &CurvePath,
        tolerance: Scalar,
    ) -> GeomResult<Vec<DistanceInterval>> {
        let joints = self
            .0
            .path_joint_distances(&self.restated(path)?, tolerance)?;
        Ok(match self.1 {
            Bounds::Shifted if has_bounds(path) => joints
                .iter()
                .map(|j| DistanceInterval::new(j.lower + 1e-3, j.upper + 1e-3))
                .collect(),
            _ => joints,
        })
    }
}

#[test]
fn the_conformance_suite_reads_bounds() {
    assert_eq!(
        conformance::check(&ReferenceCurveEvaluator::new()),
        Vec::new()
    );
    for (how, checks) in [
        (
            Bounds::Shifted,
            &[
                "a bound resolves to the distance of its parameter or point",
                "joint distances of a bounded path are certified",
            ][..],
        ),
        (
            Bounds::Invalid,
            &[
                "a bounded path is answered or refused as unsupported",
                "joint distances of a bounded path are certified or refused as unsupported",
            ][..],
        ),
        (
            Bounds::Projected,
            &["a point bound off the curve or named by several parameters is refused"][..],
        ),
    ] {
        let failures = conformance::check(&Bounded(ReferenceCurveEvaluator::new(), how));
        for check in checks {
            assert!(
                failures.iter().any(|f| f.check == *check),
                "{how:?}, {check}: {failures:?}"
            );
        }
    }
    // Refusing bounds typed, by name, is conformant.
    assert_eq!(
        conformance::check(&Bounded(
            ReferenceCurveEvaluator::new(),
            Bounds::Unsupported
        )),
        Vec::new()
    );
}
