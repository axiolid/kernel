//! Point, tangent and frame queries over a curve path on the
//! curve-evaluation contract (#290).
//!
//! A path is pieces of curves laid end to end (`axiolid_curve::CurvePath`),
//! the neutral form of a curve relation. Through the contract the
//! reference evaluator must read it exactly as `CompositeBasis::section_on`
//! reads the same pieces with the same side (#285): bitwise, since it
//! reads them by that code; and both must match closed forms at each
//! joint, on a line and arc, a trim of it, and segments placed at stations
//! of it. A provider that does not implement paths refuses them by name.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionTarget, GeomError, GeomResult, Operation,
};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Circle2, Curve2, Curve3, CurvePath, Elevated3,
    ElevationLaw, Line2, Line3, PathCurve, PathPiece, SeamSide,
};
use axiolid_curve_evaluate_contract::{
    conformance, CurveEvaluator, CurveMeasure, DistanceConvention, CURVE_PATH_UNSUPPORTED,
};
use axiolid_evaluate::station::{CompositeBasis, StationCurve, StationPiece};
use axiolid_evaluate::ReferenceCurveEvaluator;

const EPS: Scalar = 1e-9;
const PI: Scalar = core::f64::consts::PI;
const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];

fn close(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    let off = (actual - expected).length();
    assert!(
        off <= eps,
        "{what}: {actual:?} != {expected:?} (off by {off:e})"
    );
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, needle: &str) {
    match result {
        Err(error) => assert!(
            error.to_string().contains(needle),
            "expected a refusal naming {needle:?}, got {error}"
        ),
        Ok(value) => panic!("expected a refusal naming {needle:?}, got {value:?}"),
    }
}

/// `+x` from the origin.
fn line() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })
}

/// Radius 5 about `(15, 0)`: at angle `pi` it passes `(10, 0)` heading
/// `-y`.
fn circle() -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(15.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 5.0,
    })
}

/// The line over `[0, 10]`, then the quarter arc from angle `pi` to
/// `3 pi / 2`: a right-angle corner at distance 10 (incoming `+x`,
/// outgoing `-y`), length `10 + 5 pi / 2`.
fn line_and_arc(line: &Curve2, circle: &Curve2) -> CurvePath {
    CurvePath::new(vec![
        PathPiece::new(PathCurve::Two(line.clone()), 0.0, 10.0),
        StationPiece::between_parameters(StationCurve::Two(circle), PI, 1.5 * PI)
            .unwrap()
            .into(),
    ])
}

/// `d` and the distances a hair either side of it, within the seam
/// tolerance `1e-12 * max(1, d)`.
fn around(d: Scalar) -> [Scalar; 3] {
    let hair = 0.4e-12 * d.max(1.0);
    [d - hair, d, d + hair]
}

/// Every query, plain and sided, equals the composite's own reading of
/// the same pieces with the same side, bitwise.
fn agrees_with_the_composite(path: &CurvePath, distances: &[Scalar], what: &str) {
    let e = ReferenceCurveEvaluator::new();
    let composite = CompositeBasis::from_path(path).unwrap();
    assert_eq!(&composite.path(), path, "{what}: the path round-trips");
    for &d in distances {
        let at = CurveMeasure::Distance(d);
        for side in SIDES {
            let station = composite.section_on(d, side).unwrap();
            assert_eq!(
                e.path_point_at_on(path, at, side).unwrap(),
                station.point,
                "{what} {side:?} at {d}: point"
            );
            assert_eq!(
                e.path_tangent_at_on(path, at, side).unwrap(),
                station.tangent,
                "{what} {side:?} at {d}: tangent"
            );
            assert_eq!(
                e.path_frame_at_on(path, at, side).unwrap(),
                station.frame(),
                "{what} {side:?} at {d}: frame"
            );
            assert_eq!(
                e.path_frame_is_exact_at(path, at, side),
                composite.frame_is_exact_at(d, side),
                "{what} {side:?} at {d}: exactness"
            );
        }
        let outgoing = composite.section_on(d, SeamSide::Outgoing).unwrap();
        assert_eq!(e.path_point_at(path, at).unwrap(), outgoing.point);
        assert_eq!(e.path_tangent_at(path, at).unwrap(), outgoing.tangent);
        assert_eq!(e.path_frame_at(path, at).unwrap(), outgoing.frame());
    }
}

/// Each side of a joint at `d` reads its own piece: the closed-form
/// point and tangents there, a hair either side included, framed `x`
/// tangent, `y` up (`+Z`), `z` right.
fn joint_reads(path: &CurvePath, d: Scalar, point: Point3, incoming: Vec3, outgoing: Vec3) {
    let e = ReferenceCurveEvaluator::new();
    for s in around(d) {
        let at = CurveMeasure::Distance(s);
        for (side, want) in [
            (SeamSide::Incoming, incoming),
            (SeamSide::Outgoing, outgoing),
        ] {
            let what = format!("{side:?} at {s}");
            close(
                e.path_point_at_on(path, at, side).unwrap() - Point3::ZERO,
                point - Point3::ZERO,
                EPS,
                &format!("{what}: point"),
            );
            close(
                e.path_tangent_at_on(path, at, side).unwrap(),
                want,
                EPS,
                &format!("{what}: tangent"),
            );
            let frame = e.path_frame_at_on(path, at, side).unwrap();
            close(frame.x, want, EPS, &format!("{what}: frame x"));
            close(frame.y, Vec3::Z, EPS, &format!("{what}: frame y"));
            close(
                frame.z,
                want.cross(Vec3::Z),
                EPS,
                &format!("{what}: frame z"),
            );
        }
    }
}

#[test]
fn a_line_and_arc_path_reads_as_its_composite_on_both_sides_of_the_joint() {
    let (line, circle) = (line(), circle());
    let path = line_and_arc(&line, &circle);
    let length = 10.0 + 2.5 * PI;
    let mut distances = vec![0.0, 4.0, 12.0, length];
    distances.extend(around(10.0));
    agrees_with_the_composite(&path, &distances, "line and arc");
    joint_reads(&path, 10.0, Point3::new(10.0, 0.0, 0.0), Vec3::X, -Vec3::Y);

    let e = ReferenceCurveEvaluator::new();
    assert_eq!(
        e.path_distance_convention(&path),
        DistanceConvention::ArcLength3d
    );
    // Inside the arc, the closed form at angle pi + 2 / 5.
    let (s, c) = (PI + 0.4).sin_cos();
    let at = CurveMeasure::Distance(12.0);
    close(
        e.path_point_at(&path, at).unwrap() - Point3::ZERO,
        Vec3::new(15.0 + 5.0 * c, 5.0 * s, 0.0),
        EPS,
        "arc point",
    );
    close(
        e.path_tangent_at(&path, at).unwrap(),
        Vec3::new(-s, c, 0.0),
        EPS,
        "arc tangent",
    );
    // Exact on the line, not on or after the arc.
    let exact = |d: Scalar, side| e.path_frame_is_exact_at(&path, CurveMeasure::Distance(d), side);
    assert!(exact(4.0, SeamSide::Outgoing) && exact(10.0, SeamSide::Incoming));
    assert!(!exact(10.0, SeamSide::Outgoing) && !exact(12.0, SeamSide::Incoming));
}

#[test]
fn a_trim_reads_as_its_composite_forwards_and_reversed() {
    let (line, circle) = (line(), circle());
    let base = line_and_arc(&line, &circle);
    let full = CompositeBasis::from_path(&base).unwrap();
    // Trimmed to [3, 12]: the joint moves to 7.
    let trimmed: CurvePath = full
        .pieces_between(3.0, 12.0)
        .unwrap()
        .into_iter()
        .map(PathPiece::from)
        .collect();
    let mut distances = vec![0.0, 2.0, 8.0, 9.0];
    distances.extend(around(7.0));
    agrees_with_the_composite(&trimmed, &distances, "trim");
    joint_reads(
        &trimmed,
        7.0,
        Point3::new(10.0, 0.0, 0.0),
        Vec3::X,
        -Vec3::Y,
    );
    let e = ReferenceCurveEvaluator::new();
    close(
        e.path_point_at(&trimmed, CurveMeasure::Distance(0.0))
            .unwrap()
            - Point3::ZERO,
        Vec3::new(3.0, 0.0, 0.0),
        EPS,
        "the trim starts at 3",
    );

    // Reversed: the arc back to the joint at 2, then the line back to 3.
    let reversed = trimmed.reversed();
    let mut distances = vec![0.0, 1.0, 5.0, 9.0];
    distances.extend(around(2.0));
    agrees_with_the_composite(&reversed, &distances, "reversed trim");
    joint_reads(
        &reversed,
        2.0,
        Point3::new(10.0, 0.0, 0.0),
        Vec3::Y,
        -Vec3::X,
    );
    close(
        e.path_point_at(&reversed, CurveMeasure::Distance(9.0))
            .unwrap()
            - Point3::ZERO,
        Vec3::new(3.0, 0.0, 0.0),
        EPS,
        "the reversed trim ends at 3",
    );
}

#[test]
fn segments_placed_at_stations_of_a_composite_read_as_the_composite() {
    let (line, circle) = (line(), circle());
    let base_path = line_and_arc(&line, &circle);
    let base = CompositeBasis::from_path(&base_path).unwrap();
    // A local line over [0, 6] at station 4, and a local quarter arc
    // (radius 5 about (0, 5), from -pi / 2 to 0) at the joint, outgoing:
    // together they retrace the base from distance 4.
    let local_line = line.clone();
    let local_circle = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(0.0, 5.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 5.0,
    });
    let place = |d: Scalar| {
        (
            base.section_on(d, SeamSide::Outgoing).unwrap().placement(),
            base.frame_is_exact_at(d, SeamSide::Outgoing),
        )
    };
    let (at_four, four_exact) = place(4.0);
    let (at_joint, joint_exact) = place(10.0);
    assert!(four_exact && !joint_exact);
    let straight =
        PathPiece::new(PathCurve::Two(local_line.clone()), 0.0, 6.0).placed(at_four, four_exact);
    let bend = PathPiece::from(
        StationPiece::between_parameters(StationCurve::Two(&local_circle), -0.5 * PI, 0.0).unwrap(),
    )
    .placed(at_joint, joint_exact);
    let segmented = CurvePath::new(vec![straight, bend]);

    let mut distances = vec![1.0, 3.0, 7.5, 6.0 + 2.5 * PI];
    distances.extend(around(6.0));
    agrees_with_the_composite(&segmented, &distances, "placed segments");
    joint_reads(
        &segmented,
        6.0,
        Point3::new(10.0, 0.0, 0.0),
        Vec3::X,
        -Vec3::Y,
    );

    // Along the base, by a different code path: within rounding.
    let e = ReferenceCurveEvaluator::new();
    for d in [1.0, 3.0, 6.0, 7.5, 6.0 + 2.5 * PI] {
        for side in SIDES {
            let on_segments = e
                .path_frame_at_on(&segmented, CurveMeasure::Distance(d), side)
                .unwrap();
            let on_base = e
                .path_frame_at_on(&base_path, CurveMeasure::Distance(d + 4.0), side)
                .unwrap();
            let what = format!("{side:?} at {d}");
            close(
                on_segments.origin - on_base.origin,
                Vec3::ZERO,
                EPS,
                &format!("{what}: origin"),
            );
            close(on_segments.x, on_base.x, EPS, &format!("{what}: x"));
            close(on_segments.y, on_base.y, EPS, &format!("{what}: y"));
        }
    }
    // Exact on the straight segment (a line placed on the base's line),
    // not on the bend nor beyond it.
    let exact =
        |d: Scalar, side| e.path_frame_is_exact_at(&segmented, CurveMeasure::Distance(d), side);
    assert!(exact(3.0, SeamSide::Outgoing) && exact(6.0, SeamSide::Incoming));
    assert!(!exact(6.0, SeamSide::Outgoing) && !exact(8.0, SeamSide::Outgoing));
    // A line after the arc is not exact either: the arc's length makes its
    // distance.
    let line3 = Curve3::Line(Line3 {
        origin: Point3::new(15.0, -5.0, 0.0),
        direction: Vec3::X,
    });
    let mut after = line_and_arc(&line, &circle);
    after.extend([PathPiece::new(PathCurve::Three(line3.clone()), 0.0, 3.0)]);
    let tail = 10.0 + 2.5 * PI + 1.0;
    assert!(!e.path_frame_is_exact_at(&after, CurveMeasure::Distance(tail), SeamSide::Outgoing));
    close(
        e.path_tangent_at(&after, CurveMeasure::Distance(tail))
            .unwrap(),
        Vec3::X,
        EPS,
        "the line after the arc",
    );
}

#[test]
fn a_path_query_refuses_what_its_composite_refuses() {
    let (line, circle) = (line(), circle());
    let path = line_and_arc(&line, &circle);
    let e = ReferenceCurveEvaluator::new();
    for side in SIDES {
        refused(
            e.path_point_at_on(&path, CurveMeasure::Parameter(1.0), side),
            "native parameter along a curve path",
        );
        refused(
            e.path_frame_at_on(&path, CurveMeasure::Distance(Scalar::NAN), side),
            "finite",
        );
        refused(
            e.path_tangent_at_on(&path, CurveMeasure::Distance(-1.0), side),
            "before the curve's start",
        );
        refused(
            e.path_point_at_on(&path, CurveMeasure::Distance(10.0 + 2.5 * PI + 1e-3), side),
            "beyond the curve's length",
        );
        assert!(!e.path_frame_is_exact_at(&path, CurveMeasure::Parameter(1.0), side));
    }
    // Pieces that do not meet.
    let gapped = CurvePath::new(vec![
        PathPiece::new(PathCurve::Two(line.clone()), 0.0, 9.0),
        path.pieces()[1].clone(),
    ]);
    refused(e.path_point_at(&gapped, CurveMeasure::Distance(1.0)), "gap");
    assert_eq!(
        e.path_distance_convention(&gapped),
        DistanceConvention::Unsupported
    );
    assert!(!e.path_frame_is_exact_at(&gapped, CurveMeasure::Distance(1.0), SeamSide::Outgoing));
    // No pieces, and pieces measured differently.
    refused(
        e.path_point_at(&CurvePath::default(), CurveMeasure::Distance(0.0)),
        "no pieces",
    );
    let elevated = Curve3::Elevated(Elevated3::new(line.clone(), ElevationLaw::level(0.0)));
    let mixed = CurvePath::new(vec![
        PathPiece::new(PathCurve::Three(elevated.clone()), 0.0, 10.0),
        path.pieces()[1].clone(),
    ]);
    refused(
        e.path_frame_at(&mixed, CurveMeasure::Distance(1.0)),
        "measure stations differently",
    );
    let plan = CurvePath::from(PathPiece::new(
        PathCurve::Three(elevated.clone()),
        0.0,
        10.0,
    ));
    assert_eq!(
        e.path_distance_convention(&plan),
        DistanceConvention::PlanDistance
    );
    // A span outside its curve.
    let outside = CurvePath::from(PathPiece::new(PathCurve::Two(circle.clone()), 0.0, 40.0));
    refused(
        e.path_point_at(&outside, CurveMeasure::Distance(1.0)),
        "outside its curve's measure",
    );
}

/// Level and straight: 150 mm cant for 100 m, then none.
fn cant_break() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(line(), ElevationLaw::level(0.0)),
        CantLaw::new(vec![
            CantPiece::constant(100.0, 0.15),
            CantPiece::constant(50.0, 0.0),
        ]),
        CantLaw::zero(150.0),
        1.5,
        BankConvention::TangentRotation,
    ))
}

#[test]
fn another_reference_up_frames_a_path_against_it_and_refuses_a_banked_piece() {
    let (line, circle) = (line(), circle());
    let path = line_and_arc(&line, &circle);
    let up = Vec3::new(0.0, 0.6, 0.8);
    let e = ReferenceCurveEvaluator::with_up(up).unwrap();
    for d in around(10.0) {
        for side in SIDES {
            let at = CurveMeasure::Distance(d);
            let tangent = e.path_tangent_at_on(&path, at, side).unwrap();
            let frame = e.path_frame_at_on(&path, at, side).unwrap();
            let right = tangent.cross(up).normalize();
            close(frame.x, tangent, 0.0, "x is the tangent");
            close(
                frame.z,
                right,
                1e-15,
                "z is right of the tangent against up",
            );
            close(
                frame.y,
                right.cross(tangent),
                1e-15,
                "y is up made perpendicular",
            );
        }
    }
    let banked = cant_break();
    let canted = CurvePath::from(PathPiece::new(PathCurve::Three(banked.clone()), 0.0, 150.0));
    refused(
        e.path_frame_at(&canted, CurveMeasure::Distance(100.0)),
        "+Z",
    );
    // Against +Z each side of the cant break reads its own roll, as the
    // composite does.
    let z = ReferenceCurveEvaluator::new();
    agrees_with_the_composite(&canted, &around(100.0), "cant break");
    let incoming = z
        .path_frame_at_on(&canted, CurveMeasure::Distance(100.0), SeamSide::Incoming)
        .unwrap();
    let outgoing = z
        .path_frame_at_on(&canted, CurveMeasure::Distance(100.0), SeamSide::Outgoing)
        .unwrap();
    assert!((incoming.y - outgoing.y).length() > 1e-3);
    close(outgoing.y, Vec3::Z, EPS, "no cant after the break");
}

/// A provider implementing single curves only.
#[derive(Debug)]
struct SingleCurves(ReferenceCurveEvaluator);

impl Backend for SingleCurves {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("single-curves"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl CurveEvaluator for SingleCurves {
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
}

fn refused_by_name<T: core::fmt::Debug>(result: GeomResult<T>) {
    match result {
        Err(GeomError::UnsupportedInput {
            backend,
            operation,
            input,
        }) => {
            assert_eq!(backend, BackendId::new("single-curves"));
            assert_eq!(operation, Operation::CurveEvaluation);
            assert_eq!(input, CURVE_PATH_UNSUPPORTED);
        }
        other => panic!("expected a refusal naming the path, got {other:?}"),
    }
}

#[test]
fn a_provider_without_paths_refuses_them_by_name() {
    let (line, circle) = (line(), circle());
    let path = line_and_arc(&line, &circle);
    let p = SingleCurves(ReferenceCurveEvaluator::new());
    for d in [4.0, 10.0] {
        let at = CurveMeasure::Distance(d);
        refused_by_name(p.path_point_at(&path, at));
        refused_by_name(p.path_tangent_at(&path, at));
        refused_by_name(p.path_frame_at(&path, at));
        for side in SIDES {
            refused_by_name(p.path_point_at_on(&path, at, side));
            refused_by_name(p.path_tangent_at_on(&path, at, side));
            refused_by_name(p.path_frame_at_on(&path, at, side));
            assert!(!p.path_frame_is_exact_at(&path, at, side));
        }
    }
    assert_eq!(
        p.path_distance_convention(&path),
        DistanceConvention::Unsupported
    );
    // Refusing paths is conformant; so is reading them right.
    assert_eq!(conformance::check(&p), Vec::new());
    assert_eq!(
        conformance::check(&ReferenceCurveEvaluator::new()),
        Vec::new()
    );
}

/// A provider that reads every joint from its outgoing piece.
#[derive(Debug)]
struct OutgoingOnly(ReferenceCurveEvaluator);

impl Backend for OutgoingOnly {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("outgoing-only"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl CurveEvaluator for OutgoingOnly {
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
        self.0.path_distance_convention(path)
    }
    fn path_point_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        _: SeamSide,
    ) -> GeomResult<Point3> {
        self.0.path_point_at_on(path, at, SeamSide::Outgoing)
    }
    fn path_tangent_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        _: SeamSide,
    ) -> GeomResult<Vec3> {
        self.0.path_tangent_at_on(path, at, SeamSide::Outgoing)
    }
    fn path_frame_at_on(
        &self,
        path: &CurvePath,
        at: CurveMeasure,
        _: SeamSide,
    ) -> GeomResult<Frame3> {
        self.0.path_frame_at_on(path, at, SeamSide::Outgoing)
    }
}

#[test]
fn the_conformance_suite_rejects_a_joint_read_from_the_wrong_side() {
    let failures = conformance::check(&OutgoingOnly(ReferenceCurveEvaluator::new()));
    assert!(
        failures
            .iter()
            .any(|f| f.check == "a path reads each side of a joint from its own piece"),
        "{failures:?}"
    );
}

#[test]
fn a_line_is_exact_only_where_it_is_placed_exactly() {
    let line = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::new(0.0, 3.0, 4.0),
    });
    let e = ReferenceCurveEvaluator::new();
    let at = CurveMeasure::Distance(2.0);
    let piece = PathPiece::new(PathCurve::Three(line.clone()), 0.0, 5.0);
    let rigid = axiolid_core::Transform3::from_rotation_z(0.3);
    for (path, exact) in [
        (CurvePath::from(piece.clone()), true),
        (CurvePath::from(piece.clone().placed(rigid, true)), true),
        (CurvePath::from(piece.clone().placed(rigid, false)), false),
        (
            CurvePath::from(piece.placed(rigid, false).placed(rigid, true)),
            false,
        ),
    ] {
        for side in SIDES {
            assert_eq!(e.path_frame_is_exact_at(&path, at, side), exact);
        }
        // The placed line's own point, carried.
        let own = Point3::new(0.0, 1.2, 1.6);
        let want = path.pieces()[0]
            .placement
            .map_or(own, |rigid| rigid.transform_point3(own));
        close(
            e.path_point_at(&path, at).unwrap() - want,
            Vec3::ZERO,
            1e-12,
            "placed point",
        );
    }
}

/// A provider that refuses paths as bad input instead of unsupported.
#[derive(Debug)]
struct Misrefusing(ReferenceCurveEvaluator);

impl Backend for Misrefusing {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new("misrefusing"), ExecutionTarget::PortableCpu)
    }
}

fn no_paths<T>() -> GeomResult<T> {
    Err(GeomError::InvalidInput("no paths".into()))
}

impl CurveEvaluator for Misrefusing {
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
    fn path_point_at_on(&self, _: &CurvePath, _: CurveMeasure, _: SeamSide) -> GeomResult<Point3> {
        no_paths()
    }
    fn path_tangent_at_on(&self, _: &CurvePath, _: CurveMeasure, _: SeamSide) -> GeomResult<Vec3> {
        no_paths()
    }
    fn path_frame_at_on(&self, _: &CurvePath, _: CurveMeasure, _: SeamSide) -> GeomResult<Frame3> {
        no_paths()
    }
}

#[test]
fn the_conformance_suite_rejects_a_path_refused_as_bad_input() {
    let failures = conformance::check(&Misrefusing(ReferenceCurveEvaluator::new()));
    assert!(
        failures
            .iter()
            .any(|f| f.check == "a path is answered or refused as unsupported"),
        "{failures:?}"
    );
}
