//! A seam side on the curve-evaluation contract (#286), and the axes of
//! `frame_at` (#242).
//!
//! The reference evaluator's sided queries must be the station reading
//! with the same side (`station_section3_on`), at a polyline corner, a
//! grade break and a cant break, and on either side of each within the
//! seam tolerance; off a seam both sides are the side-less answer. A
//! provider that does not implement sides refuses `Incoming` by name.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionTarget, GeomError, GeomResult, Operation,
};
use axiolid_core::{Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Curve2, Curve3, Elevated3, ElevationLaw, Line2,
    Line3, Polyline3, SeamSide,
};
use axiolid_curve_evaluate_contract::{
    conformance, CurveEvaluator, CurveMeasure, DistanceConvention, SEAM_SIDE_UNSUPPORTED,
};
use axiolid_evaluate::station::{station_section3, station_section3_on};
use axiolid_evaluate::ReferenceCurveEvaluator;

const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];

fn close(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    let off = (actual - expected).length();
    assert!(
        off <= eps,
        "{what}: {actual:?} != {expected:?} (off by {off:e})"
    );
}

fn corner() -> Curve3 {
    Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 3.0),
            Point3::new(4.0, 6.0, 3.0),
        ],
        closed: false,
    })
}

fn straight_plan() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })
}

fn grade_break() -> Curve3 {
    Curve3::Elevated(Elevated3::new(
        straight_plan(),
        ElevationLaw::Piecewise {
            breaks: vec![50.0],
            laws: vec![
                ElevationLaw::constant_grade(0.0, 0.02),
                ElevationLaw::constant_grade(1.0, -0.01),
            ],
        },
    ))
}

/// Level and straight: 150 mm cant for 100 m, then none.
fn cant_break() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(straight_plan(), ElevationLaw::level(0.0)),
        CantLaw::new(vec![
            CantPiece::constant(100.0, 0.15),
            CantPiece::constant(50.0, 0.0),
        ]),
        CantLaw::zero(150.0),
        1.5,
        BankConvention::TangentRotation,
    ))
}

/// Each curve with its seam's distance.
fn seams() -> [(&'static str, Curve3, Scalar); 3] {
    [
        ("polyline corner", corner(), 5.0),
        ("grade break", grade_break(), 50.0),
        ("cant break", cant_break(), 100.0),
    ]
}

#[test]
fn a_side_on_a_seam_is_the_station_reading_with_that_side() {
    let e = ReferenceCurveEvaluator::new();
    for (what, curve, seam) in seams() {
        // On the seam, and within its tolerance (1e-12 * s) either side.
        for s in [seam, seam - 0.4e-12 * seam, seam + 0.4e-12 * seam] {
            let at = CurveMeasure::Distance(s);
            for side in SIDES {
                let station = station_section3_on(&curve, s, side).unwrap();
                assert_eq!(
                    e.point_at_on(&curve, at, side).unwrap(),
                    station.point,
                    "{what} {side:?} at {s}: point"
                );
                assert_eq!(
                    e.tangent_at_on(&curve, at, side).unwrap(),
                    station.tangent,
                    "{what} {side:?} at {s}: tangent"
                );
                assert_eq!(
                    e.frame_at_on(&curve, at, side).unwrap(),
                    station.frame(),
                    "{what} {side:?} at {s}: frame"
                );
            }
            // The two sides differ there: the tangent turns at the corner
            // and the grade break, the roll jumps at the cant break.
            let incoming = e.frame_at_on(&curve, at, SeamSide::Incoming).unwrap();
            let outgoing = e.frame_at_on(&curve, at, SeamSide::Outgoing).unwrap();
            assert!(
                (incoming.y - outgoing.y).length() > 1e-6,
                "{what} at {s}: the sides read different pieces"
            );
        }
    }
}

#[test]
fn the_sides_read_the_pieces_closed_forms() {
    let e = ReferenceCurveEvaluator::new();
    let on = |curve: &Curve3, s: Scalar, side| {
        let at = CurveMeasure::Distance(s);
        (
            e.point_at_on(curve, at, side).unwrap(),
            e.frame_at_on(curve, at, side).unwrap(),
        )
    };
    // Corner: (0,0,0) -> (4,0,3) is 5 long, then 6 along +y.
    let (point, frame) = on(&corner(), 5.0, SeamSide::Incoming);
    close(point, Point3::new(4.0, 0.0, 3.0), 0.0, "corner point");
    close(frame.x, Vec3::new(0.8, 0.0, 0.6), 1e-15, "incoming leg");
    close(frame.y, Vec3::new(-0.6, 0.0, 0.8), 1e-15, "incoming up");
    let (_, frame) = on(&corner(), 5.0, SeamSide::Outgoing);
    close(frame.x, Vec3::Y, 1e-15, "outgoing leg");
    close(frame.y, Vec3::Z, 1e-15, "outgoing up");
    // Grade break: +2% then -1%.
    for (side, grade) in [(SeamSide::Incoming, 0.02), (SeamSide::Outgoing, -0.01)] {
        let k = Scalar::hypot(grade, 1.0);
        let (point, frame) = on(&grade_break(), 50.0, side);
        close(point, Point3::new(50.0, 0.0, 1.0), 1e-12, "grade point");
        close(frame.x, Vec3::new(1.0, 0.0, grade) / k, 1e-12, "grade");
        close(frame.y, Vec3::new(-grade, 0.0, 1.0) / k, 1e-12, "up");
        close(frame.z, -Vec3::Y, 1e-12, "right");
    }
    // Cant break: rolled by 150 mm over a 1.5 m base, then level.
    let (_, frame) = on(&cant_break(), 100.0, SeamSide::Incoming);
    let lateral = Vec3::new(0.0, 0.99_f64.sqrt(), 0.1);
    close(frame.x, Vec3::X, 1e-12, "banked tangent");
    close(frame.z, -lateral, 1e-12, "banked right");
    close(frame.y, Vec3::X.cross(lateral), 1e-12, "banked up");
    let (_, frame) = on(&cant_break(), 100.0, SeamSide::Outgoing);
    close(frame.y, Vec3::Z, 1e-12, "level up");
}

#[test]
fn off_a_seam_both_sides_are_the_side_less_answer() {
    let e = ReferenceCurveEvaluator::new();
    let line = Curve3::Line(Line3 {
        origin: Point3::new(1.0, 2.0, 3.0),
        direction: Vec3::new(2.0, 1.0, 0.0),
    });
    let mut cases: Vec<(Curve3, Scalar)> = vec![(line.clone(), -3.0), (line, 7.0)];
    for (_, curve, seam) in seams() {
        // Start, inside each piece, just outside the tolerance, and end.
        for s in [0.0, 0.5 * seam, seam - 1e-9, seam + 1e-9, seam + 1.0] {
            cases.push((curve.clone(), s));
        }
    }
    for (curve, s) in cases {
        let at = CurveMeasure::Distance(s);
        let (point, tangent, frame) = (
            e.point_at(&curve, at).unwrap(),
            e.tangent_at(&curve, at).unwrap(),
            e.frame_at(&curve, at).unwrap(),
        );
        for side in SIDES {
            assert_eq!(e.point_at_on(&curve, at, side).unwrap(), point, "at {s}");
            assert_eq!(e.tangent_at_on(&curve, at, side).unwrap(), tangent);
            assert_eq!(e.frame_at_on(&curve, at, side).unwrap(), frame);
        }
    }
}

#[test]
fn outgoing_exactly_on_a_seam_is_the_side_less_reading() {
    let e = ReferenceCurveEvaluator::new();
    for (what, curve, seam) in seams() {
        let at = CurveMeasure::Distance(seam);
        let sided = e.frame_at_on(&curve, at, SeamSide::Outgoing).unwrap();
        let plain = e.frame_at(&curve, at).unwrap();
        close(sided.origin, plain.origin, 1e-12, what);
        close(sided.x, plain.x, 1e-12, what);
        close(sided.y, plain.y, 1e-12, what);
        assert_eq!(station_section3(&curve, seam).unwrap().frame(), sided);
    }
}

#[test]
fn a_sided_query_keeps_the_side_less_refusals() {
    let e = ReferenceCurveEvaluator::new();
    for (_, curve, seam) in seams() {
        for side in SIDES {
            for bad in [Scalar::NAN, Scalar::INFINITY, -1.0, 10.0 * seam] {
                let at = CurveMeasure::Distance(bad);
                assert_eq!(
                    e.point_at_on(&curve, at, side).is_err(),
                    e.point_at(&curve, at).is_err(),
                    "point at {bad}"
                );
                assert_eq!(
                    e.frame_at_on(&curve, at, side).is_err(),
                    e.frame_at(&curve, at).is_err(),
                    "frame at {bad}"
                );
            }
        }
    }
}

#[test]
fn a_native_parameter_reads_a_side_where_it_is_the_station_measure() {
    let e = ReferenceCurveEvaluator::new();
    // An elevated or banked curve's parameter is its plan distance.
    for curve in [grade_break(), cant_break()] {
        let seam = if matches!(curve, Curve3::Elevated(_)) {
            50.0
        } else {
            100.0
        };
        for side in SIDES {
            assert_eq!(
                e.frame_at_on(&curve, CurveMeasure::Parameter(seam), side)
                    .unwrap(),
                e.frame_at_on(&curve, CurveMeasure::Distance(seam), side)
                    .unwrap(),
            );
        }
    }
    // A polyline's parameter is a vertex index, not located as a seam:
    // the incoming side is refused by name, the outgoing one side-less.
    let at = CurveMeasure::Parameter(1.0);
    let error = e
        .tangent_at_on(&corner(), at, SeamSide::Incoming)
        .unwrap_err();
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains("native parameter")),
        "{error}"
    );
    assert_eq!(
        e.tangent_at_on(&corner(), at, SeamSide::Outgoing).unwrap(),
        e.tangent_at(&corner(), at).unwrap()
    );
}

#[test]
fn another_reference_up_frames_each_side_against_it() {
    let up = Vec3::new(0.0, 1.0, 1.0);
    let e = ReferenceCurveEvaluator::with_up(up).unwrap();
    let curve = corner();
    // Each side's frame is the reference-up frame of its own piece: the
    // side-less frame a metre inside that piece, moved to the seam.
    for (side, inside) in [(SeamSide::Incoming, 4.0), (SeamSide::Outgoing, 6.0)] {
        let sided = e
            .frame_at_on(&curve, CurveMeasure::Distance(5.0), side)
            .unwrap();
        let plain = e.frame_at(&curve, CurveMeasure::Distance(inside)).unwrap();
        close(sided.origin, Point3::new(4.0, 0.0, 3.0), 0.0, "seam point");
        close(sided.x, plain.x, 1e-15, "tangent");
        close(sided.y, plain.y, 1e-15, "up against (0, 1, 1)");
        close(sided.z, plain.z, 1e-15, "right");
    }
    // A banked curve's roll is measured against +Z: its frame is refused
    // on either side, as side-lessly, while its point is not.
    let banked = cant_break();
    let at = CurveMeasure::Distance(100.0);
    for side in SIDES {
        assert!(e.frame_at_on(&banked, at, side).is_err());
        assert!(e.point_at_on(&banked, at, side).is_ok());
    }
    assert!(e.frame_at(&banked, at).is_err());
}

/// A provider that implements the side-less queries only, by delegation.
#[derive(Debug)]
struct SideLess(ReferenceCurveEvaluator);

const SIDE_LESS: BackendId = BackendId::new("side-less-stub");

impl Backend for SideLess {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(SIDE_LESS, ExecutionTarget::PortableCpu)
    }
}

impl CurveEvaluator for SideLess {
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

#[test]
fn a_provider_without_sides_refuses_incoming_by_name() {
    let stub = SideLess(ReferenceCurveEvaluator::new());
    let named = GeomError::UnsupportedInput {
        backend: SIDE_LESS,
        operation: Operation::CurveEvaluation,
        input: SEAM_SIDE_UNSUPPORTED,
    };
    for (what, curve, seam) in seams() {
        // On the seam and off it: the default never answers Incoming,
        // since it cannot tell a seam it would read wrongly.
        for s in [seam, 0.5 * seam] {
            let at = CurveMeasure::Distance(s);
            assert_eq!(
                stub.point_at_on(&curve, at, SeamSide::Incoming),
                Err(named.clone()),
                "{what} at {s}"
            );
            assert_eq!(
                stub.tangent_at_on(&curve, at, SeamSide::Incoming),
                Err(named.clone())
            );
            assert_eq!(
                stub.frame_at_on(&curve, at, SeamSide::Incoming),
                Err(named.clone())
            );
            // Outgoing is the side-less method itself.
            assert_eq!(
                stub.point_at_on(&curve, at, SeamSide::Outgoing),
                stub.point_at(&curve, at)
            );
            assert_eq!(
                stub.tangent_at_on(&curve, at, SeamSide::Outgoing),
                stub.tangent_at(&curve, at)
            );
            assert_eq!(
                stub.frame_at_on(&curve, at, SeamSide::Outgoing),
                stub.frame_at(&curve, at)
            );
        }
    }
    // The defaults keep a side-less provider conformant.
    let failures = conformance::check(&stub);
    assert!(failures.is_empty(), "conformance failures: {failures:#?}");
}

/// A provider whose Incoming answer is its outgoing one fails the
/// conformance suite: the suite tells a silent outgoing frame from a
/// refusal.
#[derive(Debug)]
struct Silent(ReferenceCurveEvaluator);

impl Backend for Silent {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(SIDE_LESS, ExecutionTarget::PortableCpu)
    }
}

impl CurveEvaluator for Silent {
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
    fn point_at_on(&self, curve: &Curve3, at: CurveMeasure, _: SeamSide) -> GeomResult<Point3> {
        self.0.point_at(curve, at)
    }
    fn tangent_at_on(&self, curve: &Curve3, at: CurveMeasure, _: SeamSide) -> GeomResult<Vec3> {
        self.0.tangent_at(curve, at)
    }
    fn frame_at_on(&self, curve: &Curve3, at: CurveMeasure, _: SeamSide) -> GeomResult<Frame3> {
        self.0.frame_at(curve, at)
    }
}

#[test]
fn the_conformance_suite_rejects_a_silent_outgoing_frame() {
    let failures = conformance::check(&Silent(ReferenceCurveEvaluator::new()));
    assert!(
        failures
            .iter()
            .any(|f| f.check == "a seam side reads its own piece"),
        "{failures:#?}"
    );
}

/// A provider that refuses `Incoming` as bad input, not as unsupported.
#[derive(Debug)]
struct Misrefusing(ReferenceCurveEvaluator);

impl Backend for Misrefusing {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(SIDE_LESS, ExecutionTarget::PortableCpu)
    }
}

fn misrefused<T>(side: SeamSide, outgoing: GeomResult<T>) -> GeomResult<T> {
    match side {
        SeamSide::Outgoing => outgoing,
        _ => Err(GeomError::InvalidInput("no incoming side".into())),
    }
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
    fn point_at_on(&self, curve: &Curve3, at: CurveMeasure, side: SeamSide) -> GeomResult<Point3> {
        misrefused(side, self.0.point_at(curve, at))
    }
    fn tangent_at_on(&self, curve: &Curve3, at: CurveMeasure, side: SeamSide) -> GeomResult<Vec3> {
        misrefused(side, self.0.tangent_at(curve, at))
    }
    fn frame_at_on(&self, curve: &Curve3, at: CurveMeasure, side: SeamSide) -> GeomResult<Frame3> {
        misrefused(side, self.0.frame_at(curve, at))
    }
}

#[test]
fn the_conformance_suite_rejects_a_side_refused_as_bad_input() {
    let failures = conformance::check(&Misrefusing(ReferenceCurveEvaluator::new()));
    assert!(
        failures
            .iter()
            .any(|f| f.check == "a side is answered or refused as unsupported"),
        "{failures:#?}"
    );
}

/// #242: `frame_at` is `x` tangent, `y` up, `z` right, as the contract
/// documents, on a level line, a grade and a banked curve.
#[test]
fn frame_at_lays_out_x_tangent_y_up_z_right() {
    let e = ReferenceCurveEvaluator::new();
    let at = CurveMeasure::Distance(10.0);
    let level = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    });
    let frame = e.frame_at(&level, at).unwrap();
    close(frame.x, Vec3::X, 0.0, "level: x tangent");
    close(frame.y, Vec3::Z, 0.0, "level: y up");
    close(frame.z, -Vec3::Y, 0.0, "level: z right");

    let g: Scalar = 0.05;
    let k = g.hypot(1.0);
    let grade = Curve3::Elevated(Elevated3::new(
        straight_plan(),
        ElevationLaw::constant_grade(0.0, g),
    ));
    let frame = e.frame_at(&grade, at).unwrap();
    close(
        frame.x,
        Vec3::new(1.0, 0.0, g) / k,
        1e-15,
        "grade: x tangent",
    );
    close(frame.y, Vec3::new(-g, 0.0, 1.0) / k, 1e-15, "grade: y up");
    close(frame.z, -Vec3::Y, 1e-15, "grade: z right");

    let banked = cant_break();
    let frame = e.frame_at(&banked, at).unwrap();
    let lateral = Vec3::new(0.0, 0.99_f64.sqrt(), 0.1);
    close(frame.x, Vec3::X, 1e-12, "banked: x tangent");
    close(
        frame.y,
        Vec3::X.cross(lateral),
        1e-12,
        "banked: y section up",
    );
    close(frame.z, -lateral, 1e-12, "banked: z right");
}
