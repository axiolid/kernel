//! Cutting a composite basis (#424): the native parameters a piece's
//! measures are read at, the certified length up to a cut, and the points a
//! composite starts and ends at. Expected values are the curve evaluators'
//! own points and a Simpson quadrature of the derivative, never the
//! station code under test.

use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, Circle2, CurvatureLaw, Curve2, Curve3, Ellipse2, Intrinsic2, KnotSpec, Line2,
    Line3, Polyline2, SeamSide,
};
use axiolid_evaluate::station::{CompositeBasis, DistanceInterval, StationCurve, StationPiece};
use axiolid_evaluate::{derivative2, evaluate2};

const PI: Scalar = core::f64::consts::PI;

fn frame(origin: Point2) -> Frame2 {
    Frame2 {
        origin,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

/// Speed 5 from `(1, 2)`.
fn line() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::new(1.0, 2.0),
        direction: Vec2::new(3.0, 4.0),
    })
}

fn circle() -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: frame(Point2::new(15.0, 0.0)),
        radius: 5.0,
    })
}

fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: frame(Point2::new(30.0, 0.0)),
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    })
}

/// An L: 3 along `+x`, then 4 along `+y`.
fn polyline() -> Curve2 {
    Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 4.0),
        ],
        closed: false,
    })
}

/// A cubic Bezier from `from`, about 95 long.
fn spline(from: Point2) -> Curve2 {
    let points = [
        Point2::new(0.0, 0.0),
        Point2::new(30.0, 0.0),
        Point2::new(60.0, 20.0),
        Point2::new(90.0, 20.0),
    ];
    Curve2::BSpline(BSplineCurve {
        degree: 3,
        control_points: points
            .iter()
            .map(|p| Point2::new(p.x + from.x, p.y + from.y))
            .collect(),
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

fn intrinsic() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        frame(Point2::ZERO),
        CurvatureLaw::circular(0.25),
        6.0,
    ))
}

fn lift(p: Point2) -> Point3 {
    Point3::new(p.x, p.y, 0.0)
}

fn near(a: Point3, b: Point3, eps: Scalar, what: &str) {
    assert!((a - b).abs().max_element() <= eps, "{what}: {a:?} != {b:?}");
}

/// The piece's native span is `expected`, and the curve evaluated at its
/// ends is the point a one-piece composite reads at its ends.
fn round_trip(curve: &Curve2, start: Scalar, end: Scalar, expected: (Scalar, Scalar)) {
    let piece = StationPiece::between(StationCurve::Two(curve), start, end).unwrap();
    let (t0, t1) = piece.native_span().unwrap();
    assert!((t0 - expected.0).abs() < 1e-9, "{t0} != {}", expected.0);
    assert!((t1 - expected.1).abs() < 1e-9, "{t1} != {}", expected.1);
    let basis = CompositeBasis::new(vec![piece]).unwrap();
    for (distance, t) in [(0.0, t0), (basis.length(), t1)] {
        let read = basis
            .section_on(distance, SeamSide::Outgoing)
            .unwrap()
            .point;
        near(
            read,
            lift(evaluate2(curve, t).unwrap()),
            1e-9,
            "native point",
        );
    }
}

#[test]
fn a_lines_span_is_its_measure_over_its_speed() {
    round_trip(&line(), 2.0, 7.5, (0.4, 1.5));
}

#[test]
fn a_circles_span_runs_past_its_seam() {
    // Radius 5, from 3 pi / 2 for a half turn.
    let start = 5.0 * 1.5 * PI;
    round_trip(&circle(), start, start + 5.0 * PI, (1.5 * PI, 2.5 * PI));
}

#[test]
fn an_ellipses_span_runs_past_its_seam() {
    let curve = ellipse();
    let length = StationPiece::whole(StationCurve::Two(&curve))
        .unwrap()
        .length();
    let piece =
        StationPiece::between(StationCurve::Two(&curve), 0.8 * length, 1.2 * length).unwrap();
    let (t0, t1) = piece.native_span().unwrap();
    assert!(t0 < 2.0 * PI && t1 > 2.0 * PI, "{t0} {t1}");
    let basis = CompositeBasis::new(vec![piece]).unwrap();
    for (distance, t) in [(0.0, t0), (basis.length(), t1)] {
        let read = basis
            .section_on(distance, SeamSide::Outgoing)
            .unwrap()
            .point;
        near(
            read,
            lift(evaluate2(&curve, t).unwrap()),
            1e-9,
            "ellipse point",
        );
    }
}

#[test]
fn a_polylines_span_is_the_arc_length_inverse() {
    let curve = polyline();
    let piece = StationPiece::between(StationCurve::Two(&curve), 1.0, 5.0).unwrap();
    let (t0, t1) = piece.native_span().unwrap();
    near(
        lift(evaluate2(&curve, t0).unwrap()),
        Point3::new(1.0, 0.0, 0.0),
        1e-9,
        "start",
    );
    near(
        lift(evaluate2(&curve, t1).unwrap()),
        Point3::new(3.0, 2.0, 0.0),
        1e-9,
        "end",
    );
}

#[test]
fn a_b_splines_span_is_the_arc_length_inverse() {
    let curve = spline(Point2::ZERO);
    let piece = StationPiece::between(StationCurve::Two(&curve), 10.0, 60.0).unwrap();
    let (t0, t1) = piece.native_span().unwrap();
    let basis = CompositeBasis::new(vec![piece]).unwrap();
    near(
        basis.section_on(0.0, SeamSide::Outgoing).unwrap().point,
        lift(evaluate2(&curve, t0).unwrap()),
        1e-9,
        "start",
    );
    near(
        basis.section_on(50.0, SeamSide::Outgoing).unwrap().point,
        lift(evaluate2(&curve, t1).unwrap()),
        1e-9,
        "end",
    );
}

#[test]
fn an_intrinsic_curves_span_is_its_measure() {
    round_trip(&intrinsic(), 1.0, 4.0, (1.0, 4.0));
}

#[test]
fn a_3d_lines_span_is_its_measure() {
    let curve = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    });
    let piece = StationPiece::between(StationCurve::Three(&curve), 0.0, 2.0).unwrap();
    assert_eq!(piece.native_span().unwrap(), (0.0, 2.0));
}

/// Simpson's rule, `1e5` panels, of the speed over `[a, b]`.
fn simpson_length(curve: &Curve2, a: Scalar, b: Scalar) -> Scalar {
    let panels = 100_000;
    let h = (b - a) / Scalar::from(panels);
    let speed = |t: Scalar| derivative2(curve, t).unwrap().length();
    let mut sum = speed(a) + speed(b);
    for i in 1..panels {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * speed(a + h * Scalar::from(i));
    }
    sum * h / 3.0
}

#[test]
fn a_distance_interval_contains_the_true_length_to_the_cut() {
    // A segment of 3, then the spline from its end.
    let segment = Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    });
    let curve = spline(Point2::new(3.0, 0.0));
    let basis = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(&segment), 0.0, 3.0).unwrap(),
        StationPiece::whole(StationCurve::Two(&curve)).unwrap(),
    ])
    .unwrap();
    let tolerance = 1e-6;
    for distance in [20.0, 60.0, 90.0] {
        let interval = basis.distance_interval_at(distance, tolerance).unwrap();
        assert!(interval.width() <= tolerance, "{distance}: {interval:?}");
        assert!(!interval.is_point(), "{distance}: {interval:?}");
        // The cut's spline piece, from the curve's start to where the cut
        // reads it, plus the segment before it.
        let cut = basis.pieces_between(0.0, distance).unwrap();
        let (t0, t1) = cut[1].native_span().unwrap();
        let reference = 3.0 + simpson_length(&curve, t0, t1);
        assert!(
            interval.contains(reference),
            "{distance}: {interval:?} misses {reference}"
        );
    }
    // Inside the segment, its own length.
    assert_eq!(
        basis.distance_interval_at(1.5, tolerance).unwrap(),
        DistanceInterval::point(1.5)
    );
}

#[test]
fn a_distance_interval_on_stored_lengths_is_a_point() {
    let a = Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    });
    let b = Curve2::Line(Line2 {
        origin: Point2::new(2.0, 0.0),
        direction: Vec2::Y,
    });
    let basis = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(&a), 0.0, 2.0).unwrap(),
        StationPiece::between(StationCurve::Two(&b), 0.0, 3.0).unwrap(),
    ])
    .unwrap();
    for distance in [0.0, 0.7, 2.0, 4.25, 5.0] {
        let interval = basis.distance_interval_at(distance, 1e-3).unwrap();
        assert!(
            interval.is_point() && interval.lower == distance,
            "{interval:?}"
        );
    }
    assert!(basis.distance_interval_at(5.5, 1e-3).is_err());
    assert!(basis.distance_interval_at(-0.5, 1e-3).is_err());
    assert!(basis.distance_interval_at(1.0, 0.0).is_err());
}

#[test]
fn end_points_are_the_first_start_and_the_last_end() {
    let a = Curve2::Line(Line2 {
        origin: Point2::new(1.0, 1.0),
        direction: Vec2::X,
    });
    let b = Curve2::Line(Line2 {
        origin: Point2::new(3.0, 1.0),
        direction: Vec2::Y,
    });
    let basis = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(&a), 0.0, 2.0).unwrap(),
        StationPiece::between(StationCurve::Two(&b), 0.0, 3.0).unwrap(),
    ])
    .unwrap();
    let [first, last] = basis.end_points().unwrap();
    assert_eq!(first, Point3::new(1.0, 1.0, 0.0));
    near(last, Point3::new(3.0, 4.0, 0.0), 1e-12, "last");
}

#[test]
fn a_distance_interval_inside_a_reversed_piece_measures_from_its_end() {
    // A segment ending at the spline's end `(93, 20)`, then the spline run
    // backwards: the cut lies at 40 from its END, not from its start.
    let curve = spline(Point2::new(3.0, 0.0));
    let segment = Curve2::Line(Line2 {
        origin: Point2::new(90.0, 20.0),
        direction: Vec2::X,
    });
    let basis = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(&segment), 0.0, 3.0).unwrap(),
        StationPiece::whole(StationCurve::Two(&curve))
            .unwrap()
            .reversed(),
    ])
    .unwrap();
    let tolerance = 1e-6;
    let distance = 3.0 + 40.0;
    let interval = basis.distance_interval_at(distance, tolerance).unwrap();
    assert!(interval.width() <= tolerance && !interval.is_point());
    let cut = basis.pieces_between(0.0, distance).unwrap();
    let (t0, t1) = cut[1].native_span().unwrap();
    let reference = 3.0 + simpson_length(&curve, t0, t1);
    assert!(
        interval.contains(reference),
        "{interval:?} misses {reference}"
    );
    // The last 40 are the curve's final 40, so its parameter is past the middle.
    assert!(t0 > 0.5 && t1 == 1.0, "{t0} {t1}");
}
