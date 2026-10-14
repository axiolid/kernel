//! Certified lengths of offsets by distances along an ellipse and a
//! B-spline, and window tolerances along them (#307).
//!
//! Every reference here is independent of the code under test: the offset
//! curve is written out in components, its speed differentiated by hand,
//! and its length integrated in the base's native parameter by a classical
//! fourth-order Runge-Kutta scheme that carries the base's arc length as a
//! second state (the law runs in it). The kernel measures by chords with
//! certified remainders; nothing is shared but the curve's definition.

use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{
    BSplineCurve, Curve2, Ellipse2, KnotSpec, Line2, OffsetFrame, OffsetLaw, PathOffsets,
};
use axiolid_curve_evaluate_contract::DistanceInterval;
use axiolid_evaluate::station::{offset_pieces, CompositeBasis, StationCurve, StationPiece};

const TAU: Scalar = core::f64::consts::TAU;

/// A plane curve written out: its first and second derivatives in its
/// native parameter.
struct Written<'a> {
    first: &'a dyn Fn(Scalar) -> Vec2,
    second: &'a dyn Fn(Scalar) -> Vec2,
}

/// The base's arc length over `[t0, t1]`, by Runge-Kutta (Simpson's rule,
/// for a quadrature).
fn arc(curve: &Written<'_>, t0: Scalar, t1: Scalar) -> Scalar {
    let steps = 100_000;
    let h = (t1 - t0) / Scalar::from(steps);
    let speed = |t: Scalar| (curve.first)(t).length();
    (0..steps)
        .map(|i| {
            let t = t0 + h * Scalar::from(i);
            h / 6.0 * (speed(t) + 4.0 * speed(t + 0.5 * h) + speed(t + h))
        })
        .sum()
}

/// The native parameter where the base's arc length from `t0` is `s`, by
/// bisection on [`arc`].
fn parameter_at(curve: &Written<'_>, t0: Scalar, hi: Scalar, s: Scalar) -> Scalar {
    let (mut a, mut b) = (t0, hi);
    for _ in 0..60 {
        let mid = 0.5 * (a + b);
        if arc(curve, t0, mid) < s {
            a = mid;
        } else {
            b = mid;
        }
    }
    0.5 * (a + b)
}

/// The length of the offset of `curve` over `[t0, t1]`, travelled from
/// `t0` (or from `t1`, `backwards`), its lateral offset (to the left of the
/// travel) running from `lateral.0` to `lateral.1` and its vertical one
/// from `vertical.0` to `vertical.1`, both linear in the base's arc length
/// from the travel's start.
fn reference(
    curve: &Written<'_>,
    (t0, t1): (Scalar, Scalar),
    backwards: bool,
    lateral: (Scalar, Scalar),
    vertical: (Scalar, Scalar),
) -> Scalar {
    let total = arc(curve, t0, t1);
    let rate = (t1 - t0) * if backwards { -1.0 } else { 1.0 };
    let at = |u: Scalar| {
        if backwards {
            t1 - u * (t1 - t0)
        } else {
            t0 + u * (t1 - t0)
        }
    };
    // `d(s, L)/du` at the travel's fraction `u`, `s` the base's arc length.
    let field = |u: Scalar, s: Scalar| -> (Scalar, Scalar) {
        let t = at(u);
        let d1 = (curve.first)(t) * rate;
        let d2 = (curve.second)(t) * (rate * rate);
        let speed = d1.length();
        let tangent = d1 / speed;
        let turn = (d2 - tangent * tangent.dot(d2)) / speed;
        let normal = tangent.perp();
        let dnormal = turn.perp();
        let a = lateral.0 + (lateral.1 - lateral.0) * s / total;
        let da = (lateral.1 - lateral.0) / total * speed;
        let dc = (vertical.1 - vertical.0) / total * speed;
        let plan = d1 + normal * da + dnormal * a;
        (speed, (plan.length_squared() + dc * dc).sqrt())
    };
    let steps = 200_000;
    let h = 1.0 / Scalar::from(steps);
    let (mut s, mut length) = (0.0, 0.0);
    for i in 0..steps {
        let u = h * Scalar::from(i);
        let k1 = field(u, s);
        let k2 = field(u + 0.5 * h, s + 0.5 * h * k1.0);
        let k3 = field(u + 0.5 * h, s + 0.5 * h * k2.0);
        let k4 = field(u + h, s + h * k3.0);
        s += h / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
        length += h / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
    }
    length
}

fn contains(
    interval: DistanceInterval,
    truth: Scalar,
    tolerance: Scalar,
    slack: Scalar,
    what: &str,
) {
    assert!(
        interval.lower <= interval.upper && interval.width() <= tolerance,
        "{what}: {interval:?} is wider than {tolerance:e}"
    );
    assert!(
        interval.lower - slack <= truth && truth <= interval.upper + slack,
        "{what}: {interval:?} does not contain {truth} (width {:e})",
        interval.width()
    );
}

fn frame(origin: Point2) -> Frame2 {
    Frame2 {
        origin,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

/// The deck's ellipse (openbimrs/ifc#443): centre `(30, 0)`, semi-axes 10
/// and 6.
fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: frame(Point2::new(30.0, 0.0)),
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    })
}

fn ellipse_first(t: Scalar) -> Vec2 {
    Vec2::new(-10.0 * t.sin(), 6.0 * t.cos())
}

fn ellipse_second(t: Scalar) -> Vec2 {
    Vec2::new(-10.0 * t.cos(), -6.0 * t.sin())
}

const ELLIPSE: Written<'static> = Written {
    first: &ellipse_first,
    second: &ellipse_second,
};

/// The control points of a cubic Bezier.
const BEZIER: [Point2; 4] = [
    Point2::new(0.0, 0.0),
    Point2::new(30.0, 0.0),
    Point2::new(60.0, 20.0),
    Point2::new(90.0, 20.0),
];

fn bezier() -> Curve2 {
    Curve2::BSpline(BSplineCurve {
        degree: 3,
        control_points: BEZIER.to_vec(),
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

fn bezier_first(t: Scalar) -> Vec2 {
    let p = BEZIER;
    let u = 1.0 - t;
    3.0 * (u * u * (p[1] - p[0]) + 2.0 * u * t * (p[2] - p[1]) + t * t * (p[3] - p[2]))
}

fn bezier_second(t: Scalar) -> Vec2 {
    let p = BEZIER;
    6.0 * ((1.0 - t) * (p[2] - 2.0 * p[1] + p[0]) + t * (p[3] - 2.0 * p[2] + p[1]))
}

const BEZIER_WRITTEN: Written<'static> = Written {
    first: &bezier_first,
    second: &bezier_second,
};

fn by_distances(start: (Scalar, Scalar), end: (Scalar, Scalar)) -> OffsetLaw {
    OffsetLaw::Linear {
        start: PathOffsets::new(start.0, start.1, 0.0),
        end: PathOffsets::new(end.0, end.1, 0.0),
        frame: OffsetFrame::Section,
    }
}

/// `base` offset by distances from `start` to `end`, then a line from
/// `origin` along `direction` offset by `end`, so that the first joint is
/// the offset's length.
fn widening_then_line(
    base: StationPiece<'_>,
    start: (Scalar, Scalar),
    end: (Scalar, Scalar),
    line: &Curve2,
) -> Vec<DistanceInterval> {
    let mut pieces = offset_pieces(&[base], by_distances(start, end)).unwrap();
    let tail = StationPiece::between(StationCurve::Two(line), 0.0, 10.0).unwrap();
    pieces.extend(offset_pieces(&[tail], by_distances(end, end)).unwrap());
    let basis = CompositeBasis::new(pieces).unwrap();
    [1e-6, 1e-8, 1e-9]
        .iter()
        .map(|&tolerance| basis.joint_intervals(tolerance).unwrap()[0])
        .collect()
}

#[test]
fn an_offset_by_distances_along_an_ellipse_contains_the_reference() {
    let curve = ellipse();
    let (start, end) = ((0.5, 0.0), (0.75, 0.3));
    // Forwards, from `(40, 0)` heading `+y` and back to it; the line goes
    // on along `+y`.
    let line = Curve2::Line(Line2 {
        origin: Point2::new(40.0, 0.0),
        direction: Vec2::Y,
    });
    let whole = StationPiece::whole(StationCurve::Two(&curve)).unwrap();
    let truth = reference(&ELLIPSE, (0.0, TAU), false, (0.5, 0.75), (0.0, 0.3));
    for (interval, tolerance) in widening_then_line(whole, start, end, &line)
        .into_iter()
        .zip([1e-6, 1e-8, 1e-9])
    {
        contains(interval, truth, tolerance, 0.0, "ellipse");
    }
    // Backwards, heading `-y` from `(40, 0)`; the line goes on along `-y`.
    let line = Curve2::Line(Line2 {
        origin: Point2::new(40.0, 0.0),
        direction: -Vec2::Y,
    });
    let truth = reference(&ELLIPSE, (0.0, TAU), true, (0.5, 0.75), (0.0, 0.3));
    for (interval, tolerance) in widening_then_line(whole.reversed(), start, end, &line)
        .into_iter()
        .zip([1e-6, 1e-8, 1e-9])
    {
        contains(interval, truth, tolerance, 0.0, "reversed ellipse");
    }
}

#[test]
fn an_offset_by_distances_along_a_b_spline_contains_the_reference() {
    let curve = bezier();
    let line = Curve2::Line(Line2 {
        origin: BEZIER[3],
        direction: Vec2::X,
    });
    let whole = StationPiece::whole(StationCurve::Two(&curve)).unwrap();
    let truth = reference(&BEZIER_WRITTEN, (0.0, 1.0), false, (1.0, 2.0), (0.0, 0.3));
    for (interval, tolerance) in widening_then_line(whole, (1.0, 0.0), (2.0, 0.3), &line)
        .into_iter()
        .zip([1e-6, 1e-8, 1e-9])
    {
        contains(interval, truth, tolerance, 0.0, "B-spline");
    }
    // Backwards, from `(90, 20)` heading `-x` to the origin; the line goes
    // on along `-x`. The Bezier is not symmetric, so a law read from the
    // wrong end gives another length.
    let line = Curve2::Line(Line2 {
        origin: BEZIER[0],
        direction: -Vec2::X,
    });
    let truth = reference(&BEZIER_WRITTEN, (0.0, 1.0), true, (1.0, 2.0), (0.0, 0.3));
    for (interval, tolerance) in widening_then_line(whole.reversed(), (1.0, 0.0), (2.0, 0.3), &line)
        .into_iter()
        .zip([1e-6, 1e-8, 1e-9])
    {
        contains(interval, truth, tolerance, 0.0, "reversed B-spline");
    }
}

/// The deck's offset (openbimrs/ifc#443): beside the ellipse by `0.5` to
/// station 2, by distances from `0.5` to `0.75` to station 20, by `0.75`
/// on to the end. Its joints at a window's certification of a `1e-5`
/// window, `1e-8`, contain the reference (within the arc-length tolerance
/// the kernel reads the stations' native parameters to), and a station 30
/// along, past the varying piece, is read where it stands.
#[test]
fn the_deck_s_offset_certifies_at_a_window_s_tolerance() {
    let curve = ellipse();
    let length = StationPiece::whole(StationCurve::Two(&curve))
        .unwrap()
        .length();
    let span = |a, b| StationPiece::between(StationCurve::Two(&curve), a, b).unwrap();
    let mut pieces =
        offset_pieces(&[span(0.0, 2.0)], by_distances((0.5, 0.0), (0.5, 0.0))).unwrap();
    pieces
        .extend(offset_pieces(&[span(2.0, 20.0)], by_distances((0.5, 0.0), (0.75, 0.0))).unwrap());
    pieces.extend(
        offset_pieces(
            &[span(20.0, length)],
            by_distances((0.75, 0.0), (0.75, 0.0)),
        )
        .unwrap(),
    );
    let basis = CompositeBasis::new(pieces).unwrap();
    let (t2, t20) = (
        parameter_at(&ELLIPSE, 0.0, 1.0, 2.0),
        parameter_at(&ELLIPSE, 0.0, 3.0, 20.0),
    );
    let first = reference(&ELLIPSE, (0.0, t2), false, (0.5, 0.5), (0.0, 0.0));
    let second = reference(&ELLIPSE, (t2, t20), false, (0.5, 0.75), (0.0, 0.0));
    let window = 1e-5;
    let tolerance = 1e-8;
    let joints = basis.joint_intervals(tolerance).unwrap();
    contains(joints[0], first, tolerance, 1e-10, "station 2");
    contains(joints[1], first + second, tolerance, 1e-10, "station 20");
    // A station 30 along, and one past the joint at station 20 but within
    // its window, at the window's certification.
    assert_eq!(basis.snapped_distance(30.0, window).unwrap(), 30.0);
    let nominal = basis.pieces()[0].length() + basis.pieces()[1].length();
    assert_eq!(
        basis
            .snapped_distance(nominal + 0.4 * window, window)
            .unwrap(),
        nominal
    );
}
