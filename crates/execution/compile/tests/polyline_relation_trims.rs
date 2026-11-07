//! Trims of trimmed, offset and placed curves compiled to polylines (#424).
//!
//! A trim reads its basis's parameter (ISO 10303-42), an offset its basis's,
//! a placed curve its source's, as a sweep directrix does. Expected points
//! are the reference evaluator's own, written here, and every dense sample
//! of the exact curve must lie within the reported bound.

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve2, Circle2, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2,
    KnotSpec, Line2, Polyline2, Polyline3,
};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{PolylineCompiler, PolylineOutcome};
use axiolid_model::{
    CurveRelation, CurveSegment, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode,
    Instance, InstanceAtStation, NodeId, Station, StationFrame, Transition, TrimSelector,
    TrimmingPreference,
};
use axiolid_reference::curve::{derivative2, evaluate2};
use axiolid_reference::station::CompositeBasis;

const CHORD: Scalar = 1e-3;

/// Dense samples per curve.
const SAMPLES: usize = 20_000;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(CHORD)
        .expect("positive chord budget")
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn graph_of(node: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let root = node(&mut b);
    (b.finish(vec![root]).unwrap(), root)
}

fn polyline(graph: &GeometryGraph, root: NodeId) -> GeomResult<PolylineOutcome> {
    compiler().compile_polyline(graph, root, &options())
}

fn lift(p: Point2) -> Point3 {
    Point3::new(p.x, p.y, 0.0)
}

fn segment_distance(p: Point3, a: Point3, b: Point3) -> Scalar {
    let ab = b - a;
    let len2 = ab.length_squared();
    let t = if len2 == 0.0 {
        0.0
    } else {
        ((p - a).dot(ab) / len2).clamp(0.0, 1.0)
    };
    (p - (a + ab * t)).length()
}

fn distance_to(points: &[Point3], p: Point3) -> Scalar {
    points
        .windows(2)
        .map(|s| segment_distance(p, s[0], s[1]))
        .fold(Scalar::INFINITY, Scalar::min)
}

/// Every one of `samples` exact points within the reported bound, the
/// bound within the budget; the largest sampled distance is returned.
fn within_bound(
    outcome: &PolylineOutcome,
    exact: impl Fn(Scalar) -> Point3,
    samples: usize,
) -> Scalar {
    let bound = outcome.deviation.bound;
    assert_eq!(outcome.deviation.requested, CHORD);
    assert!(outcome.deviation.meets_requested(), "{bound} > {CHORD}");
    let scale = outcome
        .points
        .iter()
        .map(|p| p.abs().max_element())
        .fold(1.0, Scalar::max);
    let mut worst: Scalar = 0.0;
    for i in 0..=samples {
        let p = exact(i as Scalar / samples as Scalar);
        let d = distance_to(&outcome.points, p);
        assert!(
            d <= bound + 1e-12 * scale,
            "sample {i} at {p:?} is {d} from the polyline, beyond the bound {bound}"
        );
        worst = worst.max(d);
    }
    worst
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, needle: &str) {
    match result {
        Err(GeomError::UnsupportedInput {
            operation: Operation::CurveEvaluation,
            input,
            ..
        }) => assert!(input.contains(needle), "{input:?} does not name {needle:?}"),
        other => panic!("expected a refusal naming {needle:?}, got {other:?}"),
    }
}

fn plane() -> Frame2 {
    Frame2 {
        origin: Point2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

fn circle2(centre: Point2, radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: centre,
            ..plane()
        },
        radius,
    })
}

fn trim(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    start: TrimSelector,
    end: TrimSelector,
    sense: bool,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![start],
        end: vec![end],
        sense_agreement: sense,
        preference: TrimmingPreference::Parameter,
    }))
    .unwrap()
}

fn offset(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    distance: Scalar,
    reference_direction: Option<Vec3>,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Offset {
        basis,
        distance,
        reference_direction,
    }))
    .unwrap()
}

fn spline2() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 3.0),
            Point2::new(4.0, 2.5),
            Point2::new(6.0, 4.0),
            Point2::new(8.0, 1.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// No vertex repeated by its successor: every joint is one vertex.
fn no_repeated_vertex(points: &[Point3]) {
    for (i, pair) in points.windows(2).enumerate() {
        assert_ne!(pair[0], pair[1], "vertex {i} is repeated");
    }
}

fn composite(b: &mut GeometryGraphBuilder, segments: &[(NodeId, bool)]) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments: segments
            .iter()
            .map(|&(curve, same_sense)| CurveSegment {
                curve,
                same_sense,
                transition: Transition::Continuous,
            })
            .collect(),
    }))
    .unwrap()
}

fn segment2(b: &mut GeometryGraphBuilder, from: (Scalar, Scalar), to: (Scalar, Scalar)) -> NodeId {
    b.push_value(Curve2::Polyline(Polyline2 {
        points: vec![Point2::new(from.0, from.1), Point2::new(to.0, to.1)],
        closed: false,
    }))
    .unwrap()
}

/// A segment `(0, 0) -> (3, 0)`, a quarter arc about `(3, 1)` to `(4, 1)`,
/// and a segment `(4, 4) -> (4, 1)` run backwards: `3 + pi / 2 + 3` long.
fn bent(b: &mut GeometryGraphBuilder) -> NodeId {
    let first = segment2(b, (0.0, 0.0), (3.0, 0.0));
    let basis = b.push_value(circle2(Point2::new(3.0, 1.0), 1.0)).unwrap();
    let arc = trim(b, basis, p(-FRAC_PI_2), p(0.0), true);
    let last = segment2(b, (4.0, 4.0), (4.0, 1.0));
    composite(b, &[(first, true), (arc, true), (last, false)])
}

/// The point of [`bent`] at arc length `d`.
fn bent_at(d: Scalar) -> Point3 {
    if d <= 3.0 {
        Point3::new(d, 0.0, 0.0)
    } else if d <= 3.0 + FRAC_PI_2 {
        let angle = -FRAC_PI_2 + (d - 3.0);
        Point3::new(3.0 + angle.cos(), 1.0 + angle.sin(), 0.0)
    } else {
        Point3::new(4.0, 1.0 + (d - 3.0 - FRAC_PI_2), 0.0)
    }
}

fn p(value: Scalar) -> TrimSelector {
    TrimSelector::Parameter(value)
}

fn close(a: Point3, b: Point3, what: &str) {
    assert!((a - b).length() <= 1e-12, "{what}: {a:?} != {b:?}");
}

fn invalid<T: core::fmt::Debug>(result: GeomResult<T>) {
    assert!(
        matches!(result, Err(GeomError::InvalidInput(_))),
        "expected invalid input, got {result:?}"
    );
}

#[test]
fn a_trim_of_a_trim_reads_its_basis_parameter() {
    let circle = circle2(Point2::ZERO, 2.0);
    let on = |t: Scalar| lift(evaluate2(&circle, t).unwrap());
    for sense in [true, false] {
        let (graph, root) = graph_of(|b| {
            let basis = b.push_value(circle.clone()).unwrap();
            let inner = trim(b, basis, p(0.0), p(PI), true);
            trim(b, inner, p(FRAC_PI_4), p(FRAC_PI_2), sense)
        });
        let outcome = polyline(&graph, root).unwrap();
        // The arc of the circle at angles [pi / 4, pi / 2], bitwise at its
        // ends, reversed by the outer sense.
        let (first, last) = if sense {
            (FRAC_PI_4, FRAC_PI_2)
        } else {
            (FRAC_PI_2, FRAC_PI_4)
        };
        assert_eq!(outcome.points[0], on(first), "sense {sense}");
        assert_eq!(*outcome.points.last().unwrap(), on(last), "sense {sense}");
        within_bound(&outcome, |s| on(FRAC_PI_4 + FRAC_PI_4 * s), SAMPLES);
    }
}

#[test]
fn a_trim_of_a_reversed_trim_reads_the_circles_angles() {
    let circle = circle2(Point2::ZERO, 2.0);
    let on = |t: Scalar| lift(evaluate2(&circle, t).unwrap());
    // The inner trim is theta in [-pi, 0], run decreasing.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        let inner = trim(b, basis, p(0.0), p(PI), false);
        trim(b, inner, p(1.5 * PI), p(1.75 * PI), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    // The arc [-pi / 2, -pi / 4], in the inner trim's direction: from
    // circle(-pi / 4) down to circle(-pi / 2).
    close(outcome.points[0], on(-FRAC_PI_4), "first");
    close(*outcome.points.last().unwrap(), on(-FRAC_PI_2), "last");
    within_bound(&outcome, |s| on(-FRAC_PI_2 + FRAC_PI_4 * s), SAMPLES);
    // The first quarter lies on the half the inner trim leaves out.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        let inner = trim(b, basis, p(0.0), p(PI), false);
        trim(b, inner, p(FRAC_PI_4), p(FRAC_PI_2), true)
    });
    invalid(polyline(&graph, root));
}

#[test]
fn a_trim_of_a_bsplines_trim_narrows_and_refuses_outside() {
    let spline = spline2();
    let on = |t: Scalar| lift(evaluate2(&spline, t).unwrap());
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(spline.clone()).unwrap();
        let inner = trim(b, basis, p(0.1), p(0.9), true);
        trim(b, inner, p(0.3), p(0.6), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    assert_eq!(outcome.points[0], on(0.3));
    assert_eq!(*outcome.points.last().unwrap(), on(0.6));
    within_bound(&outcome, |s| on(0.3 + 0.3 * s), SAMPLES);

    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(spline.clone()).unwrap();
        let inner = trim(b, basis, p(0.1), p(0.9), true);
        trim(b, inner, p(0.05), p(0.6), true)
    });
    invalid(polyline(&graph, root));
}

#[test]
fn a_trim_of_an_offset_reads_its_basis_parameter() {
    let circle = circle2(Point2::ZERO, 2.0);
    // Anticlockwise the left is inwards: radius 1.5.
    let on = |t: Scalar| Point3::new(1.5 * t.cos(), 1.5 * t.sin(), 0.0);
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        let moved = offset(b, basis, 0.5, None);
        trim(b, moved, p(0.0), p(FRAC_PI_2), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    close(outcome.points[0], on(0.0), "first");
    close(*outcome.points.last().unwrap(), on(FRAC_PI_2), "last");
    within_bound(&outcome, |s| on(FRAC_PI_2 * s), SAMPLES);

    // Across the seam: from 7 pi / 4 to pi / 4 a turn on.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        let moved = offset(b, basis, 0.5, None);
        trim(b, moved, p(1.75 * PI), p(FRAC_PI_4), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    close(outcome.points[0], on(-FRAC_PI_4), "wrapped first");
    close(
        *outcome.points.last().unwrap(),
        on(FRAC_PI_4),
        "wrapped last",
    );
    within_bound(&outcome, |s| on(-FRAC_PI_4 + FRAC_PI_2 * s), SAMPLES);
}

#[test]
fn a_trim_of_a_placed_curve_is_the_placed_arc() {
    let circle = circle2(Point2::ZERO, 2.0);
    let rigid =
        Transform3::from_translation(Vec3::new(1.0, 2.0, 0.0)) * Transform3::from_rotation_z(0.3);
    let on = |t: Scalar| rigid.transform_point3(lift(evaluate2(&circle, t).unwrap()));
    let (graph, root) = graph_of(|b| {
        let source = b.push_value(circle.clone()).unwrap();
        let placed = b
            .push(GeometryNode::Instance(Instance {
                source,
                transform: rigid,
            }))
            .unwrap();
        trim(b, placed, p(0.0), p(FRAC_PI_2), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.deviation.meets_requested());
    close(outcome.points[0], on(0.0), "first");
    close(*outcome.points.last().unwrap(), on(FRAC_PI_2), "last");
    within_bound(&outcome, |s| on(FRAC_PI_2 * s), SAMPLES);
}

#[test]
fn a_selector_other_than_a_parameter_on_a_relation_is_refused() {
    let circle = circle2(Point2::ZERO, 2.0);
    let selectors = |start: TrimSelector, end: TrimSelector, preference| {
        graph_of(|b| {
            let basis = b.push_value(circle.clone()).unwrap();
            let inner = trim(b, basis, p(0.0), p(PI), true);
            b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
                basis: inner,
                start: vec![start],
                end: vec![end],
                sense_agreement: true,
                preference,
            }))
            .unwrap()
        })
    };
    let (graph, root) = selectors(
        TrimSelector::Point2(Point2::new(2.0, 0.0)),
        TrimSelector::Point2(Point2::new(0.0, 2.0)),
        TrimmingPreference::Cartesian,
    );
    refused(
        polyline(&graph, root),
        "a point or arc-length trim selector",
    );
    let (graph, root) = selectors(
        TrimSelector::ArcLength(1.0),
        TrimSelector::ArcLength(2.0),
        TrimmingPreference::Parameter,
    );
    refused(
        polyline(&graph, root),
        "a point or arc-length trim selector",
    );
    // A parameter beside a Cartesian preference is refused too.
    let (graph, root) = selectors(p(0.5), p(1.0), TrimmingPreference::Cartesian);
    refused(
        polyline(&graph, root),
        "a point or arc-length trim selector",
    );
}

#[test]
fn a_trim_of_a_composite_is_cut_in_its_arc_length() {
    let (from, to) = (1.5, 3.0 + FRAC_PI_4);
    let (graph, root) = graph_of(|b| {
        let whole = bent(b);
        trim(b, whole, p(from), p(to), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    close(outcome.points[0], Point3::new(1.5, 0.0, 0.0), "first");
    close(
        *outcome.points.last().unwrap(),
        Point3::new(3.0 + (-FRAC_PI_4).cos(), 1.0 + (-FRAC_PI_4).sin(), 0.0),
        "last",
    );
    within_bound(&outcome, |s| bent_at(from + (to - from) * s), SAMPLES);
    no_repeated_vertex(&outcome.points);

    // Against the sense, the same points the other way round.
    let (graph, root) = graph_of(|b| {
        let whole = bent(b);
        trim(b, whole, p(from), p(to), false)
    });
    let against = polyline(&graph, root).unwrap();
    let mut reversed = outcome.points.clone();
    reversed.reverse();
    assert_eq!(against.points, reversed);
}

#[test]
fn a_cut_at_a_joint_leaves_no_repeated_vertex() {
    let (graph, root) = graph_of(|b| {
        let whole = bent(b);
        trim(b, whole, p(3.0), p(4.0 + FRAC_PI_2), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    no_repeated_vertex(&outcome.points);
    close(outcome.points[0], Point3::new(3.0, 0.0, 0.0), "first");
    close(
        *outcome.points.last().unwrap(),
        Point3::new(4.0, 2.0, 0.0),
        "last",
    );
    // Ending on the first joint.
    let (graph, root) = graph_of(|b| {
        let whole = bent(b);
        trim(b, whole, p(1.0), p(3.0), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    no_repeated_vertex(&outcome.points);
    close(
        *outcome.points.last().unwrap(),
        Point3::new(3.0, 0.0, 0.0),
        "last",
    );
}

/// A cubic Bezier from `(3, 0)`, about 95 long, whose arc length is a
/// quadrature.
fn bezier() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(3.0, 0.0),
            Point2::new(33.0, 0.0),
            Point2::new(63.0, 20.0),
            Point2::new(93.0, 20.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// A segment `(0, 0) -> (3, 0)` and a degree-one B-spline from `(3, 0)`
/// along the x axis, a quadrature too: the exact curve is easy to write.
fn segment_then_spline(b: &mut GeometryGraphBuilder) -> NodeId {
    let first = segment2(b, (0.0, 0.0), (3.0, 0.0));
    let spline = b
        .push_value(Curve2::BSpline(BSplineCurve2 {
            degree: 1,
            control_points: vec![
                Point2::new(3.0, 0.0),
                Point2::new(4.0, 0.0),
                Point2::new(6.0, 0.0),
            ],
            knots: vec![0.0, 0.5, 1.0],
            multiplicities: vec![2, 1, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }))
        .unwrap();
    composite(b, &[(first, true), (spline, true)])
}

/// Simpson's rule, `1e5` panels, of the speed of `curve` over `[a, b]`.
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
fn a_cut_inside_a_quadrature_piece_adds_its_certified_error() {
    let curve = bezier();
    let (from, to) = (1.0, 20.0);
    let (graph, root) = graph_of(|b| {
        let first = segment2(b, (0.0, 0.0), (3.0, 0.0));
        let spline = b.push_value(curve.clone()).unwrap();
        let whole = composite(b, &[(first, true), (spline, true)]);
        trim(b, whole, p(from), p(to), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    let whole = match graph.get(root) {
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed { basis, .. })) => *basis,
        other => panic!("{other:?}"),
    };
    let path = axiolid_mesh_compile::station::curve_path(&graph, whole).unwrap();
    let basis = CompositeBasis::from_path(&path).unwrap();
    // Where the polyline ends: the curve at the native parameter the cut
    // reads, whose true length from the start is a quadrature of our own.
    let cut = basis.pieces_between(from, to).unwrap();
    let (t0, t1) = cut.last().unwrap().native_span().unwrap();
    close(
        *outcome.points.last().unwrap(),
        lift(evaluate2(&curve, t1).unwrap()),
        "last",
    );
    let reached = 3.0 + simpson_length(&curve, t0, t1);
    let off = (reached - to).abs();
    // The certified error of the cut, from the same basis.
    let tolerance = 0.125 * CHORD;
    let start = basis.distance_interval_at(from, tolerance).unwrap();
    let end = basis.distance_interval_at(to, tolerance).unwrap();
    let expected = (start.upper - from).max(to - end.lower).max(0.0);
    assert!(expected > 0.0, "the cut is exact: {start:?} {end:?}");
    assert!(end.contains(reached), "{end:?} misses the true {reached}");
    assert!(
        outcome.deviation.bound >= expected && outcome.deviation.bound >= off,
        "{} < {expected} / {off}",
        outcome.deviation.bound
    );
}

fn left_offset(b: &mut GeometryGraphBuilder, reversed_first: bool) -> NodeId {
    let first = if reversed_first {
        segment2(b, (3.0, 0.0), (0.0, 0.0))
    } else {
        segment2(b, (0.0, 0.0), (3.0, 0.0))
    };
    let basis = b.push_value(circle2(Point2::new(3.0, 1.0), 1.0)).unwrap();
    let arc = trim(b, basis, p(-FRAC_PI_2), p(0.0), true);
    let whole = composite(b, &[(first, !reversed_first), (arc, true)]);
    offset(b, whole, 0.25, None)
}

#[test]
fn a_trim_of_an_offset_composite_is_offset_on_the_side_it_is_traversed() {
    for reversed in [false, true] {
        let (graph, root) = graph_of(|b| {
            let moved = left_offset(b, reversed);
            trim(b, moved, p(1.0), p(2.0), true)
        });
        let outcome = polyline(&graph, root).unwrap();
        close(outcome.points[0], Point3::new(1.0, 0.25, 0.0), "first");
        close(
            *outcome.points.last().unwrap(),
            Point3::new(2.0, 0.25, 0.0),
            "last",
        );
        assert_eq!(outcome.deviation.bound, 0.0, "reversed {reversed}");
    }
}

#[test]
fn a_trim_of_a_trimmed_composite_is_the_trim_of_the_composite() {
    let (graph, root) = graph_of(|b| {
        let whole = bent(b);
        let inner = trim(b, whole, p(0.5), p(5.0), true);
        trim(b, inner, p(1.0), p(2.0), true)
    });
    let nested = polyline(&graph, root).unwrap();
    let (graph, root) = graph_of(|b| {
        let whole = bent(b);
        trim(b, whole, p(1.0), p(2.0), true)
    });
    let direct = polyline(&graph, root).unwrap();
    assert_eq!(nested.points, direct.points);
    assert_eq!(nested.deviation.bound, direct.deviation.bound);
}

#[test]
fn a_composite_holding_a_trim_of_a_trim_reads_the_basis_parameter() {
    // The arc of the radius-2 circle at the angles [pi / 4, pi / 2], then a
    // segment from where it ends.
    let circle = circle2(Point2::ZERO, 2.0);
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        let inner = trim(b, basis, p(0.0), p(PI), true);
        let arc = trim(b, inner, p(FRAC_PI_4), p(FRAC_PI_2), true);
        let tail = segment2(b, (0.0, 2.0), (-3.0, 2.0));
        let whole = composite(b, &[(arc, true), (tail, true)]);
        trim(b, whole, p(0.0), p(1.0 + FRAC_PI_2), true)
    });
    let outcome = polyline(&graph, root).unwrap();
    let on = |t: Scalar| lift(evaluate2(&circle, t).unwrap());
    close(outcome.points[0], on(FRAC_PI_4), "first");
    close(
        *outcome.points.last().unwrap(),
        Point3::new(-1.0, 2.0, 0.0),
        "last",
    );
}

#[test]
fn a_trimmed_curve_placed_at_a_station_in_a_composite_comes_out_placed() {
    let mut b = GeometryGraphBuilder::new();
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    // A curve placed at a station is 3D, and so is the segment before it.
    let first = b
        .push_value(Curve3::Polyline(Polyline3 {
            points: vec![Point3::ZERO, Point3::new(3.0, 0.0, 0.0)],
            closed: false,
        }))
        .unwrap();
    let source = segment2(&mut b, (0.0, 0.0), (0.0, 2.0));
    let station = CurveStation {
        basis: line,
        station: Station::at(3.0),
        frame: StationFrame::Section,
    };
    let placed = b
        .push_value(InstanceAtStation::new(source, station.into()))
        .unwrap();
    let whole = composite(&mut b, &[(first, true), (placed, true)]);
    let cut = trim(&mut b, whole, p(1.0), p(4.0), true);
    let graph = b.finish(vec![cut]).unwrap();
    let placement = axiolid_mesh_compile::station::placement(&graph, placed).unwrap();
    assert!(placement.exact);
    let outcome = polyline(&graph, cut).unwrap();
    // One metre into the placed segment, local `(0, 1)`.
    close(
        *outcome.points.last().unwrap(),
        placement
            .transform
            .transform_point3(Point3::new(0.0, 1.0, 0.0)),
        "placed last",
    );
    close(outcome.points[0], Point3::new(1.0, 0.0, 0.0), "first");
}

#[test]
fn a_composite_whose_joint_leaves_a_gap_is_refused_trimmed() {
    let hair = 2e-7;
    let (graph, root) = graph_of(|b| {
        let first = segment2(b, (0.0, 0.0), (3.0, 0.0));
        let second = segment2(b, (3.0 + hair, 0.0), (6.0, 0.0));
        let whole = composite(b, &[(first, true), (second, true)]);
        trim(b, whole, p(1.0), p(5.0), true)
    });
    let result = polyline(&graph, root);
    assert!(
        matches!(&result, Err(GeomError::InvalidInput(detail)) if detail.contains("gap")),
        "{result:?}"
    );
}

#[test]
fn the_remaining_trims_of_composites_are_refused_by_name() {
    // A composite measured in plan distance.
    let (graph, root) = graph_of(|b| {
        let plan = Intrinsic2::new(
            Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
            60.0,
        );
        let lifted = b
            .push_value(Curve3::Elevated(Elevated3::new(
                Curve2::Intrinsic(plan),
                ElevationLaw::constant_grade(60.0, 0.02),
            )))
            .unwrap();
        let whole = composite(b, &[(lifted, true)]);
        trim(b, whole, p(0.0), p(10.0), true)
    });
    refused(polyline(&graph, root), "measured in plan distance");

    // A closed composite trimmed against its sense.
    let closed = |start: Scalar, end: Scalar| {
        graph_of(|b| {
            let basis = b.push_value(circle2(Point2::ZERO, 1.0)).unwrap();
            let upper = trim(b, basis, p(0.0), p(PI), true);
            let lower = trim(b, basis, p(PI), p(2.0 * PI), true);
            let whole = composite(b, &[(upper, true), (lower, true)]);
            trim(b, whole, p(start), p(end), true)
        })
    };
    let (graph, root) = closed(3.0, 1.0);
    refused(polyline(&graph, root), "closing joint");
    let (graph, root) = closed(1.0, 3.0);
    polyline(&graph, root).unwrap();

    // An offset segment ends inside.
    let offset_segment = |start: Scalar, end: Scalar| {
        graph_of(|b| {
            let basis = segment2(b, (0.0, 0.0), (3.0, 0.0));
            let moved = offset(b, basis, 0.25, None);
            let whole = composite(b, &[(moved, true)]);
            trim(b, whole, p(start), p(end), true)
        })
    };
    let (graph, root) = offset_segment(1.0, 2.0);
    refused(polyline(&graph, root), "inside an offset segment");
    let (graph, root) = offset_segment(0.0, 3.0);
    let whole = polyline(&graph, root).unwrap();
    close(whole.points[0], Point3::new(0.0, 0.25, 0.0), "offset first");

    // An offset of a trimmed composite cut inside a quadrature.
    let (graph, root) = graph_of(|b| {
        let whole = segment_then_spline(b);
        let moved = offset(b, whole, 0.1, None);
        trim(b, moved, p(0.5), p(4.2), true)
    });
    refused(
        polyline(&graph, root),
        "cut inside a piece whose length is a quadrature",
    );

    // A composite no arc length runs along.
    let (graph, root) = graph_of(|b| {
        let source = b.push_value(circle2(Point2::ZERO, 1.0)).unwrap();
        let placed = b
            .push(GeometryNode::Instance(Instance {
                source,
                transform: Transform3::from_translation(Vec3::X),
            }))
            .unwrap();
        let whole = composite(b, &[(placed, true)]);
        trim(b, whole, p(0.0), p(1.0), true)
    });
    refused(polyline(&graph, root), "no arc length is measured along it");
}

#[test]
fn a_trim_of_a_full_turn_keeps_the_arc_it_names() {
    let circle = circle2(Point2::ZERO, 2.0);
    let on = |t: Scalar| lift(evaluate2(&circle, t).unwrap());
    for (from, to) in [(PI, 2.0 * PI), (1.5 * PI, 2.0 * PI), (0.0, PI)] {
        let (graph, root) = graph_of(|b| {
            let basis = b.push_value(circle.clone()).unwrap();
            let whole = trim(b, basis, p(0.0), p(2.0 * PI), true);
            trim(b, whole, p(from), p(to), true)
        });
        let outcome = polyline(&graph, root).unwrap();
        close(outcome.points[0], on(from), "first");
        close(*outcome.points.last().unwrap(), on(to), "last");
        within_bound(&outcome, |s| on(from + (to - from) * s), SAMPLES);
    }
}

#[test]
fn a_deep_chain_of_trims_is_compiled_in_linear_time() {
    let spline = spline2();
    let (graph, root) = graph_of(|b| {
        let mut node = b.push_value(spline.clone()).unwrap();
        for level in 0..60 {
            let margin = 0.001 * Scalar::from(level);
            node = trim(b, node, p(margin), p(1.0 - margin), true);
        }
        node
    });
    let started = std::time::Instant::now();
    let outcome = polyline(&graph, root).unwrap();
    assert!(started.elapsed().as_secs() < 20, "{:?}", started.elapsed());
    let on = |t: Scalar| lift(evaluate2(&spline, t).unwrap());
    close(outcome.points[0], on(0.059), "first");
}

#[test]
fn a_window_past_the_inner_span_is_refused_in_parameter_units() {
    let spline = spline2();
    let nested = |from: Scalar| {
        graph_of(|b| {
            let basis = b.push_value(spline.clone()).unwrap();
            let inner = trim(b, basis, p(0.1), p(0.9), true);
            trim(b, inner, p(from), p(0.6), true)
        })
    };
    // Half a thousandth of a parameter past the inner span is not rounding.
    let (graph, root) = nested(0.0995);
    invalid(polyline(&graph, root));
    // A rounding residue is.
    let (graph, root) = nested(0.1 - 1e-15);
    let outcome = polyline(&graph, root).unwrap();
    close(
        outcome.points[0],
        lift(evaluate2(&spline, 0.1).unwrap()),
        "first",
    );
}
