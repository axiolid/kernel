//! Trims of trimmed, offset and placed curves compiled to polylines (#424).
//!
//! A trim reads its basis's parameter (ISO 10303-42), an offset its basis's,
//! a placed curve its source's, as a sweep directrix does. Expected points
//! are the reference evaluator's own, written here, and every dense sample
//! of the exact curve must lie within the reported bound.

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{BSplineCurve2, Circle2, Curve2, KnotSpec};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{PolylineCompiler, PolylineOutcome};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId,
    TrimSelector, TrimmingPreference,
};
use axiolid_reference::curve::evaluate2;

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
