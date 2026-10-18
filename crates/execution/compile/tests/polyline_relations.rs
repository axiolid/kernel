//! Trimmed, composite, offset and intrinsic curves compiled to polylines
//! (#317).
//!
//! Each reported chord bound is checked against the exact curve: the tests
//! write the curve here from the reference evaluator (`axiolid_reference::
//! curve`: a basis between its trim parameters, each segment of a
//! composite, an offset as `B + d u` from the basis's point and tangent, a
//! circle's offset as the concentric circle) and measure every dense
//! sample's distance to the polyline. Every sample must lie within the
//! bound, and the bound within the budget. Joints of a composite must be
//! one vertex, the earlier segment's end, bitwise.

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve2, BSplineCurve3, BankConvention, Banked3, CantLaw, CantPiece, Circle2,
    CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2, Intrinsic3, KnotSpec, Line3,
    Polyline2,
};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{PolylineCompiler, PolylineOutcome};
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance,
    NodeId, OpenProfile, Station, StationFrame, StationOffsets, Transition, TrimSelector,
    TrimmingPreference,
};
use axiolid_reference::curve::{derivative2, derivative3, evaluate2, evaluate3};

const CHORD: Scalar = 1e-3;

/// A joint's offset, small enough to stay within the tolerance.
const HAIR: Scalar = 2e-7;

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

/// No vertex repeated by its successor: every joint is one vertex.
fn no_repeated_vertex(points: &[Point3]) {
    for (i, pair) in points.windows(2).enumerate() {
        assert_ne!(pair[0], pair[1], "vertex {i} is repeated");
    }
}

/// How many vertices equal `p` bitwise.
fn count(points: &[Point3], p: Point3) -> usize {
    points.iter().filter(|q| **q == p).count()
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

fn tilted_frame() -> Frame3 {
    let z = Vec3::new(1.0, -2.0, 2.0).normalize();
    let x = Vec3::new(2.0, 1.0, 0.0).normalize();
    Frame3 {
        origin: Point3::new(5.0, -3.0, 1.5),
        x,
        y: z.cross(x),
        z,
    }
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

fn segment2(b: &mut GeometryGraphBuilder, from: Point2, to: Point2) -> NodeId {
    b.push_value(Curve2::Polyline(Polyline2 {
        points: vec![from, to],
        closed: false,
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

fn spline3() -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 3.0, 0.5),
            Point3::new(3.0, -2.0, 1.0),
            Point3::new(5.0, 4.0, -1.0),
            Point3::new(7.0, 0.0, 2.0),
            Point3::new(9.0, 1.0, 0.0),
        ],
        knots: vec![0.0, 1.0, 2.5, 4.0],
        multiplicities: vec![4, 1, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// The exact left offset of a 2D curve at its parameter `t`.
fn left_offset2(curve: &Curve2, distance: Scalar, t: Scalar) -> Point3 {
    let p = evaluate2(curve, t).unwrap();
    let d = derivative2(curve, t).unwrap();
    let left = Vec2::new(-d.y, d.x) / d.length();
    lift(p + left * distance)
}

/// The exact offset of a 3D curve along `V x B'` at its parameter `t`.
fn directed_offset3(curve: &Curve3, distance: Scalar, v: Vec3, t: Scalar) -> Point3 {
    let p = evaluate3(curve, t).unwrap();
    let w = v.normalize().cross(derivative3(curve, t).unwrap());
    p + w / w.length() * distance
}

fn clothoid(length: Scalar) -> Intrinsic2 {
    Intrinsic2::new(
        Frame2 {
            origin: Point2::new(3.0, -1.0),
            x: Vec2::new(0.6, 0.8),
            y: Vec2::new(-0.8, 0.6),
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 30.0, length),
        length,
    )
}

/// A clothoid into an arc, then an arc the other way: seams at 30 and 50.
fn joined_law() -> CurvatureLaw {
    CurvatureLaw::Piecewise {
        breaks: vec![30.0, 50.0],
        laws: vec![
            CurvatureLaw::clothoid(0.0, 1.0 / 25.0, 30.0),
            CurvatureLaw::Constant {
                curvature: 1.0 / 25.0,
            },
            CurvatureLaw::Constant {
                curvature: -1.0 / 10.0,
            },
        ],
    }
}

// --- trims ----------------------------------------------------------------------

#[test]
fn a_trim_of_a_circle_wraps_its_seam_by_its_sense() {
    let circle = circle2(Point2::new(1.0, 2.0), 2.0);
    let on = |t: Scalar| lift(evaluate2(&circle, t).unwrap());
    for sense in [true, false] {
        let (graph, root) = graph_of(|b| {
            let basis = b.push_value(circle.clone()).unwrap();
            trim(
                b,
                basis,
                TrimSelector::Parameter(5.5),
                TrimSelector::Parameter(1.0),
                sense,
            )
        });
        let outcome = polyline(&graph, root).unwrap();
        // With the sense: 5.5 forward across the seam to 1.0 + 2 pi;
        // against it: 5.5 backward to 1.0, the complementary arc.
        let (from, travel) = if sense {
            (5.5, 1.0 + TAU - 5.5)
        } else {
            (5.5, 1.0 - 5.5)
        };
        assert!((outcome.points[0] - on(5.5)).length() <= 1e-12);
        assert!((*outcome.points.last().unwrap() - on(1.0)).length() <= 1e-12);
        within_bound(&outcome, |s| on(from + travel * s), SAMPLES);
        // The other arc is not on the polyline.
        let elsewhere = on(from - 0.5 * (TAU - travel.abs()) * travel.signum());
        assert!(
            distance_to(&outcome.points, elsewhere) > 1.0,
            "sense {sense}"
        );
        let centre = Point3::new(1.0, 2.0, 0.0);
        for p in &outcome.points {
            assert!(((*p - centre).length() - 2.0).abs() <= 1e-12);
        }
    }
}

#[test]
fn a_trim_of_a_bspline_or_a_line_flattens_between_its_parameters() {
    let spline = spline3();
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(spline.clone()).unwrap();
        trim(
            b,
            basis,
            TrimSelector::Parameter(0.7),
            TrimSelector::Parameter(3.1),
            false,
        )
    });
    let outcome = polyline(&graph, root).unwrap();
    // Against the sense it runs from the end parameter back.
    assert_eq!(outcome.points[0], evaluate3(&spline, 3.1).unwrap());
    assert_eq!(
        *outcome.points.last().unwrap(),
        evaluate3(&spline, 0.7).unwrap()
    );
    let worst = within_bound(
        &outcome,
        |s| evaluate3(&spline, 0.7 + 2.4 * s).unwrap(),
        SAMPLES,
    );
    assert!(worst > 0.0);
    // The trimmed-off ends are not on it.
    assert!(distance_to(&outcome.points, evaluate3(&spline, 0.0).unwrap()) > 0.1);
    assert!(distance_to(&outcome.points, evaluate3(&spline, 4.0).unwrap()) > 0.1);

    // A trim past the basis's domain names no part of it.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(spline.clone()).unwrap();
        trim(
            b,
            basis,
            TrimSelector::Parameter(1.0),
            TrimSelector::Parameter(5.0),
            true,
        )
    });
    assert!(matches!(
        polyline(&graph, root),
        Err(GeomError::InvalidInput(detail)) if detail.contains("outside its basis")
    ));

    // A trim of a line is a segment, its own curve.
    let line = Curve3::Line(Line3 {
        origin: Point3::new(1.0, 2.0, 3.0),
        direction: Vec3::new(0.0, 2.0, 1.0),
    });
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(line.clone()).unwrap();
        trim(
            b,
            basis,
            TrimSelector::Parameter(-1.0),
            TrimSelector::Parameter(2.0),
            true,
        )
    });
    let outcome = polyline(&graph, root).unwrap();
    assert_eq!(
        outcome.points,
        vec![
            evaluate3(&line, -1.0).unwrap(),
            evaluate3(&line, 2.0).unwrap()
        ]
    );
    assert_eq!(outcome.deviation.bound, 0.0);
}

// --- composites -----------------------------------------------------------------

#[test]
fn a_composite_shares_its_joints_and_reports_its_worst_segment() {
    // A segment, a quarter arc tangent to it, and a segment traversed
    // against its own direction. The first segment ends a hair (2e-7)
    // above the arc's start, so the joint can be told apart.
    let arc = circle2(Point2::new(4.0, 2.0), 2.0);
    let mut arc_node = None;
    let (graph, root) = graph_of(|b| {
        let first = segment2(b, Point2::new(0.0, 0.0), Point2::new(4.0, HAIR));
        let basis = b.push_value(arc.clone()).unwrap();
        let bend = trim(
            b,
            basis,
            TrimSelector::Parameter(-FRAC_PI_2),
            TrimSelector::Parameter(0.0),
            true,
        );
        arc_node = Some(bend);
        let last = segment2(b, Point2::new(6.0, 5.0), Point2::new(6.0, 2.0));
        composite(b, &[(first, true), (bend, true), (last, false)])
    });
    let outcome = polyline(&graph, root).unwrap();
    no_repeated_vertex(&outcome.points);
    // The first joint is the first segment's end, once; the arc's own
    // start is not a vertex.
    let joint = Point3::new(4.0, HAIR, 0.0);
    assert_eq!(outcome.points[1], joint);
    assert_eq!(count(&outcome.points, joint), 1);
    let arc_start = lift(evaluate2(&arc, -FRAC_PI_2).unwrap());
    assert_ne!(arc_start, joint);
    assert_eq!(count(&outcome.points, arc_start), 0);
    // The second joint is the arc's end, once, then the reversed segment.
    let second = lift(evaluate2(&arc, 0.0).unwrap());
    assert_eq!(count(&outcome.points, second), 1);
    assert_eq!(*outcome.points.last().unwrap(), Point3::new(6.0, 5.0, 0.0));

    // The bound is the worst segment's (the arc's), plus the joint gap.
    let bend = polyline(&graph, arc_node.unwrap()).unwrap();
    let gap = (arc_start - joint).length();
    assert!(outcome.deviation.bound >= bend.deviation.bound + gap);
    assert!(outcome.deviation.bound <= bend.deviation.bound + gap + 1e-15);
    within_bound(
        &outcome,
        |s| {
            let t = 3.0 * s;
            if t <= 1.0 {
                Point3::new(4.0 * t, HAIR * t, 0.0)
            } else if t <= 2.0 {
                lift(evaluate2(&arc, -FRAC_PI_2 * (2.0 - t)).unwrap())
            } else {
                Point3::new(6.0, 2.0 + 3.0 * (t - 2.0), 0.0)
            }
        },
        SAMPLES,
    );
}

#[test]
fn a_closed_composite_ends_on_its_first_vertex() {
    // A stadium: two segments and two half circles.
    let (left, right) = (
        circle2(Point2::new(0.0, 1.0), 1.0),
        circle2(Point2::new(4.0, 1.0), 1.0),
    );
    let (graph, root) = graph_of(|b| {
        let bottom = segment2(b, Point2::new(0.0, 0.0), Point2::new(4.0, 0.0));
        let r = b.push_value(right.clone()).unwrap();
        let r = trim(
            b,
            r,
            TrimSelector::Parameter(-FRAC_PI_2),
            TrimSelector::Parameter(FRAC_PI_2),
            true,
        );
        let top = segment2(b, Point2::new(0.0, 2.0), Point2::new(4.0, 2.0));
        let l = b.push_value(left.clone()).unwrap();
        let l = trim(
            b,
            l,
            TrimSelector::Parameter(FRAC_PI_2),
            TrimSelector::Parameter(3.0 * FRAC_PI_2),
            true,
        );
        composite(b, &[(bottom, true), (r, true), (top, false), (l, true)])
    });
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.closed);
    assert_eq!(outcome.points.first(), outcome.points.last());
    no_repeated_vertex(&outcome.points);
    within_bound(
        &outcome,
        |s| {
            let t = 4.0 * s;
            if t <= 1.0 {
                Point3::new(4.0 * t, 0.0, 0.0)
            } else if t <= 2.0 {
                lift(evaluate2(&right, -FRAC_PI_2 + PI * (t - 1.0)).unwrap())
            } else if t <= 3.0 {
                Point3::new(4.0 * (3.0 - t), 2.0, 0.0)
            } else {
                lift(evaluate2(&left, FRAC_PI_2 + PI * (t - 3.0)).unwrap())
            }
        },
        SAMPLES,
    );
}

#[test]
fn a_composite_gap_within_tolerance_is_bounded_and_beyond_it_refused() {
    let gapped = |gap: Scalar| {
        graph_of(|b| {
            let first = segment2(b, Point2::new(0.0, 0.0), Point2::new(3.0, 0.0));
            let second = segment2(b, Point2::new(3.0, gap), Point2::new(3.0, 4.0));
            composite(b, &[(first, true), (second, true)])
        })
    };
    // Half a millimetre: the joint is the first segment's end, and the
    // second segment's first chord moved by the gap is within it.
    let (graph, root) = gapped(5e-4);
    let outcome = polyline(&graph, root).unwrap();
    assert_eq!(
        outcome.points,
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 0.0)
        ]
    );
    assert!(outcome.deviation.bound >= 5e-4, "{:?}", outcome.deviation);
    within_bound(
        &outcome,
        |s| {
            if s <= 0.5 {
                Point3::new(6.0 * s, 0.0, 0.0)
            } else {
                Point3::new(3.0, 5e-4 + (4.0 - 5e-4) * (2.0 * s - 1.0), 0.0)
            }
        },
        SAMPLES,
    );
    // Two millimetres, beyond the tolerance: no single polyline.
    let (graph, root) = gapped(2e-3);
    refused(polyline(&graph, root), "do not meet");
}

#[test]
fn a_composite_gap_leaves_room_in_the_segments_budget() {
    // A segment ending half a millimetre from an arc whose own chords, at
    // the full budget, leave less than that: the arc is compiled again
    // within the budget less the gap, so the whole stays within it.
    let arc = circle2(Point2::new(4.0, 2.0), 2.0);
    let (graph, root) = graph_of(|b| {
        let first = segment2(b, Point2::new(0.0, 0.0), Point2::new(4.0, 5e-4));
        let basis = b.push_value(arc.clone()).unwrap();
        let bend = trim(
            b,
            basis,
            TrimSelector::Parameter(-FRAC_PI_2),
            TrimSelector::Parameter(1.3),
            true,
        );
        composite(b, &[(first, true), (bend, true)])
    });
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.deviation.bound >= 5e-4);
    within_bound(
        &outcome,
        |s| {
            if s <= 0.5 {
                Point3::new(8.0 * s, 1e-3 * s, 0.0)
            } else {
                lift(evaluate2(&arc, -FRAC_PI_2 + (1.3 + FRAC_PI_2) * (2.0 * s - 1.0)).unwrap())
            }
        },
        SAMPLES,
    );
}

#[test]
fn an_open_profile_whose_path_is_a_composite_compiles_it() {
    let spline = spline2();
    let (graph, root) = graph_of(|b| {
        let lead = segment2(b, Point2::new(-2.0, 0.0), Point2::new(0.0, 0.0));
        let basis = b.push_value(spline.clone()).unwrap();
        let path = composite(b, &[(lead, true), (basis, true)]);
        b.push_value(OpenProfile::new(path)).unwrap()
    });
    let outcome = polyline(&graph, root).unwrap();
    no_repeated_vertex(&outcome.points);
    assert_eq!(count(&outcome.points, Point3::ZERO), 1);
    within_bound(
        &outcome,
        |s| {
            if s <= 0.5 {
                Point3::new(-2.0 + 4.0 * s, 0.0, 0.0)
            } else {
                lift(evaluate2(&spline, 2.0 * s - 1.0).unwrap())
            }
        },
        SAMPLES,
    );
}

// --- offsets --------------------------------------------------------------------

#[test]
fn an_offset_of_a_circle_is_the_concentric_circle_on_its_left() {
    let circle = circle2(Point2::new(1.0, -1.0), 2.0);
    // Anticlockwise, the left is inwards: radius 1.5.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        offset(b, basis, 0.5, None)
    });
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.closed);
    let centre = Point3::new(1.0, -1.0, 0.0);
    for p in &outcome.points {
        assert!(((*p - centre).length() - 1.5).abs() <= 1e-12);
    }
    within_bound(
        &outcome,
        |s| centre + Vec3::new((TAU * s).cos(), (TAU * s).sin(), 0.0) * 1.5,
        SAMPLES,
    );

    // A half circle trimmed against its sense runs clockwise, so its
    // left is outwards: radius 2.5.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        let half = trim(
            b,
            basis,
            TrimSelector::Parameter(PI),
            TrimSelector::Parameter(0.0),
            false,
        );
        offset(b, half, 0.5, None)
    });
    let outcome = polyline(&graph, root).unwrap();
    for p in &outcome.points {
        assert!(((*p - centre).length() - 2.5).abs() <= 1e-12);
    }
    assert!((outcome.points[0] - (centre + Vec3::new(-2.5, 0.0, 0.0))).length() <= 1e-12);
    within_bound(
        &outcome,
        |s| centre + Vec3::new((PI * s).cos(), (PI * s).sin(), 0.0) * 2.5,
        SAMPLES,
    );
}

#[test]
fn an_offset_of_a_bspline_holds_its_bound_against_the_exact_offset() {
    let spline = spline2();
    for distance in [0.4, -0.25] {
        let (graph, root) = graph_of(|b| {
            let basis = b.push_value(spline.clone()).unwrap();
            offset(b, basis, distance, None)
        });
        let outcome = polyline(&graph, root).unwrap();
        assert_eq!(outcome.points[0], left_offset2(&spline, distance, 0.0));
        let worst = within_bound(&outcome, |s| left_offset2(&spline, distance, s), SAMPLES);
        assert!(worst > 0.0);
    }

    // A double knot leaves a cubic only C^1, so its offset only C^0 in
    // the parameter: it is cut there, and each side bounded on its own.
    let kinked = Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(3.0, 3.0),
            Point2::new(5.0, 2.0),
            Point2::new(6.0, 4.0),
            Point2::new(8.0, 3.0),
        ],
        knots: vec![0.0, 1.3, 4.0],
        multiplicities: vec![4, 2, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(kinked.clone()).unwrap();
        offset(b, basis, 0.3, None)
    });
    let outcome = polyline(&graph, root).unwrap();
    assert_eq!(count(&outcome.points, left_offset2(&kinked, 0.3, 1.3)), 1);
    within_bound(&outcome, |s| left_offset2(&kinked, 0.3, 4.0 * s), SAMPLES);

    // In 3D, along V x B'.
    let spline = spline3();
    let v = Vec3::new(0.2, -0.1, 1.0);
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(spline.clone()).unwrap();
        offset(b, basis, 0.3, Some(v))
    });
    let outcome = polyline(&graph, root).unwrap();
    within_bound(
        &outcome,
        |s| directed_offset3(&spline, 0.3, v, 4.0 * s),
        SAMPLES,
    );
}

#[test]
fn an_offset_of_a_composite_offsets_each_piece_and_shares_its_joints() {
    // A segment and a tangent quarter arc; offset half a unit to the left.
    let arc = circle2(Point2::new(4.0, 2.0), 2.0);
    let (graph, root) = graph_of(|b| {
        let first = segment2(b, Point2::new(0.0, 0.0), Point2::new(4.0, 0.0));
        let basis = b.push_value(arc.clone()).unwrap();
        let bend = trim(
            b,
            basis,
            TrimSelector::Parameter(-FRAC_PI_2),
            TrimSelector::Parameter(0.0),
            true,
        );
        // Traversed upwards, against its own direction: its left is -x.
        let up = segment2(b, Point2::new(6.0, 5.0), Point2::new(6.0, 2.0));
        let path = composite(b, &[(first, true), (bend, true), (up, false)]);
        offset(b, path, 0.5, None)
    });
    let outcome = polyline(&graph, root).unwrap();
    no_repeated_vertex(&outcome.points);
    let joint = Point3::new(4.0, 0.5, 0.0);
    assert_eq!(outcome.points[1], joint);
    assert_eq!(count(&outcome.points, joint), 1);
    assert_eq!(*outcome.points.last().unwrap(), Point3::new(5.5, 5.0, 0.0));
    within_bound(
        &outcome,
        |s| {
            let t = 3.0 * s;
            if t <= 1.0 {
                Point3::new(4.0 * t, 0.5, 0.0)
            } else if t <= 2.0 {
                let angle = -FRAC_PI_2 * (2.0 - t);
                Point3::new(4.0, 2.0, 0.0) + Vec3::new(angle.cos(), angle.sin(), 0.0) * 1.5
            } else {
                Point3::new(5.5, 2.0 + 3.0 * (t - 2.0), 0.0)
            }
        },
        SAMPLES,
    );
}

// --- intrinsic curves -----------------------------------------------------------

/// Nodes and weights of the 8-point Gauss-Legendre rule on `[-1, 1]`.
const GAUSS: [(Scalar, Scalar); 8] = [
    (-0.960_289_856_497_536_2, 0.101_228_536_290_376_3),
    (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (-0.525_532_409_916_328_9, 0.313_706_645_877_887_3),
    (-0.183_434_642_495_649_8, 0.362_683_783_378_361_9),
    (0.183_434_642_495_649_8, 0.362_683_783_378_361_9),
    (0.525_532_409_916_328_9, 0.313_706_645_877_887_3),
    (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
    (0.960_289_856_497_536_2, 0.101_228_536_290_376_3),
];

/// An independent high-precision reference for an intrinsic curve, in
/// its start frame's coordinates: the unit tangent of the law's
/// closed-form heading integrated cumulatively by 8-point Gauss-Legendre
/// on panels a few millimetres wide, split at the law's seams (error far
/// below 1e-12), at `n + 1` arc lengths evenly from `from` to `to`. Each
/// entry is the point and the unit tangent there.
fn reference(curve: &Intrinsic2, from: Scalar, to: Scalar, n: usize) -> Vec<(Vec2, Vec2)> {
    let seams = curve.curvature.seams_within(curve.length);
    let heading = |s: Scalar| curve.heading_at(s).unwrap();
    let integrate = |a: Scalar, b: Scalar| -> Vec2 {
        let mut cuts = vec![a];
        cuts.extend(seams.iter().copied().filter(|&t| t > a && t < b));
        cuts.push(b);
        let mut out = Vec2::ZERO;
        for pair in cuts.windows(2) {
            let panels = ((pair[1] - pair[0]) / 0.005).ceil().max(1.0) as usize;
            let w = (pair[1] - pair[0]) / panels as Scalar;
            for k in 0..panels {
                let mid = pair[0] + (k as Scalar + 0.5) * w;
                for (x, weight) in GAUSS {
                    let theta = heading(mid + 0.5 * w * x);
                    out += Vec2::new(theta.cos(), theta.sin()) * (0.5 * w * weight);
                }
            }
        }
        out
    };
    let mut at = integrate(0.0, from);
    let mut out = Vec::with_capacity(n + 1);
    let mut last = from;
    for i in 0..=n {
        let s = from + (to - from) * i as Scalar / n as Scalar;
        at += integrate(last, s);
        last = s;
        let theta = heading(s);
        out.push((at, Vec2::new(theta.cos(), theta.sin())));
    }
    out
}

/// A local point or vector of a 2D start frame, placed.
fn place(frame: &Frame2, local: Vec2) -> Point2 {
    frame.origin + frame.x * local.x + frame.y * local.y
}

/// The reference sample nearest a fraction `s` of `n` samples.
fn nth(s: Scalar, n: usize) -> usize {
    (s * n as Scalar).round() as usize
}

/// The reference curve of `curve` from `from` to `to`, offset by
/// `distance` to its left, as `within_bound` samples it.
fn reference2(
    curve: &Intrinsic2,
    from: Scalar,
    to: Scalar,
    distance: Scalar,
) -> impl Fn(Scalar) -> Point3 {
    let table = reference(curve, from, to, SAMPLES);
    let frame = curve.start;
    move |s| {
        let (p, t) = table[nth(s, SAMPLES)];
        let tangent = frame.x * t.x + frame.y * t.y;
        let left = Vec2::new(-tangent.y, tangent.x) / tangent.length();
        lift(place(&frame, p) + left * distance)
    }
}

#[test]
fn an_intrinsic_curve_is_bounded_by_its_curvature_law() {
    // Constant curvature is an arc, in closed form.
    let arc = Intrinsic2::new(
        clothoid(40.0).start,
        CurvatureLaw::Constant {
            curvature: 1.0 / 12.0,
        },
        40.0,
    );
    let (graph, root) = graph_of(|b| b.push_value(Curve2::Intrinsic(arc.clone())).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    within_bound(
        &outcome,
        |s| {
            let a = 40.0 * s / 12.0;
            lift(place(
                &arc.start,
                Vec2::new(12.0 * a.sin(), 12.0 * (1.0 - a.cos())),
            ))
        },
        SAMPLES,
    );

    // A straight one has no chord error, but its vertices are integrated:
    // the bound carries their certified error, however small.
    let straight = Intrinsic2::new(
        clothoid(40.0).start,
        CurvatureLaw::Constant { curvature: 0.0 },
        40.0,
    );
    let (graph, root) = graph_of(|b| b.push_value(Curve2::Intrinsic(straight.clone())).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.deviation.bound > 0.0 && outcome.deviation.bound < 1e-9);
    within_bound(
        &outcome,
        |s| lift(place(&straight.start, Vec2::new(40.0 * s, 0.0))),
        SAMPLES,
    );

    let spiral = clothoid(60.0);
    let (graph, root) = graph_of(|b| b.push_value(Curve2::Intrinsic(spiral.clone())).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    let worst = within_bound(&outcome, reference2(&spiral, 0.0, 60.0, 0.0), SAMPLES);
    assert!(worst > 0.0);

    // Placed by a frame that doubles it, its chords bend twice as much.
    let scaled = Intrinsic2 {
        start: Frame2 {
            origin: Point2::new(3.0, -1.0),
            x: Vec2::new(1.2, 1.6),
            y: Vec2::new(-1.6, 1.2),
        },
        ..clothoid(60.0)
    };
    let (graph, root) = graph_of(|b| b.push_value(Curve2::Intrinsic(scaled.clone())).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    within_bound(&outcome, reference2(&scaled, 0.0, 60.0, 0.0), SAMPLES);

    // A piecewise law: a clothoid into an arc, its curvature continuous
    // but its rate jumping at the seam, and an arc against it.
    let joined = Intrinsic2::new(plane(), joined_law(), 70.0);
    let (graph, root) = graph_of(|b| b.push_value(Curve2::Intrinsic(joined.clone())).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    within_bound(&outcome, reference2(&joined, 0.0, 70.0, 0.0), SAMPLES);

    // A planar space curve, in a tilted plane.
    let frame = tilted_frame();
    let space = Curve3::Intrinsic(Intrinsic3::new(
        frame,
        joined_law(),
        CurvatureLaw::Constant { curvature: 0.0 },
        70.0,
    ));
    let (graph, root) = graph_of(|b| b.push_value(space.clone()).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    let table = reference(&joined, 0.0, 70.0, SAMPLES);
    within_bound(
        &outcome,
        |s| {
            let (p, _) = table[nth(s, SAMPLES)];
            frame.origin + frame.x * p.x + frame.y * p.y
        },
        SAMPLES,
    );
}

#[test]
fn a_trim_and_an_offset_of_a_clothoid_are_bounded() {
    let spiral = clothoid(60.0);
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(Curve2::Intrinsic(spiral.clone())).unwrap();
        let part = trim(
            b,
            basis,
            TrimSelector::Parameter(10.0),
            TrimSelector::Parameter(55.0),
            true,
        );
        offset(b, part, 1.5, None)
    });
    let outcome = polyline(&graph, root).unwrap();
    let exact = reference2(&spiral, 10.0, 55.0, 1.5);
    assert!((outcome.points[0] - exact(0.0)).length() <= outcome.deviation.bound);
    within_bound(&outcome, exact, SAMPLES);

    // Beside a law whose curvature jumps at its seams the offset's speed
    // jumps: it is cut at each seam.
    let joined = Intrinsic2::new(plane(), joined_law(), 70.0);
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(Curve2::Intrinsic(joined.clone())).unwrap();
        offset(b, basis, 1.5, None)
    });
    let outcome = polyline(&graph, root).unwrap();
    let exact = reference2(&joined, 0.0, 70.0, 1.5);
    for seam in [30.0, 50.0] {
        let (p, t) = reference(&joined, seam, seam, 1)[0];
        let at = lift(place(&joined.start, p) + Vec2::new(-t.y, t.x) * 1.5);
        assert!(
            // Within the share of the budget its position may be off.
            outcome
                .points
                .iter()
                .any(|p| (*p - at).length() <= CHORD / 8.0),
            "the seam at {seam} is a vertex"
        );
    }
    within_bound(&outcome, exact, SAMPLES);
}

// --- refused by name ------------------------------------------------------------

#[test]
fn the_remaining_relations_and_families_are_refused_by_name() {
    let circle = circle2(Point2::ZERO, 1.0);

    // An offset of an offset.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(circle.clone()).unwrap();
        let once = offset(b, basis, 0.1, None);
        offset(b, once, 0.1, None)
    });
    refused(polyline(&graph, root), "an offset of an offset");

    // An offset of a placed curve.
    let (graph, root) = graph_of(|b| {
        let source = b.push_value(circle.clone()).unwrap();
        let placed = b
            .push(GeometryNode::Instance(Instance {
                source,
                transform: Transform3::from_translation(Vec3::X),
            }))
            .unwrap();
        offset(b, placed, 0.1, None)
    });
    refused(polyline(&graph, root), "an offset of a placed curve");

    // A trim of a placed curve.
    let (graph, root) = graph_of(|b| {
        let source = b.push_value(circle.clone()).unwrap();
        let placed = b
            .push(GeometryNode::Instance(Instance {
                source,
                transform: Transform3::from_translation(Vec3::X),
            }))
            .unwrap();
        trim(
            b,
            placed,
            TrimSelector::Parameter(0.0),
            TrimSelector::Parameter(1.0),
            true,
        )
    });
    refused(
        polyline(&graph, root),
        "a trim of a curve relation or of a placed curve",
    );

    // An offset across a polyline's corner.
    let (graph, root) = graph_of(|b| {
        let basis = b
            .push_value(Curve2::Polyline(Polyline2 {
                points: vec![Point2::ZERO, Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)],
                closed: false,
            }))
            .unwrap();
        offset(b, basis, 0.1, None)
    });
    refused(polyline(&graph, root), "corner of its basis");

    // A 3D offset needs its reference direction.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(spline3()).unwrap();
        offset(b, basis, 0.1, None)
    });
    refused(polyline(&graph, root), "without a reference direction");

    // An offset by stations.
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(spline3()).unwrap();
        b.push(GeometryNode::CurveRelation(
            CurveRelation::OffsetByStations {
                basis,
                stations: vec![
                    Station::new(0.0, StationOffsets::new(0.5, 0.0, 0.0)),
                    Station::new(2.0, StationOffsets::new(0.5, 0.0, 0.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap()
    });
    refused(polyline(&graph, root), "an offset by stations");

    // A banked curve, alone or offset.
    let banked = Curve3::Banked(Banked3::new(
        Elevated3::new(
            Curve2::Intrinsic(clothoid(60.0)),
            ElevationLaw::constant_grade(60.0, 0.02),
        ),
        CantLaw::new(vec![CantPiece::constant(60.0, 0.1)]),
        CantLaw::new(vec![CantPiece::constant(60.0, 0.0)]),
        1.5,
        BankConvention::TangentRotation,
    ));
    let (graph, root) = graph_of(|b| b.push_value(banked.clone()).unwrap());
    refused(polyline(&graph, root), "banked 3D curve");
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(banked.clone()).unwrap();
        offset(b, basis, 0.1, Some(Vec3::Z))
    });
    refused(polyline(&graph, root), "banked 3D curve");

    // A space curve with torsion, alone or trimmed.
    let twisted = Curve3::Intrinsic(Intrinsic3::new(
        Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        CurvatureLaw::Constant { curvature: 0.5 },
        CurvatureLaw::Constant { curvature: 0.1 },
        10.0,
    ));
    let (graph, root) = graph_of(|b| {
        let basis = b.push_value(twisted.clone()).unwrap();
        trim(
            b,
            basis,
            TrimSelector::Parameter(1.0),
            TrimSelector::Parameter(4.0),
            true,
        )
    });
    refused(polyline(&graph, root), "intrinsic 3D curve with torsion");
}
