//! Curve nodes compiled to polylines (#311).
//!
//! The reported chord bound is checked against the exact curve: each test
//! evaluates the curve densely (through `axiolid_reference::curve`, the
//! exact evaluation, and the placement computed here) and measures every
//! sample's distance to the polyline's segments. Every sample must lie
//! within the bound, and the bound within the requested budget. A circle's
//! bound is its exact sagitta, so there the densest sample must also come
//! close to it: a bound that is merely large would not pass.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve2, BSplineCurve3, Circle2, Circle3, CurvatureLaw, Curve2, Curve3, Ellipse2,
    Intrinsic2, KnotSpec, Line2, Line3, Polyline2, Polyline3,
};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{MeshCompiler, PolylineCompiler, PolylineOutcome};
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance,
    InstanceAtStation, NodeId, OpenProfile, Station, StationFrame, TrimSelector,
    TrimmingPreference,
};

const CHORD: Scalar = 1e-3;

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

// --- bounded by the exact curve ----------------------------------------------

#[test]
fn a_polyline_is_its_own_curve() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 1.0),
        Point3::new(4.0, 3.0, 1.0),
    ];
    let (graph, root) = graph_of(|b| {
        b.push_value(Curve3::Polyline(Polyline3 {
            points: points.clone(),
            closed: false,
        }))
        .unwrap()
    });
    let outcome = polyline(&graph, root).unwrap();
    assert_eq!(outcome.points, points);
    assert!(!outcome.closed);
    assert_eq!(outcome.deviation.bound, 0.0);
    // Along each straight leg, every point is on the polyline.
    within_bound(
        &outcome,
        |s| {
            let t = 2.0 * s;
            let (a, b, f) = if t <= 1.0 {
                (points[0], points[1], t)
            } else {
                (points[1], points[2], t - 1.0)
            };
            a + (b - a) * f
        },
        400,
    );

    // A closed 2D polyline ends on its first vertex, lifted to z = 0.
    let (graph, root) = graph_of(|b| {
        b.push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::ZERO, Point2::new(2.0, 0.0), Point2::new(0.0, 1.0)],
            closed: true,
        }))
        .unwrap()
    });
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.closed);
    assert_eq!(outcome.points.len(), 4);
    assert_eq!(outcome.points[3], Point3::ZERO);
    assert!(outcome.points.iter().all(|p| p.z == 0.0));
}

#[test]
fn a_circle_is_chorded_to_its_exact_sagitta() {
    let circle = Curve3::Circle(Circle3 {
        frame: tilted_frame(),
        radius: 2.5,
    });
    let (graph, root) = graph_of(|b| b.push_value(circle.clone()).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.closed, "a full circle closes");
    assert_eq!(outcome.points.first(), outcome.points.last());
    let worst = within_bound(
        &outcome,
        |s| axiolid_reference::curve::evaluate3(&circle, s * core::f64::consts::TAU).unwrap(),
        20_000,
    );
    // Equal spans of a circle: the sagitta of one is the bound, and the
    // densest sample lands within a hair of it.
    assert!(
        worst >= 0.99 * outcome.deviation.bound,
        "{worst} vs {}",
        outcome.deviation.bound
    );
    // The fewest chords that fit: half as many would not.
    let n = (outcome.points.len() - 1) as Scalar;
    let sagitta = |chords: Scalar| 2.5 * (1.0 - (core::f64::consts::PI / chords).cos());
    assert!(
        sagitta(n) <= CHORD && sagitta(n / 2.0) > CHORD,
        "{n} chords"
    );
}

fn cubic(control: Vec<Point3>, knots: Vec<Scalar>, multiplicities: Vec<u32>) -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree: 3,
        control_points: control,
        knots,
        multiplicities,
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

#[test]
fn a_bspline_is_chorded_within_its_certified_bound() {
    let control = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 3.0, 0.5),
        Point3::new(3.0, -2.0, 1.0),
        Point3::new(5.0, 4.0, -1.0),
        Point3::new(7.0, 0.0, 2.0),
        Point3::new(9.0, 1.0, 0.0),
    ];
    let spline = cubic(control.clone(), vec![0.0, 1.0, 2.5, 4.0], vec![4, 1, 1, 4]);
    let (graph, root) = graph_of(|b| b.push_value(spline.clone()).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    assert!(!outcome.closed);
    assert_eq!(outcome.points[0], control[0]);
    assert_eq!(*outcome.points.last().unwrap(), control[5]);
    let worst = within_bound(
        &outcome,
        |s| axiolid_reference::curve::evaluate3(&spline, 4.0 * s).unwrap(),
        20_000,
    );
    assert!(worst > 0.0);

    // A knot of full multiplicity is a corner: it is a vertex, and the
    // chords on both sides are bounded on their own.
    let mut seven = control;
    seven.push(Point3::new(10.0, -2.0, 1.0));
    let cornered = cubic(seven, vec![0.0, 1.0, 2.0], vec![4, 3, 4]);
    let (graph, root) = graph_of(|b| b.push_value(cornered.clone()).unwrap());
    let outcome = polyline(&graph, root).unwrap();
    let corner = axiolid_reference::curve::evaluate3(&cornered, 1.0).unwrap();
    assert!(
        outcome
            .points
            .iter()
            .any(|p| (*p - corner).length() <= 1e-12),
        "the corner {corner:?} is a vertex"
    );
    within_bound(
        &outcome,
        |s| axiolid_reference::curve::evaluate3(&cornered, 2.0 * s).unwrap(),
        20_000,
    );
}

#[test]
fn a_placed_2d_curve_is_moved_and_its_bound_grows_by_the_stretch() {
    let ellipse = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(1.0, -0.5),
            x: Vec2::new(0.6, 0.8),
            y: Vec2::new(-0.8, 0.6),
        },
        semi_axis_x: 1.5,
        semi_axis_y: 0.4,
    });
    let scale = 3.0;
    let transform = Transform3::from_translation(Vec3::new(10.0, -4.0, 2.0))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
        * Transform3::from_scale(Vec3::splat(scale));
    let (graph, root) = graph_of(|b| {
        let source = b.push_value(ellipse.clone()).unwrap();
        b.push(GeometryNode::Instance(Instance { source, transform }))
            .unwrap()
    });
    let outcome = polyline(&graph, root).unwrap();
    assert!(outcome.closed);
    let worst = within_bound(
        &outcome,
        |s| {
            let p =
                axiolid_reference::curve::evaluate2(&ellipse, s * core::f64::consts::TAU).unwrap();
            transform.transform_point3(Point3::new(p.x, p.y, 0.0))
        },
        20_000,
    );
    assert!(worst > 0.0);

    // Unplaced, the same ellipse lies in z = 0 and needs fewer chords:
    // its budget is not shrunk by the scale.
    let (graph, root) = graph_of(|b| b.push_value(ellipse.clone()).unwrap());
    let local = polyline(&graph, root).unwrap();
    assert!(local.points.iter().all(|p| p.z == 0.0));
    assert!(local.points.len() < outcome.points.len());
}

#[test]
fn an_open_profile_compiles_its_path() {
    let spline = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 3.0),
            Point2::new(4.0, -1.0),
            Point2::new(6.0, 2.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 1, 3],
        weights: Some(vec![1.0, 2.0, 0.5, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let (graph, root) = graph_of(|b| {
        let path = b.push_value(spline.clone()).unwrap();
        b.push_value(OpenProfile::new(path)).unwrap()
    });
    let outcome = polyline(&graph, root).unwrap();
    assert!(!outcome.closed);
    assert!(outcome.points.iter().all(|p| p.z == 0.0));
    within_bound(
        &outcome,
        |s| {
            let p = axiolid_reference::curve::evaluate2(&spline, s).unwrap();
            Point3::new(p.x, p.y, 0.0)
        },
        20_000,
    );
}

#[test]
fn a_curve_at_a_station_on_a_line_is_placed_exactly() {
    let circle = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 0.75,
    });
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::new(2.0, 1.0),
            direction: Vec2::new(0.0, 2.0),
        }))
        .unwrap();
    let source = b.push_value(circle.clone()).unwrap();
    let station = CurveStation {
        basis,
        station: Station::at(3.0),
        frame: StationFrame::Section,
    };
    let placed = b
        .push_value(InstanceAtStation::new(source, station.into()))
        .unwrap();
    let graph = b.finish(vec![placed]).unwrap();
    let placement = axiolid_mesh_compile::station::placement(&graph, placed).unwrap();
    assert!(placement.exact);
    let outcome = polyline(&graph, placed).unwrap();
    within_bound(
        &outcome,
        |s| {
            let p =
                axiolid_reference::curve::evaluate2(&circle, s * core::f64::consts::TAU).unwrap();
            placement
                .transform
                .transform_point3(Point3::new(p.x, p.y, 0.0))
        },
        20_000,
    );
}

// --- refused by name ----------------------------------------------------------

#[test]
fn a_curve_with_no_certified_chord_bound_is_refused_by_name() {
    let (graph, root) = graph_of(|b| {
        b.push_value(Curve2::Intrinsic(Intrinsic2::new(
            Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
            60.0,
        )))
        .unwrap()
    });
    refused(polyline(&graph, root), "intrinsic 2D curve");
}

#[test]
fn unbounded_relations_collections_and_solids_are_refused() {
    let (graph, root) = graph_of(|b| {
        b.push_value(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }))
        .unwrap()
    });
    refused(polyline(&graph, root), "unbounded 3D line");

    let (graph, root) = graph_of(|b| {
        let basis = b
            .push_value(Curve3::Circle(Circle3 {
                frame: tilted_frame(),
                radius: 1.0,
            }))
            .unwrap();
        b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start: vec![TrimSelector::Parameter(0.0)],
            end: vec![TrimSelector::Parameter(1.0)],
            sense_agreement: true,
            preference: TrimmingPreference::Parameter,
        }))
        .unwrap()
    });
    refused(polyline(&graph, root), "curve relation");

    let (graph, root) = graph_of(|b| {
        let a = b
            .push_value(Curve3::Circle(Circle3 {
                frame: tilted_frame(),
                radius: 1.0,
            }))
            .unwrap();
        b.push(GeometryNode::Collection(vec![a])).unwrap()
    });
    refused(polyline(&graph, root), "collection");

    let (graph, root) = graph_of(|b| {
        b.push(GeometryNode::Primitive(
            axiolid_primitive::Primitive::Sphere { radius: 1.0 },
        ))
        .unwrap()
    });
    assert!(matches!(
        polyline(&graph, root),
        Err(GeomError::InvalidInput(detail)) if detail.contains("not a curve")
    ));
}

#[test]
fn a_curve_at_a_numerically_resolved_station_is_refused() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 10.0,
        }))
        .unwrap();
    let source = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::ZERO, Point2::X],
            closed: false,
        }))
        .unwrap();
    let station = CurveStation {
        basis,
        station: Station::at(3.0),
        frame: StationFrame::Section,
    };
    let placed = b
        .push_value(InstanceAtStation::new(source, station.into()))
        .unwrap();
    let graph = b.finish(vec![placed]).unwrap();
    refused(polyline(&graph, placed), "numerically");
}

#[test]
fn compile_mesh_still_refuses_curves() {
    let (graph, root) = graph_of(|b| {
        b.push_value(Curve3::Circle(Circle3 {
            frame: tilted_frame(),
            radius: 1.0,
        }))
        .unwrap()
    });
    assert!(matches!(
        compiler().compile_mesh(&graph, root, &options()),
        Err(GeomError::Unsupported {
            operation: Operation::CurveEvaluation,
            ..
        })
    ));
}

// --- placed far from the origin (#295) ----------------------------------------

/// The issue's georeferenced placement.
const SITE: Vec3 = Vec3::new(6e6, 6e6, 100.0);

fn options_at(chord: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(chord, 1e-9).unwrap())
        .with_chord_error(chord)
        .expect("positive chord budget")
}

fn refused_below_rounding<T: core::fmt::Debug>(result: GeomResult<T>) {
    match result {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert_eq!(
                input,
                axiolid_mesh_compile::deviation::BUDGET_BELOW_PLACEMENT_ROUNDING
            );
        }
        other => panic!("expected the placement rounding refusal, got {other:?}"),
    }
}

/// A polyline is its own curve, so placed far out its whole bound is the
/// rounding of applying the transform: `4 eps` of the largest placed
/// coordinate's summed magnitude, never `0`.
#[test]
fn a_polyline_placed_far_out_reports_the_rounding_of_placing_it() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 1.0),
        Point3::new(4.0, 3.0, 1.0),
    ];
    let (graph, root) = graph_of(|b| {
        let source = b
            .push_value(Curve3::Polyline(Polyline3 {
                points: points.clone(),
                closed: false,
            }))
            .unwrap();
        b.push(GeometryNode::Instance(Instance {
            source,
            transform: Transform3::from_translation(SITE),
        }))
        .unwrap()
    });
    let outcome = compiler()
        .compile_polyline(&graph, root, &options_at(1e-7))
        .unwrap();
    let expected = 4.0 * Scalar::EPSILON * (6e6 + 4.0);
    assert_eq!(outcome.deviation.bound, expected);
    assert!(outcome.deviation.meets_requested());
    // Below the floor at the site (about 5.3e-9), refused by name.
    refused_below_rounding(compiler().compile_polyline(&graph, root, &options_at(1e-9)));
}

/// A small circle placed at the site with a budget a few times the
/// rounding: within the budget, and dense samples of the exact circle
/// against the f64 polyline moved back exactly (Sterbenz) within the bound.
#[test]
fn a_circle_placed_far_out_holds_its_bound_in_f64() {
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: 1e-3,
    });
    let (graph, root) = graph_of(|b| {
        let source = b.push_value(circle.clone()).unwrap();
        b.push(GeometryNode::Instance(Instance {
            source,
            transform: Transform3::from_translation(SITE),
        }))
        .unwrap()
    });
    let chord = 2e-8;
    let outcome = compiler()
        .compile_polyline(&graph, root, &options_at(chord))
        .unwrap();
    let bound = outcome.deviation.bound;
    assert!(outcome.deviation.meets_requested(), "{bound:e}");
    assert!(bound >= 4.0 * Scalar::EPSILON * 6e6, "{bound:e}");
    let back: Vec<Point3> = outcome
        .points
        .iter()
        .map(|&p| {
            let q = p - SITE;
            assert_eq!(q + SITE, p, "the subtraction is exact");
            q
        })
        .collect();
    let mut worst: Scalar = 0.0;
    for i in 0..20_000 {
        let s = (i as Scalar + 0.5) / 20_000.0 * core::f64::consts::TAU;
        let p = axiolid_reference::curve::evaluate3(&circle, s).unwrap();
        worst = worst.max(distance_to(&back, p));
    }
    assert!(worst <= bound, "{worst:e} over {bound:e}");
    assert!(worst > 0.0);
    refused_below_rounding(compiler().compile_polyline(&graph, root, &options_at(1e-9)));
}
