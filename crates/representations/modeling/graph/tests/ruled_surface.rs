//! A ruled surface joins two curve sides, or a curve and an apex (#402).

use axiolid_core::{Interval, Point3, Transform3, Vec2, Vec3};
use axiolid_curve::{Curve2, Curve3, Line2, Line3};
use axiolid_model::{
    CurveRelation, GeometryGraphBuilder, GeometryNode, GraphError, Instance, NodeId, RuledSide,
    SurfaceRelation, TrimSelector, TrimmingPreference,
};
use axiolid_surface::{Plane, Surface};

fn line(builder: &mut GeometryGraphBuilder, z: f64) -> NodeId {
    builder
        .push(GeometryNode::Curve3(Curve3::Line(Line3 {
            origin: Vec3::new(0.0, 0.0, z),
            direction: Vec3::X,
        })))
        .expect("a 3d line is a valid node")
}

fn side(curve: NodeId, start: f64, end: f64) -> RuledSide {
    RuledSide::Curve {
        curve,
        interval: Interval { start, end },
    }
}

fn ruled(
    builder: &mut GeometryGraphBuilder,
    first: RuledSide,
    second: RuledSide,
) -> Result<NodeId, GraphError> {
    builder.push(GeometryNode::SurfaceRelation(SurfaceRelation::Ruled {
        first,
        second,
    }))
}

fn detail(result: Result<NodeId, GraphError>) -> &'static str {
    match result {
        Err(GraphError::InvalidRuledSurface { detail }) => detail,
        other => panic!("expected an InvalidRuledSurface refusal, got {other:?}"),
    }
}

#[test]
fn two_curves_a_curve_and_an_apex_are_accepted() {
    let mut builder = GeometryGraphBuilder::new();
    let a = line(&mut builder, 0.0);
    let b = line(&mut builder, 1.0);
    let apex = RuledSide::Apex(Point3::new(0.0, 1.0, 0.0));
    ruled(&mut builder, side(a, 0.0, 1.0), side(b, 0.0, 1.0)).unwrap();
    ruled(&mut builder, side(a, 0.0, 1.0), apex).unwrap();
    ruled(&mut builder, apex, side(b, 0.0, 1.0)).unwrap();
}

#[test]
fn one_curve_run_over_different_or_reversed_intervals_is_accepted() {
    let mut builder = GeometryGraphBuilder::new();
    let a = line(&mut builder, 0.0);
    ruled(&mut builder, side(a, 0.0, 1.0), side(a, 1.0, 2.0)).unwrap();
    ruled(&mut builder, side(a, 0.0, 1.0), side(a, 1.0, 0.0)).unwrap();
}

#[test]
fn every_refusal_names_its_side() {
    let mut builder = GeometryGraphBuilder::new();
    let a = line(&mut builder, 0.0);
    let b = line(&mut builder, 1.0);
    let good = side(b, 0.0, 1.0);
    for (bad, first, second) in [
        (
            side(a, f64::NAN, 1.0),
            "the first side's interval is not finite",
            "the second side's interval is not finite",
        ),
        (
            side(a, 0.0, f64::INFINITY),
            "the first side's interval is not finite",
            "the second side's interval is not finite",
        ),
        (
            side(a, 0.5, 0.5),
            "the first side's interval is empty",
            "the second side's interval is empty",
        ),
        (
            RuledSide::Apex(Point3::new(f64::NAN, 0.0, 0.0)),
            "the first apex is not finite",
            "the second apex is not finite",
        ),
    ] {
        assert_eq!(detail(ruled(&mut builder, bad, good)), first);
        assert_eq!(detail(ruled(&mut builder, good, bad)), second);
    }
}

#[test]
fn the_first_side_is_judged_before_the_second() {
    let mut builder = GeometryGraphBuilder::new();
    let a = line(&mut builder, 0.0);
    let result = ruled(&mut builder, side(a, 0.5, 0.5), side(a, f64::NAN, 1.0));
    assert_eq!(detail(result), "the first side's interval is empty");
}

#[test]
fn two_apexes_and_one_curve_twice_are_refused() {
    let mut builder = GeometryGraphBuilder::new();
    let a = line(&mut builder, 0.0);
    let twin = line(&mut builder, 0.0);
    let apex = RuledSide::Apex(Point3::ZERO);
    assert_eq!(
        detail(ruled(&mut builder, apex, apex)),
        "both sides are apexes"
    );
    let same = "both sides run the same curve over the same interval";
    assert_eq!(
        detail(ruled(&mut builder, side(a, 0.0, 1.0), side(a, 0.0, 1.0))),
        same
    );
    // Equal values under two nodes are the same curve too.
    assert_eq!(
        detail(ruled(&mut builder, side(a, 0.0, 1.0), side(twin, 0.0, 1.0))),
        same
    );
    ruled(&mut builder, side(a, 0.0, 1.0), side(twin, 0.0, 2.0)).unwrap();
}

#[test]
fn a_side_must_be_an_atomic_curve3() {
    let mut builder = GeometryGraphBuilder::new();
    let a = line(&mut builder, 0.0);
    let flat = builder
        .push(GeometryNode::Curve2(Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::X,
        })))
        .unwrap();
    let trimmed = builder
        .push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis: a,
            start: vec![TrimSelector::Parameter(0.0)],
            end: vec![TrimSelector::Parameter(1.0)],
            sense_agreement: true,
            preference: TrimmingPreference::Parameter,
        }))
        .unwrap();
    let instance = builder
        .push(GeometryNode::Instance(Instance {
            source: a,
            transform: Transform3::IDENTITY,
        }))
        .unwrap();
    let plane = builder
        .push(GeometryNode::Surface(Surface::Plane(Plane {
            frame: axiolid_core::Frame3 {
                origin: Vec3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        })))
        .unwrap();
    for (bad, actual) in [
        (flat, "curve2"),
        (trimmed, "curve-relation"),
        (instance, "instance"),
        (plane, "surface"),
    ] {
        for (first, second) in [
            (side(bad, 0.0, 1.0), side(a, 0.0, 1.0)),
            (side(a, 0.0, 1.0), side(bad, 0.0, 1.0)),
        ] {
            assert_eq!(
                ruled(&mut builder, first, second),
                Err(GraphError::InvalidReferenceType {
                    reference: bad,
                    expected: "curve3 value",
                    actual,
                })
            );
        }
    }
}

#[test]
fn references_list_the_curve_sides_first_then_second() {
    let mut builder = GeometryGraphBuilder::new();
    let a = line(&mut builder, 0.0);
    let b = line(&mut builder, 1.0);
    let both = ruled(&mut builder, side(b, 0.0, 1.0), side(a, 0.0, 1.0)).unwrap();
    let apex_first = ruled(
        &mut builder,
        RuledSide::Apex(Point3::ZERO),
        side(b, 0.0, 1.0),
    )
    .unwrap();
    let apex_second = ruled(
        &mut builder,
        side(a, 0.0, 1.0),
        RuledSide::Apex(Point3::ZERO),
    )
    .unwrap();
    let graph = builder.finish(vec![both, apex_first, apex_second]).unwrap();
    assert_eq!(graph.get(both).unwrap().references(), vec![b, a]);
    assert_eq!(graph.get(apex_first).unwrap().references(), vec![b]);
    assert_eq!(graph.get(apex_second).unwrap().references(), vec![a]);
}

#[test]
fn the_refusal_reads_in_the_callers_terms() {
    let error = GraphError::InvalidRuledSurface {
        detail: "both sides are apexes",
    };
    assert_eq!(
        error.to_string(),
        "invalid ruled surface: both sides are apexes"
    );
}
