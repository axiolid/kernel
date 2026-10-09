//! A bounded half-space's boundary is a 2D curve (#162).
//!
//! The boundary is authored in the placement frame's XY plane, and the
//! compilers read it as a `Curve2` or a profile. The graph used to accept any curve,
//! so a graph with a `Curve3` boundary validated and then could never
//! compile. It is refused when the graph is built.
//!
//! A boundary of line and circular-arc segments is a profile node (#277):
//! accepted, and stored unchanged.

use axiolid_core::{Frame2, Interval, Plane3, Point2, Point3, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Line2, Polyline2, Polyline3};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, GraphError, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

fn square() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
}

fn bounded(boundary: GeometryNode, wrap: bool) -> Result<(), GraphError> {
    built(boundary, wrap).map(|_| ())
}

fn built(boundary: GeometryNode, wrap: bool) -> Result<(GeometryGraph, NodeId), GraphError> {
    let mut builder = GeometryGraphBuilder::new();
    let half_space = builder.push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: Vec3::ZERO,
            normal: Vec3::Z,
        },
        agreement: true,
    }))?;
    let mut curve = builder.push(boundary)?;
    if wrap {
        curve = builder.push(GeometryNode::Instance(Instance {
            source: curve,
            transform: Transform3::IDENTITY,
        }))?;
    }
    let op = builder.push(GeometryNode::SolidOperation(
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary: curve,
            placement: Transform3::IDENTITY,
        },
    ))?;
    builder.finish(vec![op]).map(|graph| (graph, op))
}

fn curve2() -> GeometryNode {
    GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
        points: square()
            .into_iter()
            .map(|(x, y)| Point2::new(x, y))
            .collect(),
        closed: true,
    }))
}

fn curve3() -> GeometryNode {
    GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
        points: square()
            .into_iter()
            .map(|(x, y)| Point3::new(x, y, 0.0))
            .collect(),
        closed: true,
    }))
}

#[test]
fn a_2d_boundary_is_accepted() {
    bounded(curve2(), false).expect("a Curve2 boundary");
    bounded(curve2(), true).expect("an instanced Curve2 boundary");
}

#[test]
fn a_3d_boundary_is_refused_at_graph_construction() {
    for wrap in [false, true] {
        match bounded(curve3(), wrap) {
            Err(GraphError::InvalidReferenceType { expected, .. }) => {
                assert_eq!(expected, "curve2 or profile");
            }
            other => panic!("wrap {wrap}: {other:?}"),
        }
    }
}

/// A unit square with its top edge bulged by a half circle: three lines
/// and one exact arc.
fn arched() -> GeometryNode {
    let line = |a: Point2, b: Point2| ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    };
    let arc = ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.5, 1.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 0.5,
        }),
        domain: Interval::new(0.0, std::f64::consts::PI),
        same_sense: true,
    };
    GeometryNode::Profile(Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            line(Point2::new(0.0, 1.0), Point2::new(0.0, 0.0)),
            line(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
            line(Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)),
            arc,
        ]),
        holes: Vec::new(),
    }))
}

#[test]
fn a_profile_boundary_with_arcs_is_accepted_and_round_trips() {
    let node = arched();
    for wrap in [false, true] {
        let (graph, op) = built(node.clone(), wrap).expect("a profile boundary");
        let Some(GeometryNode::SolidOperation(SolidOperation::BoundedHalfSpace {
            boundary, ..
        })) = graph.get(op)
        else {
            panic!("the bounded half-space comes back as itself");
        };
        let mut stored = graph.get(*boundary).expect("the boundary is present");
        if let GeometryNode::Instance(instance) = stored {
            assert!(wrap);
            stored = graph.get(instance.source).expect("the instanced profile");
        }
        assert_eq!(stored, &node, "the arcs survive storage unchanged");
    }
}

#[test]
fn a_node_that_is_neither_a_curve2_nor_a_profile_is_refused() {
    let point = GeometryNode::Point2(Point2::new(0.0, 0.0));
    match bounded(point, false) {
        Err(GraphError::InvalidReferenceType {
            expected, actual, ..
        }) => {
            assert_eq!(expected, "curve2 or profile");
            assert_eq!(actual, "point2");
        }
        other => panic!("{other:?}"),
    }
}
