//! A bounded half-space's boundary is a 2D curve (#162).
//!
//! The boundary is authored in the placement frame's XY plane, and the
//! compiler reads it only as a `Curve2`. The graph used to accept any curve,
//! so a graph with a `Curve3` boundary validated and then could never
//! compile. It is refused when the graph is built.

use axiolid_core::{Plane3, Point2, Point3, Transform3, Vec3};
use axiolid_curve::{Curve2, Curve3, Polyline2, Polyline3};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, GraphError, Instance, SolidOperation};
use axiolid_primitive::HalfSpace;

fn square() -> Vec<(f64, f64)> {
    vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
}

fn bounded(boundary: GeometryNode, wrap: bool) -> Result<(), GraphError> {
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
    builder.finish(vec![op]).map(|_| ())
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
                assert_eq!(expected, "curve2");
            }
            other => panic!("wrap {wrap}: {other:?}"),
        }
    }
}
