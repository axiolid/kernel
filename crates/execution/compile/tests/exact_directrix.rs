//! `exact_directrix` reads a swept disk's directrix as the exact compiler
//! does (#230): one segment or one arc, everything else refused by name.

use std::f64::consts::PI;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame3, Point3, Tolerance, Vec3};
use axiolid_curve::{Circle3, Curve3, Line3, Polyline3};
use axiolid_mesh_compile::{exact_directrix, ExactDirectrix};
use axiolid_model::{GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId};

fn graph_of(curve: Curve3) -> (GeometryGraph, NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let id = builder
        .push(GeometryNode::Curve3(curve))
        .expect("a valid curve");
    (builder.finish(vec![id]).expect("a valid graph"), id)
}

fn read(curve: Curve3, range: Option<(f64, f64)>) -> Result<ExactDirectrix, GeomError> {
    let (graph, id) = graph_of(curve);
    exact_directrix(&graph, id, range, &ExecutionOptions::new(Tolerance::METRE))
}

fn polyline(points: Vec<Point3>) -> Curve3 {
    Curve3::Polyline(Polyline3 {
        points,
        closed: false,
    })
}

#[test]
fn a_line_over_a_range_is_the_segment_between_its_ends() {
    let line = Curve3::Line(Line3 {
        origin: Point3::new(1.0, 2.0, 3.0),
        direction: Vec3::new(0.0, 0.0, 2.0),
    });
    let ExactDirectrix::Segment(start, end) = read(line, Some((0.5, 1.5))).expect("a segment")
    else {
        panic!("a line reads as a segment");
    };
    assert!((start - Point3::new(1.0, 2.0, 4.0)).length() < 1e-12);
    assert!((end - Point3::new(1.0, 2.0, 6.0)).length() < 1e-12);
}

#[test]
fn a_two_point_polyline_is_a_segment() {
    let (a, b) = (Point3::ZERO, Point3::new(0.0, 0.0, 2.0));
    assert_eq!(
        read(polyline(vec![a, b]), None).expect("a segment"),
        ExactDirectrix::Segment(a, b)
    );
}

#[test]
fn a_circle_over_a_range_is_an_arc_with_that_span() {
    let circle = Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.2, 0.5, 3.0),
            x: Vec3::X,
            y: Vec3::Z,
            z: -Vec3::Y,
        },
        radius: 2.0,
    };
    let ExactDirectrix::Arc(read_circle, span) =
        read(Curve3::Circle(circle), Some((PI, 2.0 * PI - 0.2))).expect("an arc")
    else {
        panic!("a circle reads as an arc");
    };
    assert_eq!(read_circle.radius, circle.radius);
    assert!((span.length() - (PI - 0.2)).abs() < 1e-12);
}

#[test]
fn a_polyline_with_a_corner_is_refused_by_name() {
    let corner = polyline(vec![
        Point3::ZERO,
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ]);
    assert!(matches!(
        read(corner, None),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

#[test]
fn an_unbounded_line_without_a_range_is_refused_by_name() {
    let line = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    });
    assert!(matches!(
        read(line, None),
        Err(GeomError::UnsupportedInput { .. })
    ));
}

#[test]
fn an_empty_range_is_degenerate() {
    let line = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    });
    assert!(matches!(
        read(line, Some((1.0, 1.0))),
        Err(GeomError::Degenerate(_))
    ));
}
