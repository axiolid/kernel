//! Section edges between exact B-reps (#167, ADR 0075 stage 1).
//!
//! Oracles:
//! - every edge lies on BOTH faces' support surfaces (inverted and
//!   re-evaluated, sharing nothing with the trimming);
//! - the edges close into loops: the section of two closed boundaries has
//!   no loose ends;
//! - lengths and counts are closed forms from the inputs.

use axiolid_brep::ExactBRep;
use axiolid_brep_boolean::{section_edges, SectionEdge};
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, clip_arc_prism_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Plane3, Point2, Point3, Tolerance, Vec3};
use axiolid_curve::Curve3;
use axiolid_evaluate::evaluate3;
use axiolid_evaluate::surface::{evaluate, invert};
use axiolid_overlay::ArcRing;
use axiolid_primitive::HalfSpace;

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn prism(section: ArcRing, bottom: f64, top: f64) -> ArcPrism {
    ArcPrism {
        section,
        bottom,
        top,
    }
}

/// An exact solid over `section` between two heights.
fn solid(section: ArcRing, bottom: f64, top: f64) -> ExactBRep {
    // Intersecting a prism with itself returns it through the exact
    // column builder, faces and pcurves included.
    boolean_arc_prisms_exact(
        &prism(section.clone(), bottom, top),
        &prism(section, bottom, top),
        BooleanOperator::Intersection,
        tol(),
    )
    .expect("a solid")
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

/// Every sampled point of every edge lies on both faces' supports.
fn on_both(a: &ExactBRep, b: &ExactBRep, edges: &[SectionEdge]) {
    for edge in edges {
        for (brep, face) in [(a, edge.face_a), (b, edge.face_b)] {
            let surface = &brep.surfaces()[brep.topology().faces()[face.index()]
                .surface
                .unwrap()
                .index()];
            for i in 0..=20 {
                let t = edge.span.start + (edge.span.end - edge.span.start) * i as f64 / 20.0;
                let p = evaluate3(&edge.curve, t).unwrap();
                let (u, v) = invert(surface, p, tol()).unwrap();
                let back = evaluate(surface, u, v).unwrap();
                assert!((p - back).length() < 1e-9, "edge point {p:?} off its face");
            }
        }
    }
}

/// Every edge end meets another edge end (or the edge closes on itself).
fn closed(edges: &[SectionEdge]) {
    let ends: Vec<Point3> = edges.iter().flat_map(|e| [e.start, e.end]).collect();
    for (i, p) in ends.iter().enumerate() {
        let partner = ends
            .iter()
            .enumerate()
            .any(|(j, q)| j != i && (*p - *q).length() < 1e-9);
        assert!(partner, "a loose section end at {p:?}");
    }
}

fn length(edge: &SectionEdge) -> f64 {
    let n = 2000;
    (0..n)
        .map(|i| {
            let t0 = edge.span.start + (edge.span.end - edge.span.start) * i as f64 / n as f64;
            let t1 =
                edge.span.start + (edge.span.end - edge.span.start) * (i + 1) as f64 / n as f64;
            (evaluate3(&edge.curve, t1).unwrap() - evaluate3(&edge.curve, t0).unwrap()).length()
        })
        .sum()
}

#[test]
fn a_pipe_through_a_box_cuts_two_circles_split_at_the_seams() {
    let block = solid(square(-1.0, -1.0, 1.0, 1.0), 0.0, 2.0);
    let pipe = solid(ArcRing::circle(Point2::new(0.2, -0.1), 0.5), -1.0, 3.0);
    let edges = section_edges(&block, &pipe, tol()).expect("sections");
    on_both(&block, &pipe, &edges);
    closed(&edges);
    // The pipe's column is two half-walls with seams, so each circle is cut
    // there: one half-circle per (floor or roof, half-wall) pair.
    assert_eq!(edges.len(), 4, "{edges:#?}");
    for edge in &edges {
        assert!(matches!(edge.curve, Curve3::Circle(_)));
        assert!((length(edge) - PI * 0.5).abs() < 1e-5);
    }
}

#[test]
fn a_pipe_at_a_corner_cuts_one_loop_of_arcs_and_rulings() {
    // A pipe of radius 0.5 about (1, 0.9), next to the box's corner (1, 1).
    // Its column splits the circle into half-walls with seams at y = 0.9,
    // clear of the box's faces. Section: an arc on the floor and on the
    // roof where x <= 1 and y <= 1, and a ruling on each of the two side
    // faces near the corner.
    let block = solid(square(-1.0, -1.0, 1.0, 1.0), 0.0, 2.0);
    let (cx, cy, r) = (1.0, 0.9, 0.5);
    let pipe = solid(ArcRing::circle(Point2::new(cx, cy), r), -1.0, 3.0);
    let edges = section_edges(&block, &pipe, tol()).expect("sections");
    on_both(&block, &pipe, &edges);
    closed(&edges);
    let arcs: Vec<&SectionEdge> = edges
        .iter()
        .filter(|e| matches!(e.curve, Curve3::Circle(_)))
        .collect();
    let lines: Vec<&SectionEdge> = edges
        .iter()
        .filter(|e| matches!(e.curve, Curve3::Line(_)))
        .collect();
    // Over the box the circle runs from angle pi - asin(0.2) (where it
    // leaves y <= 1) to 3 pi / 2 (where it reaches x = 1); the column's
    // seam at angle pi splits that into two edges on each of floor and
    // roof.
    let arc_total: f64 = arcs.iter().map(|e| length(e)).sum();
    let expected = 2.0 * r * (PI / 2.0 + ((1.0 - cy) / r).asin());
    assert!(
        (arc_total - expected).abs() < 1e-5,
        "{arc_total} vs {expected}"
    );
    assert_eq!(lines.len(), 2, "{edges:#?}");
    for line in lines {
        assert!((length(line) - 2.0).abs() < 1e-9, "{}", length(line));
    }
}

#[test]
fn a_pipe_through_a_sloped_roof_cuts_an_ellipse() {
    // A 4 x 4 column under the roof z = 3 + 0.4 x - 0.2 y; a pipe of radius
    // 0.6 through it: an ellipse on the roof, a circle on the floor.
    let column = prism(square(-2.0, -2.0, 2.0, 2.0), 0.0, 10.0);
    let roof = HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, 3.0),
            normal: Vec3::new(-0.4, 0.2, 1.0),
        },
        agreement: false,
    };
    let block = clip_arc_prism_exact(&column, &roof, tol()).expect("a sloped block");
    let pipe = solid(ArcRing::circle(Point2::new(0.3, 0.2), 0.6), -1.0, 12.0);
    let edges = section_edges(&block, &pipe, tol()).expect("sections");
    on_both(&block, &pipe, &edges);
    closed(&edges);
    // Each is split at the pipe's two seams.
    assert_eq!(edges.len(), 4, "{edges:#?}");
    let ellipses = edges
        .iter()
        .filter(|e| matches!(e.curve, Curve3::Ellipse(_)))
        .count();
    let circles = edges
        .iter()
        .filter(|e| matches!(e.curve, Curve3::Circle(_)))
        .count();
    assert_eq!((ellipses, circles), (2, 2));
    let floor: f64 = edges
        .iter()
        .filter(|e| matches!(e.curve, Curve3::Circle(_)))
        .map(length)
        .sum();
    assert!((floor - 2.0 * PI * 0.6).abs() < 1e-5, "{floor}");
}

#[test]
fn apart_operands_have_no_section() {
    let a = solid(square(-1.0, -1.0, 1.0, 1.0), 0.0, 1.0);
    let b = solid(ArcRing::circle(Point2::new(5.0, 5.0), 0.5), 0.0, 1.0);
    // Their floors and roofs are coplanar: refused in this stage, by name.
    assert!(section_edges(&a, &b, tol()).is_err());
    let c = solid(ArcRing::circle(Point2::new(5.0, 5.0), 0.5), 2.0, 3.0);
    assert!(section_edges(&a, &c, tol()).expect("apart").is_empty());
}
