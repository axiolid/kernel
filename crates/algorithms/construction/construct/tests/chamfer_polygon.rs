//! Chamfering arbitrary polygon corners (ADR 0053).
//!
//! A chamfer replaces a corner with a flat. Geometry is checked against the
//! closed form -- the cut length for a setback `d` at interior angle `theta`
//! is `2 d sin(theta/2)` -- rather than against another call of the solver.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::feature::chamfer_polygon_corners;
use axiolid_core::{Point2, Tolerance};
use axiolid_surface::Surface;

fn square(half: f64) -> Vec<Point2> {
    vec![
        Point2::new(-half, -half),
        Point2::new(half, -half),
        Point2::new(half, half),
        Point2::new(-half, half),
    ]
}

/// A regular polygon, counter-clockwise.
fn regular(sides: usize, radius: f64) -> Vec<Point2> {
    (0..sides)
        .map(|index| {
            let angle = core::f64::consts::TAU * (index as f64) / (sides as f64);
            Point2::new(radius * angle.cos(), radius * angle.sin())
        })
        .collect()
}

#[test]
fn a_chamfer_replaces_the_corner_with_a_flat_of_the_closed_form_length() {
    // Hexagon interior angle is 120 deg, so a setback of d gives a cut of
    // 2 d sin(60 deg) = d*sqrt(3) -- a right-angle assumption would give
    // d*sqrt(2) instead.
    let sides = 6;
    let ring = regular(sides, 2.0);
    let distance = 0.3;
    let solid = chamfer_polygon_corners(&ring, &[(0, distance)], 1.0).expect("hexagon chamfer");

    let theta = (sides as f64 - 2.0) * core::f64::consts::PI / sides as f64;
    let expected = 2.0 * distance * (theta / 2.0).sin();
    assert!(
        (expected - distance * 3.0_f64.sqrt()).abs() < 1e-12,
        "closed form sanity"
    );

    // The two new vertices are the chamfer's ends; they must be `expected`
    // apart and each `distance` from the original corner.
    let corner = ring[0];
    let mut on_base: Vec<Point2> = Vec::new();
    for vertex in solid.topology().vertices() {
        if vertex.position.z.abs() > 1e-12 {
            continue;
        }
        let planar = Point2::new(vertex.position.x, vertex.position.y);
        if ((planar - corner).length() - distance).abs() < 1e-9 {
            on_base.push(planar);
        }
    }
    assert_eq!(
        on_base.len(),
        2,
        "a chamfer introduces exactly two vertices"
    );
    let cut = (on_base[0] - on_base[1]).length();
    assert!(
        (cut - expected).abs() < 1e-9,
        "chamfer flat should measure {expected}, got {cut}"
    );

    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
}

#[test]
fn every_chamfered_corner_adds_one_planar_wall() {
    let ring = square(2.0);
    let solid = chamfer_polygon_corners(&ring, &[(0, 0.4), (2, 0.7)], 1.0).expect("two chamfers");
    let planes = solid
        .surfaces()
        .iter()
        .filter(|surface| matches!(surface, Surface::Plane(_)))
        .count();
    // Four original walls plus two chamfer flats plus two caps.
    assert_eq!(planes, 8, "got {planes}");
    assert!(geometric_audit(&solid, Tolerance::METRE).is_consistent());
}

#[test]
fn distances_that_individually_fit_but_collide_together_are_refused() {
    // Edge length 4 on this square. 2.5 + 1.0 = 3.5 fits; 2.5 + 1.6 = 4.1
    // overruns, and no per-corner check can see it.
    let ring = square(2.0);
    assert!(chamfer_polygon_corners(&ring, &[(0, 2.5), (1, 1.0)], 1.0).is_ok());
    let error = chamfer_polygon_corners(&ring, &[(0, 2.5), (1, 1.6)], 1.0)
        .expect_err("the two chamfers cross on their shared edge");
    let text = format!("{error:?}");
    assert!(text.contains("too large for the edge"), "got {text}");
}

#[test]
fn a_reflex_corner_is_refused() {
    // L-shape with the inner corner at index 3.
    let ring = vec![
        Point2::new(0.0, 0.0),
        Point2::new(3.0, 0.0),
        Point2::new(3.0, 1.0),
        Point2::new(1.0, 1.0),
        Point2::new(1.0, 3.0),
        Point2::new(0.0, 3.0),
    ];
    assert!(chamfer_polygon_corners(&ring, &[(1, 0.3)], 1.0).is_ok());
    let error = chamfer_polygon_corners(&ring, &[(3, 0.3)], 1.0)
        .expect_err("a reflex corner cuts from outside the material");
    assert!(format!("{error:?}").contains("reflex"));
}
