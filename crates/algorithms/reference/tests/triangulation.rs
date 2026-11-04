//! Gate: a triangulation must cover exactly the input area, no more, no less.
//!
//! Area conservation is the triangulation-invariant check. Counting triangles
//! only proves a fan was emitted; summing their areas proves it covers the
//! polygon, and for holes it proves the holes were actually subtracted.

use axiolid_core::Point2;
use axiolid_reference::{signed_area2, triangulate_simple};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn tri_area(a: Point2, b: Point2, c: Point2) -> f64 {
    ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) / 2.0
}

fn covered_area(verts: &[Point2], tris: &[[u32; 3]]) -> f64 {
    tris.iter()
        .map(|t| {
            tri_area(
                verts[t[0] as usize],
                verts[t[1] as usize],
                verts[t[2] as usize],
            )
            .abs()
        })
        .sum()
}

#[test]
fn a_square_triangulates_into_two_triangles_covering_its_area() {
    let square = vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)];
    let tris = triangulate_simple(&square).expect("triangulate");
    assert_eq!(tris.len(), 2, "n-gon yields n-2 triangles");
    assert!((covered_area(&square, &tris) - 16.0).abs() < 1e-12);
}

/// An L-shape has a reflex vertex, which naive fan triangulation gets wrong.
#[test]
fn a_reflex_vertex_is_handled_without_covering_outside_area() {
    let l = vec![
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 4.0),
        p(0.0, 4.0),
    ];
    let tris = triangulate_simple(&l).expect("triangulate");
    assert_eq!(tris.len(), 4);
    let expected = signed_area2(&l) / 2.0;
    assert!((covered_area(&l, &tris) - expected).abs() < 1e-12);
    assert!((expected - 12.0).abs() < 1e-12, "L area is 12");
}

/// A reflex vertex can lie exactly on an otherwise plausible candidate ear's
/// diagonal. It must block that ear; treating only *strictly* interior points
/// as blockers emits outside area for this ordinary L profile.
#[test]
fn a_reflex_vertex_on_an_ear_diagonal_does_not_leak_outside_area() {
    let l = vec![
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(2.0, 1.0),
        p(1.0, 1.0),
        p(1.0, 2.0),
        p(0.0, 2.0),
    ];
    let tris = triangulate_simple(&l).expect("triangulate");
    assert_eq!(tris.len(), 4);
    assert!(tris.iter().all(|triangle| {
        tri_area(
            l[triangle[0] as usize],
            l[triangle[1] as usize],
            l[triangle[2] as usize],
        ) > 0.0
    }));
    assert!((covered_area(&l, &tris) - 3.0).abs() < 1e-12);
}

#[test]
fn a_degenerate_ring_is_refused_not_silently_empty() {
    assert!(triangulate_simple(&[p(0.0, 0.0), p(1.0, 1.0)]).is_err());
    let collinear = vec![p(0.0, 0.0), p(1.0, 0.0), p(2.0, 0.0)];
    assert!(triangulate_simple(&collinear).is_err());
}

/// A NaN or infinite coordinate is refused with a typed error: no orientation
/// predicate can decide it, and the triangulation must not panic or guess.
#[test]
fn non_finite_coordinates_are_refused_not_panicked() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let ring = [p(0.0, 0.0), p(4.0, 0.0), p(4.0, bad), p(0.0, 4.0)];
        match triangulate_simple(&ring) {
            Err(axiolid_contracts::GeomError::InvalidInput(_)) => {}
            other => panic!("expected InvalidInput for {bad}, got {other:?}"),
        }
    }
}
