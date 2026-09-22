// SPDX-License-Identifier: MPL-2.0

//! Behavioural tests for the editable planar arrangement.

use axiolid_arrangement::{Arrangement, BuildError, EditError, FaceId};
use axiolid_core::Point2;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn unit_square() -> Arrangement {
    Arrangement::from_polygon(&[p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)])
        .expect("a unit square is a valid polygon")
}

#[test]
fn a_square_builds_one_bounded_face_plus_the_outer_one() {
    let arrangement = unit_square();
    assert_eq!(arrangement.face_count(), 2);
    assert_eq!(arrangement.vertex_count(), 4);
    assert_eq!(arrangement.halfedge_count(), 8);
    assert_eq!(arrangement.bounded_faces().count(), 1);
    assert!(arrangement.audit().is_sound(), "{:?}", arrangement.audit());
}

#[test]
fn face_area_matches_the_polygon() {
    let arrangement = unit_square();
    let face = arrangement.bounded_faces().next().expect("one face");
    assert!((arrangement.face_area(face) - 1.0).abs() < 1e-12);
}

/// The unbounded face reports zero area rather than a negative number that
/// would poison a sum over "all faces".
#[test]
fn the_outer_face_has_no_area() {
    let arrangement = unit_square();
    assert_eq!(arrangement.face_area(arrangement.outer_face()), 0.0);
}

#[test]
fn clockwise_input_is_normalised_not_rejected() {
    let clockwise =
        Arrangement::from_polygon(&[p(0.0, 0.0), p(0.0, 1.0), p(1.0, 1.0), p(1.0, 0.0)])
            .expect("a clockwise square is still a square");
    let face = clockwise.bounded_faces().next().expect("one face");
    assert!(
        clockwise.face_area(face) > 0.0,
        "the bounded face should end up counter-clockwise regardless of input winding"
    );
}

#[test]
fn degenerate_boundaries_are_refused() {
    assert_eq!(
        Arrangement::from_polygon(&[p(0.0, 0.0), p(1.0, 0.0)]).unwrap_err(),
        BuildError::DegenerateBoundary
    );
    assert_eq!(
        Arrangement::from_polygon(&[p(0.0, 0.0), p(1.0, 0.0), p(2.0, 0.0)]).unwrap_err(),
        BuildError::ZeroArea
    );
}

/// The point of the crate: an edit does not invalidate handles.
///
/// `overlay` cannot do this. Its output polygons are fresh values with no
/// relationship to the previous call's, so "the same face" is not expressible.
#[test]
fn handles_survive_an_edit() {
    let mut arrangement = unit_square();
    let face = arrangement.bounded_faces().next().expect("one face");
    let area_before = arrangement.face_area(face);

    let edge = arrangement.face_halfedges(face)[0];
    let new_vertex = arrangement
        .split_edge(edge, p(0.5, 0.0))
        .expect("the midpoint lies on the bottom edge");

    // Same handle, still valid, still the same region.
    let area_after = arrangement.face_area(face);
    assert!(
        (area_before - area_after).abs() < 1e-12,
        "splitting an edge must not change the area it bounds"
    );
    assert_eq!(arrangement.vertex_count(), 5);
    assert_eq!(arrangement.degree(new_vertex), 2);
    assert!(arrangement.audit().is_sound(), "{:?}", arrangement.audit());
}

#[test]
fn splitting_off_the_edge_is_refused() {
    let mut arrangement = unit_square();
    let face = arrangement.bounded_faces().next().expect("one face");
    let edge = arrangement.face_halfedges(face)[0];
    assert_eq!(
        arrangement.split_edge(edge, p(0.5, 0.5)),
        Err(EditError::PointNotOnEdge),
        "a point off the segment is not a split point"
    );
    // And a point collinear but beyond the endpoint.
    assert_eq!(
        arrangement.split_edge(edge, p(2.0, 0.0)),
        Err(EditError::PointNotOnEdge)
    );
}

#[test]
fn a_harmless_drag_is_applied() {
    let mut arrangement = unit_square();
    let face = arrangement.bounded_faces().next().expect("one face");
    let corner = arrangement.face_halfedges(face)[2];
    let vertex = arrangement.halfedge_origin(corner);

    arrangement
        .drag_vertex(vertex, p(1.5, 1.5))
        .expect("moving a corner outward keeps the face simple");
    assert!(arrangement.face_area(face) > 1.0);
    assert!(arrangement.audit().is_sound(), "{:?}", arrangement.audit());
}

/// A drag that would turn the face inside out is refused, and leaves the
/// arrangement exactly as it was.
#[test]
fn a_corrupting_drag_is_refused_and_rolled_back() {
    let mut arrangement = unit_square();
    let face = arrangement.bounded_faces().next().expect("one face");
    let area_before = arrangement.face_area(face);
    let corner = arrangement.face_halfedges(face)[2];
    let vertex = arrangement.halfedge_origin(corner);
    let position_before = arrangement.position(vertex);

    let result = arrangement.drag_vertex(vertex, p(-5.0, -5.0));
    assert_eq!(result, Err(EditError::WouldSelfIntersect));
    assert_eq!(
        arrangement.position(vertex),
        position_before,
        "a refused drag must not leave the vertex moved"
    );
    assert!((arrangement.face_area(face) - area_before).abs() < 1e-12);
    assert!(arrangement.audit().is_sound(), "{:?}", arrangement.audit());
}

#[test]
fn neighbour_across_an_outer_edge_is_the_unbounded_face() {
    let arrangement = unit_square();
    let face = arrangement.bounded_faces().next().expect("one face");
    let edge = arrangement.face_halfedges(face)[0];
    assert_eq!(arrangement.neighbour_across(edge), FaceId::OUTER);
}

#[test]
fn an_empty_arrangement_is_still_valid() {
    let arrangement = Arrangement::new();
    assert_eq!(arrangement.face_count(), 1);
    assert_eq!(arrangement.bounded_faces().count(), 0);
    assert!(arrangement
        .face_halfedges(arrangement.outer_face())
        .is_empty());
    assert!(arrangement.audit().is_sound());
}
