//! Retriangulation tests for the scalar reference.
//!
//! Expectations come from the geometry of each fixture, not from a previous
//! run: a test that records whatever the code did cannot detect the code
//! being wrong.

use std::collections::BTreeMap;

use axiolid_core::Point3;
use axiolid_reference::intersection::{EdgeKey, IntersectionSegment, NodeKey, Operand, PointKey};
use axiolid_reference::retriangulate_face;

/// A node name for a test point. The specific edge is irrelevant here --
/// only that distinct points get distinct names.
fn node(seed: u32, point: Point3) -> NodeKey {
    NodeKey::EdgeSurface {
        edge: EdgeKey::new(Operand::Subject, seed, seed + 1),
        at: PointKey::new(point),
    }
}

/// Twice the area of a 3D triangle, used to compare areas exactly.
fn double_area(a: Point3, b: Point3, c: Point3) -> f64 {
    (b - a).cross(c - a).length()
}

/// An uncut face must come back unchanged.
#[test]
fn a_face_with_no_segments_is_returned_as_itself() {
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
    ];
    let patch = retriangulate_face(corners, &[], &BTreeMap::new()).expect("trivial patch");

    assert_eq!(
        patch.triangles,
        vec![[0, 1, 2]],
        "the face is its own triangulation"
    );
    assert_eq!(patch.points, corners.to_vec(), "no points are added");
    assert!(
        patch.sources.iter().all(Option::is_none),
        "no curve nodes exist"
    );
}

/// A cut spanning two edges splits the face into exactly three triangles.
///
/// Cutting a triangle with a chord that meets two different edges leaves a
/// triangle on one side and a quadrilateral on the other. The quad needs two
/// triangles, so three is the only possible count.
#[test]
fn a_chord_across_two_edges_splits_the_face_in_three() {
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
    ];
    // Midpoints of the two edges meeting at corner 0.
    let first = Point3::new(2.0, 0.0, 0.0);
    let second = Point3::new(0.0, 2.0, 0.0);
    let start = node(10, first);
    let end = node(20, second);

    let mut positions = BTreeMap::new();
    positions.insert(start, first);
    positions.insert(end, second);
    let segment = IntersectionSegment::between(start, end).expect("distinct nodes");

    let patch = retriangulate_face(corners, &[segment], &positions).expect("patch");

    assert_eq!(
        patch.triangles.len(),
        3,
        "a chord across two edges gives 3 triangles"
    );

    // Area is conserved: the pieces must tile the original face exactly.
    let whole = double_area(corners[0], corners[1], corners[2]);
    let parts: f64 = patch
        .triangles
        .iter()
        .map(|tri| {
            let [a, b, c] = tri.map(|i| patch.points[i as usize]);
            double_area(a, b, c)
        })
        .sum();
    assert!(
        (whole - parts).abs() < 1e-9,
        "pieces must tile the face: {whole} vs {parts}"
    );
}

/// The cut edge must actually appear in the output triangulation.
///
/// This is the whole point of a CONSTRAINED triangulation: a valid-looking
/// mesh that ignores the constraint would leave the two operands' surfaces
/// disagreeing about where they meet.
#[test]
fn the_cut_edge_survives_in_the_output() {
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
    ];
    let first = Point3::new(2.0, 0.0, 0.0);
    let second = Point3::new(0.0, 2.0, 0.0);
    let start = node(10, first);
    let end = node(20, second);
    let mut positions = BTreeMap::new();
    positions.insert(start, first);
    positions.insert(end, second);
    let segment = IntersectionSegment::between(start, end).expect("distinct");

    let patch = retriangulate_face(corners, &[segment], &positions).expect("patch");

    let start_index = patch
        .sources
        .iter()
        .position(|s| *s == Some(start))
        .expect("the start node is in the patch") as u32;
    let end_index = patch
        .sources
        .iter()
        .position(|s| *s == Some(end))
        .expect("the end node is in the patch") as u32;

    let present = patch.triangles.iter().any(|tri| {
        (0..3).any(|i| {
            let (a, b) = (tri[i], tri[(i + 1) % 3]);
            (a == start_index && b == end_index) || (a == end_index && b == start_index)
        })
    });
    assert!(
        present,
        "the constraint edge must be an edge of some triangle"
    );
}

/// Winding must be preserved: an inside-out patch inverts the solid.
///
/// Every output triangle must wind the same way as the face it replaces,
/// otherwise the rebuilt mesh has inconsistent normals and any downstream
/// volume or containment test silently flips sign.
#[test]
fn every_output_triangle_keeps_the_source_winding() {
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(0.0, 4.0, 0.0),
    ];
    let expected = (corners[1] - corners[0]).cross(corners[2] - corners[0]);

    let first = Point3::new(2.0, 0.0, 0.0);
    let second = Point3::new(0.0, 2.0, 0.0);
    let start = node(10, first);
    let end = node(20, second);
    let mut positions = BTreeMap::new();
    positions.insert(start, first);
    positions.insert(end, second);
    let segment = IntersectionSegment::between(start, end).expect("distinct");

    let patch = retriangulate_face(corners, &[segment], &positions).expect("patch");

    for tri in &patch.triangles {
        let [a, b, c] = tri.map(|i| patch.points[i as usize]);
        let normal = (b - a).cross(c - a);
        assert!(
            normal.dot(expected) > 0.0,
            "triangle {tri:?} winds against the source face"
        );
    }
}

/// A cut running edge-to-edge across a face far from any corner.
///
/// Area conservation is the load-bearing check: it fails for a triangulation
/// that overlaps itself, leaves a hole, or spills outside the face, none of
/// which a triangle-count assertion would catch on its own.
#[test]
fn a_cut_between_two_edges_tiles_the_face_exactly() {
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(6.0, 0.0, 0.0),
        Point3::new(0.0, 6.0, 0.0),
    ];
    // On edge 0->1 and on edge 1->2 respectively.
    let first = Point3::new(4.0, 0.0, 0.0);
    let second = Point3::new(2.0, 4.0, 0.0);
    let start = node(30, first);
    let end = node(40, second);
    let mut positions = BTreeMap::new();
    positions.insert(start, first);
    positions.insert(end, second);
    let segment = IntersectionSegment::between(start, end).expect("distinct");

    let patch = retriangulate_face(corners, &[segment], &positions).expect("patch");

    let whole = double_area(corners[0], corners[1], corners[2]);
    let parts: f64 = patch
        .triangles
        .iter()
        .map(|tri| {
            let [a, b, c] = tri.map(|i| patch.points[i as usize]);
            double_area(a, b, c)
        })
        .sum();
    assert!(
        (whole - parts).abs() < 1e-9,
        "pieces must tile the face exactly: {whole} vs {parts}"
    );
}

/// A cut ending strictly inside the face is tiled around that point.
///
/// I first expected this to be refused. It is not, and the output is
/// correct: the interior endpoint becomes a fan centre, the cut survives as
/// an edge, and the pieces still tile the face exactly. The expectation was
/// wrong, not the code -- recorded here so it is not 'fixed' back.
#[test]
fn a_cut_dangling_inside_the_face_fans_around_the_endpoint() {
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(6.0, 0.0, 0.0),
        Point3::new(0.0, 6.0, 0.0),
    ];
    let on_edge = Point3::new(3.0, 0.0, 0.0);
    let inside = Point3::new(2.0, 2.0, 0.0);
    let start = node(50, on_edge);
    let end = node(60, inside);
    let mut positions = BTreeMap::new();
    positions.insert(start, on_edge);
    positions.insert(end, inside);
    let segment = IntersectionSegment::between(start, end).expect("distinct");

    let patch = retriangulate_face(corners, &[segment], &positions).expect("patch");

    let whole = double_area(corners[0], corners[1], corners[2]);
    let parts: f64 = patch
        .triangles
        .iter()
        .map(|tri| {
            let [a, b, c] = tri.map(|i| patch.points[i as usize]);
            double_area(a, b, c)
        })
        .sum();
    assert!((whole - parts).abs() < 1e-9, "the fan must tile the face");

    let cut_present = patch.triangles.iter().any(|tri| {
        [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
            .iter()
            .any(|&(u, v)| (u.min(v), u.max(v)) == (3, 4))
    });
    assert!(cut_present, "the dangling cut must still appear as an edge");
}

/// Two crossing-adjacent cuts on one face must both survive.
///
/// With several constraint points on the boundary, a candidate triangle can
/// span across a cut. This is the case that exercises the
/// constraint-crossing filter and the coverage refusal: with two cuts, a
/// naive triangulation can satisfy one while erasing the other.
#[test]
fn two_cuts_on_one_face_both_survive() {
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(8.0, 0.0, 0.0),
        Point3::new(0.0, 8.0, 0.0),
    ];
    // Two chords, each spanning edge 0->1 to edge 0->2, nested.
    let a1 = Point3::new(2.0, 0.0, 0.0);
    let a2 = Point3::new(0.0, 2.0, 0.0);
    let b1 = Point3::new(5.0, 0.0, 0.0);
    let b2 = Point3::new(0.0, 5.0, 0.0);
    let (n1, n2, n3, n4) = (node(1, a1), node(2, a2), node(3, b1), node(4, b2));
    let mut positions = BTreeMap::new();
    for (k, v) in [(n1, a1), (n2, a2), (n3, b1), (n4, b2)] {
        positions.insert(k, v);
    }
    let cuts = [
        IntersectionSegment::between(n1, n2).expect("distinct"),
        IntersectionSegment::between(n3, n4).expect("distinct"),
    ];

    let patch = retriangulate_face(corners, &cuts, &positions).expect("patch");

    // Both cuts must appear as edges.
    for (u, v) in [(a1, a2), (b1, b2)] {
        let iu = patch.points.iter().position(|p| *p == u).expect("u kept") as u32;
        let iv = patch.points.iter().position(|p| *p == v).expect("v kept") as u32;
        let present = patch.triangles.iter().any(|tri| {
            [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
                .iter()
                .any(|&(x, y)| (x.min(y), x.max(y)) == (iu.min(iv), iu.max(iv)))
        });
        assert!(present, "cut {u:?}->{v:?} must survive");
    }

    let whole = double_area(corners[0], corners[1], corners[2]);
    let parts: f64 = patch
        .triangles
        .iter()
        .map(|tri| {
            let [x, y, z] = tri.map(|i| patch.points[i as usize]);
            double_area(x, y, z)
        })
        .sum();
    assert!((whole - parts).abs() < 1e-9, "pieces must tile the face");
}
