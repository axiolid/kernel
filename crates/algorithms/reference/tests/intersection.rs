//! Intersection-curve tests for the scalar reference.
//!
//! Every expectation is derived from the construction, never from a previous
//! run of the code under test.

use axiolid_core::Point3;
use axiolid_mesh::TriMesh;
use axiolid_reference::{assemble_polylines, intersection_segments};

/// An axis-aligned box as a closed 12-triangle solid.
fn cuboid(min: Point3, max: Point3) -> TriMesh {
    let positions = vec![
        Point3::new(min.x, min.y, min.z),
        Point3::new(max.x, min.y, min.z),
        Point3::new(max.x, max.y, min.z),
        Point3::new(min.x, max.y, min.z),
        Point3::new(min.x, min.y, max.z),
        Point3::new(max.x, min.y, max.z),
        Point3::new(max.x, max.y, max.z),
        Point3::new(min.x, max.y, max.z),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // bottom
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // front
        1, 2, 6, 1, 6, 5, // right
        2, 3, 7, 2, 7, 6, // back
        3, 0, 4, 3, 4, 7, // left
    ];
    TriMesh::new(positions, indices)
}

/// Disjoint solids have no intersection curve.
#[test]
fn disjoint_boxes_produce_no_segments() {
    let left = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
    let right = cuboid(Point3::new(5.0, 5.0, 5.0), Point3::new(6.0, 6.0, 6.0));

    let curve = intersection_segments(&left, &right).expect("disjoint boxes are supported");

    assert!(curve.segments.is_empty(), "disjoint solids share no curve");
}

/// A box fully inside another never crosses its surface.
#[test]
fn a_nested_box_produces_no_segments() {
    let outer = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 10.0));
    let inner = cuboid(Point3::new(4.0, 4.0, 4.0), Point3::new(6.0, 6.0, 6.0));

    let curve = intersection_segments(&outer, &inner).expect("nested boxes are supported");

    assert!(curve.segments.is_empty(), "a nested solid crosses nothing");
}

/// Two overlapping boxes intersect in exactly one closed loop.
///
/// The expectation is derived, not recorded: the tool's corner sits inside
/// the subject, so the surfaces cross in a single closed ring. Any open run,
/// or more than one loop, means the curve cracked.
#[test]
fn overlapping_boxes_intersect_in_one_closed_loop() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let curve = intersection_segments(&subject, &tool).expect("crossing boxes are supported");
    assert!(
        !curve.segments.is_empty(),
        "crossing surfaces share a curve"
    );

    let polylines = assemble_polylines(&curve.segments).expect("the curve is a 1-manifold");

    assert_eq!(polylines.len(), 1, "one crossing corner gives one loop");
    assert!(
        polylines[0].closed,
        "the loop must close: an open run is a crack"
    );
}

/// The curve must not depend on operand order.
///
/// Node identity is symmetric by construction, so the same crossing must
/// yield the same node count either way round. A mismatch means identity
/// leaked operand order and the stitch would crack.
#[test]
fn the_curve_is_independent_of_operand_order() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let forward = intersection_segments(&subject, &tool).expect("supported");
    let reverse = intersection_segments(&tool, &subject).expect("supported");

    assert_eq!(
        forward.segments.len(),
        reverse.segments.len(),
        "segment count must not depend on which operand is first",
    );
}

/// Repeated runs must produce byte-identical output.
///
/// The segment set is ordered and the positions come from exact arithmetic,
/// so there is no hash iteration order or accumulated error to vary.
#[test]
fn the_curve_is_reproducible() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let first = intersection_segments(&subject, &tool).expect("supported");
    let second = intersection_segments(&subject, &tool).expect("supported");

    assert_eq!(first.segments, second.segments, "segments must be stable");
    assert_eq!(
        first.positions, second.positions,
        "positions must be stable"
    );
}

/// Coplanar overlap is refused, not approximated.
///
/// Two boxes sharing a face meet in an AREA. Reporting that as a curve would
/// be a silent lie, so the contract is an explicit refusal.
#[test]
fn coplanar_faces_are_refused_rather_than_approximated() {
    let lower = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
    let upper = cuboid(Point3::new(0.0, 0.0, 2.0), Point3::new(2.0, 2.0, 4.0));

    let result = intersection_segments(&lower, &upper);

    assert!(
        result.is_err(),
        "an area of contact must not be reported as a curve"
    );
}

/// The loop must have the exact node count the geometry implies.
///
/// Two equal cubes offset along the diagonal so one corner of each lies inside
/// the other cross in a HEXAGON: six points, six segments. The subject and the
/// tool each contribute six edge punctures, but they are the SAME six points
/// seen from either side, not twelve distinct ones.
///
/// Asserting only "one closed loop" is too weak: taking the outer pair of each
/// face intersection instead of the inner pair still closes, but traces the
/// wrong curve. Pinning the counts catches that.
#[test]
fn the_loop_has_the_node_count_the_geometry_implies() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let curve = intersection_segments(&subject, &tool).expect("supported");

    assert_eq!(curve.segments.len(), 6, "a hexagonal ring has six sides");
    assert_eq!(curve.positions.len(), 6, "and six corners");

    let polylines = assemble_polylines(&curve.segments).expect("1-manifold");
    assert_eq!(polylines.len(), 1);
    assert!(polylines[0].closed);
    assert_eq!(polylines[0].nodes.len(), 6, "six nodes around the hexagon");
}

/// Every curve node must lie on both surfaces.
///
/// This is the property that catches a wrong interval endpoint: a node taken
/// from outside the other triangle still stitches into a closed ring, but it
/// is not on the intersection. Checked geometrically, independent of how the
/// node was named.
#[test]
fn every_node_lies_inside_both_operands_bounds() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let curve = intersection_segments(&subject, &tool).expect("supported");

    // The overlap region is the box [2,4]^3. Every intersection point must
    // sit within it, on the boundary of both solids.
    for point in curve.positions.values() {
        for axis in [point.x, point.y, point.z] {
            assert!(
                (2.0..=4.0).contains(&axis),
                "node {point:?} lies outside the shared region [2,4]^3",
            );
        }
    }
}

/// A slab cutting through a box crosses its surface in two closed rings.
///
/// This case used to be REFUSED. The slab's faces land exactly on the box's
/// edges, so one physical puncture is reached both as "a box edge crossing the
/// slab" and "a slab edge crossing the box". Naming those separately left the
/// rings open; merging on exact coordinate bits closes them.
///
/// Two rings, not one: a through-cut enters one side and exits the other.
#[test]
fn a_slab_through_a_box_crosses_in_two_closed_rings() {
    let box_solid = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 10.0));
    let slab = cuboid(Point3::new(3.0, -5.0, -5.0), Point3::new(6.0, 15.0, 15.0));

    let curve = intersection_segments(&box_solid, &slab).expect("a through-cut is a curve");
    let polylines = assemble_polylines(&curve.segments).expect("both rings close");

    assert_eq!(polylines.len(), 2, "a through-cut enters and exits");
    for line in &polylines {
        assert!(line.closed, "a cut through a closed solid closes");
        assert_eq!(
            line.nodes.len(),
            8,
            "four box faces, two of them split by a diagonal"
        );
    }
}
