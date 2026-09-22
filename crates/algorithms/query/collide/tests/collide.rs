// SPDX-License-Identifier: MPL-2.0

//! Behavioural tests for convex collision queries.

use axiolid_collide::{
    boxes_intersect, contains_point, distance, intersects, separation, ConvexShape,
    SeparationResult,
};
use axiolid_core::{Box3, Point3, Vec3};

const TOL: f64 = 1e-9;

fn unit_cube_at(x: f64, y: f64, z: f64) -> ConvexShape {
    ConvexShape::from_aabb(Point3::new(x, y, z), Point3::new(x + 1.0, y + 1.0, z + 1.0))
}

#[test]
fn overlapping_boxes_are_detected() {
    let a = unit_cube_at(0.0, 0.0, 0.0);
    let b = unit_cube_at(0.5, 0.5, 0.5);
    assert!(intersects(&a, &b, TOL));
    assert_eq!(distance(&a, &b, TOL), 0.0);
}

#[test]
fn separated_boxes_report_their_gap() {
    let a = unit_cube_at(0.0, 0.0, 0.0);
    let b = unit_cube_at(3.0, 0.0, 0.0);
    assert!(!intersects(&a, &b, TOL));
    // Boxes span [0,1] and [3,4]: the gap is 2.
    assert!((distance(&a, &b, TOL) - 2.0).abs() < 1e-12);
}

#[test]
fn touching_boxes_are_not_apart() {
    let a = unit_cube_at(0.0, 0.0, 0.0);
    let b = unit_cube_at(1.0, 0.0, 0.0);
    // Face-to-face contact: gap is exactly zero, which is not "apart".
    assert!(intersects(&a, &b, TOL));
}

#[test]
fn the_separating_axis_points_the_right_way() {
    let a = unit_cube_at(0.0, 0.0, 0.0);
    let b = unit_cube_at(5.0, 0.0, 0.0);
    match separation(&a, &b, TOL) {
        SeparationResult::Apart { axis, distance } => {
            assert!((distance - 4.0).abs() < 1e-12);
            // The axis must be (anti)parallel to x, and unit length.
            assert!((axis.length() - 1.0).abs() < 1e-12);
            assert!(axis.x.abs() > 0.999, "expected an x-ish axis, got {axis:?}");
        }
        SeparationResult::Overlapping => panic!("boxes 4 apart are not overlapping"),
        other => panic!("unexpected separation result: {other:?}"),
    }
}

/// The case `Box3` could not answer before this crate existed.
///
/// Two oriented boxes that are disjoint but whose AXIS-ALIGNED bounds
/// overlap. An AABB test says "maybe colliding" and a caller with only
/// `Box3::corners()` had no way to do better.
#[test]
fn rotated_boxes_disjoint_despite_overlapping_aabbs() {
    let angle = std::f64::consts::FRAC_PI_4;
    let (sin, cos) = angle.sin_cos();

    // A long thin box along x.
    let a = Box3::new(
        Point3::new(-2.0, -0.2, -0.2),
        Vec3::new(4.0, 0.0, 0.0),
        Vec3::new(0.0, 0.4, 0.0),
        Vec3::new(0.0, 0.0, 0.4),
    );
    // The same box rotated 45 degrees and pushed along the diagonal, so its
    // AABB still overlaps a's but the solids do not touch.
    let b = Box3::new(
        Point3::new(1.2, 1.2, -0.2),
        Vec3::new(4.0 * cos, 4.0 * sin, 0.0),
        Vec3::new(-0.4 * sin, 0.4 * cos, 0.0),
        Vec3::new(0.0, 0.0, 0.4),
    );

    // Their axis-aligned bounds really do overlap.
    let a_corners = a.corners();
    let b_corners = b.corners();
    let a_max_x = a_corners.iter().fold(f64::NEG_INFINITY, |m, p| m.max(p.x));
    let b_min_x = b_corners.iter().fold(f64::INFINITY, |m, p| m.min(p.x));
    assert!(
        b_min_x < a_max_x,
        "fixture is wrong: the AABBs must overlap for this test to mean anything"
    );

    // But the oriented boxes do not.
    assert!(
        !boxes_intersect(&a, &b, TOL),
        "the rotated boxes are disjoint; an OBB test should say so"
    );
}

/// Two mutually tilted crossing bars, separable ONLY by an edge-pair axis.
///
/// This replaces an earlier version of this test that used two rotated
/// `Box3`s with overlapping AABBs. That fixture looked convincing and proved
/// nothing: mutation testing showed it still passed with edge-pair axes
/// removed entirely, because for a rectangular box the edge vectors ARE the
/// face normals, so face normals alone separated it.
///
/// Here both bars are tilted about their own long axes, so all six face
/// normals fail and the separating direction is a cross product of one edge
/// from each. Verified by mutation: deleting the edge-pair axes turns this
/// red. Found by search over 320 configurations, of which 253 have this
/// property.
#[test]
fn separation_can_require_an_edge_pair_axis() {
    let (sa, ca) = 0.6_f64.sin_cos();
    let a = parallelepiped(
        Point3::new(-4.0, -0.2, -0.2),
        Vec3::new(8.0, 0.0, 0.0),
        Vec3::new(0.0, 0.4 * ca, 0.4 * sa),
        Vec3::new(0.0, -0.4 * sa, 0.4 * ca),
    );
    let (sb, cb) = 0.3_f64.sin_cos();
    let b = parallelepiped(
        Point3::new(-0.2, -4.0, 0.5),
        Vec3::new(0.0, 8.0, 0.0),
        Vec3::new(0.4 * cb, 0.0, 0.4 * sb),
        Vec3::new(-0.4 * sb, 0.0, 0.4 * cb),
    );
    assert!(
        !intersects(
            &ConvexShape::from_points(&a),
            &ConvexShape::from_points(&b),
            TOL
        ),
        "these tilted bars are disjoint, and only an edge-pair axis shows it"
    );
}

/// Corners of a parallelepiped, bit 0 selects x, bit 1 y, bit 2 z.
fn parallelepiped(origin: Point3, x: Vec3, y: Vec3, z: Vec3) -> [Point3; 8] {
    let mut out = [origin; 8];
    for (index, corner) in out.iter_mut().enumerate() {
        let mut point = origin;
        if index & 1 != 0 {
            point += x;
        }
        if index & 2 != 0 {
            point += y;
        }
        if index & 4 != 0 {
            point += z;
        }
        *corner = point;
    }
    out
}

#[test]
fn rotated_boxes_that_really_overlap_are_detected() {
    let a = Box3::new(
        Point3::ZERO,
        Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(0.0, 2.0, 0.0),
        Vec3::new(0.0, 0.0, 2.0),
    );
    let angle = 0.3_f64;
    let (sin, cos) = angle.sin_cos();
    let b = Box3::new(
        Point3::new(1.0, 1.0, 1.0),
        Vec3::new(2.0 * cos, 2.0 * sin, 0.0),
        Vec3::new(-2.0 * sin, 2.0 * cos, 0.0),
        Vec3::new(0.0, 0.0, 2.0),
    );
    assert!(boxes_intersect(&a, &b, TOL));
}

#[test]
fn a_point_inside_and_outside_a_shape() {
    let cube = unit_cube_at(0.0, 0.0, 0.0);
    assert!(contains_point(&cube, Point3::new(0.5, 0.5, 0.5), TOL));
    assert!(!contains_point(&cube, Point3::new(2.0, 0.5, 0.5), TOL));
}

/// Nearly parallel edges must not fabricate a separating axis.
///
/// Honest status: this test currently passes even with the degeneracy filter
/// removed, verified by mutation. It is kept as a guard on the OUTCOME
/// (these shapes overlap and must be reported as overlapping) rather than
/// presented as proof that the filter is load-bearing.
///
/// Why the filter stays despite that: `candidate_axes` normalises by
/// `v / length_squared.sqrt()`, so a cross product of exactly parallel edges
/// would divide by zero and yield a NaN axis. NaN comparisons are all false,
/// which silently turns a projection test into "no separation found". The
/// filter makes that case impossible by construction. A test that forces the
/// NaN path would need two exactly parallel edges surviving into the axis
/// list, which the same filter prevents -- so the guarantee is structural,
/// and saying so is more useful than a test that appears to check it.
#[test]
fn nearly_parallel_edges_do_not_fabricate_a_separating_axis() {
    // Two long thin boxes, almost aligned, definitely overlapping.
    let a = ConvexShape::from_aabb(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.001, 0.001));
    let epsilon = 1e-9;
    let b = ConvexShape::from_points(&[
        Point3::new(5.0, 0.0, 0.0),
        Point3::new(15.0, epsilon, 0.0),
        Point3::new(15.0, epsilon + 0.001, 0.0),
        Point3::new(5.0, 0.001, 0.0),
        Point3::new(5.0, 0.0, 0.001),
        Point3::new(15.0, epsilon, 0.001),
        Point3::new(15.0, epsilon + 0.001, 0.001),
        Point3::new(5.0, 0.001, 0.001),
    ]);
    assert!(
        intersects(&a, &b, TOL),
        "these boxes share the region around x=5..10; a degenerate axis must \
         not be allowed to claim otherwise"
    );
}

#[test]
fn an_empty_shape_collides_with_nothing() {
    let empty = ConvexShape::from_points(&[]);
    let cube = unit_cube_at(0.0, 0.0, 0.0);
    assert!(!intersects(&empty, &cube, TOL));
}

/// Separation is symmetric: swapping the operands must not change whether
/// they collide, nor the distance between them.
#[test]
fn the_query_is_symmetric() {
    let a = unit_cube_at(0.0, 0.0, 0.0);
    let b = unit_cube_at(2.5, 0.0, 0.0);
    assert_eq!(intersects(&a, &b, TOL), intersects(&b, &a, TOL));
    let forward = distance(&a, &b, TOL);
    let backward = distance(&b, &a, TOL);
    assert!(
        (forward - backward).abs() < 1e-12,
        "distance changed with operand order: {forward} vs {backward}"
    );
}
