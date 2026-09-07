//! Exact boolean assembly tests.
//!
//! The load-bearing assertion is VOLUME against a hand-computed value. A
//! boolean can produce a plausible triangle soup that looks correct and
//! encloses the wrong space; only measuring the enclosed volume catches that.

use axiolid_core::{BooleanOperator, Point3};
use axiolid_mesh::TriMesh;
use axiolid_reference::exact_boolean;

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
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Signed volume via the divergence theorem: sum of tetrahedra on the origin.
///
/// Correct only for a CLOSED, consistently-wound surface, which is exactly
/// what makes it a good check: an open or inside-out result gives a wrong
/// number rather than a plausible one.
fn volume(mesh: &TriMesh) -> f64 {
    let mut total = 0.0;
    for triangle in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[triangle[0] as usize];
        let b = mesh.positions[triangle[1] as usize];
        let c = mesh.positions[triangle[2] as usize];
        total += a.dot(b.cross(c)) / 6.0;
    }
    total
}

/// Subtracting a corner-overlapping cube removes exactly the shared volume.
///
/// Subject is the cube [0,4]^3, volume 64. The tool overlaps it in [2,4]^3,
/// a 2x2x2 cube of volume 8. The difference must therefore be 64 - 8 = 56.
/// That number comes from the geometry, not from running the code.
#[test]
fn subtracting_an_overlapping_corner_removes_exactly_that_volume() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let result = exact_boolean(&subject, &tool, BooleanOperator::Difference)
        .expect("interpenetrating solids are the case this handles");

    assert!(
        (volume(&result) - 56.0).abs() < 1e-9,
        "expected 64 - 8 = 56, got {}",
        volume(&result)
    );
}

/// Intersection of the same pair is the shared 2x2x2 cube.
#[test]
fn intersecting_two_cubes_keeps_only_the_shared_corner() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let result = exact_boolean(&subject, &tool, BooleanOperator::Intersection).expect("supported");

    assert!(
        (volume(&result) - 8.0).abs() < 1e-9,
        "the shared corner is 2x2x2 = 8, got {}",
        volume(&result)
    );
}

/// Union of the same pair is both cubes minus the double-counted overlap.
#[test]
fn uniting_two_cubes_counts_the_overlap_once() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));

    let result = exact_boolean(&subject, &tool, BooleanOperator::Union).expect("supported");

    // 64 + 64 - 8, the inclusion-exclusion identity.
    assert!(
        (volume(&result) - 120.0).abs() < 1e-9,
        "expected 64 + 64 - 8 = 120, got {}",
        volume(&result)
    );
}

/// A slab cut clean through a box removes a rectangular tunnel.
///
/// Box [0,10]^3 has volume 1000. The slab spans x in [3,6] and covers the
/// box entirely in y and z, so it removes a 3x10x10 = 300 slice, leaving
/// 700. This is the case that needed the coincident-puncture fix, so it also
/// guards that the fix survives.
#[test]
fn cutting_a_slab_through_a_box_removes_the_slice() {
    let box_solid = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 10.0));
    let slab = cuboid(Point3::new(3.0, -5.0, -5.0), Point3::new(6.0, 15.0, 15.0));

    let result = exact_boolean(&box_solid, &slab, BooleanOperator::Difference).expect("supported");

    assert!(
        (volume(&result) - 700.0).abs() < 1e-9,
        "expected 1000 - 300 = 700, got {}",
        volume(&result)
    );
}

/// The result must be a closed surface: every edge shared by two triangles.
///
/// A boolean that leaves a crack still reports a plausible volume, because
/// the divergence sum does not check closure. This does.
#[test]
fn the_result_is_a_closed_surface() {
    use std::collections::BTreeMap;

    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(2.0, 2.0, 2.0), Point3::new(6.0, 6.0, 6.0));
    let result = exact_boolean(&subject, &tool, BooleanOperator::Difference).expect("supported");

    let mut edge_use: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for triangle in result.indices.chunks_exact(3) {
        for i in 0..3 {
            let (u, v) = (triangle[i], triangle[(i + 1) % 3]);
            *edge_use.entry((u.min(v), u.max(v))).or_insert(0) += 1;
        }
    }

    let open: Vec<_> = edge_use.iter().filter(|(_, n)| **n != 2).collect();
    assert!(
        open.is_empty(),
        "a closed surface uses every edge exactly twice; {} edges do not",
        open.len()
    );
}

/// Closure must hold for the through-cut too, not just the corner overlap.
///
/// The slab case exercises a different topology: two separate rings rather
/// than one, and T-junctions on four faces instead of one.
#[test]
fn the_through_cut_result_is_also_closed() {
    use std::collections::BTreeMap;

    let box_solid = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 10.0, 10.0));
    let slab = cuboid(Point3::new(3.0, -5.0, -5.0), Point3::new(6.0, 15.0, 15.0));
    let result = exact_boolean(&box_solid, &slab, BooleanOperator::Difference).expect("supported");

    let mut edge_use: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for triangle in result.indices.chunks_exact(3) {
        for i in 0..3 {
            let (u, v) = (triangle[i], triangle[(i + 1) % 3]);
            *edge_use.entry((u.min(v), u.max(v))).or_insert(0) += 1;
        }
    }
    let open = edge_use.values().filter(|n| **n != 2).count();
    assert_eq!(open, 0, "{open} edges are not shared by exactly two faces");
}

/// Subtracting a disjoint tool leaves the subject's volume untouched.
///
/// Guards the boundary between this path and `ScalarBoolean`: a solid that
/// does not interpenetrate must still come back whole, not mangled.
#[test]
fn subtracting_a_disjoint_tool_changes_nothing() {
    let subject = cuboid(Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 4.0));
    let tool = cuboid(Point3::new(10.0, 10.0, 10.0), Point3::new(12.0, 12.0, 12.0));

    let result = exact_boolean(&subject, &tool, BooleanOperator::Difference).expect("supported");

    assert!(
        (volume(&result) - 64.0).abs() < 1e-9,
        "a disjoint subtraction keeps all 64, got {}",
        volume(&result)
    );
}
