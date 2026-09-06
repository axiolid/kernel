//! Characterisation: `UnifyOrientation` where a degenerate triangle bridges faces.
//!
//! `unify_orientation` walks face adjacency. `EdgeAdjacency` skips triangles
//! with a repeated corner, so a degenerate triangle cannot act as a bridge
//! between two otherwise separate faces. This pins the CURRENT behaviour
//! before that structure is shared, so a change of result is visible.

use axiolid_core::{Point3, Tolerance};
use axiolid_heal::mesh::MeshHealer;
use axiolid_heal::{Repair, RepairAction, RepairPlan};
use axiolid_mesh::TriMesh;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

/// Two triangles sharing edge 1-2, the second wound the same way round it
/// (inconsistent), plus a degenerate triangle touching both.
///
/// Vertices 0..3 form the two real faces. Triangle 2 is degenerate: it
/// repeats corner 1, so it has no meaningful adjacency.
fn bridged_by_degenerate() -> TriMesh {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
    ];
    // 0,1,2 then 1,2,3 traverses shared edge 1->2 the SAME way: inconsistent.
    // 1,1,3 is degenerate.
    let indices = vec![0, 1, 2, 1, 2, 3, 1, 1, 3];
    TriMesh::new(positions, indices)
}

/// The degenerate triangle does not stop the two real faces being unified.
#[test]
fn unify_orientation_flips_across_a_degenerate_neighbour() {
    let mesh = bridged_by_degenerate();
    let plan = RepairPlan {
        actions: vec![RepairAction::UnifyOrientation],
    };
    let (out, report) = MeshHealer.repair(&mesh, &plan, tol()).expect("repair runs");

    assert!(
        report.applied.contains(&RepairAction::UnifyOrientation),
        "the inconsistent pair must be reported as fixed, got {report:?}"
    );
    // Triangle count is untouched: unify only rewinds, never deletes.
    assert_eq!(out.indices.len(), mesh.indices.len());
}

/// A consistent pair is left alone and reported as skipped.
#[test]
fn unify_orientation_leaves_a_consistent_pair_untouched() {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
    ];
    // 0,1,2 and 2,1,3 traverse shared edge 1-2 in OPPOSITE directions.
    let mesh = TriMesh::new(positions, vec![0, 1, 2, 2, 1, 3]);
    let plan = RepairPlan {
        actions: vec![RepairAction::UnifyOrientation],
    };
    let (out, report) = MeshHealer.repair(&mesh, &plan, tol()).expect("repair runs");

    assert!(
        report.skipped.contains(&RepairAction::UnifyOrientation),
        "nothing to fix must be reported as skipped, got {report:?}"
    );
    assert_eq!(out.indices, mesh.indices, "a consistent mesh must not move");
}
