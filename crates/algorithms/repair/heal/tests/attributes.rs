//! Repairs keep attribute channels and normals in step with geometry (#114).
//!
//! Every case ends in `validate_structure`, because the original defect was a
//! mesh that the healer returned as "repaired" and that failed it.

use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_heal::mesh::MeshHealer;
use axiolid_heal::{Repair, RepairAction, RepairPlan};
use axiolid_mesh::{AttributeChannel, AttributeFate, Blend, DropReason, NormalAttribute, TriMesh};

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("valid tolerance")
}

fn run(mesh: &TriMesh, actions: Vec<RepairAction>) -> (TriMesh, axiolid_heal::RepairReport) {
    let (out, report) = MeshHealer
        .repair(mesh, &RepairPlan { actions }, tol())
        .expect("repair");
    out.validate_structure()
        .expect("a repaired mesh must stay structurally valid");
    (out, report)
}

/// Two triangles sharing an edge, with the shared edge's vertices
/// duplicated: 0..3 is the first triangle, 3..6 the second, and 3 and 4
/// sit on top of 1 and 2.
fn split_quad() -> TriMesh {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    TriMesh::new(
        vec![
            p(0.0, 0.0),
            p(1.0, 0.0),
            p(0.0, 1.0),
            p(1.0, 0.0),
            p(0.0, 1.0),
            p(1.0, 1.0),
        ],
        vec![0, 1, 2, 3, 5, 4],
    )
}

#[test]
fn weld_compacts_a_channel_with_the_positions_it_describes() {
    let mut mesh = split_quad();
    // Duplicates carry the SAME id as their twin, so nothing conflicts.
    mesh.attributes.push(AttributeChannel::new(
        "id",
        vec![10.0, 11.0, 12.0, 11.0, 12.0, 15.0],
        1,
        Blend::Nearest,
    ));
    let (out, report) = run(&mesh, vec![RepairAction::WeldVertices]);
    assert_eq!(out.positions.len(), 4);
    // Each surviving position keeps the value it had, including the one
    // after the merged pair -- the original bug shifted it.
    for (v, want) in [(0, 10.0), (1, 11.0), (2, 12.0), (3, 15.0)] {
        assert_eq!(out.attributes[0].get(v), Some(&[want][..]), "vertex {v}");
    }
    assert_eq!(
        report.attribute_fates,
        vec![("id".to_owned(), AttributeFate::Preserved)]
    );
}

#[test]
fn weld_across_a_seam_drops_the_channel_by_name() {
    let mut mesh = split_quad();
    // Vertex 3 sits on vertex 1 but carries a different UV: a seam.
    mesh.attributes.push(AttributeChannel::new(
        "uv",
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.0, 1.0, 1.0, 1.0],
        2,
        Blend::Linear,
    ));
    let (out, report) = run(&mesh, vec![RepairAction::WeldVertices]);
    assert!(
        out.attributes.is_empty(),
        "a smeared channel must not be returned"
    );
    assert_eq!(
        report.attribute_fates,
        vec![(
            "uv".to_owned(),
            AttributeFate::Dropped(DropReason::ConflictingValues)
        )]
    );
}

#[test]
fn weld_keeps_per_vertex_normals_in_step() {
    let mut mesh = split_quad();
    mesh.normals = Some(NormalAttribute {
        values: vec![Vec3::Z; 6],
        indices: None,
    });
    let (out, _) = run(&mesh, vec![RepairAction::WeldVertices]);
    assert_eq!(out.normals.expect("normals kept").values.len(), 4);
}

#[test]
fn weld_across_a_hard_edge_switches_normals_to_corner_indexed() {
    let mut mesh = split_quad();
    // The duplicates carry a different normal from their twins. One normal
    // per position cannot hold both, so the weld must not pick one.
    let values = vec![Vec3::Z, Vec3::Z, Vec3::Z, Vec3::X, Vec3::X, Vec3::X];
    mesh.normals = Some(NormalAttribute {
        values: values.clone(),
        indices: None,
    });
    let (out, _) = run(&mesh, vec![RepairAction::WeldVertices]);
    let normals = out.normals.expect("normals kept");
    let corners = normals.indices.as_ref().expect("corner-indexed after weld");
    // Every corner still resolves to the normal it had before the weld.
    for (corner, &was) in mesh.indices.iter().enumerate() {
        assert_eq!(
            normals.values[corners[corner] as usize], values[was as usize],
            "corner {corner}"
        );
    }
}

#[test]
fn dropping_a_degenerate_triangle_drops_its_corner_normals() {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let mut mesh = TriMesh::new(
        vec![p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0), p(2.0, 0.0)],
        // Middle triangle is collinear: zero area.
        vec![0, 1, 2, 0, 1, 3, 1, 3, 2],
    );
    mesh.normals = Some(NormalAttribute {
        values: vec![Vec3::X, Vec3::Y, Vec3::Z],
        // Corner normals tagged by triangle: 0 -> X, 1 -> Y, 2 -> Z.
        indices: Some(vec![0, 0, 0, 1, 1, 1, 2, 2, 2]),
    });
    let (out, _) = run(&mesh, vec![RepairAction::DropDegenerateElements]);
    let corners = out
        .normals
        .expect("normals kept")
        .indices
        .expect("still indexed");
    // The survivors keep THEIR normals, not the dropped triangle's.
    assert_eq!(corners, vec![0, 0, 0, 2, 2, 2]);
}

#[test]
fn flipping_a_triangle_flips_its_corner_normals_with_it() {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    // A closed tetrahedron wound inward, so OrientOutward flips everything.
    let mut mesh = TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ],
        vec![0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2],
    );
    let corner_ids: Vec<u32> = (0..12).collect();
    mesh.normals = Some(NormalAttribute {
        values: vec![Vec3::Z; 12],
        indices: Some(corner_ids),
    });
    let (out, report) = run(&mesh, vec![RepairAction::OrientOutward]);
    assert_eq!(report.applied, vec![RepairAction::OrientOutward]);
    let corners = out.normals.expect("normals kept").indices.expect("indexed");
    // Each corner entry moved with its corner: position index and normal
    // index still pair up exactly as they did before the flip.
    for (c, &original) in corners.iter().enumerate() {
        assert_eq!(
            out.indices[c], mesh.indices[original as usize],
            "corner {c}"
        );
    }
}

#[test]
fn unifying_orientation_flips_corner_normals_with_the_triangle() {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    // Two triangles sharing edge 1-2, both traversing it 1 -> 2: inconsistent.
    // Whichever one the fill flips, its corner normals must move with it.
    let mut mesh = TriMesh::new(
        vec![p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0), p(1.0, 1.0)],
        vec![0, 1, 2, 3, 1, 2],
    );
    mesh.normals = Some(NormalAttribute {
        values: vec![Vec3::Z; 6],
        indices: Some((0..6).collect()),
    });
    let (out, report) = run(&mesh, vec![RepairAction::UnifyOrientation]);
    assert_eq!(report.applied, vec![RepairAction::UnifyOrientation]);
    let corners = out.normals.expect("normals kept").indices.expect("indexed");
    for (c, &original) in corners.iter().enumerate() {
        assert_eq!(
            out.indices[c], mesh.indices[original as usize],
            "corner {c}"
        );
    }
}
