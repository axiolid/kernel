//! Flipping an inside-out closed shell the right way round (#74).
//!
//! `UnifyOrientation` makes neighbouring faces agree but leaves the
//! absolute sense as its seed found it. An inside-out shell is closed,
//! two-manifold and consistently wound, so no topological audit finds it.
//! Only enclosed volume can tell, which is why this repair needs the
//! measurement provider.

use axiolid_core::Tolerance;
use axiolid_heal::mesh::MeshHealer;
use axiolid_heal::{Repair, RepairAction, RepairPlan};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;

fn tol() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).expect("valid tolerance")
}

/// Unit cube wound so its faces point OUTWARD (positive volume).
fn outward_cube() -> TriMesh {
    let positions = vec![
        [0.0, 0.0, 0.0].into(),
        [1.0, 0.0, 0.0].into(),
        [1.0, 1.0, 0.0].into(),
        [0.0, 1.0, 0.0].into(),
        [0.0, 0.0, 1.0].into(),
        [1.0, 0.0, 1.0].into(),
        [1.0, 1.0, 1.0].into(),
        [0.0, 1.0, 1.0].into(),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// The same cube with every triangle reversed: inside-out.
fn inward_cube() -> TriMesh {
    let mut mesh = outward_cube();
    for triangle in mesh.indices.as_chunks_mut::<3>().0 {
        triangle.swap(1, 2);
    }
    mesh
}

fn signed_volume(mesh: &TriMesh) -> f64 {
    volume_properties(mesh, tol())
        .expect("closed shell measures")
        .signed_volume
}

fn plan() -> RepairPlan {
    RepairPlan {
        actions: vec![RepairAction::OrientOutward],
    }
}

#[test]
fn an_inside_out_cube_is_flipped_to_enclose_positive_volume() {
    let broken = inward_cube();
    assert!(
        signed_volume(&broken) < 0.0,
        "fixture must start inside-out, else the test proves nothing"
    );

    let (mesh, report) = MeshHealer
        .repair(&broken, &plan(), tol())
        .expect("a closed shell is repairable");

    assert_eq!(
        report.applied,
        vec![RepairAction::OrientOutward],
        "the flip must be reported, not applied silently"
    );
    let volume = signed_volume(&mesh);
    assert!(
        (volume - 1.0).abs() < 1e-12,
        "a unit cube must enclose +1, got {volume}"
    );
}

#[test]
fn an_already_outward_cube_is_left_exactly_alone() {
    // The issue's own requirement: a repair that always finds something
    // to fix is not trustworthy.
    let clean = outward_cube();
    let (mesh, report) = MeshHealer
        .repair(&clean, &plan(), tol())
        .expect("a clean mesh is repairable");

    assert!(
        report.applied.is_empty(),
        "nothing was wrong, so nothing may be reported as fixed"
    );
    assert_eq!(
        report.skipped,
        vec![RepairAction::OrientOutward],
        "an inapplicable repair is skipped, not silently dropped"
    );
    assert_eq!(
        mesh.indices, clean.indices,
        "a no-op repair must return byte-identical topology"
    );
    assert_eq!(mesh.positions, clean.positions);
}

#[test]
fn an_open_surface_is_not_flipped() {
    // A single triangle encloses no volume, so `inward` is undefined for
    // it. Flipping on a guess would corrupt a caller's surface normals.
    let sheet = TriMesh::new(
        vec![
            [0.0, 0.0, 0.0].into(),
            [1.0, 0.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
        ],
        vec![0, 1, 2],
    );
    let (mesh, report) = MeshHealer
        .repair(&sheet, &plan(), tol())
        .expect("an open surface is still a valid input");
    assert!(
        report.applied.is_empty(),
        "an open surface has no enclosed volume to correct"
    );
    assert_eq!(mesh.indices, sheet.indices);
}
