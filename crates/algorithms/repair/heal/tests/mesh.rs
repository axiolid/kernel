//! Mesh diagnosis and repair against known-defective inputs.
//!
//! Every fixture is defective on purpose and is checked for its SPECIFIC
//! defect, not merely for "some" defect. The strongest case is the
//! inverted cube: closed, two-manifold, zero boundary edges, and still
//! wrong. That is the defect class that cost two debugging rounds in the
//! sweep work, and a count-based audit cannot locate it.

use axiolid_core::{Point3, Tolerance};
use axiolid_heal::mesh::MeshHealer;
use axiolid_heal::{DefectKind, Diagnose, Repair, RepairAction, RepairPlan};
use axiolid_mesh::TriMesh;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

/// A closed unit cube wound consistently outward.
fn cube() -> TriMesh {
    let p = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
    ];
    let f: [[u32; 3]; 12] = [
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    TriMesh {
        positions: p.to_vec(),
        indices: f.concat(),
        normals: None,
        attributes: Vec::new(),
    }
}

#[test]
fn a_clean_cube_is_diagnosed_clean() {
    let d = MeshHealer.diagnose(&cube(), tol()).expect("diagnose");
    assert!(d.is_clean(), "clean cube reported {:?}", d.defects);
    assert!(!d.blocks_boolean());
}

#[test]
fn an_inverted_face_is_located_not_just_counted() {
    // Flip ONE triangle. The mesh stays closed and two-manifold: only the
    // shared-edge traversal reveals it.
    let mut mesh = cube();
    mesh.indices.swap(1, 2);
    let d = MeshHealer.diagnose(&mesh, tol()).expect("diagnose");
    let located: Vec<_> = d
        .defects
        .iter()
        .filter(|x| x.kind == DefectKind::InconsistentOrientation)
        .collect();
    assert!(!located.is_empty(), "inconsistent winding not detected");
    assert!(
        located.iter().all(|x| x.element.is_some()),
        "a defect without an element is not actionable"
    );
    assert!(d.blocks_boolean());
}

/// Signed volume by the divergence theorem, computed here rather than
/// imported: this test must not depend on the crate it is judging.
fn signed_volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            a.cross(b).dot(c) / 6.0
        })
        .sum()
}

#[test]
fn unifying_orientation_restores_the_signed_volume() {
    // The cube's volume is 1. Flipping half its faces leaves it closed and
    // manifold while the signed volume collapses; a real repair must put
    // the magnitude back, which reporting "applied" alone cannot fake.
    let good = cube();
    assert!((signed_volume(&good) - 1.0).abs() < 1e-12);

    let mut broken = cube();
    for t in [0usize, 2, 4, 6] {
        broken.indices.swap(t * 3 + 1, t * 3 + 2);
    }
    assert!(
        (signed_volume(&broken) - 1.0).abs() > 0.1,
        "fixture must actually be broken"
    );

    let plan = RepairPlan {
        actions: vec![RepairAction::UnifyOrientation],
    };
    let (fixed, report) = MeshHealer.repair(&broken, &plan, tol()).expect("repair");
    assert_eq!(report.applied, vec![RepairAction::UnifyOrientation]);
    assert!(
        (signed_volume(&fixed).abs() - 1.0).abs() < 1e-12,
        "volume after repair was {}",
        signed_volume(&fixed)
    );
    assert!(
        MeshHealer
            .diagnose(&fixed, tol())
            .expect("re-diagnose")
            .is_clean(),
        "repaired mesh must diagnose clean"
    );
}

#[test]
fn welding_merges_split_vertices_and_closes_the_shell() {
    // Split every vertex per-face: 36 coincident positions, a mesh that
    // renders identically to the cube but has 36 boundary edges.
    let cube = cube();
    let mut split = TriMesh {
        positions: Vec::new(),
        indices: Vec::new(),
        normals: None,
        attributes: Vec::new(),
    };
    for &i in &cube.indices {
        split.indices.push(split.positions.len() as u32);
        split.positions.push(cube.positions[i as usize]);
    }
    assert_eq!(split.positions.len(), 36);

    let before = MeshHealer.diagnose(&split, tol()).expect("diagnose");
    assert!(
        before
            .defects
            .iter()
            .any(|d| d.kind == DefectKind::DuplicateVertex),
        "split vertices not detected"
    );
    assert!(
        before
            .defects
            .iter()
            .any(|d| d.kind == DefectKind::OpenShell),
        "split cube must read as open"
    );

    let plan = RepairPlan {
        actions: vec![RepairAction::WeldVertices],
    };
    let (welded, report) = MeshHealer.repair(&split, &plan, tol()).expect("repair");
    assert_eq!(report.applied, vec![RepairAction::WeldVertices]);
    assert_eq!(welded.positions.len(), 8, "weld must compact positions");
    assert!(
        MeshHealer
            .diagnose(&welded, tol())
            .expect("re-diagnose")
            .is_clean(),
        "welding must close the shell"
    );
}

#[test]
fn repair_does_only_what_the_plan_asked() {
    // A mesh with BOTH a degenerate triangle and split vertices, repaired
    // with a plan naming only one action. There is deliberately no `All`
    // variant in RepairAction, and this is what enforces that intent.
    let mut mesh = cube();
    let dup = mesh.positions[0];
    mesh.positions.push(dup);
    let n = mesh.positions.len() as u32 - 1;
    mesh.indices.extend_from_slice(&[0, n, 0]);

    let plan = RepairPlan {
        actions: vec![RepairAction::DropDegenerateElements],
    };
    let (fixed, report) = MeshHealer.repair(&mesh, &plan, tol()).expect("repair");
    assert_eq!(report.applied, vec![RepairAction::DropDegenerateElements]);
    assert_eq!(fixed.indices.len(), 36, "degenerate triangle not removed");
    assert_eq!(fixed.positions.len(), 9, "an unrequested weld must not run");
}

#[test]
fn an_inapplicable_repair_is_reported_as_skipped() {
    // A clean cube has nothing to weld. Reporting that honestly is the
    // difference between an audit trail and a rubber stamp.
    let plan = RepairPlan {
        actions: vec![RepairAction::WeldVertices],
    };
    let (_, report) = MeshHealer.repair(&cube(), &plan, tol()).expect("repair");
    assert!(report.applied.is_empty());
    assert_eq!(report.skipped, vec![RepairAction::WeldVertices]);
}

#[test]
fn a_ragged_index_buffer_is_refused() {
    let mut mesh = cube();
    mesh.indices.pop();
    assert!(MeshHealer.diagnose(&mesh, tol()).is_err());
    assert!(MeshHealer
        .repair(&mesh, &RepairPlan::default(), tol())
        .is_err());
}

#[test]
fn a_sphere_with_split_poles_welds_to_a_closed_shell() {
    // The pole defect from the primitive work: every triangle around a
    // pole carries its own coincident pole vertex, so the shell reads as
    // open even though it looks correct. This is a real defect class from
    // this codebase, not an invented one.
    let bands = 8usize;
    let mut mesh = TriMesh {
        positions: Vec::new(),
        indices: Vec::new(),
        normals: None,
        attributes: Vec::new(),
    };
    let ring = |k: usize, lat: f64| {
        let phi = std::f64::consts::TAU * (k % bands) as f64 / bands as f64;
        Point3::new(lat.sin() * phi.cos(), lat.sin() * phi.sin(), lat.cos())
    };
    let mid = std::f64::consts::FRAC_PI_2;
    for k in 0..bands {
        // North fan: a fresh pole vertex per triangle, which is the bug.
        let mut push = |p: Point3| {
            mesh.indices.push(mesh.positions.len() as u32);
            mesh.positions.push(p);
        };
        push(Point3::new(0.0, 0.0, 1.0));
        push(ring(k, mid));
        push(ring(k + 1, mid));
        push(Point3::new(0.0, 0.0, -1.0));
        push(ring(k + 1, mid));
        push(ring(k, mid));
    }

    let before = MeshHealer.diagnose(&mesh, tol()).expect("diagnose");
    assert!(
        before
            .defects
            .iter()
            .any(|d| d.kind == DefectKind::OpenShell),
        "split poles must read as an open shell"
    );

    let plan = RepairPlan {
        actions: vec![RepairAction::WeldVertices],
    };
    let (welded, report) = MeshHealer.repair(&mesh, &plan, tol()).expect("repair");
    assert_eq!(report.applied, vec![RepairAction::WeldVertices]);
    assert_eq!(
        welded.positions.len(),
        bands + 2,
        "poles and ring must collapse to one vertex each"
    );
    let after = MeshHealer.diagnose(&welded, tol()).expect("re-diagnose");
    assert!(
        !after
            .defects
            .iter()
            .any(|d| d.kind == DefectKind::OpenShell),
        "welding must close the sphere: {:?}",
        after.defects
    );
}
