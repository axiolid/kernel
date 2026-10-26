use axiolid_core::{Point3, Tolerance};
use axiolid_measure::{surface_properties, volume_properties, MeshMeasureError};
use axiolid_mesh::TriMesh;

fn tetra() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3],
    )
}

#[test]
fn tetra_surface_and_volume_properties_are_raw_and_deterministic() {
    let mesh = tetra();
    let surface = surface_properties(&mesh, Tolerance::ZERO).unwrap();
    assert!((surface.area - (1.5 + 3.0_f64.sqrt() / 2.0)).abs() < 1e-12);
    let expected_surface_coordinate = (2.0 + 3.0_f64.sqrt()) / (3.0 * (3.0 + 3.0_f64.sqrt()));
    assert!((surface.centroid.x - expected_surface_coordinate).abs() < 1e-12);
    assert!((surface.centroid.y - expected_surface_coordinate).abs() < 1e-12);
    assert!((surface.centroid.z - expected_surface_coordinate).abs() < 1e-12);
    let volume = volume_properties(&mesh, Tolerance::ZERO).unwrap();
    assert!((volume.signed_volume - 1.0 / 6.0).abs() < 1e-12);
    assert_eq!(volume.centroid, Point3::splat(0.25));
    for _ in 0..16 {
        assert_eq!(surface_properties(&mesh, Tolerance::ZERO).unwrap(), surface);
        assert_eq!(volume_properties(&mesh, Tolerance::ZERO).unwrap(), volume);
    }
}

#[test]
fn volume_rejects_open_mesh() {
    let mut mesh = tetra();
    mesh.indices.truncate(9);
    assert!(volume_properties(&mesh, Tolerance::ZERO).is_err());
}

#[test]
fn volume_rejects_closed_mesh_with_inconsistent_winding() {
    let mut mesh = tetra();
    mesh.indices[..3].reverse();
    assert!(volume_properties(&mesh, Tolerance::ZERO).is_err());
}

#[test]
fn reversed_orientation_flips_only_signed_volume() {
    let mut mesh = tetra();
    for face in mesh.indices.as_chunks_mut::<3>().0 {
        face.swap(1, 2);
    }
    let forward = volume_properties(&tetra(), Tolerance::ZERO).unwrap();
    let reversed = volume_properties(&mesh, Tolerance::ZERO).unwrap();
    assert_eq!(forward.signed_volume, -reversed.signed_volume);
    assert_eq!(forward.centroid, reversed.centroid);
}

#[test]
fn translation_preserves_area_volume_and_translates_centroids() {
    let original = tetra();
    let mut moved = tetra();
    let offset = Point3::new(7.0, -3.0, 11.0);
    for p in &mut moved.positions {
        *p += offset;
    }
    let a = surface_properties(&original, Tolerance::ZERO).unwrap();
    let b = surface_properties(&moved, Tolerance::ZERO).unwrap();
    let v = volume_properties(&moved, Tolerance::ZERO).unwrap();
    assert_eq!(a.area, b.area);
    assert!((b.centroid - (a.centroid + offset)).length() < 1e-12);
    assert!((v.signed_volume - 1.0 / 6.0).abs() < 1e-12);
    assert!((v.centroid - (Point3::splat(0.25) + offset)).length() < 1e-12);
}

#[test]
fn zero_area_and_non_manifold_meshes_have_typed_outcomes() {
    let flat = TriMesh::new(vec![Point3::ZERO; 3], vec![0, 1, 2]);
    assert!(matches!(
        surface_properties(&flat, Tolerance::ZERO),
        Err(MeshMeasureError::MeshNotSurfaceUsable(_))
    ));
    let mut non_manifold = tetra();
    non_manifold.indices.extend_from_slice(&[0, 2, 1]);
    assert!(matches!(
        volume_properties(&non_manifold, Tolerance::ZERO),
        Err(MeshMeasureError::MeshNotVolumeUsable(_))
    ));
}
