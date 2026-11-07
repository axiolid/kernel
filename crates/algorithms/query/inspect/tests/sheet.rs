//! Faces of one mesh that face each other within a linear distance (#400).

use axiolid_core::Point3;
use axiolid_inspect::thin_sheet;
use axiolid_mesh::TriMesh;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// The corners and outward-facing triangles of the box `lo..hi`, appended to
/// `positions` and `indices`.
fn add_box(positions: &mut Vec<Point3>, indices: &mut Vec<u32>, lo: [f64; 3], hi: [f64; 3]) {
    let base = positions.len() as u32;
    for corner in 0..8u32 {
        positions.push(p(
            if corner & 1 == 0 { lo[0] } else { hi[0] },
            if corner & 2 == 0 { lo[1] } else { hi[1] },
            if corner & 4 == 0 { lo[2] } else { hi[2] },
        ));
    }
    const FACES: [[u32; 3]; 12] = [
        [0, 2, 3],
        [0, 3, 1],
        [4, 5, 7],
        [4, 7, 6],
        [0, 1, 5],
        [0, 5, 4],
        [2, 6, 7],
        [2, 7, 3],
        [0, 4, 6],
        [0, 6, 2],
        [1, 3, 7],
        [1, 7, 5],
    ];
    for face in FACES {
        indices.extend(face.map(|corner| base + corner));
    }
}

fn boxes(parts: &[([f64; 3], [f64; 3])]) -> TriMesh {
    let (mut positions, mut indices) = (Vec::new(), Vec::new());
    for (lo, hi) in parts {
        add_box(&mut positions, &mut indices, *lo, *hi);
    }
    TriMesh::new(positions, indices)
}

#[test]
fn a_box_has_no_sheet() {
    let mesh = boxes(&[([0.0; 3], [1.0; 3])]);
    assert_eq!(thin_sheet(&mesh, 1e-9), None);
}

#[test]
fn two_boxes_stacked_face_to_face_are_a_sheet_of_zero_offset() {
    let mesh = boxes(&[([0.0; 3], [1.0; 3]), ([0.0, 0.0, 1.0], [1.0, 1.0, 2.0])]);
    let found = thin_sheet(&mesh, 1e-9).expect("the stacked faces face each other");
    assert!(found.first < found.second);
    assert_eq!(found.offset, 0.0);
}

#[test]
fn a_slab_thinner_than_rounding_is_a_sheet() {
    let mesh = boxes(&[([0.0; 3], [1.0, 1.0, 1e-16])]);
    let found = thin_sheet(&mesh, 1e-9).expect("top and bottom face each other");
    assert!(found.offset <= 1e-9);
}

#[test]
fn a_slab_is_a_sheet_only_within_the_linear_distance() {
    let mesh = boxes(&[([0.0; 3], [1.0, 1.0, 1e-6])]);
    assert_eq!(thin_sheet(&mesh, 1e-9), None);
    assert!(thin_sheet(&mesh, 1e-5).is_some());
}

#[test]
fn boxes_touching_along_an_edge_are_no_sheet() {
    let mesh = boxes(&[([0.0; 3], [1.0; 3]), ([1.0, 0.0, -1.0], [2.0, 1.0, 0.0])]);
    assert_eq!(thin_sheet(&mesh, 1e-9), None);
}

#[test]
fn a_fin_folded_back_about_a_shared_edge_is_a_sheet() {
    let mesh = TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 1.0, 1e-12),
        ],
        vec![0, 1, 2, 1, 0, 3],
    );
    let found = thin_sheet(&mesh, 1e-9).expect("the fin lies on its face");
    assert_eq!((found.first, found.second), (0, 1));
    assert!(found.offset > 0.0 && found.offset <= 1e-9);
}

#[test]
fn a_linear_distance_that_is_negative_or_not_finite_finds_nothing() {
    let mesh = boxes(&[([0.0; 3], [1.0, 1.0, 1e-16])]);
    for linear in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(thin_sheet(&mesh, linear), None, "{linear}");
    }
}

#[test]
fn a_malformed_index_buffer_finds_nothing() {
    let corners = vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)];
    for indices in [vec![0, 1, 7], vec![0, 1, 2, 2, 1], vec![0, 1]] {
        let mesh = TriMesh::new(corners.clone(), indices);
        assert_eq!(thin_sheet(&mesh, 1e-9), None);
    }
}
