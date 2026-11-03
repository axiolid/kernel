//! Volume over closed oriented cycles: closed up to shells that touch along
//! an edge or a T-junction already split, refusing everything open.

use axiolid_core::{Point3, Tolerance};
use axiolid_measure::{
    cycle_volume_properties, volume_properties, CycleVolumeError, MeshMeasureError,
};
use axiolid_mesh::TriMesh;

const TOL: Tolerance = Tolerance::MILLIMETRE;

/// Outward-wound faces of a box: bottom, top, front, right, back, left.
const BOX: [u32; 36] = [
    0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0,
    4, 3, 4, 7,
];

fn box_positions(origin: [f64; 3], size: [f64; 3]) -> Vec<Point3> {
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    corners
        .iter()
        .map(|c| {
            Point3::new(
                origin[0] + c[0] * size[0],
                origin[1] + c[1] * size[1],
                origin[2] + c[2] * size[2],
            )
        })
        .collect()
}

fn cube_at(origin: [f64; 3]) -> TriMesh {
    TriMesh::new(box_positions(origin, [1.0; 3]), BOX.to_vec())
}

/// Cube A plus a second box welded to it through `weld` (second-box corner
/// to first-cube corner); `flip` reverses the second box's winding.
fn pair(second: Vec<Point3>, weld: &[(u32, u32)], flip: bool) -> TriMesh {
    let mut positions = box_positions([0.0; 3], [1.0; 3]);
    let mut map = [0u32; 8];
    for (corner, slot) in map.iter_mut().enumerate() {
        *slot = match weld.iter().find(|(from, _)| *from as usize == corner) {
            Some(&(_, to)) => to,
            None => {
                positions.push(second[corner]);
                positions.len() as u32 - 1
            }
        };
    }
    let mut indices = BOX.to_vec();
    for triangle in BOX.chunks(3) {
        let [a, b, c] = [triangle[0], triangle[1], triangle[2]].map(|i| map[i as usize]);
        indices.extend(if flip { [a, c, b] } else { [a, b, c] });
    }
    TriMesh::new(positions, indices)
}

fn edge_touching() -> TriMesh {
    pair(
        box_positions([1.0, 1.0, 0.0], [1.0; 3]),
        &[(0, 2), (4, 6)],
        false,
    )
}

#[test]
fn two_cubes_sharing_an_edge_sum_and_volume_properties_refuses_them() {
    let mesh = edge_touching();
    assert!(matches!(
        volume_properties(&mesh, TOL),
        Err(MeshMeasureError::MeshNotVolumeUsable(_))
    ));
    let volume = cycle_volume_properties(&mesh, TOL).unwrap();
    assert!((volume.signed_volume - 2.0).abs() < 1e-12);
    assert!((volume.centroid - Point3::new(1.0, 1.0, 0.5)).length() < 1e-12);
}

#[test]
fn two_cubes_sharing_a_face_sum() {
    let mesh = pair(
        box_positions([1.0, 0.0, 0.0], [1.0; 3]),
        &[(0, 1), (3, 2), (4, 5), (7, 6)],
        false,
    );
    let volume = cycle_volume_properties(&mesh, TOL).unwrap();
    assert!((volume.signed_volume - 2.0).abs() < 1e-12);
    assert!((volume.centroid - Point3::new(1.0, 0.5, 0.5)).length() < 1e-12);
}

/// A vertex mid-edge on the top-front edge 4-5. The front face is always
/// split at it; the top face keeps the long edge unless `split_top`.
fn t_junction(split_top: bool) -> TriMesh {
    let mut positions = box_positions([0.0; 3], [1.0; 3]);
    positions.push(Point3::new(0.5, 0.0, 1.0));
    let mut indices = BOX.to_vec();
    indices.splice(15..18, [0, 5, 8, 0, 8, 4]);
    if split_top {
        indices.splice(6..9, [4, 8, 6, 8, 5, 6]);
    }
    TriMesh::new(positions, indices)
}

#[test]
fn a_t_junction_is_refused_until_both_faces_meet_at_its_vertex() {
    let audit = match cycle_volume_properties(&t_junction(false), TOL) {
        Ok(volume) => panic!("a T-junction was measured: {volume:?}"),
        Err(CycleVolumeError::NotClosedOrientedCycle(audit)) => audit,
        Err(other) => panic!("unexpected refusal {other:?}"),
    };
    assert_eq!(audit.unbalanced_edges, 3);
    let volume = cycle_volume_properties(&t_junction(true), TOL).unwrap();
    assert!((volume.signed_volume - 1.0).abs() < 1e-12);
    assert!((volume.centroid - Point3::new(0.5, 0.5, 0.5)).length() < 1e-12);
}

#[test]
fn a_box_without_top_and_bottom_is_refused() {
    // Its vector area vanishes, so a translation-invariance test alone
    // would wrongly accept it.
    let sides = BOX[12..].to_vec();
    let mesh = TriMesh::new(box_positions([0.0; 3], [1.0; 3]), sides);
    match cycle_volume_properties(&mesh, TOL) {
        Err(CycleVolumeError::NotClosedOrientedCycle(audit)) => {
            assert_eq!(audit.unbalanced_edges, 8);
        }
        other => panic!("expected the open box refused, got {other:?}"),
    }
}

fn tetra() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 0, 3, 2],
    )
}

#[test]
fn on_two_manifold_input_it_is_bit_identical_to_volume_properties() {
    for mesh in [tetra(), cube_at([0.0; 3]), cube_at([6.0e5, 5.6e6, 204.5])] {
        let want = volume_properties(&mesh, TOL).unwrap();
        let got = cycle_volume_properties(&mesh, TOL).unwrap();
        assert_eq!(got.signed_volume.to_bits(), want.signed_volume.to_bits());
        assert_eq!(got.centroid, want.centroid);
    }
}

#[test]
fn edge_touching_cubes_on_a_national_grid_match_the_closed_form() {
    let origin = [6.0e5, 5.6e6, 204.5];
    let mesh = {
        let mut m = edge_touching();
        for p in &mut m.positions {
            *p += Point3::new(origin[0], origin[1], origin[2]) - Point3::ZERO;
        }
        m
    };
    let volume = cycle_volume_properties(&mesh, TOL).unwrap();
    assert!((volume.signed_volume - 2.0).abs() < 1e-9);
    let want = Point3::new(origin[0] + 1.0, origin[1] + 1.0, origin[2] + 0.5);
    assert!((volume.centroid - want).length() < 1e-6);
}

/// The result is the integral of the winding number, not the enclosed
/// volume: it counts overlap twice and subtracts an inverted shell.
#[test]
fn it_integrates_the_winding_number_not_the_enclosed_volume() {
    // B = [1,3] x [1,2] x [0,1] touches A = [0,1]^3 along an edge.
    let inverted = pair(
        box_positions([1.0, 1.0, 0.0], [2.0, 1.0, 1.0]),
        &[(0, 2), (4, 6)],
        true,
    );
    let volume = cycle_volume_properties(&inverted, TOL).unwrap();
    assert!((volume.signed_volume + 1.0).abs() < 1e-12);

    let mut doubled_indices = BOX.to_vec();
    doubled_indices.extend(BOX);
    let doubled = TriMesh::new(box_positions([0.0; 3], [1.0; 3]), doubled_indices);
    let volume = cycle_volume_properties(&doubled, TOL).unwrap();
    assert!((volume.signed_volume - 2.0).abs() < 1e-12);
}
