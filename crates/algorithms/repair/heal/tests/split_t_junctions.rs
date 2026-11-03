//! T-junction splitting: restores the edge pairing a cut opening broke,
//! and declines rather than guess when a split would be wrong.

use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_heal::mesh::{MeshHealError, MeshHealer};
use axiolid_heal::{split_t_junctions, Repair, RepairAction, RepairPlan};
use axiolid_measure::cycle_volume_properties;
use axiolid_mesh::{
    audit_mesh, audit_mesh_cycle, AttributeChannel, AttributeFate, Blend, NormalAttribute, TriMesh,
};

const TOL: f64 = 1e-6;

fn tol() -> Tolerance {
    Tolerance::new(TOL, 1e-9).expect("tolerance")
}

const CORNERS: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];

/// Outward-wound faces: bottom, top (triangles 2, 3), front (4, 5), right,
/// back, left (10, 11).
const CUBE: [u32; 36] = [
    0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0,
    4, 3, 4, 7,
];

fn cube_positions() -> Vec<Point3> {
    CORNERS
        .iter()
        .map(|c| Point3::new(c[0], c[1], c[2]))
        .collect()
}

/// A unit cube with extra vertices on the top-front edge 4-5 (`on_edge`,
/// ascending along it). The top face is split at them; the front face keeps
/// the long edge, which is the T-junction.
fn t_junction(on_edge: &[Point3]) -> TriMesh {
    let mut positions = cube_positions();
    let first = positions.len() as u32;
    positions.extend_from_slice(on_edge);
    let mut indices = CUBE.to_vec();
    let mut top = Vec::new();
    // Fan the top triangle (4, 5, 6) from its third corner over the edge.
    let mut chain = vec![4u32];
    chain.extend((0..on_edge.len() as u32).map(|k| first + k));
    chain.push(5);
    for pair in chain.windows(2) {
        top.extend([pair[0], pair[1], 6]);
    }
    indices.splice(6..9, top);
    TriMesh::new(positions, indices)
}

fn volume(mesh: &TriMesh) -> f64 {
    cycle_volume_properties(mesh, tol())
        .expect("cycle")
        .signed_volume
}

#[test]
fn an_exact_t_junction_is_split_into_a_closed_two_manifold() {
    let mesh = t_junction(&[Point3::new(0.5, 0.0, 1.0)]);
    assert!(!audit_mesh_cycle(&mesh, tol()).is_closed_oriented_cycle());
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.t_junctions, 1);
    assert_eq!(report.split_triangles, 1);
    assert_eq!(report.triangles_added, 1);
    assert_eq!(report.declined_triangles, 0);
    assert_eq!(report.max_deviation, 0.0);
    assert!(audit_mesh(&out, tol()).is_closed_two_manifold());
    assert!((volume(&out) - 1.0).abs() < 1e-15);
    assert_eq!(out.positions, mesh.positions);
}

#[test]
fn a_vertex_off_the_edge_by_half_the_tolerance_is_split_and_reported() {
    let delta = TOL / 2.0;
    let mesh = t_junction(&[Point3::new(0.5, -delta, 1.0)]);
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert!(audit_mesh(&out, tol()).is_closed_two_manifold());
    assert_eq!(report.max_deviation, delta);
    assert!((volume(&out) - (1.0 + delta / 6.0)).abs() < 1e-12);
}

#[test]
fn two_vertices_on_one_edge_are_both_used() {
    let mesh = t_junction(&[
        Point3::new(1.0 / 3.0, 0.0, 1.0),
        Point3::new(2.0 / 3.0, 0.0, 1.0),
    ]);
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.t_junctions, 2);
    assert_eq!(report.triangles_added, 2);
    assert!(audit_mesh(&out, tol()).is_closed_two_manifold());
    assert!((volume(&out) - 1.0).abs() < 1e-15);
}

#[test]
fn vertices_on_two_sides_of_one_triangle_split_it_into_three() {
    // m on 5-4 (top-front), n on 0-4 (front-left vertical edge); the top
    // and left faces are split at them, the front triangle (0, 5, 4) is not.
    let mut positions = cube_positions();
    positions.push(Point3::new(0.5, 0.0, 1.0));
    positions.push(Point3::new(0.0, 0.0, 0.5));
    let mut indices = CUBE.to_vec();
    // Left (3, 0, 4) -> (3, 0, n), (3, n, 4); do it first, it is later.
    indices.splice(30..33, [3, 0, 9, 3, 9, 4]);
    indices.splice(6..9, [4, 8, 6, 8, 5, 6]);
    let mesh = TriMesh::new(positions, indices);
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.t_junctions, 2);
    assert_eq!(report.split_triangles, 1);
    assert_eq!(report.triangles_added, 2);
    assert!(audit_mesh(&out, tol()).is_closed_two_manifold());
    assert!((volume(&out) - 1.0).abs() < 1e-15);
}

#[test]
fn a_vertex_too_far_or_too_near_an_endpoint_is_not_split() {
    for m in [
        Point3::new(0.5, 2.0 * TOL, 1.0),
        Point3::new(TOL / 2.0, 0.0, 1.0),
        Point3::new(1.0 - TOL / 2.0, 0.0, 1.0),
    ] {
        let mesh = t_junction(&[m]);
        let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
        assert_eq!(report.t_junctions, 0, "{m:?}");
        assert_eq!(report.split_triangles, 0);
        assert_eq!(out, mesh);
    }
}

#[test]
fn a_balanced_edge_is_never_split_even_with_a_vertex_on_it() {
    // A second cube's corner sits on the first cube's edge; both are closed
    // two-manifolds, so no edge is unbalanced.
    let mut positions = cube_positions();
    positions.extend(
        CORNERS
            .iter()
            .map(|c| Point3::new(c[0] + 0.5, c[1] - 1.0, c[2] + 1.0)),
    );
    let mut indices = CUBE.to_vec();
    indices.extend(CUBE.iter().map(|i| i + 8));
    let mesh = TriMesh::new(positions, indices);
    assert!(audit_mesh(&mesh, tol()).is_closed_two_manifold());
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.t_junctions, 0);
    assert_eq!(out, mesh);
    let (_, applied) = MeshHealer
        .repair(
            &mesh,
            &RepairPlan {
                actions: vec![RepairAction::SplitTJunctions],
            },
            tol(),
        )
        .unwrap();
    assert_eq!(applied.skipped, vec![RepairAction::SplitTJunctions]);
}

/// Two vertices a hair apart on one edge, the left one just below it and
/// the right one just above: the piece between them and the apex would wind
/// the wrong way, though it is not degenerate.
#[test]
fn a_split_that_inverts_a_piece_is_declined() {
    let tolerance = Tolerance::new(1e-3, 1e-9).unwrap();
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let mesh = TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.5, 1.0, 0.0),
            p(0.3, -9e-4, 0.0),
            p(-5.0, -3.0, 0.0),
            p(-5.0, -4.0, 0.0),
            p(0.3001, 9e-4, 0.0),
            p(5.0, -3.0, 0.0),
            p(5.0, -4.0, 0.0),
        ],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
    );
    let (out, report) = split_t_junctions(&mesh, tolerance).unwrap();
    assert_eq!(report.t_junctions, 2);
    assert_eq!(report.split_triangles, 0);
    assert_eq!(report.declined_triangles, 1);
    assert_eq!(out, mesh);
}

/// Two vertices a hair apart on one edge would leave a sliver between them.
#[test]
fn a_split_that_leaves_a_degenerate_piece_is_declined() {
    let tolerance = Tolerance::new(1e-3, 1e-9).unwrap();
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let mesh = TriMesh::new(
        vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.5, 1.0, 0.0),
            p(0.5, 0.0, 0.0),
            p(0.0, 5.0, 0.0),
            p(1.0, 5.0, 0.0),
            p(0.5 + 1e-7, 0.0, 0.0),
            p(0.0, 6.0, 0.0),
            p(1.0, 6.0, 0.0),
        ],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
    );
    let (out, report) = split_t_junctions(&mesh, tolerance).unwrap();
    assert_eq!(report.t_junctions, 2);
    assert_eq!(report.split_triangles, 0);
    assert_eq!(report.declined_triangles, 1);
    assert_eq!(out, mesh);
}

fn channel_with(values_at_corner: impl Fn(usize) -> u32) -> AttributeChannel {
    let corners: Vec<u32> = (0..39).map(&values_at_corner).collect();
    let values: Vec<f64> = (0..60).flat_map(|i| [f64::from(i), 0.0]).collect();
    AttributeChannel::corner_indexed("uv", values, 2, Blend::Linear, corners)
}

/// The front triangle is index 5, so its corners are slots 15, 16, 17; the
/// long side 5 -> 4 runs from slot 16 to slot 17.
fn with_channel(channel: AttributeChannel) -> TriMesh {
    let mut mesh = t_junction(&[Point3::new(0.5, 0.0, 1.0)]);
    // The top split added a triangle at index 2, so the front moved on by 1.
    assert_eq!(mesh.indices.len(), 39);
    mesh.attributes.push(channel);
    mesh
}

#[test]
fn a_corner_channel_with_equal_ends_is_carried_and_preserved() {
    // Front triangle is now index 6: slots 18..21, long side slots 19, 20.
    let mesh = with_channel(channel_with(|c| if c == 20 { 19 } else { c as u32 }));
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.split_triangles, 1);
    let corners = out.attributes[0].corner_indices.as_ref().unwrap();
    assert_eq!(corners.len(), out.indices.len());
    // Every new corner at the split vertex 8 takes the start corner's value.
    for (slot, &vertex) in out.indices.iter().enumerate() {
        if vertex == 8 && slot >= 18 {
            assert_eq!(corners[slot], 19, "slot {slot}");
        }
    }
    let (_, repaired) = MeshHealer
        .repair(
            &mesh,
            &RepairPlan {
                actions: vec![RepairAction::SplitTJunctions],
            },
            tol(),
        )
        .unwrap();
    assert_eq!(repaired.applied, vec![RepairAction::SplitTJunctions]);
    assert_eq!(
        repaired.attribute_fates,
        vec![("uv".to_owned(), AttributeFate::Preserved)]
    );
}

#[test]
fn bit_equal_values_under_different_indices_count_as_equal() {
    let mut channel = channel_with(|c| c as u32);
    // Pool entry 20 repeats entry 19.
    channel.values[40] = channel.values[38];
    let (_, report) = split_t_junctions(&with_channel(channel), tol()).unwrap();
    assert_eq!(report.split_triangles, 1);
}

#[test]
fn a_corner_channel_with_different_ends_declines_the_split() {
    let mesh = with_channel(channel_with(|c| c as u32));
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.split_triangles, 0);
    assert_eq!(report.declined_triangles, 1);
    assert_eq!(out, mesh);
}

#[test]
fn corner_normals_follow_the_same_rule() {
    let mut mesh = t_junction(&[Point3::new(0.5, 0.0, 1.0)]);
    let values: Vec<Vec3> = (0..39).map(|i| Vec3::new(f64::from(i), 0.0, 1.0)).collect();
    let mut indices: Vec<u32> = (0..39).collect();
    mesh.normals = Some(NormalAttribute {
        values: values.clone(),
        indices: Some(indices.clone()),
    });
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.declined_triangles, 1);
    assert_eq!(out, mesh);

    indices[20] = 19;
    mesh.normals = Some(NormalAttribute {
        values,
        indices: Some(indices),
    });
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.split_triangles, 1);
    let carried = out.normals.unwrap().indices.unwrap();
    assert_eq!(carried.len(), out.indices.len());
}

#[test]
fn triangle_indices_are_stable_and_new_pieces_are_appended() {
    let mesh = t_junction(&[Point3::new(0.5, 0.0, 1.0)]);
    let (out, _) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(out.indices.len(), mesh.indices.len() + 3);
    // Untouched triangles keep their slots; the split one keeps its first
    // piece's slot and the rest comes last.
    assert_eq!(out.indices[..18], mesh.indices[..18]);
    assert_eq!(out.indices[21..mesh.indices.len()], mesh.indices[21..]);
    let first_piece = &out.indices[18..21];
    assert!(first_piece.iter().all(|v| [0, 5, 8, 4].contains(v)));
}

#[test]
fn the_split_is_deterministic() {
    let mesh = t_junction(&[Point3::new(0.25, 0.0, 1.0), Point3::new(0.75, 0.0, 1.0)]);
    assert_eq!(
        split_t_junctions(&mesh, tol()).unwrap(),
        split_t_junctions(&mesh, tol()).unwrap()
    );
}

#[test]
fn a_ragged_index_buffer_is_refused() {
    let mesh = TriMesh::new(cube_positions(), vec![0, 1, 2, 3]);
    assert_eq!(
        split_t_junctions(&mesh, tol()).unwrap_err(),
        MeshHealError::RaggedIndices { indices: 4 }
    );
}

#[test]
fn through_the_trait_the_action_is_applied_or_skipped() {
    let plan = RepairPlan {
        actions: vec![RepairAction::SplitTJunctions],
    };
    let (_, applied) = MeshHealer
        .repair(&t_junction(&[Point3::new(0.5, 0.0, 1.0)]), &plan, tol())
        .unwrap();
    assert_eq!(applied.applied, vec![RepairAction::SplitTJunctions]);
    let clean = TriMesh::new(cube_positions(), CUBE.to_vec());
    let (same, skipped) = MeshHealer.repair(&clean, &plan, tol()).unwrap();
    assert_eq!(skipped.skipped, vec![RepairAction::SplitTJunctions]);
    assert_eq!(same, clean);
}

#[test]
fn the_reported_deviation_is_the_largest_one_used() {
    let delta = TOL / 2.0;
    let mesh = t_junction(&[
        Point3::new(1.0 / 3.0, 0.0, 1.0),
        Point3::new(2.0 / 3.0, -delta, 1.0),
    ]);
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.split_triangles, 1);
    assert_eq!(report.max_deviation, delta);
    assert!(audit_mesh(&out, tol()).is_closed_two_manifold());
}

/// Per-vertex channel over the T-junction cube: vertex 4 holds 0, vertex 5
/// holds 1, the inserted vertex 8 (halfway between them) holds `middle`.
fn with_vertex_channel(blend: Blend, ends: (f64, f64), middle: f64) -> TriMesh {
    let mut mesh = t_junction(&[Point3::new(0.5, 0.0, 1.0)]);
    let mut values = vec![9.0; 9];
    values[4] = ends.0;
    values[5] = ends.1;
    values[8] = middle;
    mesh.attributes
        .push(AttributeChannel::new("w", values, 1, blend));
    mesh
}

#[test]
fn a_per_vertex_value_must_be_the_interpolation_of_its_side() {
    let (_, report) =
        split_t_junctions(&with_vertex_channel(Blend::Linear, (0.0, 1.0), 7.0), tol()).unwrap();
    assert_eq!((report.split_triangles, report.declined_triangles), (0, 1));
    let mesh = with_vertex_channel(Blend::Linear, (0.0, 1.0), 0.5);
    let (out, report) = split_t_junctions(&mesh, tol()).unwrap();
    assert_eq!(report.split_triangles, 1);
    assert_eq!(out.attributes, mesh.attributes);
}

#[test]
fn a_non_linear_per_vertex_value_must_equal_both_ends() {
    for blend in [Blend::Nearest, Blend::None] {
        let (_, report) =
            split_t_junctions(&with_vertex_channel(blend, (3.0, 3.0), 3.0), tol()).unwrap();
        assert_eq!(report.split_triangles, 1);
        let (_, report) =
            split_t_junctions(&with_vertex_channel(blend, (0.0, 1.0), 0.5), tol()).unwrap();
        assert_eq!(report.declined_triangles, 1);
    }
}

#[test]
fn a_per_vertex_normal_must_be_the_interpolation_of_its_side() {
    let mut mesh = t_junction(&[Point3::new(0.5, 0.0, 1.0)]);
    let mut values = vec![Vec3::new(0.0, 0.0, 1.0); 9];
    values[4] = Vec3::new(0.0, 0.0, 0.0);
    values[5] = Vec3::new(1.0, 0.0, 0.0);
    values[8] = Vec3::new(0.5, 0.0, 0.0);
    mesh.normals = Some(NormalAttribute {
        values: values.clone(),
        indices: None,
    });
    assert_eq!(
        split_t_junctions(&mesh, tol()).unwrap().1.split_triangles,
        1
    );
    values[8] = Vec3::new(0.0, 1.0, 0.0);
    mesh.normals = Some(NormalAttribute {
        values,
        indices: None,
    });
    assert_eq!(
        split_t_junctions(&mesh, tol())
            .unwrap()
            .1
            .declined_triangles,
        1
    );
}

/// The triangle's own apex lies within the tolerance of its base, so
/// splitting the base at it would reuse a corner of the triangle itself.
#[test]
fn a_triangle_whose_own_corner_lies_on_its_edge_is_declined() {
    let tolerance = Tolerance::new(1e-3, 1e-9).unwrap();
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 9e-4, 0.0),
        ],
        vec![0, 1, 2],
    );
    let (out, report) = split_t_junctions(&mesh, tolerance).unwrap();
    assert_eq!(report.t_junctions, 1);
    assert_eq!(report.declined_triangles, 1);
    assert_eq!(out, mesh);
}
