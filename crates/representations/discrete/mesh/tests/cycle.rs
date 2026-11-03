//! Oriented-cycle audit: closed up to edge-touching shells, refusing open
//! surfaces, flipped faces, T-junctions and degenerate triangles.

use axiolid_core::{Point3, Tolerance};
use axiolid_mesh::{audit_mesh, audit_mesh_cycle, try_audit_mesh_cycle, TriMesh};

const TOL: Tolerance = Tolerance::MILLIMETRE;

fn cube_positions(offset: [f64; 3]) -> Vec<Point3> {
    let [ox, oy, oz] = offset;
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
    .iter()
    .map(|p| Point3::new(p[0] + ox, p[1] + oy, p[2] + oz))
    .collect()
}

/// Outward-wound faces: bottom, top, front, right, back, left.
const CUBE: [u32; 36] = [
    0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0,
    4, 3, 4, 7,
];

fn cube() -> TriMesh {
    TriMesh::new(cube_positions([0.0; 3]), CUBE.to_vec())
}

#[test]
fn a_closed_box_is_a_cycle_with_no_unbalanced_edge() {
    let audit = audit_mesh_cycle(&cube(), TOL);
    assert_eq!(audit.unbalanced_edges, 0);
    assert_eq!(audit.first_unbalanced_edge, None);
    assert!(audit.is_closed_oriented_cycle());
    assert_eq!(audit.health, audit_mesh(&cube(), TOL));
    assert_eq!(try_audit_mesh_cycle(&cube(), TOL).unwrap(), audit);
}

#[test]
fn two_cubes_welded_along_an_edge_are_a_cycle_but_not_two_manifold() {
    // The second cube is shifted by (1, 1, 0): its corner 0 is A's corner 2
    // and its corner 4 is A's corner 6, so both share the edge x = y = 1.
    let mut positions = cube_positions([0.0; 3]);
    positions.extend(cube_positions([1.0, 1.0, 0.0]));
    let weld = |i: u32| match i {
        0 => 2,
        4 => 6,
        other => other + 8,
    };
    let mut indices = CUBE.to_vec();
    indices.extend(CUBE.iter().map(|&i| weld(i)));
    let audit = audit_mesh_cycle(&TriMesh::new(positions, indices), TOL);
    assert_eq!(audit.health.non_manifold_edges, 1);
    assert!(!audit.health.is_closed_two_manifold());
    assert_eq!(audit.unbalanced_edges, 0);
    assert!(audit.is_closed_oriented_cycle());
}

#[test]
fn a_box_without_top_and_bottom_is_unbalanced_on_its_eight_rim_edges() {
    let sides = CUBE[12..].to_vec();
    let audit = audit_mesh_cycle(&TriMesh::new(cube_positions([0.0; 3]), sides), TOL);
    assert_eq!(audit.unbalanced_edges, 8);
    assert_eq!(audit.first_unbalanced_edge, Some((0, 1)));
    assert!(!audit.is_closed_oriented_cycle());
}

#[test]
fn one_flipped_triangle_unbalances_its_edges() {
    let mut indices = CUBE.to_vec();
    indices.swap(1, 2);
    let audit = audit_mesh_cycle(&TriMesh::new(cube_positions([0.0; 3]), indices), TOL);
    assert!(audit.unbalanced_edges > 0);
    assert!(!audit.is_closed_oriented_cycle());
}

#[test]
fn an_unsplit_t_junction_is_three_unbalanced_edges() {
    // Vertex 8 sits mid-edge on the top-front edge 4-5; the top face keeps
    // the long edge, the front face is split at 8.
    let mut positions = cube_positions([0.0; 3]);
    positions.push(Point3::new(0.5, 0.0, 1.0));
    let mut indices = CUBE.to_vec();
    indices.splice(15..18, [0, 5, 8, 0, 8, 4]);
    let audit = audit_mesh_cycle(&TriMesh::new(positions, indices), TOL);
    assert_eq!(audit.unbalanced_edges, 3);
    assert!(!audit.is_closed_oriented_cycle());
}

/// Edge (0, 1) is shared by four triangles; it balances when two use each
/// direction, not merely when the count is even.
#[test]
fn four_use_edges_balance_by_direction_not_by_count() {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, -1.0),
    ];
    let balanced = TriMesh::new(positions.clone(), vec![0, 1, 2, 0, 1, 3, 1, 0, 4, 1, 0, 5]);
    let audit = audit_mesh_cycle(&balanced, TOL);
    assert_eq!(audit.health.non_manifold_edges, 1);
    // The other edges are open, so (0, 1) is not the first unbalanced one.
    assert_ne!(audit.first_unbalanced_edge, Some((0, 1)));
    assert_eq!(audit.unbalanced_edges, 8);
    let skewed = TriMesh::new(positions, vec![0, 1, 2, 0, 1, 3, 0, 1, 4, 1, 0, 5]);
    let audit = audit_mesh_cycle(&skewed, TOL);
    assert_eq!(audit.unbalanced_edges, 9);
    assert_eq!(audit.first_unbalanced_edge, Some((0, 1)));
}

#[test]
fn a_degenerate_triangle_refuses_even_when_balanced() {
    // A sliver and its reverse cancel but are degenerate.
    let mut positions = cube_positions([0.0; 3]);
    positions.extend([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ]);
    let mut indices = CUBE.to_vec();
    indices.extend([8, 9, 10, 8, 10, 9]);
    let audit = audit_mesh_cycle(&TriMesh::new(positions, indices), TOL);
    assert_eq!(audit.unbalanced_edges, 0);
    assert_eq!(audit.health.degenerate_triangles, 2);
    assert!(!audit.is_closed_oriented_cycle());
}

#[test]
fn a_two_manifold_mesh_is_always_a_cycle() {
    let tetra = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 0, 3, 2],
    );
    for mesh in [cube(), tetra] {
        let audit = audit_mesh_cycle(&mesh, TOL);
        assert!(audit.health.is_closed_two_manifold());
        assert!(audit.is_closed_oriented_cycle());
    }
}

#[test]
fn an_invalid_index_is_not_a_cycle() {
    let mesh = TriMesh::new(cube_positions([0.0; 3]), vec![0, 1, 99]);
    assert!(!audit_mesh_cycle(&mesh, TOL).is_closed_oriented_cycle());
}
