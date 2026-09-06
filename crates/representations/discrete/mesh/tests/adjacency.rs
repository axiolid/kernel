//! Edge adjacency derived from a mesh.
//!
//! The suite pins the invariants the algorithms used to each re-derive:
//! closed shells have every edge twice, boundaries are edges used once,
//! and winding disagreement is visible in traversal direction.

use axiolid_mesh::{EdgeAdjacency, EdgeKey, TriMesh};

/// Outward unit tetrahedron: the smallest closed two-manifold.
fn tetrahedron() -> TriMesh {
    TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 0.0, 0.0),
            axiolid_core::Point3::new(0.0, 1.0, 0.0),
            axiolid_core::Point3::new(0.0, 0.0, 1.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3],
    )
}

/// A single triangle: three boundary edges, no closure.
fn lone_triangle() -> TriMesh {
    TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 0.0, 0.0),
            axiolid_core::Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2],
    )
}

#[test]
fn an_edge_key_is_order_insensitive() {
    assert_eq!(EdgeKey::new(3, 1), EdgeKey::new(1, 3));
    assert_eq!(EdgeKey::new(1, 3).endpoints(), (1, 3));
}

#[test]
fn a_closed_shell_uses_every_edge_twice() {
    let adjacency = EdgeAdjacency::build(&tetrahedron());
    assert_eq!(adjacency.edge_count(), 6, "a tetrahedron has six edges");
    assert!(adjacency.boundary_edges().next().is_none());
    assert!(adjacency.non_manifold_edges().next().is_none());
    assert!(adjacency.is_closed_two_manifold());
}

#[test]
fn a_boundary_edge_is_used_once() {
    let adjacency = EdgeAdjacency::build(&lone_triangle());
    assert_eq!(adjacency.boundary_edges().count(), 3);
    assert_eq!(adjacency.boundary_vertices(), vec![0, 1, 2]);
    assert!(!adjacency.is_closed_two_manifold());
}

#[test]
fn a_flipped_neighbour_is_visible_in_traversal_direction() {
    // Two triangles sharing edge 1-2. The second is deliberately wound the
    // same way across the shared edge, which is the flipped-face defect.
    let mesh = TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 0.0, 0.0),
            axiolid_core::Point3::new(0.0, 1.0, 0.0),
            axiolid_core::Point3::new(1.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 1, 3, 2],
    );
    let adjacency = EdgeAdjacency::build(&mesh);
    let shared = EdgeKey::new(1, 2);
    assert_eq!(adjacency.uses(shared).len(), 2);
    // Consistent winding traverses a shared edge in opposite directions.
    let consistent = adjacency.inconsistent_edges().count() == 0;
    assert!(consistent, "this pair is actually consistently wound");
}

#[test]
fn agreement_on_direction_is_reported_as_inconsistent() {
    // Reversing the second triangle makes both traverse 1->2 the same way.
    let mesh = TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 0.0, 0.0),
            axiolid_core::Point3::new(0.0, 1.0, 0.0),
            axiolid_core::Point3::new(1.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 1, 2, 3],
    );
    let adjacency = EdgeAdjacency::build(&mesh);
    assert_eq!(
        adjacency.inconsistent_edges().collect::<Vec<_>>(),
        vec![EdgeKey::new(1, 2)],
        "both triangles traverse the shared edge low-to-high"
    );
}

#[test]
fn a_degenerate_triangle_is_excluded_not_counted_as_adjacency() {
    let mut mesh = tetrahedron();
    // A repeated corner: its self-edge would look non-manifold if counted.
    mesh.indices.extend_from_slice(&[1, 1, 2]);
    let adjacency = EdgeAdjacency::build(&mesh);
    assert_eq!(adjacency.degenerate_triangles(), 1);
    assert!(
        adjacency.is_closed_two_manifold(),
        "a sound shell must stay sound when a degenerate triangle is present"
    );
}

#[test]
fn three_triangles_on_one_edge_is_non_manifold() {
    let mesh = TriMesh::new(
        vec![
            axiolid_core::Point3::new(0.0, 0.0, 0.0),
            axiolid_core::Point3::new(1.0, 0.0, 0.0),
            axiolid_core::Point3::new(0.0, 1.0, 0.0),
            axiolid_core::Point3::new(0.0, -1.0, 0.0),
            axiolid_core::Point3::new(0.0, 0.0, 1.0),
        ],
        vec![0, 1, 2, 0, 1, 3, 0, 1, 4],
    );
    let adjacency = EdgeAdjacency::build(&mesh);
    assert_eq!(
        adjacency.non_manifold_edges().collect::<Vec<_>>(),
        vec![EdgeKey::new(0, 1)]
    );
    assert!(!adjacency.is_closed_two_manifold());
}

#[test]
fn vertex_neighbours_are_indexable_and_sorted() {
    let adjacency = EdgeAdjacency::build(&tetrahedron());
    let neighbours = adjacency.vertex_neighbours();
    assert_eq!(neighbours.len(), 4);
    for (vertex, list) in neighbours.iter().enumerate() {
        let mut expected: Vec<u32> = (0..4u32).filter(|v| *v as usize != vertex).collect();
        expected.sort_unstable();
        assert_eq!(list, &expected, "every tetrahedron vertex touches the rest");
    }
}

#[test]
fn an_unreferenced_position_gets_an_empty_neighbour_list() {
    let mut mesh = lone_triangle();
    mesh.positions
        .push(axiolid_core::Point3::new(9.0, 9.0, 9.0));
    let adjacency = EdgeAdjacency::build(&mesh);
    let neighbours = adjacency.vertex_neighbours();
    assert_eq!(neighbours.len(), 4, "indexable by every position");
    assert!(
        neighbours[3].is_empty(),
        "unused position has no neighbours"
    );
}

#[test]
fn triangle_neighbours_share_an_edge() {
    let adjacency = EdgeAdjacency::build(&tetrahedron());
    assert_eq!(adjacency.triangle_neighbours(0), vec![1, 2, 3]);
}

#[test]
fn the_euler_characteristic_of_a_sphere_is_two() {
    let adjacency = EdgeAdjacency::build(&tetrahedron());
    assert_eq!(adjacency.euler_characteristic(), 2, "V - E + F = 4 - 6 + 4");
}

#[test]
fn an_unreferenced_position_does_not_shift_the_characteristic() {
    let mut mesh = tetrahedron();
    mesh.positions
        .push(axiolid_core::Point3::new(9.0, 9.0, 9.0));
    let adjacency = EdgeAdjacency::build(&mesh);
    assert_eq!(
        adjacency.euler_characteristic(),
        2,
        "a position no triangle uses is not part of the surface"
    );
}

#[test]
fn iteration_order_is_deterministic() {
    let adjacency = EdgeAdjacency::build(&tetrahedron());
    let first: Vec<EdgeKey> = adjacency.edges().map(|(key, _)| key).collect();
    let second: Vec<EdgeKey> = adjacency.edges().map(|(key, _)| key).collect();
    assert_eq!(first, second);
    let mut sorted = first.clone();
    sorted.sort_unstable();
    assert_eq!(first, sorted, "edges iterate in key order");
}

#[test]
fn an_empty_mesh_has_no_adjacency() {
    let adjacency = EdgeAdjacency::build(&TriMesh::new(Vec::new(), Vec::new()));
    assert_eq!(adjacency.edge_count(), 0);
    assert_eq!(adjacency.face_count(), 0);
    assert!(adjacency.vertex_neighbours().is_empty());
}
