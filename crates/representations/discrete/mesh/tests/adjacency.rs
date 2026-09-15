//! Edge adjacency derived from a mesh.
//!
//! The suite pins the invariants the algorithms used to each re-derive:
//! closed shells have every edge twice, boundaries are edges used once,
//! and winding disagreement is visible in traversal direction.

use axiolid_core::Point3;
use axiolid_mesh::{EdgeAdjacency, EdgeKey, EdgeUse, TriMesh};

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

/// `face_count` counts triangles that contributed adjacency. The O(1)
/// form (total minus degenerate) must agree with counting distinct
/// triangle indices in the edge map: degenerate, duplicated, and
/// fully-shared triangles all have to land the same way.
#[test]
fn face_count_matches_counting_distinct_triangles() {
    let quad = || {
        (
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![0u32, 1, 2, 0, 2, 3],
        )
    };

    // Two clean triangles.
    let (positions, indices) = quad();
    check_face_count("quad", &positions, &indices, 2);

    // One degenerate triangle appended: skipped by build, so not a face.
    let (positions, mut indices) = quad();
    indices.extend_from_slice(&[1, 1, 2]);
    check_face_count("one degenerate", &positions, &indices, 2);

    // Every triangle degenerate: no edges at all.
    let (positions, _) = quad();
    check_face_count("all degenerate", &positions, &[0, 0, 0, 1, 1, 1], 0);

    // The same triangle twice: two distinct indices, both real faces.
    let (positions, _) = quad();
    check_face_count("duplicated", &positions, &[0, 1, 2, 0, 1, 2], 2);

    // No triangles.
    check_face_count("empty", &[], &[], 0);
}

/// Compare the stored count against recounting distinct triangle indices
/// from the edge map, which is what `face_count` used to do.
fn check_face_count(label: &str, positions: &[Point3], indices: &[u32], want: usize) {
    let mesh = TriMesh::new(positions.to_vec(), indices.to_vec());
    let adjacency = EdgeAdjacency::build(&mesh);

    let mut distinct = std::collections::BTreeSet::new();
    for (_, uses) in adjacency.edges() {
        for use_ in uses {
            distinct.insert(use_.triangle);
        }
    }

    assert_eq!(
        adjacency.face_count(),
        distinct.len(),
        "{label}: stored count disagrees with recounting the edge map"
    );
    assert_eq!(
        adjacency.face_count(),
        want,
        "{label}: unexpected face count"
    );
}

/// Every observable property of the adjacency, over meshes that stress
/// the awkward cases: shared edges, boundary, non-manifold fans,
/// flipped winding, degenerate and duplicated triangles.
///
/// This exists to pin the CSR rewrite of the edge store. It was written
/// and run against the BTreeMap implementation first, so it is an oracle
/// rather than a restatement of whatever the code currently does.
#[test]
fn the_adjacency_surface_is_stable_across_representations() {
    for (label, mesh) in adjacency_corpus() {
        let adjacency = EdgeAdjacency::build(&mesh);

        // Iteration is ordered by EdgeKey, and `uses` agrees with it.
        let listed: Vec<(EdgeKey, Vec<EdgeUse>)> = adjacency
            .edges()
            .map(|(key, uses)| (key, uses.to_vec()))
            .collect();

        let mut sorted = listed.clone();
        sorted.sort_by_key(|(key, _)| *key);
        assert_eq!(listed, sorted, "{label}: edges() is not ordered by EdgeKey");

        for (key, uses) in &listed {
            assert_eq!(adjacency.uses(*key), &uses[..], "{label}: uses() disagrees");
            // Within an edge, uses are ascending by triangle.
            let mut triangles: Vec<usize> = uses.iter().map(|u| u.triangle).collect();
            let original = triangles.clone();
            triangles.sort_unstable();
            assert_eq!(
                original, triangles,
                "{label}: uses are not in triangle order"
            );
        }

        // An absent edge yields an empty slice, not a panic.
        assert!(
            adjacency
                .uses(EdgeKey::new(u32::MAX - 1, u32::MAX))
                .is_empty(),
            "{label}: absent edge should be empty"
        );

        // Counts and derived queries.
        assert_eq!(adjacency.edge_count(), listed.len(), "{label}: edge_count");
        let boundary: Vec<EdgeKey> = adjacency.boundary_edges().collect();
        let non_manifold: Vec<EdgeKey> = adjacency.non_manifold_edges().collect();
        let inconsistent: Vec<EdgeKey> = adjacency.inconsistent_edges().collect();
        assert_eq!(
            boundary,
            listed
                .iter()
                .filter(|(_, u)| u.len() == 1)
                .map(|(k, _)| *k)
                .collect::<Vec<_>>(),
            "{label}: boundary_edges"
        );
        assert_eq!(
            non_manifold,
            listed
                .iter()
                .filter(|(_, u)| u.len() > 2)
                .map(|(k, _)| *k)
                .collect::<Vec<_>>(),
            "{label}: non_manifold_edges"
        );
        assert_eq!(
            inconsistent,
            listed
                .iter()
                .filter(|(_, u)| u.len() == 2 && u[0].forward == u[1].forward)
                .map(|(k, _)| *k)
                .collect::<Vec<_>>(),
            "{label}: inconsistent_edges"
        );

        // Everything else that callers can observe.
        assert_eq!(
            adjacency.is_closed_two_manifold(),
            listed.iter().all(|(_, u)| u.len() == 2) && inconsistent.is_empty(),
            "{label}: is_closed_two_manifold"
        );
        assert_eq!(
            adjacency.face_count(),
            mesh.indices
                .chunks_exact(3)
                .filter(|c| c[0] != c[1] && c[1] != c[2] && c[2] != c[0])
                .count(),
            "{label}: face_count"
        );
        for triangle in 0..mesh.indices.len() / 3 {
            let neighbours = adjacency.triangle_neighbours(triangle);
            let mut expected: Vec<usize> = listed
                .iter()
                .filter(|(_, u)| u.iter().any(|x| x.triangle == triangle))
                .flat_map(|(_, u)| u.iter().map(|x| x.triangle))
                .filter(|t| *t != triangle)
                .collect();
            expected.sort_unstable();
            expected.dedup();
            assert_eq!(
                neighbours, expected,
                "{label}: triangle_neighbours({triangle})"
            );
        }

        // Cloning and equality must survive the representation change.
        assert_eq!(adjacency.clone(), adjacency, "{label}: Clone/PartialEq");
    }
}

/// Meshes chosen so every branch of the edge store is exercised.
fn adjacency_corpus() -> Vec<(&'static str, TriMesh)> {
    let square = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.5, 0.5, 1.0),
    ];
    vec![
        ("empty", TriMesh::new(Vec::new(), Vec::new())),
        ("lone triangle", TriMesh::new(square.clone(), vec![0, 1, 2])),
        (
            "two triangles sharing an edge",
            TriMesh::new(square.clone(), vec![0, 1, 2, 0, 2, 3]),
        ),
        ("tetrahedron", tetrahedron()),
        // Three triangles on one edge: a non-manifold fan.
        (
            "non-manifold fan",
            TriMesh::new(square.clone(), vec![0, 1, 2, 0, 1, 3, 0, 1, 4]),
        ),
        // Second triangle traverses the shared edge the same way.
        (
            "flipped neighbour",
            TriMesh::new(square.clone(), vec![0, 1, 2, 0, 2, 1]),
        ),
        (
            "degenerate mixed in",
            TriMesh::new(square.clone(), vec![0, 1, 2, 1, 1, 2, 0, 2, 3]),
        ),
        (
            "all degenerate",
            TriMesh::new(square.clone(), vec![0, 0, 0, 1, 1, 1]),
        ),
        (
            "duplicated triangle",
            TriMesh::new(square.clone(), vec![0, 1, 2, 0, 1, 2]),
        ),
        // Indices past the position array: build takes corners as given.
        (
            "index beyond positions",
            TriMesh::new(square, vec![0, 1, 9]),
        ),
    ]
}

/// Uses within one edge must be ascending by triangle.
///
/// The small fixtures above cannot catch a violation: with one or two
/// records per edge there is nothing for an unstable sort to reorder.
/// This uses a mesh large enough that `sort_unstable_by_key` on the key
/// alone genuinely permutes equal keys.
///
/// The order is observable, not incidental: `inconsistent_edges` pairs
/// `uses[0]` against `uses[1]`, and callers index the slice directly.
#[test]
fn uses_within_an_edge_are_ascending_by_triangle() {
    // A closed fan: many triangles, every interior edge shared by two,
    // and triangle indices deliberately not in key order.
    let mut positions = vec![Point3::new(0.0, 0.0, 1.0)];
    let ring = 64;
    for step in 0..ring {
        let angle = f64::from(step) / f64::from(ring) * std::f64::consts::TAU;
        positions.push(Point3::new(angle.cos(), angle.sin(), 0.0));
    }
    positions.push(Point3::new(0.0, 0.0, -1.0));
    let apex = 0u32;
    let base = (ring + 1) as u32;

    // Emit the lower cap first so triangle indices run opposite to the
    // edge-key order for the shared rim edges.
    let mut indices = Vec::new();
    for step in 0..ring {
        let a = 1 + step as u32;
        let b = 1 + ((step as u32 + 1) % ring as u32);
        indices.extend_from_slice(&[base, b, a]);
    }
    for step in 0..ring {
        let a = 1 + step as u32;
        let b = 1 + ((step as u32 + 1) % ring as u32);
        indices.extend_from_slice(&[apex, a, b]);
    }

    let mesh = TriMesh::new(positions, indices);
    let adjacency = EdgeAdjacency::build(&mesh);

    let mut shared = 0;
    for (key, uses) in adjacency.edges() {
        if uses.len() > 1 {
            shared += 1;
        }
        let triangles: Vec<usize> = uses.iter().map(|use_| use_.triangle).collect();
        let mut sorted = triangles.clone();
        sorted.sort_unstable();
        assert_eq!(
            triangles,
            sorted,
            "edge {:?}: uses are not ascending by triangle",
            key.endpoints()
        );
    }
    assert!(
        shared > ring,
        "fixture is too weak: only {shared} shared edges, nothing to reorder"
    );
}

/// A few triangles over a huge index space takes the comparison-sort
/// fallback, not the counting sort. Without this every test above
/// exercises only the counted path, so a broken fallback would ship
/// unnoticed.
#[test]
fn a_sparse_index_space_still_sorts_correctly() {
    // Indices far apart: buckets would dwarf the record count.
    let far = 5_000_000u32;
    let positions = vec![Point3::new(0.0, 0.0, 0.0); 4];
    // Two triangles sharing edge (0, far): the shared edge must still
    // report both, in ascending triangle order.
    let indices = vec![0, far, far + 1, 0, far, far + 2];
    let mesh = TriMesh::new(positions, indices);
    let adjacency = EdgeAdjacency::build(&mesh);

    let keys: Vec<(u32, u32)> = adjacency.edges().map(|(k, _)| k.endpoints()).collect();
    let mut sorted = keys.clone();
    sorted.sort_unstable();
    assert_eq!(keys, sorted, "keys are not ascending on the sparse path");

    let shared = adjacency.uses(EdgeKey::new(0, far));
    let triangles: Vec<usize> = shared.iter().map(|u| u.triangle).collect();
    assert_eq!(triangles, vec![0, 1], "shared edge lost a use or its order");
}
