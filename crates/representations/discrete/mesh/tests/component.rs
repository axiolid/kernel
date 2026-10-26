//! Connected-component decomposition and recomposition (#85).
//!
//! The load-bearing case is a boolean that cuts a part in two: the caller
//! needs the pieces, not the number 2.

use axiolid_mesh::{component_count, compose, decompose, TriMesh};

/// An axis-aligned box as a closed, outward-wound triangle mesh.
fn box_mesh(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions = vec![
        [x0, y0, z0].into(),
        [x1, y0, z0].into(),
        [x1, y1, z0].into(),
        [x0, y1, z0].into(),
        [x0, y0, z1].into(),
        [x1, y0, z1].into(),
        [x1, y1, z1].into(),
        [x0, y1, z1].into(),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // bottom
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // -y
        1, 2, 6, 1, 6, 5, // +x
        2, 3, 7, 2, 7, 6, // +y
        3, 0, 4, 3, 4, 7, // -x
    ];
    TriMesh::new(positions, indices)
}

#[test]
fn two_disjoint_boxes_decompose_into_two_closed_bodies() {
    // The shape a difference leaves when it cuts a bar in two.
    let cut = compose(&[
        box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        box_mesh([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]),
    ]);
    assert_eq!(component_count(&cut), 2);

    let parts = decompose(&cut);
    assert_eq!(parts.len(), 2);
    for part in &parts {
        // Each piece is a whole box, not a fragment of one.
        assert_eq!(part.positions.len(), 8);
        assert_eq!(part.indices.len(), 36);
        let health = axiolid_mesh::audit_mesh(part, axiolid_core::Tolerance::METRE);
        assert!(
            health.is_closed_two_manifold(),
            "component is not a closed solid: {health:?}"
        );
    }

    // The pieces are the two originals, distinguished by position.
    let mut x_mins: Vec<f64> = parts
        .iter()
        .map(|p| p.positions.iter().map(|v| v.x).fold(f64::MAX, f64::min))
        .collect();
    x_mins.sort_by(f64::total_cmp);
    assert_eq!(x_mins, vec![0.0, 3.0]);
}

#[test]
fn a_single_body_decomposes_to_itself_unchanged() {
    // Decomposing a connected mesh must not be a reindexing hazard: a
    // caller that always decomposes should not pay for it on single bodies.
    let solid = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let parts = decompose(&solid);
    assert_eq!(parts.len(), 1);
    assert_eq!(parts[0].indices, solid.indices);
    assert_eq!(parts[0].positions, solid.positions);
}

#[test]
fn compose_of_decompose_preserves_every_triangle() {
    let original = compose(&[
        box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        box_mesh([3.0, 0.0, 0.0], [4.0, 1.0, 1.0]),
        box_mesh([6.0, 0.0, 0.0], [7.0, 1.0, 1.0]),
    ]);
    let round_tripped = compose(&decompose(&original));

    assert_eq!(round_tripped.positions.len(), original.positions.len());
    assert_eq!(round_tripped.indices.len(), original.indices.len());

    // Compare as resolved coordinate triples: index numbering may differ,
    // the geometry may not.
    let resolve = |m: &TriMesh| {
        let mut tris: Vec<[[u64; 3]; 3]> = m
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .map(|t| {
                let mut corners: Vec<[u64; 3]> = t
                    .iter()
                    .map(|&i| {
                        let p = m.positions[i as usize];
                        [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
                    })
                    .collect();
                corners.sort();
                [corners[0], corners[1], corners[2]]
            })
            .collect();
        tris.sort();
        tris
    };
    assert_eq!(resolve(&round_tripped), resolve(&original));
}

#[test]
fn component_order_is_stable_across_runs() {
    // #85 requires determinism. Hash-map iteration order would satisfy the
    // count tests and fail this one.
    let mesh = compose(&[
        box_mesh([9.0, 0.0, 0.0], [10.0, 1.0, 1.0]),
        box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        box_mesh([5.0, 0.0, 0.0], [6.0, 1.0, 1.0]),
    ]);

    let signature = |m: &TriMesh| -> Vec<u64> {
        decompose(m)
            .iter()
            .map(|p| p.positions[0].x.to_bits())
            .collect()
    };
    let first = signature(&mesh);
    for _ in 0..8 {
        assert_eq!(signature(&mesh), first, "component order drifted");
    }

    // Order follows first appearance in the input, so it is predictable
    // rather than merely repeatable.
    assert_eq!(
        first,
        vec![9.0f64.to_bits(), 0.0f64.to_bits(), 5.0f64.to_bits()]
    );
}

#[test]
fn an_empty_mesh_has_no_components() {
    let empty = TriMesh::new(Vec::new(), Vec::new());
    assert_eq!(component_count(&empty), 0);
    assert!(decompose(&empty).is_empty());
}

#[test]
fn unreferenced_positions_are_not_components() {
    // A stray position is unused data, not a body. Counting it would make
    // a mesh with slack storage look fragmented.
    let mut solid = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    solid.positions.push([50.0, 50.0, 50.0].into());
    assert_eq!(component_count(&solid), 1);
    assert_eq!(decompose(&solid).len(), 1);
    // The component drops the unreferenced position.
    assert_eq!(decompose(&solid)[0].positions.len(), 8);
}
