//! `SymmetricDifference` composed through a real backend (ADR 0017 section 1).
//!
//! The contract exposes four operations; `boolmesh` implements three. The
//! fourth is composed as `(A ∪ B) \ (A ∩ B)`. Composition is only acceptable
//! if it is *correct*, so this measures volume against the set identity rather
//! than asserting the call merely returned.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// Axis-aligned box as an outward-oriented closed mesh.
fn box_at(min: [f64; 3], max: [f64; 3]) -> TriMesh {
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
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Enclosed volume by the divergence theorem.
fn volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

#[test]
fn symmetric_difference_matches_the_set_identity() {
    let provider = BoolmeshBoolean::new();
    // Two unit cubes overlapping in an eighth of their volume.
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    let xor = provider
        .boolean(&a, &b, BooleanOperator::SymmetricDifference, &options())
        .expect("composed symmetric difference");

    let union = provider
        .boolean(&a, &b, BooleanOperator::Union, &options())
        .expect("union")
        .mesh;
    let intersection = provider
        .boolean(&a, &b, BooleanOperator::Intersection, &options())
        .expect("intersection")
        .mesh;

    // vol(A △ B) == vol(A ∪ B) - vol(A ∩ B)
    let expected = volume(&union) - volume(&intersection);
    let actual = volume(&xor.mesh);
    assert!(
        (actual - expected).abs() < 1e-9,
        "symmetric difference volume {actual} != {expected}"
    );

    // Independently: each cube is 1.0, overlap is 0.5^3 = 0.125, so the
    // symmetric difference is 2*1.0 - 2*0.125 = 1.75.
    assert!(
        (actual - 1.75).abs() < 1e-9,
        "expected analytic 1.75, got {actual}"
    );
}

#[test]
fn composition_is_visible_in_the_evidence() {
    let provider = BoolmeshBoolean::new();
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    let primitive = provider
        .boolean(&a, &b, BooleanOperator::Union, &options())
        .unwrap();
    assert_eq!(
        primitive.evidence.sub_operations, 1,
        "a native operation is one pass"
    );

    let composed = provider
        .boolean(&a, &b, BooleanOperator::SymmetricDifference, &options())
        .unwrap();
    assert_eq!(
        composed.evidence.sub_operations, 3,
        "union + intersection + difference"
    );
}

#[test]
fn disjoint_operands_shortcut_to_the_union() {
    let provider = BoolmeshBoolean::new();
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let far = box_at([10.0, 10.0, 10.0], [11.0, 11.0, 11.0]);

    let xor = provider
        .boolean(&a, &far, BooleanOperator::SymmetricDifference, &options())
        .expect("disjoint symmetric difference");

    // A ∩ B is empty, so A △ B == A ∪ B and the final difference is skipped.
    assert_eq!(xor.evidence.sub_operations, 2);
    assert!((volume(&xor.mesh) - 2.0).abs() < 1e-9);
    assert_eq!(xor.evidence.output_components, 2, "two separate solids");
}

#[test]
fn a_difference_that_severs_the_subject_reports_two_components() {
    let provider = BoolmeshBoolean::new();
    // A bar cut clean through the middle by a taller, wider blade.
    let bar = box_at([0.0, 0.0, 0.0], [3.0, 1.0, 1.0]);
    let blade = box_at([1.0, -1.0, -1.0], [2.0, 2.0, 2.0]);

    let outcome = provider
        .boolean(&bar, &blade, BooleanOperator::Difference, &options())
        .expect("severing difference");

    assert_eq!(
        outcome.evidence.output_components, 2,
        "the cut splits the bar, and the caller learns it here"
    );
    assert!((volume(&outcome.mesh) - 2.0).abs() < 1e-9);
}
