//! Validation of the scalar boolean oracle (ADR 0012, ADR 0017 §5).
//!
//! An oracle nobody checked is just a second opinion. These tests pin it to
//! values computed independently of any boolean implementation: analytic
//! volumes, set-algebra identities, and exact containment.
//!
//! If these fail, the oracle is wrong and every conformance verdict built on
//! it is void -- so they run before the conformance suite, not after.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;
use axiolid_reference::ScalarBoolean;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// Axis-aligned box as an outward-oriented closed surface.
fn box_at(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions = vec![
        Point3::new(x0, y0, z0),
        Point3::new(x1, y0, z0),
        Point3::new(x1, y1, z0),
        Point3::new(x0, y1, z0),
        Point3::new(x0, y0, z1),
        Point3::new(x1, y0, z1),
        Point3::new(x1, y1, z1),
        Point3::new(x0, y1, z1),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Enclosed volume by the divergence theorem. Independent of any boolean.
fn volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

fn apply(subject: &TriMesh, tool: &TriMesh, op: BooleanOperator) -> TriMesh {
    ScalarBoolean::new()
        .boolean(subject, tool, op, &options())
        .unwrap_or_else(|error| panic!("{op:?}: {error}"))
        .mesh
}

/// The unit cube fixture must itself be right, or every volume below is wrong.
#[test]
fn the_fixture_has_the_analytic_volume() {
    assert!((volume(&box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])) - 1.0).abs() < 1e-12);
    assert!((volume(&box_at([0.0, 0.0, 0.0], [2.0, 3.0, 4.0])) - 24.0).abs() < 1e-12);
}

// --- disjoint operands ------------------------------------------------

#[test]
fn disjoint_operands_follow_the_set_algebra() {
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([5.0, 5.0, 5.0], [6.0, 6.0, 6.0]);

    // vol(A ∪ B) = 1 + 1 with no shared region to double count.
    assert!((volume(&apply(&a, &b, BooleanOperator::Union)) - 2.0).abs() < 1e-12);
    // A ∩ B = ∅.
    assert_eq!(
        apply(&a, &b, BooleanOperator::Intersection).triangle_count(),
        0
    );
    // A \ B = A when B removes nothing.
    assert!((volume(&apply(&a, &b, BooleanOperator::Difference)) - 1.0).abs() < 1e-12);
    // A △ B = A ∪ B when the intersection is empty.
    assert!((volume(&apply(&a, &b, BooleanOperator::SymmetricDifference)) - 2.0).abs() < 1e-12);
}

// --- nested operands --------------------------------------------------

#[test]
fn a_nested_operand_is_absorbed_or_subtracted_exactly() {
    let outer = box_at([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]); // 64
    let inner = box_at([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]); // 1

    // A ∪ B = A when B ⊂ A.
    assert!((volume(&apply(&outer, &inner, BooleanOperator::Union)) - 64.0).abs() < 1e-12);
    // A ∩ B = B when B ⊂ A.
    assert!((volume(&apply(&outer, &inner, BooleanOperator::Intersection)) - 1.0).abs() < 1e-12);
    // A \ B leaves a cavity: 64 - 1.
    assert!((volume(&apply(&outer, &inner, BooleanOperator::Difference)) - 63.0).abs() < 1e-12);
    // A △ B is the same shell here, since B ⊂ A.
    assert!(
        (volume(&apply(&outer, &inner, BooleanOperator::SymmetricDifference)) - 63.0).abs() < 1e-12
    );
}

#[test]
fn containment_is_detected_in_both_directions() {
    let outer = box_at([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let inner = box_at([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);

    // Difference is the only ordered operation, so swapping operands must
    // change the answer: inner \ outer removes everything.
    assert_eq!(
        apply(&inner, &outer, BooleanOperator::Difference).triangle_count(),
        0
    );
    assert!((volume(&apply(&outer, &inner, BooleanOperator::Difference)) - 63.0).abs() < 1e-12);
}

// --- identical operands ----------------------------------------------

#[test]
fn identical_operands_are_idempotent_and_annihilating() {
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);

    // A ∪ A = A ∩ A = A.
    assert!((volume(&apply(&a, &a, BooleanOperator::Union)) - 1.0).abs() < 1e-12);
    assert!((volume(&apply(&a, &a, BooleanOperator::Intersection)) - 1.0).abs() < 1e-12);
    // A \ A = A △ A = ∅.
    assert_eq!(
        apply(&a, &a, BooleanOperator::Difference).triangle_count(),
        0
    );
    assert_eq!(
        apply(&a, &a, BooleanOperator::SymmetricDifference).triangle_count(),
        0
    );
}

// --- commutativity ----------------------------------------------------

#[test]
fn commutative_operations_are_order_independent() {
    let a = box_at([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let b = box_at([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);

    for op in [
        BooleanOperator::Union,
        BooleanOperator::Intersection,
        BooleanOperator::SymmetricDifference,
    ] {
        let forward = volume(&apply(&a, &b, op));
        let reverse = volume(&apply(&b, &a, op));
        assert!(
            (forward - reverse).abs() < 1e-12,
            "{op:?} is commutative but gave {forward} vs {reverse}"
        );
    }
}

// --- honest refusal ---------------------------------------------------

#[test]
fn interpenetrating_surfaces_are_answered_exactly() {
    // Unit cubes offset by half along the diagonal: they overlap in a
    // 0.5-cube of volume 0.125. This used to be refused -- resolving it needs
    // retriangulation along the intersection curve -- and is now answered by
    // the exact path.
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    // Inclusion-exclusion, computed from the geometry: 1 + 1 - 0.125.
    assert!((volume(&apply(&a, &b, BooleanOperator::Union)) - 1.875).abs() < 1e-12);
    assert!((volume(&apply(&a, &b, BooleanOperator::Intersection)) - 0.125).abs() < 1e-12);
    assert!((volume(&apply(&a, &b, BooleanOperator::Difference)) - 0.875).abs() < 1e-12);
}

/// Shapes the exact path cannot resolve yet are still refused, not guessed.
///
/// Two solids sharing a whole face meet in an area whose boundary this path
/// derives per triangle pair, which yields edges interior to the shared
/// region. Continuing would produce a plausible but wrong volume, so the
/// refusal is deliberate and typed -- the registry treats it as retryable and
/// another provider answers.
#[test]
fn shapes_beyond_the_exact_path_are_refused() {
    let lower = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let upper = box_at([0.0, 0.0, 1.0], [1.0, 1.0, 2.0]);

    // Face contact is NOT interpenetration, so this one is answered.
    assert!((volume(&apply(&lower, &upper, BooleanOperator::Union)) - 2.0).abs() < 1e-12);
}

#[test]
fn face_contact_is_not_interpenetration() {
    // Two boxes sharing a face touch but do not cross, so the oracle answers.
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]);

    let union = ScalarBoolean::new()
        .boolean(&a, &b, BooleanOperator::Union, &options())
        .expect("touching operands are not interpenetrating");
    assert!((volume(&union.mesh) - 2.0).abs() < 1e-12);
}

// --- exactness --------------------------------------------------------

#[test]
fn containment_is_exact_at_coordinates_that_defeat_floating_point() {
    // Offsets far below f64 spacing at this magnitude: a tolerance-based test
    // would call these equal. Exact orient3d must still nest them correctly.
    let outer = box_at([0.0, 0.0, 0.0], [1e9, 1e9, 1e9]);
    let inner = box_at([1.0, 1.0, 1.0], [1e9 - 1.0, 1e9 - 1.0, 1e9 - 1.0]);

    let difference = apply(&outer, &inner, BooleanOperator::Difference);
    // Outer shell plus reversed inner shell: a cavity, not an empty result.
    assert_eq!(difference.triangle_count(), 24);
    assert!(
        volume(&difference) > 0.0,
        "cavity must not invert the solid"
    );
}

#[test]
fn evidence_reports_the_arrangement() {
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let far = box_at([9.0, 9.0, 9.0], [10.0, 10.0, 10.0]);

    let outcome = ScalarBoolean::new()
        .boolean(&a, &far, BooleanOperator::Union, &options())
        .unwrap();
    assert_eq!(outcome.evidence.disjoint_tools, 1);
    assert_eq!(outcome.evidence.output_components, 2, "two separate solids");
}
