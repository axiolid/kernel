//! A failure inside the solve is the provider's defect, not the caller's.
//!
//! Every operand here passes the input gates: each box is a closed,
//! outward, manifold solid, and so is every partial union. When the solve
//! still refuses (an odd edge-point count in `pair_up`, #101), the error
//! must say so as `BackendContractViolation`, never `Degenerate`, which
//! would blame operands that were admissible.

mod support;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::boxx;

#[test]
fn a_refusal_inside_the_solve_is_a_backend_contract_violation() {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    // A 5 x 5 x 5 grid of overlapping unit boxes at pitch 0.8: the
    // sequential fold reaches the `pair_up` refusal (#101) at its sixth
    // union (#203). Pitches 0.65 and 0.7 fail too; 0.6 no longer does.
    let pitch = 0.8;
    let mut current: Option<axiolid_mesh::TriMesh> = None;
    let mut refusal = None;
    'fold: for i in 0..5 {
        for j in 0..5 {
            for l in 0..5 {
                let (x, y, z) = (i as f64 * pitch, j as f64 * pitch, l as f64 * pitch);
                let solid = boxx(x, y, z - 0.5, 1.0, 1.0, 1.0, 0.0);
                let Some(acc) = current.take() else {
                    current = Some(solid);
                    continue;
                };
                match provider.boolean(&acc, &solid, BooleanOperator::Union, &options) {
                    Ok(outcome) => current = Some(outcome.mesh),
                    Err(error) => {
                        refusal = Some(error);
                        break 'fold;
                    }
                }
            }
        }
    }
    let error = refusal.expect(
        "the overlapping grid no longer reaches a solve refusal; \
         pick a new reproducer so this mapping stays tested",
    );
    assert!(
        matches!(
            &error,
            GeomError::BackendContractViolation { backend, detail }
                if *backend == BoolmeshBoolean::ID && detail.contains("inside the solve")
        ),
        "{error:?}"
    );
}
