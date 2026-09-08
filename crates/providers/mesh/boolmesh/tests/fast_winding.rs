//! The opt-in `boolean_fast` path.
//!
//! Sphere fixtures already cover the winding classification directly in
//! kernel03's unit tests. This file checks the public entry point end to
//! end: boolean_fast agrees with boolean on real corpus fixtures, and it
//! refuses SymmetricDifference instead of silently composing three slow
//! calls.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_fixtures::corpus;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::{boxx, volume};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

#[test]
fn agrees_with_the_default_path_on_the_fixture_corpus() {
    let tool = boxx(0.3, 0.3, 0.1, 1.2, 1.2, 1.2, 0.0);
    let provider = BoolmeshBoolean::new();
    for fixture in corpus() {
        for op in [
            BooleanOperator::Union,
            BooleanOperator::Intersection,
            BooleanOperator::Difference,
        ] {
            let slow = provider.boolean(&fixture.mesh, &tool, op, &options());
            let fast = provider.boolean_fast(&fixture.mesh, &tool, op, &options());
            match (slow, fast) {
                (Ok(s), Ok(f)) => {
                    let vs = volume(&s.mesh);
                    let vf = volume(&f.mesh);
                    assert!(
                        (vs - vf).abs() < 1e-9 * vs.abs().max(1.0),
                        "{}: {op:?} slow={vs} fast={vf}",
                        fixture.name
                    );
                }
                (Err(_), Err(_)) => {}
                (s, f) => panic!(
                    "{}: {op:?} disagreed on refusal: slow={s:?} fast={f:?}",
                    fixture.name
                ),
            }
        }
    }
}

#[test]
fn refuses_symmetric_difference_rather_than_composing_slow_calls() {
    let subject = boxx(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 0.0);
    let tool = boxx(0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 0.0);
    let provider = BoolmeshBoolean::new();
    let result = provider.boolean_fast(
        &subject,
        &tool,
        BooleanOperator::SymmetricDifference,
        &options(),
    );
    assert!(
        result.is_err(),
        "boolean_fast must refuse SymmetricDifference, not silently compose it"
    );
}
