//! Gates for the `union_many` batch override.
//!
//! The override's entire justification is that it is FASTER while producing
//! the SAME solid. These tests hold the second half: the tree reduction must
//! agree with the sequential fold on every layout it will meet.
//!
//! Union is associative and commutative, so unlike `subtract_many` there is
//! no disjointness precondition and no fusing -- but that makes the gates
//! MORE important, not less: an order-independent operation that turns out to
//! be order-dependent in practice is exactly the bug worth catching.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use support::{boxx, volume};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

/// The trait-default behaviour, kept as the reference implementation.
fn sequential(solids: &[TriMesh]) -> TriMesh {
    let provider = BoolmeshBoolean::new();
    let mut current = solids[0].clone();
    for solid in &solids[1..] {
        current = provider
            .boolean(&current, solid, BooleanOperator::Union, &options())
            .expect("sequential union")
            .mesh;
    }
    current
}

/// Volumes must agree to a relative tolerance, not bitwise: the two paths
/// sum a differently ordered triangle list, so the last bits legitimately
/// differ. Bitwise equality here would be testing the float summation order,
/// not the geometry.
fn assert_same_volume(left: &TriMesh, right: &TriMesh, what: &str) {
    let (a, b) = (volume(left), volume(right));
    assert!(
        (a - b).abs() <= 1e-9 * a.abs().max(1.0),
        "{what}: tree volume {b} disagrees with sequential {a}"
    );
}

/// A row of `n` boxes, each overlapping its neighbour.
fn overlapping_row(n: usize) -> Vec<TriMesh> {
    (0..n)
        .map(|i| boxx(i as f64 * 0.6, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0))
        .collect()
}

/// A row of `n` boxes with clear air between them.
fn disjoint_row(n: usize) -> Vec<TriMesh> {
    (0..n)
        .map(|i| boxx(i as f64 * 3.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0))
        .collect()
}

/// Overlapping operands: every union does real cutting work.
#[test]
fn overlapping_solids_agree_with_the_sequential_path() {
    let solids = overlapping_row(7);
    let expected = sequential(&solids);
    let actual = BoolmeshBoolean::new()
        .union_many(&solids, &options())
        .expect("tree union")
        .mesh;

    assert_same_volume(&expected, &actual, "overlapping");
    assert!(
        volume(&actual) > volume(&solids[0]),
        "a union of overlapping solids must exceed any single one"
    );
}

/// Disjoint operands: the result is multi-component, which is the shape that
/// has historically been mishandled -- see the determinism notes in the
/// `benchmarks` sibling repo.
#[test]
fn disjoint_solids_agree_with_the_sequential_path() {
    let solids = disjoint_row(6);
    let expected = sequential(&solids);
    let actual = BoolmeshBoolean::new()
        .union_many(&solids, &options())
        .expect("tree union")
        .mesh;

    assert_same_volume(&expected, &actual, "disjoint");
    // Six unit cubes that never touch: the union is exactly their total.
    let want: f64 = solids.iter().map(volume).sum();
    assert!(
        (volume(&actual) - want).abs() <= 1e-9 * want,
        "disjoint union must preserve total volume: got {} want {want}",
        volume(&actual)
    );
}

/// An ODD count exercises the remainder branch: the trailing solid must ride
/// to the next level uncombined, not be dropped and not be double-counted.
/// An off-by-one there is invisible at even counts.
#[test]
fn odd_counts_do_not_drop_the_trailing_solid() {
    for n in [1usize, 3, 5, 7, 9, 11] {
        let solids = disjoint_row(n);
        let actual = BoolmeshBoolean::new()
            .union_many(&solids, &options())
            .expect("tree union")
            .mesh;
        let want: f64 = solids.iter().map(volume).sum();
        assert!(
            (volume(&actual) - want).abs() <= 1e-9 * want,
            "n={n}: got {} want {want}",
            volume(&actual)
        );
    }
}

/// The union of nothing is nothing -- a legitimate answer, not an error.
#[test]
fn an_empty_batch_yields_an_empty_solid() {
    let outcome = BoolmeshBoolean::new()
        .union_many(&[], &options())
        .expect("empty batch is a valid request");
    assert!(outcome.mesh.indices.is_empty(), "must be empty");
}

/// A single solid must come back unchanged, having done no boolean at all.
#[test]
fn a_single_solid_is_returned_without_a_boolean() {
    let solid = boxx(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0);
    let outcome = BoolmeshBoolean::new()
        .union_many(std::slice::from_ref(&solid), &options())
        .expect("single");
    assert_eq!(
        outcome.evidence.sub_operations, 0,
        "one solid needs no boolean"
    );
    let want = volume(&solid);
    assert!((volume(&outcome.mesh) - want).abs() <= 1e-12 * want);
}

/// Evidence must report the booleans actually performed. A tree over n
/// solids issues exactly n-1 unions, same as the fold -- the win is operand
/// SIZE, not call count, and evidence that claimed otherwise would hide that.
#[test]
fn evidence_reports_one_boolean_per_merge() {
    for n in [2usize, 3, 4, 8, 9] {
        let outcome = BoolmeshBoolean::new()
            .union_many(&disjoint_row(n), &options())
            .expect("tree union");
        assert_eq!(
            outcome.evidence.sub_operations,
            n - 1,
            "n={n}: a union of n solids takes n-1 booleans"
        );
    }
}

/// Operand ORDER must not change the answer. Union is commutative, so a
/// reversed input list must give the same solid -- this is the property the
/// tree reduction relies on to be free to regroup at all.
#[test]
fn reversing_the_operands_gives_the_same_solid() {
    let solids = overlapping_row(6);
    let mut reversed = solids.clone();
    reversed.reverse();

    let provider = BoolmeshBoolean::new();
    let forward = provider.union_many(&solids, &options()).expect("fwd").mesh;
    let backward = provider
        .union_many(&reversed, &options())
        .expect("rev")
        .mesh;
    assert_same_volume(&forward, &backward, "reversed operands");
}
