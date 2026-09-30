//! Unions of overlapping axis-aligned boxes complete exactly (#203).
//!
//! A k x k x k grid of unit boxes at a pitch below one overlaps into a
//! single cube of side `(k - 1) * pitch + 1`. Every partial union shares
//! whole planes with the next box, so these folds put operand vertices
//! exactly on the other operand's extreme planes over and over. Until
//! #203 the winding number's broad phase dropped such a vertex whenever
//! its rounded grid bound fell one ulp short of the true box, and the
//! solve refused with an odd edge-point count at pitches 0.65, 0.7 and
//! 0.8.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use std::collections::BTreeMap;
use support::{boxx, volume};

const PITCHES: [f64; 5] = [0.6, 0.65, 0.7, 0.75, 0.8];

fn grid(k: usize, pitch: f64) -> Vec<TriMesh> {
    let mut boxes = Vec::with_capacity(k * k * k);
    for i in 0..k {
        for j in 0..k {
            for l in 0..k {
                let (x, y, z) = (i as f64 * pitch, j as f64 * pitch, l as f64 * pitch);
                boxes.push(boxx(x, y, z - 0.5, 1.0, 1.0, 1.0, 0.0));
            }
        }
    }
    boxes
}

/// The union is the cube the grid spans, closed and consistently wound.
fn assert_is_the_spanned_cube(mesh: &TriMesh, k: usize, pitch: f64, how: &str) {
    let side = (k - 1) as f64 * pitch + 1.0;
    let expected = side * side * side;
    let got = volume(mesh);
    assert!(
        (got - expected).abs() <= 1e-9 * expected,
        "{how} k={k} pitch={pitch}: volume {got}, expected {expected}"
    );
    // Topology over every triangle, zero-area ones included: each directed
    // edge occurs once and is matched by its reverse. `audit_mesh` would
    // skip the collinear slivers the triangulator leaves on subdivided
    // faces and report their edges as open, which is not what this checks.
    let mut directed: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((a, b)).or_default() += 1;
        }
    }
    for (&(a, b), &count) in &directed {
        assert!(
            count == 1 && directed.get(&(b, a)) == Some(&1),
            "{how} k={k} pitch={pitch}: edge {a}-{b} is not closed and consistently wound"
        );
    }
}

#[test]
fn a_sequential_fold_of_overlapping_boxes_is_the_spanned_cube() {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    for k in [5, 6] {
        for pitch in PITCHES {
            let mut boxes = grid(k, pitch).into_iter();
            let mut acc = boxes.next().expect("non-empty grid");
            for (step, solid) in boxes.enumerate() {
                acc = provider
                    .boolean(&acc, &solid, BooleanOperator::Union, &options)
                    .unwrap_or_else(|e| panic!("k={k} pitch={pitch} step {}: {e}", step + 1))
                    .mesh;
            }
            assert_is_the_spanned_cube(&acc, k, pitch, "fold");
        }
    }
}

#[test]
fn union_many_of_overlapping_boxes_is_the_spanned_cube() {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    for k in [5, 6] {
        for pitch in PITCHES {
            let outcome = provider
                .union_many(&grid(k, pitch), &options)
                .unwrap_or_else(|e| panic!("k={k} pitch={pitch}: {e}"));
            assert_is_the_spanned_cube(&outcome.mesh, k, pitch, "union_many");
        }
    }
}

/// The opt-in fast winding path shares the broad phase, so it gets the
/// same check.
#[test]
fn the_fast_winding_fold_is_the_spanned_cube() {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    for pitch in PITCHES {
        let mut boxes = grid(5, pitch).into_iter();
        let mut acc = boxes.next().expect("non-empty grid");
        for (step, solid) in boxes.enumerate() {
            acc = provider
                .boolean_fast(&acc, &solid, BooleanOperator::Union, &options)
                .unwrap_or_else(|e| panic!("pitch={pitch} step {}: {e}", step + 1))
                .mesh;
        }
        assert_is_the_spanned_cube(&acc, 5, pitch, "fast fold");
    }
}
