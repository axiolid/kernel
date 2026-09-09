//! `union_many` throughput against the sequential-fold baseline.
//!
//! The override claims a COMPLEXITY win, not a constant factor: the fold
//! makes step `i` union an accumulator already holding `i` solids, while the
//! tree keeps operands small until the final levels. A complexity claim is
//! only credible if the ratio GROWS with n, so this sweeps n and prints the
//! speedup per size rather than quoting one number.
//!
//! Volumes are compared at every size: a faster path that unions to a
//! different solid has not won anything.

use std::time::Instant;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;

/// Axis-aligned box as a closed, outward-oriented triangle mesh.
fn boxx(cx: f64, cy: f64, cz: f64, s: f64) -> TriMesh {
    let h = s / 2.0;
    let mut positions = Vec::with_capacity(8);
    for &(dx, dy) in &[(-h, -h), (h, -h), (h, h), (-h, h)] {
        for &dz in &[-h, h] {
            positions.push(Point3::new(cx + dx, cy + dy, cz + dz));
        }
    }
    let indices = vec![
        0, 4, 2, 0, 6, 4, 1, 3, 5, 1, 5, 7, 0, 3, 1, 0, 2, 3, 2, 5, 3, 2, 4, 5, 4, 7, 5, 4, 6, 7,
        6, 1, 7, 6, 0, 1,
    ];
    TriMesh::new(positions, indices)
}

/// Enclosed volume by the divergence theorem.
fn volume(mesh: &TriMesh) -> f64 {
    let mut total = 0.0;
    for corner in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[corner[0] as usize];
        let b = mesh.positions[corner[1] as usize];
        let c = mesh.positions[corner[2] as usize];
        total += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (total / 6.0).abs()
}

/// A `k` x `k` x `k` grid of boxes, spaced by `pitch` box widths.
///
/// Boxes rather than icospheres: the provider crate has no sphere builder,
/// and the reduction-order effect is a property of operand SIZE growth, not
/// of the operand's curvature.
fn grid(k: usize, pitch: f64) -> Vec<TriMesh> {
    let mut out = Vec::with_capacity(k * k * k);
    for i in 0..k {
        for j in 0..k {
            for l in 0..k {
                out.push(boxx(
                    i as f64 * pitch,
                    j as f64 * pitch,
                    l as f64 * pitch,
                    1.0,
                ));
            }
        }
    }
    out
}

/// Best-of-3 wall clock, in milliseconds, plus the result volume.
fn timed(mut run: impl FnMut() -> TriMesh) -> (f64, f64) {
    let mut best = f64::MAX;
    let mut result = run();
    for _ in 0..3 {
        let start = Instant::now();
        result = run();
        best = best.min(start.elapsed().as_secs_f64() * 1e3);
    }
    (best, volume(&result))
}

/// The trait default: a serial left fold.
fn sequential(
    provider: &BoolmeshBoolean,
    solids: &[TriMesh],
    options: &ExecutionOptions,
) -> TriMesh {
    let mut current = solids[0].clone();
    for solid in &solids[1..] {
        current = provider
            .boolean(&current, solid, BooleanOperator::Union, options)
            .expect("sequential union")
            .mesh;
    }
    current
}

fn main() {
    let provider = BoolmeshBoolean::new();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);

    // 3.0 leaves clear air between boxes; 0.6 makes every neighbour
    // overlap, so each union does real cutting work.
    // Sizes are per-case. The overlapping sweep stops at 64: a 125-box
    // overlapping grid trips an assert inside the absorbed kernel
    // (`boolean45.rs`'s `pair_up`, odd edge-point count). That fault is
    // PRE-EXISTING and unrelated to reduction order -- verified by
    // reproducing it on the sequential fold, which this override does not
    // touch. Benching up to the cliff measures what is measurable without
    // pretending the cliff is not there.
    let cases: [(&str, f64, &[usize]); 2] = [
        ("disjoint grid (multi-component result)", 3.0, &[2, 3, 4, 5]),
        ("overlapping grid (fuses into one solid)", 0.6, &[2, 3, 4]),
    ];

    for (label, pitch, sizes) in cases {
        println!("\n{label}");
        println!(
            "{:>5}  {:>13}  {:>13}  {:>9}  {:>10}",
            "n", "sequential", "tree", "speedup", "volumes"
        );
        for &k in sizes {
            let solids = grid(k, pitch);
            let n = solids.len();
            let (seq_ms, seq_vol) = timed(|| sequential(&provider, &solids, &options));
            let (tree_ms, tree_vol) = timed(|| {
                provider
                    .union_many(&solids, &options)
                    .expect("tree union")
                    .mesh
            });
            // The claim is an optimisation, not a behaviour change: if the
            // volumes disagree the speedup is meaningless.
            let agree = (seq_vol - tree_vol).abs() <= 1e-9 * seq_vol.abs().max(1.0);
            println!(
                "{n:>5}  {seq_ms:>11.2}ms  {tree_ms:>11.2}ms  {:>8.1}x  {:>10}",
                seq_ms / tree_ms,
                if agree { "agree" } else { "DISAGREE" }
            );
            assert!(agree, "n={n}: seq {seq_vol} vs tree {tree_vol}");
        }
    }
}
