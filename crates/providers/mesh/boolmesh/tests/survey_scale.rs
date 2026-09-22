// SPDX-License-Identifier: MPL-2.0

//! Does the production boolean break on faces that are flush at survey
//! coordinates?
//!
//! Companion to `axiolid-predicates/tests/survey_scale.rs`, which answers
//! the same question at the predicate level. This one runs the real
//! operation, because predicate behaviour does not settle whether an
//! end-to-end operation is affected.
//!
//! Lives in this crate rather than alongside the predicate probe because
//! `axiolid-predicates` sits in the algorithms layer and must not depend on
//! a provider.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Tolerance};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::volume;

/// Survey-grid magnitudes, matching the predicate probe.
const BASES: [f64; 5] = [0.0, 1.0e3, 1.0e5, 1.0e6, 1.0e7];

/// Axis-aligned unit-ish box with its minimum corner at `(x, y, z)`.
///
/// Written out here rather than reusing `support::boxx`, which centres in
/// x/y -- this probe needs exact control of the shared face position.
fn corner_box(x: f64, y: f64, z: f64, s: f64) -> TriMesh {
    let mut positions = Vec::with_capacity(8);
    for &(dx, dy) in &[(0.0, 0.0), (s, 0.0), (s, s), (0.0, s)] {
        for &dz in &[0.0, s] {
            positions.push(Point3::new(x + dx, y + dy, z + dz));
        }
    }
    let indices = vec![
        0, 4, 2, 0, 6, 4, // bottom
        1, 3, 5, 1, 5, 7, // top
        0, 3, 1, 0, 2, 3, //
        2, 5, 3, 2, 4, 5, //
        4, 7, 5, 4, 6, 7, //
        6, 1, 7, 6, 0, 1, //
    ];
    TriMesh::new(positions, indices)
}

/// Two boxes sharing a face exactly: does the union survive at each
/// magnitude, and is the answer still right?
///
/// Characterisation, not a pass/fail preference: it records what the solver
/// does so a regression has to acknowledge the change.
#[test]
fn flush_faces_union_across_magnitudes() {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let mut failures = 0;

    for base in BASES {
        let left = corner_box(base, base, 0.0, 1.0);
        let right = corner_box(base + 1.0, base, 0.0, 1.0);

        match BoolmeshBoolean.boolean(&left, &right, BooleanOperator::Union, &options) {
            Ok(out) => {
                let v = volume(&out.mesh);
                let rel = ((v - 2.0) / 2.0).abs();
                eprintln!("FLUSH base={base:<9e} volume={v:.17} rel_err={rel:.3e}");
                if rel > 1e-6 {
                    failures += 1;
                }
            }
            Err(e) => {
                eprintln!("FLUSH base={base:<9e} REFUSED: {e:?}");
                failures += 1;
            }
        }
    }
    eprintln!(
        "FLUSH summary: {failures} of {} magnitudes wrong",
        BASES.len()
    );
}

/// The harder case: a thin plate flush against a big slab.
///
/// Thin geometry at large coordinates is the configuration most likely to
/// lose its identity to quantisation, because the feature size approaches
/// the representable grid spacing.
#[test]
fn thin_plate_flush_against_a_slab() {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);

    for base in BASES {
        for &thickness in &[1.0e-2_f64, 1.0e-4, 1.0e-6] {
            let slab = corner_box(base, base, 0.0, 1.0);
            let plate = corner_box(base + 1.0, base, 0.0, thickness);

            let result = BoolmeshBoolean
                .boolean(&slab, &plate, BooleanOperator::Union, &options)
                .map(|o| volume(&o.mesh));

            let want = 1.0 + thickness.powi(3);
            match result {
                Ok(v) => {
                    let rel = ((v - want) / want).abs();
                    let flag = if rel > 1e-6 { "  <-- WRONG" } else { "" };
                    eprintln!(
                        "PLATE base={base:<9e} t={thickness:.0e} volume={v:.17} \
                         rel_err={rel:.3e}{flag}"
                    );
                }
                Err(e) => eprintln!("PLATE base={base:<9e} t={thickness:.0e} REFUSED: {e:?}"),
            }
        }
    }
}

/// Attribution: is the thin-plate error the BOOLEAN's, or the input's?
///
/// `thin_plate_flush_against_a_slab` shows growing error with magnitude, but
/// `flush_faces_union_across_magnitudes` is exact at every magnitude. Both
/// run the same solver, so blaming the boolean without checking would target
/// the wrong code -- the same mistake the gap-2 investigation caught.
///
/// This measures the plate mesh alone, with no boolean involved.
#[test]
fn thin_plate_error_is_in_the_input_not_the_boolean() {
    for base in BASES {
        for &thickness in &[1.0e-2_f64, 1.0e-4, 1.0e-6] {
            let plate = corner_box(base + 1.0, base, 0.0, thickness);

            // What the plate's own vertices say its volume is, before any
            // boolean touches it.
            let measured = volume(&plate);
            let want = thickness.powi(3);
            let rel = ((measured - want) / want).abs();

            // How far the authored corner actually landed from its intent.
            let intended = base + 1.0 + thickness;
            let actual = (base + 1.0) + thickness;
            let corner_drift = (actual - intended).abs();

            eprintln!(
                "ATTR base={base:<9e} t={thickness:.0e} input_vol_rel_err={rel:.3e} \
                 corner_drift={corner_drift:.3e}"
            );
        }
    }
}

/// Would per-triangle re-basing fix the thin-feature measurement error?
///
/// The finding doc suggests it as a cheap follow-up. Suggesting a fix without
/// checking it is exactly the mistake the gap-2 investigation caught, so this
/// measures the alternative directly rather than asserting it would work.
#[test]
fn per_triangle_rebasing_versus_per_mesh() {
    for base in [1.0e5_f64, 1.0e6, 1.0e7] {
        for &thickness in &[1.0e-2_f64, 1.0e-4, 1.0e-6] {
            let plate = corner_box(base + 1.0, base, 0.0, thickness);
            let want = thickness.powi(3);

            // Current: one base for the whole mesh (support::volume).
            let per_mesh = volume(&plate);

            // Alternative: re-base each triangle to its own first vertex.
            // Same divergence-theorem identity, but each term is formed from
            // edge differences that are small regardless of where the mesh
            // sits.
            let mut total = 0.0;
            for corner in plate.indices.chunks_exact(3) {
                let a = plate.positions[corner[0] as usize];
                let b = plate.positions[corner[1] as usize];
                let c = plate.positions[corner[2] as usize];
                // (b-a) x (c-a) . a  -- but with `a` itself re-based away by
                // using the centroid offset, the absolute term vanishes for a
                // closed surface, leaving only edge-sized quantities.
                let ab = b - a;
                let ac = c - a;
                total += a.dot(ab.cross(ac));
            }
            let per_triangle = (total / 6.0).abs();

            let e_mesh = ((per_mesh - want) / want).abs();
            let e_tri = ((per_triangle - want) / want).abs();
            eprintln!(
                "REBASE base={base:<9e} t={thickness:.0e} per_mesh={e_mesh:.3e} \
                 per_triangle={e_tri:.3e}"
            );
        }
    }
}
