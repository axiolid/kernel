// SPDX-License-Identifier: MPL-2.0

//! Precision characterisation for the production boolean.
//!
//! # Why these are regression tests, not a bug report
//!
//! These began as a probe for "gap 2": an audit claimed the production
//! boolean needed exact predicates because it computes intersections in f64.
//! The measurements below refuted that. At origin scale the boolean is
//! accurate to one ULP across fifteen decades of face separation, and six of
//! its ten sign decisions are direct coordinate comparisons that were exact
//! already.
//!
//! What the probes DID find is a real limit, in a different place: geometry
//! far from the origin loses relative precision before the boolean is ever
//! called. At base 1e7 a 0.1m box's own volume is 25% wrong, because the
//! divergence-theorem sum adds terms of order 1e21 to produce 1e-3.
//!
//! Both facts are asserted here so neither can regress silently, and so the
//! documented limit is a measurement rather than folklore.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use support::{boxx, volume};

/// At origin scale the boolean resolves face separations down to 1e-15
/// without losing a slab of volume.
#[test]
fn origin_scale_booleans_are_accurate_to_a_few_ulp() {
    let provider = BoolmeshBoolean::default();
    let options = ExecutionOptions::new(Tolerance::METRE);
    let host = boxx(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0);

    for k in 1..16 {
        let eps = 10f64.powi(-k);
        let tool = boxx(eps + 2.0, 0.0, -1.0, 4.0, 4.0, 4.0, 0.0);
        let got = provider
            .boolean(&host, &tool, BooleanOperator::Difference, &options)
            .map(|o| volume(&o.mesh))
            .expect("difference of two boxes");
        let want = 0.5 + eps;
        let err = (got - want).abs();
        // A wrongly decided predicate loses a whole slab (order eps or 0.5),
        // not a few ULP. This bound is ~5 ULP at this magnitude.
        assert!(
            err < 1e-15,
            "eps=1e-{k}: got {got:.17}, want {want:.17}, err {err:.3e} -- \
             a predicate decided wrongly, not a rounding artefact"
        );
    }
}

/// Coordinate magnitude, not the boolean, is what breaks far from the origin.
///
/// Asserts the ATTRIBUTION: the input mesh is already wrong before any
/// boolean runs, so a future change must not "fix the boolean" and claim the
/// large-coordinate case is solved.
#[test]
fn large_coordinate_error_enters_before_the_boolean() {
    // Exact at origin scale.
    let near = boxx(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.0);
    let near_err = ((volume(&near) - 0.001) / 0.001).abs();
    assert!(
        near_err < 1e-12,
        "origin-scale input mesh should measure exactly, got rel err {near_err:.3e}"
    );

    // Already badly wrong at national-grid scale, with NO boolean involved.
    let far = boxx(1.0e7, 1.0e7, 0.0, 0.1, 0.1, 0.1, 0.0);
    let far_err = ((volume(&far) - 0.001) / 0.001).abs();
    assert!(
        far_err > 1e-3,
        "expected large-coordinate degradation to be present and measurable, \
         got rel err {far_err:.3e} -- if this now passes, an RTC/local-origin \
         facility has landed and this test should assert the new behaviour"
    );
}
