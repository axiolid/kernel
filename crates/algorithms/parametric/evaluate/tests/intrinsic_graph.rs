//! `Curve3::Intrinsic` resolves through the generic curve API (ADR 0062).
//!
//! Wiring the Frenet integrator into `evaluate3`/`derivative3` is what
//! lets EXISTING generic machinery -- trimming by parameter range,
//! composite stitching, sweep directrix sampling -- work on a torsion
//! curve without any of it special-casing the family. These tests pin
//! that dispatch, since a silent `Unsupported` would look like a missing
//! feature rather than a regression.

use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{CurvatureLaw, Curve3, Intrinsic3};
use axiolid_evaluate::{derivative3, evaluate3, frenet_point};

fn start_frame() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: Vec3::new(1.0, 0.0, 0.0),
        y: Vec3::new(0.0, 1.0, 0.0),
        z: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn helix(length: Scalar) -> Intrinsic3 {
    let (a, b): (Scalar, Scalar) = (3.0, 1.5);
    let c = a.hypot(b);
    Intrinsic3::new(
        start_frame(),
        CurvatureLaw::circular(a / (c * c)),
        CurvatureLaw::circular(b / (c * c)),
        length,
    )
}

#[test]
fn the_generic_curve_api_dispatches_a_torsion_curve() {
    // Before this wiring, evaluate3 refused Intrinsic by name and every
    // generic consumer inherited that refusal.
    let curve = Curve3::Intrinsic(helix(12.0));
    let Curve3::Intrinsic(inner) = &curve else {
        unreachable!()
    };

    for step in 0..=6 {
        let s = 12.0 * Scalar::from(step) / 6.0;
        let generic = evaluate3(&curve, s).expect("generic evaluation");
        let direct = frenet_point(inner, s).expect("direct evaluation");
        assert!(
            (generic - direct).length() < 1e-12,
            "generic API must be the same value as the integrator at s={s}"
        );
    }
}

#[test]
fn an_arc_length_parameterised_curve_reports_its_own_domain() {
    // The domain is [0, length] because the parameter IS arc length -- not
    // the unit interval a generic consumer would otherwise assume. A
    // consumer that trims by parameter range depends on this being right.
    let curve = Curve3::Intrinsic(helix(12.0));
    let domain = axiolid_evaluate::curve::domain3(&curve);
    assert!((domain.start - 0.0).abs() < 1e-12);
    assert!((domain.end - 12.0).abs() < 1e-12);
}

#[test]
fn the_derivative_of_an_arc_length_curve_is_a_unit_tangent() {
    // Arc-length parameterisation means |dP/ds| = 1 exactly. A consumer
    // normalising the derivative would hide a scaling bug; asserting the
    // norm catches it.
    let curve = Curve3::Intrinsic(helix(12.0));
    for step in 0..=6 {
        let s = 12.0 * Scalar::from(step) / 6.0;
        let tangent = derivative3(&curve, s).expect("generic derivative");
        let norm = tangent.length();
        assert!((norm - 1.0).abs() < 1e-9, "|T|={norm} at s={s}");
    }
}

#[test]
fn a_length_that_claims_no_domain_is_reported_as_empty() {
    // A non-positive or non-finite length has no knowable domain. Claiming
    // an empty one is honest; inventing [0, 1] would invite a consumer to
    // sample a curve that does not exist.
    for length in [0.0, -3.0, Scalar::INFINITY, Scalar::NAN] {
        let curve = Curve3::Intrinsic(helix(length));
        let domain = axiolid_evaluate::curve::domain3(&curve);
        assert!(
            domain.start == 0.0 && domain.end == 0.0,
            "length {length} must claim no domain"
        );
    }
}
