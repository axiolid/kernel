//! The geometric audit must accept sound solids and reject inconsistent ones.
//!
//! The motivating defect is real: an earlier arc overlay built cap loops
//! whose pcurves were straight chords across arc edges. That solid closed,
//! validated and audited clean topologically. These tests pin the behaviour
//! that catches it.

use axiolid_brep_audit::{geometric_audit, GeometricDefect};
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::feature::fillet_polygon_corners;
use axiolid_core::{BooleanOperator, Point2, Tolerance, Vec3};
use axiolid_overlay::ArcRing;
use axiolid_profile::{Profile, RectangleProfile};

fn rect() -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 3.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

#[test]
fn a_plain_extrusion_is_geometrically_consistent() {
    let solid = extrude_profile_exact(&rect(), Vec3::Z, 4.0, Tolerance::METRE).expect("extrusion");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "a box must audit clean, found {:?}",
        health.defects()
    );
}

#[test]
fn a_curved_boolean_result_is_geometrically_consistent() {
    // The arc path is where a pcurve/3D disagreement actually appeared
    // before, so it is the case most worth auditing.
    let cylinder = ArcPrism {
        section: ArcRing::circle(Point2::new(0.0, 0.0), 1.0),
        bottom: 0.0,
        top: 2.0,
    };
    let box_solid = ArcPrism {
        section: ArcRing {
            vertices: vec![
                axiolid_overlay::ArcVertex::straight(Point2::new(-0.8, -0.8)),
                axiolid_overlay::ArcVertex::straight(Point2::new(0.8, -0.8)),
                axiolid_overlay::ArcVertex::straight(Point2::new(0.8, 0.8)),
                axiolid_overlay::ArcVertex::straight(Point2::new(-0.8, 0.8)),
            ],
        },
        bottom: 0.0,
        top: 2.0,
    };
    let solid = boolean_arc_prisms_exact(
        &cylinder,
        &box_solid,
        BooleanOperator::Intersection,
        Tolerance::METRE,
    )
    .expect("curved boolean");

    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "the curved boolean result must audit clean, found {:?} worst {:?}",
        health.defects(),
        health.worst_error()
    );
}

#[test]
fn a_filleted_solid_is_geometrically_consistent() {
    let ring = vec![
        Point2::new(-1.0, -1.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 1.0),
    ];
    let solid = fillet_polygon_corners(&ring, &[(0, 0.3), (2, 0.4)], 1.0).expect("fillets");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "cylindrical blends must audit clean, found {:?}",
        health.defects()
    );
}

#[test]
fn the_audit_reports_a_worst_error_only_when_something_is_wrong() {
    let solid = extrude_profile_exact(&rect(), Vec3::Z, 4.0, Tolerance::METRE).expect("extrusion");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert_eq!(health.worst_error(), None);
    assert!(health.defects().is_empty());
    let _ = GeometricDefect::VertexOffCurve {
        edge: 0,
        error: 0.0,
    };
}

#[test]
fn the_boolean_gate_is_wired_and_can_reject() {
    // Proving the gate is live: run a boolean at a tolerance so tight that
    // the exact-but-finite-precision result cannot satisfy it. If the gate
    // were disconnected this would return Ok regardless.
    //
    // A gate nobody can trip is not a gate, and "all booleans pass" is
    // equally consistent with "the gate was never called".
    use axiolid_core::Tolerance;

    let tight = Tolerance::new(1e-18, 1e-18).expect("a valid, very tight tolerance");
    let cylinder = ArcPrism {
        section: ArcRing::circle(Point2::new(0.0, 0.0), 1.0),
        bottom: 0.0,
        top: 2.0,
    };
    let box_solid = ArcPrism {
        section: ArcRing {
            vertices: vec![
                axiolid_overlay::ArcVertex::straight(Point2::new(-0.8, -0.8)),
                axiolid_overlay::ArcVertex::straight(Point2::new(0.8, -0.8)),
                axiolid_overlay::ArcVertex::straight(Point2::new(0.8, 0.8)),
                axiolid_overlay::ArcVertex::straight(Point2::new(-0.8, 0.8)),
            ],
        },
        bottom: 0.0,
        top: 2.0,
    };

    let result =
        boolean_arc_prisms_exact(&cylinder, &box_solid, BooleanOperator::Intersection, tight);

    // Measured: the gate rejects this with 24 defects, worst deviation
    // 3.14e-16 -- real floating-point residue that 1e-18 cannot admit.
    // Asserting the rejection (rather than tolerating Ok) is what makes
    // this test fail if the gate is ever disconnected.
    let error = result.expect_err("the gate must reject at an impossible tolerance");
    let text = format!("{error:?}");
    assert!(
        text.contains("failed its geometric audit"),
        "the refusal must name the audit, got {text}"
    );
    // The reported deviation is evidence the audit measured rather than
    // guessed, so it must be present and finite.
    assert!(
        text.contains("worst deviation"),
        "the refusal must report how far off it was, got {text}"
    );
}
