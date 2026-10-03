//! A plane touching a cylinder along a ruling, read once for every face
//! pair (#243).
//!
//! A wall `[-L/2, L/2] x [-T/2, T/2] x [0, H]` minus a round hole across it
//! (its axis along `-y`) touching the wall's top face: from inside (a hole
//! whose top ruling lies in the top face) or from outside (a tool resting
//! on it). Exactly tangent on the numbers given, the exact predicates decide
//! and nothing is read. Placed by a general rotation, or a fraction of the
//! tolerance into or off the face, the plane is read as touching
//! (`PlaneTouchesCylinder`, the distance it moves), and the hole's circles
//! in the wall's faces, its caps' chords and the top face's edges all meet
//! that one contact ruling: the result sews, audits clean and measures
//! the closed form.

use std::f64::consts::PI;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{boolean_with_report, BooleanReport, ToleranceDecisionKind};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Mat3, Tolerance, Transform3, Vec3};
use axiolid_measure::exact_properties;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

const L: f64 = 6.0;
const T: f64 = 0.25;
const H: f64 = 3.0;
/// The hole's radius and length (it overshoots both faces).
const R: f64 = 0.5;
const DEPTH: f64 = 1.0;
const EPS: f64 = 1e-6;

fn wall(p: Transform3) -> ExactBRep {
    let profile = Profile::Rectangle(RectangleProfile {
        x: L,
        y: T,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    });
    extrude_profile_exact(&profile, Vec3::Z, H, Tolerance::ZERO)
        .expect("a wall")
        .transformed(&p)
        .expect("rigid")
}

/// A round hole across the wall at `x = 0.5`, its axis at height `z`.
fn hole(z: f64, p: Transform3) -> ExactBRep {
    let profile = Profile::Circle(CircleProfile {
        radius: R,
        thickness: None,
    });
    let across = Transform3::from_mat3(Mat3::from_cols(Vec3::X, Vec3::Z, -Vec3::Y));
    extrude_profile_exact(&profile, Vec3::Z, DEPTH, Tolerance::ZERO)
        .expect("a hole")
        .transformed(&(p * Transform3::from_translation(Vec3::new(0.5, DEPTH / 2.0, z)) * across))
        .expect("rigid")
}

/// A rigid placement with no axis left aligned.
fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

/// The wall minus the hole touching its top face from inside or outside,
/// `gap` higher: audited, measured against `expected` within `slack`.
fn cut(
    inside: bool,
    gap: f64,
    p: Transform3,
    tolerance: Tolerance,
    (expected, slack): (f64, f64),
) -> BooleanReport {
    let z = if inside { H - R } else { H + R } + gap;
    let (result, report) = boolean_with_report(
        &wall(p),
        &hole(z, p),
        BooleanOperator::Difference,
        tolerance,
    )
    .unwrap_or_else(|e| panic!("inside {inside} gap {gap}: {e}"));
    let health = geometric_audit(&result, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    let measured = exact_properties(&result, Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    assert!(
        (measured - expected).abs() <= slack + 1e-12 * expected,
        "inside {inside} gap {gap}: volume {measured}, expected {expected}"
    );
    report
}

fn expected(inside: bool) -> f64 {
    if inside {
        L * T * H - PI * R * R * T
    } else {
        L * T * H
    }
}

#[test]
fn an_exactly_tangent_hole_is_decided_exactly_at_zero_tolerance() {
    for inside in [true, false] {
        let report = cut(
            inside,
            0.0,
            Transform3::IDENTITY,
            Tolerance::ZERO,
            (expected(inside), 0.0),
        );
        assert!(report.is_exact(), "inside {inside}: {report:?}");
    }
}

#[test]
fn a_hole_tangent_under_a_general_placement_is_read_as_touching() {
    for inside in [true, false] {
        let report = cut(
            inside,
            0.0,
            general(),
            Tolerance::METRE,
            (expected(inside), 1e-12),
        );
        assert!(
            report.contains(ToleranceDecisionKind::PlaneTouchesCylinder),
            "inside {inside}: {report:?}"
        );
        assert!(report.linear() <= EPS);
    }
}

#[test]
fn a_hole_a_fraction_of_the_tolerance_off_touching_moves_the_plane_by_that_much() {
    // Into the top face (the plane cuts two rulings `2 sqrt(2 R gap)`
    // apart, and a chord as long across each cap) or short of it: one
    // contact ruling for every face pair, the plane moved by `gap`.
    let surface = 2.0 * (L * T + L * H + T * H);
    for inside in [true, false] {
        for p in [Transform3::IDENTITY, general()] {
            for fraction in [0.1, 0.5, 0.9, -0.1, -0.5, -0.9] {
                let gap = fraction * EPS;
                let report = cut(
                    inside,
                    gap,
                    p,
                    Tolerance::METRE,
                    (expected(inside), gap.abs() * surface),
                );
                let moved = report
                    .decisions()
                    .iter()
                    .find(|d| d.kind == ToleranceDecisionKind::PlaneTouchesCylinder)
                    .map(|d| d.linear)
                    .unwrap_or_else(|| panic!("inside {inside} {fraction}: {report:?}"));
                assert!(
                    (moved - gap.abs()).abs() <= 1e-9 * EPS + 1e-15,
                    "inside {inside} {fraction}: moved {moved}"
                );
                assert!(report.linear() <= EPS, "{report:?}");
            }
        }
    }
}

#[test]
fn ten_tolerances_off_touching_is_decided_exactly() {
    // Ten tolerances into the top face, the hole cuts a slot across the
    // wall's top `2 sqrt(2 R d)` wide; ten tolerances short, it misses it.
    let d = 10.0 * EPS;
    let segment = R * R * ((R - d) / R).acos() - (R - d) * (2.0 * R * d - d * d).sqrt();
    for (gap, removed) in [(d, segment * T), (-d, 0.0)] {
        let report = cut(
            false,
            -gap,
            Transform3::IDENTITY,
            Tolerance::METRE,
            (L * T * H - removed, 1e-3 * segment * T),
        );
        assert!(
            !report.contains(ToleranceDecisionKind::PlaneTouchesCylinder),
            "{report:?}"
        );
    }
}
