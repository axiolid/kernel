//! Composite profiles with curved members, member holes and disjoint
//! members (#111).
//!
//! Members are unioned exactly over one arc arrangement, so arcs stay arcs
//! and a member's own opening stays open unless another member fills it.
//! Areas are closed forms from the inputs; every solid audits clean and is
//! measured signed with `exact_properties`.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_core::{Point3, Tolerance, Transform2, Vec2, Vec3};
use axiolid_measure::exact_properties;
use axiolid_profile::{CircleProfile, EllipseProfile, Profile, RectangleProfile};
use axiolid_surface::Surface;

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

fn rect(x: f64, y: f64, thickness: Option<f64>) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness,
        outer_radius: None,
        inner_radius: None,
    })
}

fn circle(radius: f64) -> Profile {
    Profile::Circle(CircleProfile {
        radius,
        thickness: None,
    })
}

fn at(basis: Profile, x: f64, y: f64) -> Profile {
    Profile::Derived {
        basis: Box::new(basis),
        transform: Transform2::from_translation(Vec2::new(x, y)),
    }
}

fn checked(solid: ExactBRep) -> ExactBRep {
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    solid
}

fn volume(solid: &ExactBRep) -> f64 {
    exact_properties(solid, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

fn close(got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-10 * want.abs(),
        "expected {want}, got {got}"
    );
}

fn extrude(profile: &Profile, depth: f64) -> ExactBRep {
    checked(extrude_profile_exact(profile, Vec3::Z, depth, Tolerance::METRE).expect("extrudes"))
}

#[test]
fn a_bar_with_a_round_end_keeps_its_arc() {
    // A 4 x 2 bar whose right end is capped by a disc of radius 1 centred
    // on the end: the disc's outer half adds pi/2 to the bar's 8.
    let profile = Profile::Composite(vec![rect(4.0, 2.0, None), at(circle(1.0), 2.0, 0.0)]);
    let solid = extrude(&profile, 3.0);
    // The round end is the disc's two outer quarters: two cylinder walls,
    // both on the disc itself.
    let cylinders: Vec<_> = solid
        .surfaces()
        .iter()
        .filter_map(|surface| match surface {
            Surface::Cylinder(cylinder) => Some(*cylinder),
            _ => None,
        })
        .collect();
    assert_eq!(cylinders.len(), 2, "two quarter walls: {cylinders:?}");
    for cylinder in &cylinders {
        assert!((cylinder.radius - 1.0).abs() < 1e-12, "{cylinder:?}");
        let axis = cylinder.frame.origin;
        assert!(
            (axis.x - 2.0).abs() < 1e-12 && axis.y.abs() < 1e-12,
            "centred on the disc: {cylinder:?}"
        );
    }
    close(volume(&solid), (8.0 + PI / 2.0) * 3.0);
}

#[test]
fn a_member_opening_stays_open_unless_another_member_fills_it() {
    // A 6 x 6 hollow square (wall 1, opening 4 x 4) and a 1 x 6 bar across
    // its middle: the bar splits the opening into two 1.5 x 4 openings.
    let profile = Profile::Composite(vec![rect(6.0, 6.0, Some(1.0)), rect(1.0, 6.0, None)]);
    let solid = extrude(&profile, 2.0);
    close(volume(&solid), (36.0 - 2.0 * 1.5 * 4.0) * 2.0);
}

#[test]
fn a_round_member_inside_another_members_opening_is_a_separate_solid() {
    // A disc floating in the opening of a hollow square touches nothing:
    // two solids.
    let profile = Profile::Composite(vec![rect(6.0, 6.0, Some(1.0)), circle(1.0)]);
    let solid = extrude(&profile, 1.0);
    assert_eq!(solid.topology().solids().len(), 2);
    close(volume(&solid), 36.0 - 16.0 + PI);
}

#[test]
fn a_disjoint_composite_revolves_into_separate_solids() {
    // Two squares 4 and 8 from the axis, one with a round hole: Pappus per
    // piece.
    let hollow_disc = Profile::Circle(CircleProfile {
        radius: 0.75,
        thickness: Some(0.25),
    });
    let profile = Profile::Composite(vec![
        at(rect(1.0, 1.0, None), 4.0, 0.0),
        at(hollow_disc, 8.0, 0.0),
    ]);
    let solid = checked(
        revolve_profile_exact(&profile, Point3::ZERO, Vec3::Y, TAU, Tolerance::METRE)
            .expect("revolves"),
    );
    let solids = solid.topology().solids();
    assert_eq!(solids.len(), 2);
    assert_eq!(solids.iter().map(|s| s.voids.len()).sum::<usize>(), 1);
    let ring = PI * (0.75 * 0.75 - 0.5 * 0.5);
    close(volume(&solid), TAU * (4.0 * 1.0 + 8.0 * ring));
}

#[test]
fn a_member_that_is_not_a_contour_is_refused_by_name() {
    let profile = Profile::Composite(vec![
        rect(2.0, 2.0, None),
        Profile::Ellipse(EllipseProfile {
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        }),
    ]);
    let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("an ellipse does not lower to a contour");
    assert!(format!("{error:?}").contains("contour"), "got {error:?}");
}
