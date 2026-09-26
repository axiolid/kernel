//! Ellipse profile extrusion (ADR 0055).
//!
//! The claim: an ellipse extrudes to a solid carrying a genuine
//! `EllipticalCylinder` wall, exactly -- not a spline fit and not a polygon.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Tolerance, Vec3};
use axiolid_evaluate::surface;
use axiolid_measure::exact_properties;
use axiolid_profile::{EllipseProfile, Profile};
use axiolid_surface::Surface;

fn ellipse(a: f64, b: f64) -> Profile {
    Profile::Ellipse(EllipseProfile {
        semi_axis_x: a,
        semi_axis_y: b,
    })
}

#[test]
fn an_ellipse_extrudes_to_an_elliptical_cylinder_wall() {
    let (a, b) = (3.0, 1.0);
    let solid = extrude_profile_exact(&ellipse(a, b), Vec3::Z, 2.0, Tolerance::METRE)
        .expect("an ellipse extrudes");

    let walls: Vec<(f64, f64)> = solid
        .surfaces()
        .iter()
        .filter_map(|surface| match surface {
            Surface::EllipticalCylinder(value) => Some((value.semi_axis_x, value.semi_axis_y)),
            _ => None,
        })
        .collect();
    assert_eq!(walls.len(), 1, "one elliptical wall, got {}", walls.len());
    assert!(
        (walls[0].0 - a).abs() < 1e-15 && (walls[0].1 - b).abs() < 1e-15,
        "the wall must carry the profile's own semi-axes, got {:?}",
        walls[0]
    );

    // Not a spline fit and not tessellated into planes.
    let splines = solid
        .surfaces()
        .iter()
        .filter(|surface| matches!(surface, Surface::BSpline(_)))
        .count();
    assert_eq!(splines, 0, "the wall must be exact, not a spline fit");
    let planes = solid
        .surfaces()
        .iter()
        .filter(|surface| matches!(surface, Surface::Plane(_)))
        .count();
    assert_eq!(planes, 2, "only the two caps are planar, got {planes}");
}

#[test]
fn the_elliptical_wall_is_measured_exactly() {
    // `exact_properties` integrates the wall over its own parameters (#125):
    // the volume is pi a b h, which a wall mistaken for a circular cylinder
    // (radial normals) or for a plane would not give.
    let (a, b, h) = (3.0, 1.0, 2.0);
    let solid = extrude_profile_exact(&ellipse(a, b), Vec3::Z, h, Tolerance::METRE)
        .expect("an ellipse extrudes");
    let volume = exact_properties(&solid, Tolerance::METRE)
        .expect("an elliptical cylinder is measurable")
        .signed_volume;
    let expected = core::f64::consts::PI * a * b * h;
    assert!(
        (volume - expected).abs() < 1e-11 * expected,
        "expected {expected}, got {volume}"
    );
}

#[test]
fn the_cap_boundary_agrees_with_the_wall() {
    // This is what the geometric audit exists for: the cap pcurve is an
    // `Ellipse2` and the wall is an `EllipticalCylinder`, and both must
    // parameterise u the same way or the shared edge disagrees.
    let solid = extrude_profile_exact(&ellipse(2.5, 0.8), Vec3::Z, 1.5, Tolerance::METRE)
        .expect("an ellipse extrudes");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "cap and wall must agree, found {:?} worst {:?}",
        health.defects(),
        health.worst_error()
    );
}

#[test]
fn degenerate_semi_axes_are_refused() {
    for (a, b) in [(0.0, 1.0), (1.0, -1.0), (f64::NAN, 1.0)] {
        assert!(
            extrude_profile_exact(&ellipse(a, b), Vec3::Z, 1.0, Tolerance::METRE).is_err(),
            "semi-axes ({a}, {b}) must be refused"
        );
    }
}

#[test]
fn the_wall_surface_matches_the_closed_form_everywhere() {
    // The surface itself is checked against its closed form: every sampled point must satisfy (x/a)^2 + (y/b)^2 = 1.
    let (a, b, depth) = (3.0, 1.25, 2.0);
    let solid = extrude_profile_exact(&ellipse(a, b), Vec3::Z, depth, Tolerance::METRE)
        .expect("an ellipse extrudes");
    let wall = solid
        .surfaces()
        .iter()
        .find_map(|s| match s {
            Surface::EllipticalCylinder(c) => Some(*c),
            _ => None,
        })
        .expect("an elliptical wall");

    let mut worst: f64 = 0.0;
    for i in 0..64 {
        let u = core::f64::consts::TAU * f64::from(i) / 64.0;
        for j in 0..5 {
            let v = depth * f64::from(j) / 4.0;
            let p = surface::evaluate(&Surface::EllipticalCylinder(wall), u, v)
                .expect("wall evaluates");
            worst = worst.max(((p.x / a).powi(2) + (p.y / b).powi(2) - 1.0).abs());
            worst = worst.max((p.z - v).abs());
        }
    }
    assert!(worst < 1e-15, "wall is not the exact ellipse: {worst:e}");
}
