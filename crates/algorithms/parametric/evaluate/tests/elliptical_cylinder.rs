//! Elliptical cylinder surface (ADR 0055).
//!
//! The claim that matters: the normal is derived from the partials, not
//! copied from the circular case. For a circular cylinder the outward normal
//! IS the radial direction; for an ellipse that holds only at the four axis
//! points, so a test that samples only those would pass on wrong code.

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_evaluate::surface;
use axiolid_surface::{EllipticalCylinder, Surface};

fn elliptical(a: f64, b: f64) -> Surface {
    Surface::EllipticalCylinder(EllipticalCylinder {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        semi_axis_x: a,
        semi_axis_y: b,
    })
}

#[test]
fn points_satisfy_the_implicit_ellipse_equation() {
    let (a, b) = (3.0, 1.0);
    let value = elliptical(a, b);
    let mut worst: f64 = 0.0;
    for step in 0..72 {
        let u = core::f64::consts::TAU * f64::from(step) / 72.0;
        for v in [-2.0, 0.0, 1.7] {
            let point = surface::evaluate(&value, u, v).expect("evaluate");
            // (X/a)^2 + (Y/b)^2 = 1 everywhere on the surface.
            let residual = (point.x / a).powi(2) + (point.y / b).powi(2) - 1.0;
            worst = worst.max(residual.abs());
            assert!(
                (point.z - v).abs() < 1e-15,
                "the axis parameter must be the height"
            );
        }
    }
    assert!(worst < 1e-15, "implicit residual {worst:e}");
}

#[test]
fn the_normal_is_not_the_radial_direction() {
    // A 3:1 ellipse. Measured in Python: the true normal and the naive
    // radial direction disagree by up to 53.13 degrees, and agree exactly
    // at u = 0 and u = 90 degrees.
    let (a, b) = (3.0, 1.0);
    let value = elliptical(a, b);

    let disagreement = |u: f64| {
        let normal = surface::normal(&value, u, 0.4).expect("normal");
        let point = surface::evaluate(&value, u, 0.4).expect("point");
        let radial = Vec3::new(point.x, point.y, 0.0).normalize();
        normal
            .dot(radial)
            .abs()
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    };

    // At the axis points the two coincide -- which is why sampling only
    // there would not catch a radial-normal implementation.
    assert!(disagreement(0.0) < 1e-9, "u=0 must agree");
    assert!(
        disagreement(core::f64::consts::FRAC_PI_2) < 1e-9,
        "u=90deg must agree"
    );

    // Off-axis they genuinely differ, by the measured amount.
    let at_45 = disagreement(core::f64::consts::FRAC_PI_4);
    assert!(
        (at_45 - 53.130_102_354_156).abs() < 1e-9,
        "u=45deg disagreement should be 53.1301 degrees, got {at_45}"
    );
}

#[test]
fn the_normal_is_perpendicular_to_both_partials() {
    // The defining property, checked independently of any formula.
    let value = elliptical(2.5, 0.75);
    for step in 0..48 {
        let u = core::f64::consts::TAU * f64::from(step) / 48.0;
        let normal = surface::normal(&value, u, 0.2).expect("normal");
        let (du, dv) = surface::partials(&value, u, 0.2).expect("partials");
        assert!(
            normal.dot(du).abs() < 1e-12 && normal.dot(dv).abs() < 1e-12,
            "normal must be perpendicular to both partials at u={u}"
        );
        assert!(
            (normal.length() - 1.0).abs() < 1e-12,
            "the normal must be a unit vector"
        );
    }
}

#[test]
fn equal_semi_axes_reduce_to_the_circular_cylinder() {
    // A degenerate-case check: with a = b the surface must agree with
    // `Cylinder` of that radius, point for point.
    let radius = 1.75;
    let ellipse = elliptical(radius, radius);
    let circle = Surface::Cylinder(axiolid_surface::Cylinder {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    });
    for step in 0..36 {
        let u = core::f64::consts::TAU * f64::from(step) / 36.0;
        let a = surface::evaluate(&ellipse, u, 0.9).expect("ellipse");
        let b = surface::evaluate(&circle, u, 0.9).expect("circle");
        assert!((a - b).length() < 1e-15, "at u={u}: {a:?} vs {b:?}");
        let na = surface::normal(&ellipse, u, 0.9).expect("ellipse normal");
        let nb = surface::normal(&circle, u, 0.9).expect("circle normal");
        assert!((na - nb).length() < 1e-12, "normals must agree too");
    }
}

#[test]
fn a_non_positive_semi_axis_is_refused() {
    for (a, b) in [(0.0, 1.0), (1.0, 0.0), (-2.0, 1.0)] {
        assert!(
            surface::evaluate(&elliptical(a, b), 0.3, 0.0).is_err(),
            "semi-axes ({a}, {b}) must be refused"
        );
    }
}
