//! `Curve2::Sinusoid` (ADR 0071): the exact pcurve of a plane's cut across a
//! cylinder.
//!
//! The defining property is geometric, so it is tested geometrically: lift
//! the wave through the cylinder and every point must lie on the plane the
//! wave was derived from. The coefficients are derived here from the plane
//! by hand, not read back from the kernel.

use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Curve2, Sinusoid2};
use axiolid_evaluate::curve::{derivative2, domain2, evaluate2, invert2, second_derivative2};
use axiolid_evaluate::surface::evaluate;
use axiolid_surface::{Cylinder, EllipticalCylinder, Surface};

const TAU: f64 = std::f64::consts::TAU;

fn frame(origin: Point3) -> Frame3 {
    Frame3 {
        origin,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

#[test]
fn a_lifted_wave_lies_on_the_plane_it_was_cut_by() {
    // Plane z = 2 + 0.3 x - 0.2 y over a cylinder of radius 1.5 centred at
    // (4, -1). On the cylinder x = 4 + 1.5 cos u and y = -1 + 1.5 sin u, so
    // z - base = 2 + 0.3*4 - 0.2*(-1) - base + 0.45 cos u - 0.3 sin u.
    let base = 0.5;
    let wave = Curve2::Sinusoid(Sinusoid2 {
        mean: 2.0 + 1.2 + 0.2 - base,
        cosine: 0.3 * 1.5,
        sine: -0.2 * 1.5,
    });
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(4.0, -1.0, base)),
        radius: 1.5,
    });
    let mut worst: f64 = 0.0;
    for i in 0..=64 {
        let t = TAU * f64::from(i) / 64.0;
        let uv = evaluate2(&wave, t).expect("finite wave");
        let p = evaluate(&cylinder, uv.x, uv.y).expect("cylinder point");
        let plane_z = 2.0 + 0.3 * p.x - 0.2 * p.y;
        worst = worst.max((p.z - plane_z).abs());
    }
    assert!(worst < 1e-14, "lifted wave leaves the plane by {worst:e}");
}

#[test]
fn the_same_holds_on_an_elliptical_cylinder() {
    // On an elliptical cylinder x = a cos u and y = b sin u, so the plane's
    // slopes scale by the two semi-axes separately.
    let (a, b) = (2.0, 0.75);
    let wave = Curve2::Sinusoid(Sinusoid2 {
        mean: 1.0,
        cosine: 0.4 * a,
        sine: 0.25 * b,
    });
    let cylinder = Surface::EllipticalCylinder(EllipticalCylinder {
        frame: frame(Point3::ZERO),
        semi_axis_x: a,
        semi_axis_y: b,
    });
    for i in 0..=32 {
        let t = TAU * f64::from(i) / 32.0;
        let uv = evaluate2(&wave, t).unwrap();
        let p = evaluate(&cylinder, uv.x, uv.y).unwrap();
        let plane_z = 1.0 + 0.4 * p.x + 0.25 * p.y;
        assert!((p.z - plane_z).abs() < 1e-14, "off plane at t = {t}");
    }
}

#[test]
fn the_parameter_is_the_first_coordinate_and_the_domain_one_turn() {
    let wave = Curve2::Sinusoid(Sinusoid2 {
        mean: 3.0,
        cosine: 1.0,
        sine: 2.0,
    });
    let d = domain2(&wave);
    assert_eq!((d.start, d.end), (0.0, TAU));
    // Quadrant values by hand: cos/sin are exact at 0 and within an ulp at
    // pi/2, pi.
    assert_eq!(evaluate2(&wave, 0.0).unwrap(), Point2::new(0.0, 4.0));
    let quarter = evaluate2(&wave, TAU / 4.0).unwrap();
    assert!((quarter.x - TAU / 4.0).abs() < 1e-15 && (quarter.y - 5.0).abs() < 1e-15);
    let half = evaluate2(&wave, TAU / 2.0).unwrap();
    assert!((half.y - 2.0).abs() < 1e-15);
    // Beyond one turn it wraps rather than refusing: a pcurve may state an
    // angle span that crosses the seam.
    let wrapped = evaluate2(&wave, TAU).unwrap();
    assert!((wrapped.y - 4.0).abs() < 1e-14 && (wrapped.x - TAU).abs() < 1e-15);
}

#[test]
fn derivatives_match_central_differences() {
    let wave = Curve2::Sinusoid(Sinusoid2 {
        mean: -1.0,
        cosine: 0.7,
        sine: -1.3,
    });
    let h = 1e-5;
    for i in 0..16 {
        let t = 0.37 + 0.4 * f64::from(i);
        let d1 = derivative2(&wave, t).unwrap();
        let d2 = second_derivative2(&wave, t).unwrap();
        let (m, p) = (
            evaluate2(&wave, t - h).unwrap(),
            evaluate2(&wave, t + h).unwrap(),
        );
        let c = evaluate2(&wave, t).unwrap();
        let fd1 = (p - m) / (2.0 * h);
        // Height only: the first coordinate is linear in t.
        let fd2 = Vec2::new(0.0, (p.y - 2.0 * c.y + m.y) / (h * h));
        assert!((d1 - fd1).length() < 1e-9, "first derivative at {t}");
        assert_eq!(d1.x, 1.0, "the parameter advances one-for-one");
        assert!((d2 - fd2).length() < 1e-4, "second at {t}");
        assert_eq!(d2.x, 0.0);
    }
}

#[test]
fn inversion_returns_the_first_coordinate_and_refuses_off_curve_points() {
    let wave = Curve2::Sinusoid(Sinusoid2 {
        mean: 0.0,
        cosine: 1.0,
        sine: 0.0,
    });
    let t = 1.1;
    let on = evaluate2(&wave, t).unwrap();
    let tol = Tolerance::METRE;
    assert_eq!(invert2(&wave, on, tol).unwrap(), t);
    // Same angle, wrong height: not on the curve.
    let off = Point2::new(t, on.y + 0.5);
    assert!(invert2(&wave, off, tol).is_err());
}

#[test]
fn a_zero_amplitude_wave_is_a_horizontal_line() {
    let wave = Sinusoid2 {
        mean: 2.5,
        cosine: 0.0,
        sine: 0.0,
    };
    assert_eq!(wave.amplitude(), 0.0);
    for i in 0..8 {
        assert_eq!(wave.height(f64::from(i)), 2.5);
    }
    // Amplitude is the peak deviation, not the sum of the terms.
    let tilted = Sinusoid2 {
        mean: 0.0,
        cosine: 3.0,
        sine: 4.0,
    };
    assert_eq!(tilted.amplitude(), 5.0);
    assert!(!Sinusoid2 {
        mean: f64::NAN,
        ..tilted
    }
    .is_finite());
}
