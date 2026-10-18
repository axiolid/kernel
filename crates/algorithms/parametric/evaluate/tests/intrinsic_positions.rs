//! Certified positions of intrinsic curves (#317).
//!
//! Every expectation is independent of the midpoint rule under test: an
//! arc of constant curvature in closed form, a clothoid from the Fresnel
//! integrals' power series. Each returned point must lie within its
//! reported error of the exact one, at loose and tight targets, and the
//! error stay within the target.

use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{CurvatureLaw, Intrinsic2, Intrinsic3};
use axiolid_evaluate::bound::IntrinsicPositions;

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    }
}

fn place(frame: &Frame2, local: Vec2) -> Point2 {
    frame.origin + frame.x * local.x + frame.y * local.y
}

/// `int_0^s (cos(a u^2), sin(a u^2)) du` by its power series.
fn fresnel(a: Scalar, s: Scalar) -> Vec2 {
    let (mut x, mut y) = (0.0, 0.0);
    let mut factorial = 1.0;
    for n in 0..40 {
        if n > 0 {
            factorial *= n as Scalar;
        }
        // The series of exp(i a u^2) integrated term by term: (i a)^n
        // s^(2n+1) / (n! (2n+1)), real and imaginary parts alternating.
        let value = a.powi(n) * s.powi(2 * n + 1) / (factorial * (2 * n + 1) as Scalar);
        match n % 4 {
            0 => x += value,
            1 => y += value,
            2 => x -= value,
            _ => y -= value,
        }
    }
    Vec2::new(x, y)
}

#[test]
fn an_arc_is_within_its_certified_error() {
    let radius = 12.0;
    let arc = Intrinsic2::new(
        frame(),
        CurvatureLaw::Constant {
            curvature: 1.0 / radius,
        },
        40.0,
    );
    for target in [1e-2, 1e-5, 1e-9] {
        let mut positions = IntrinsicPositions::new2(&arc, target).unwrap();
        for i in 0..=97 {
            let s = 40.0 * i as Scalar / 97.0;
            let (p, error) = positions.point2(s).unwrap();
            let a = s / radius;
            let exact = place(
                &frame(),
                Vec2::new(radius * a.sin(), radius * (1.0 - a.cos())),
            );
            assert!((p - exact).length() <= error, "{target}: {s}");
            assert!(error <= target, "{error} > {target}");
            // Its left offset is the concentric arc of radius r - d.
            let (q, error) = positions.offset_point2(s, 2.0).unwrap();
            let centre = place(&frame(), Vec2::new(0.0, radius));
            let exact = centre + (exact - centre) * ((radius - 2.0) / radius);
            assert!((q - exact).length() <= error, "offset {target}: {s}");
        }
    }
}

#[test]
fn a_clothoid_is_within_its_certified_error_read_in_any_order() {
    // k = c s, theta = c s^2 / 2.
    let c = 1.0 / 1800.0;
    let spiral = Intrinsic2::new(frame(), CurvatureLaw::clothoid(0.0, 60.0 * c, 60.0), 60.0);
    for target in [1e-2, 1e-6] {
        let mut positions = IntrinsicPositions::new2(&spiral, target).unwrap();
        // Forward, then backward (each read restarts), then forward again.
        let forward: Vec<Scalar> = (0..=61).map(|i| 60.0 * i as Scalar / 61.0).collect();
        let backward: Vec<Scalar> = forward.iter().rev().copied().collect();
        for s in forward.iter().chain(&backward).chain(&forward) {
            let (p, error) = positions.point2(*s).unwrap();
            let exact = place(&frame(), fresnel(0.5 * c, *s));
            let off = (p - exact).length();
            assert!(off <= error, "{target}: {s}: {off} > {error}");
            assert!(error <= target);
        }
    }
}

#[test]
fn a_planar_space_curve_is_its_plane_curve_and_torsion_is_refused() {
    let c = 1.0 / 1800.0;
    let z = Vec3::new(1.0, -2.0, 2.0).normalize();
    let x = Vec3::new(2.0, 1.0, 0.0).normalize();
    let start = Frame3 {
        origin: Point3::new(5.0, -3.0, 1.5),
        x,
        y: z.cross(x),
        z,
    };
    let law = CurvatureLaw::clothoid(0.0, 60.0 * c, 60.0);
    let space = Intrinsic3::new(
        start,
        law.clone(),
        CurvatureLaw::Constant { curvature: 0.0 },
        60.0,
    );
    let mut positions = IntrinsicPositions::new3(&space, 1e-4).unwrap();
    for i in 0..=40 {
        let s = 1.5 * i as Scalar;
        let (p, error) = positions.point3(s).unwrap();
        let local = fresnel(0.5 * c, s);
        let exact = start.origin + start.x * local.x + start.y * local.y;
        assert!((p - exact).length() <= error);
        assert!(error <= 1e-4);
    }
    let twisted = Intrinsic3::new(start, law, CurvatureLaw::Constant { curvature: 0.01 }, 60.0);
    assert!(IntrinsicPositions::new3(&twisted, 1e-4).is_none());
}
