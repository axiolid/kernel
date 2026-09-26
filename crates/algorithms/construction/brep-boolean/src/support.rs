//! Whether two faces lie on one surface, whatever their frames.

use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_surface::Surface;

/// Whether two supports are the same surface point set: a column's
/// half-walls lie on one cylinder parameterised from opposite sides, and two
/// operands' coplanar faces on one plane with unrelated frames.
pub(crate) fn same_support(a: &Surface, b: &Surface, tolerance: Tolerance) -> bool {
    let eps = tolerance.linear().max(1e-9);
    let parallel = |x: Vec3, y: Vec3| x.normalize().cross(y.normalize()).length() <= 1e-9;
    let same_way = |x: Vec3, y: Vec3| parallel(x, y) && x.dot(y) > 0.0;
    let on_axis = |o1: Point3, o2: Point3, z: Vec3| (o2 - o1).cross(z.normalize()).length() <= eps;
    let close = |x: Scalar, y: Scalar| (x - y).abs() <= eps;
    match (a, b) {
        (Surface::Plane(p), Surface::Plane(q)) => {
            parallel(p.frame.z, q.frame.z)
                && (q.frame.origin - p.frame.origin)
                    .dot(p.frame.z.normalize())
                    .abs()
                    <= eps
        }
        (Surface::Cylinder(p), Surface::Cylinder(q)) => {
            close(p.radius, q.radius)
                && parallel(p.frame.z, q.frame.z)
                && on_axis(p.frame.origin, q.frame.origin, p.frame.z)
        }
        (Surface::EllipticalCylinder(p), Surface::EllipticalCylinder(q)) => {
            if !(parallel(p.frame.z, q.frame.z)
                && on_axis(p.frame.origin, q.frame.origin, p.frame.z))
            {
                return false;
            }
            // The same ellipse, with the axes possibly named the other way.
            (close(p.semi_axis_x, q.semi_axis_x)
                && close(p.semi_axis_y, q.semi_axis_y)
                && parallel(p.frame.x, q.frame.x))
                || (close(p.semi_axis_x, q.semi_axis_y)
                    && close(p.semi_axis_y, q.semi_axis_x)
                    && parallel(p.frame.x, q.frame.y))
        }
        (Surface::Cone(p), Surface::Cone(q)) => {
            // One nappe: the same apex, opening the same way at the same
            // angle.
            let (tp, tq) = (p.semi_angle.tan(), q.semi_angle.tan());
            if tp == 0.0 || tq == 0.0 {
                return false;
            }
            let apex = |c: &axiolid_surface::Cone, t: Scalar| {
                c.frame.origin - c.frame.z.normalize() * (c.radius / t)
            };
            (apex(p, tp) - apex(q, tq)).length() <= eps
                && same_way(p.frame.z * tp.signum(), q.frame.z * tq.signum())
                && (tp.abs() - tq.abs()).abs() <= 1e-12 * (1.0 + tp.abs())
        }
        (Surface::Sphere(p), Surface::Sphere(q)) => {
            close(p.radius, q.radius) && (p.frame.origin - q.frame.origin).length() <= eps
        }
        (Surface::Torus(p), Surface::Torus(q)) => {
            close(p.major_radius, q.major_radius)
                && close(p.minor_radius, q.minor_radius)
                && (p.frame.origin - q.frame.origin).length() <= eps
                && parallel(p.frame.z, q.frame.z)
        }
        _ => a == b,
    }
}
