//! Whether two faces lie on one surface, whatever their frames.

use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_surface::Surface;
use core::f64::consts::TAU;

use crate::predicate;
use crate::report::{self, ToleranceDecisionKind};

/// Whether two supports are the same surface point set: a column's
/// half-walls lie on one cylinder parameterised from opposite sides, and two
/// operands' coplanar faces on one plane with unrelated frames.
///
/// Within tolerance in `f64`, the exact predicate decides whether the two
/// are exactly one surface; a reading it does not confirm is reported
/// ([`ToleranceDecisionKind::CoincidentSupports`], #236).
pub(crate) fn same_support(a: &Surface, b: &Surface, tolerance: Tolerance) -> bool {
    let Some((linear, angular)) = support_offsets(a, b) else {
        return false;
    };
    report::support(
        ToleranceDecisionKind::CoincidentSupports,
        linear,
        angular,
        tolerance,
        || predicate::same_support(a, b),
    )
    .holds()
}

/// How far apart two supports of one family are, in `f64`: the largest
/// offset (linear) and turn (angular) a reading of them as one would need,
/// or `None` when no tolerance makes them one (other families, opposite
/// cones).
fn support_offsets(a: &Surface, b: &Surface) -> Option<(Scalar, Scalar)> {
    let turn = |x: Vec3, y: Vec3| x.normalize().cross(y.normalize()).length();
    let off_axis = |o1: Point3, o2: Point3, z: Vec3| (o2 - o1).cross(z.normalize()).length();
    let gap = |x: Scalar, y: Scalar| (x - y).abs();
    let finite = |(l, a): (Scalar, Scalar)| (l.is_finite() && a.is_finite()).then_some((l, a));
    match (a, b) {
        (Surface::Plane(p), Surface::Plane(q)) => finite((
            (q.frame.origin - p.frame.origin)
                .dot(p.frame.z.normalize())
                .abs(),
            turn(p.frame.z, q.frame.z),
        )),
        (Surface::Cylinder(p), Surface::Cylinder(q)) => finite((
            gap(p.radius, q.radius).max(off_axis(p.frame.origin, q.frame.origin, p.frame.z)),
            turn(p.frame.z, q.frame.z),
        )),
        (Surface::EllipticalCylinder(p), Surface::EllipticalCylinder(q)) => {
            let axis = off_axis(p.frame.origin, q.frame.origin, p.frame.z);
            let tilt = turn(p.frame.z, q.frame.z);
            // The same ellipse, with the axes possibly named the other way.
            let named = |sx: Scalar, sy: Scalar, x: Vec3| {
                (
                    axis.max(gap(p.semi_axis_x, sx)).max(gap(p.semi_axis_y, sy)),
                    tilt.max(turn(p.frame.x, x)),
                )
            };
            let (l1, a1) = named(q.semi_axis_x, q.semi_axis_y, q.frame.x);
            let (l2, a2) = named(q.semi_axis_y, q.semi_axis_x, q.frame.y);
            finite(if l1.max(a1) <= l2.max(a2) {
                (l1, a1)
            } else {
                (l2, a2)
            })
        }
        (Surface::Cone(p), Surface::Cone(q)) => {
            // One nappe: the same apex, opening the same way at the same
            // angle.
            let (tp, tq) = (p.semi_angle.tan(), q.semi_angle.tan());
            if tp == 0.0 || tq == 0.0 {
                return None;
            }
            let apex = |c: &axiolid_surface::Cone, t: Scalar| {
                c.frame.origin - c.frame.z.normalize() * (c.radius / t)
            };
            let (zp, zq) = (p.frame.z * tp.signum(), q.frame.z * tq.signum());
            if zp.dot(zq) <= 0.0 || (tp.abs() - tq.abs()).abs() > 1e-12 * (1.0 + tp.abs()) {
                return None;
            }
            finite(((apex(p, tp) - apex(q, tq)).length(), turn(zp, zq)))
        }
        (Surface::Sphere(p), Surface::Sphere(q)) => finite((
            gap(p.radius, q.radius).max((p.frame.origin - q.frame.origin).length()),
            0.0,
        )),
        (Surface::Torus(p), Surface::Torus(q)) => finite((
            gap(p.major_radius, q.major_radius)
                .max(gap(p.minor_radius, q.minor_radius))
                .max((p.frame.origin - q.frame.origin).length()),
            turn(p.frame.z, q.frame.z),
        )),
        _ => (a == b).then_some((0.0, 0.0)),
    }
}

/// Which parameters of a surface are angles, periodic with `2 pi`.
pub(crate) fn periods(surface: &Surface) -> (bool, bool) {
    match surface {
        Surface::Cylinder(_)
        | Surface::EllipticalCylinder(_)
        | Surface::Cone(_)
        | Surface::Sphere(_) => (true, false),
        Surface::Torus(_) => (true, true),
        _ => (false, false),
    }
}

/// The window a face's sections are traced in: its parameter box, a little
/// larger so edges on the box are inside it; a whole turn exactly where the
/// face goes all the way round, and a sphere's latitudes no further than
/// its poles.
pub(crate) fn window(surface: &Surface, lo: Point2, hi: Point2) -> (Point2, Point2) {
    let (pu, pv) = periods(surface);
    let pad = 1e-6 * (1.0 + (hi - lo).length());
    let (mut a, mut b) = (lo - Point2::splat(pad), hi + Point2::splat(pad));
    // A whole turn starts off the face's own seam angle by an irrational
    // fraction of a radian, so the window's edge is never where a symmetric
    // section crosses or turns; pieces are moved into the face's range
    // afterwards.
    let offset = 0.123_456_789_012_345_67;
    if pu && hi.x - lo.x >= TAU - 1e-9 {
        (a.x, b.x) = (lo.x + offset, lo.x + offset + TAU);
    }
    if pv && hi.y - lo.y >= TAU - 1e-9 {
        (a.y, b.y) = (lo.y + offset, lo.y + offset + TAU);
    }
    if matches!(surface, Surface::Sphere(_)) {
        a.y = a.y.max(-core::f64::consts::FRAC_PI_2);
        b.y = b.y.min(core::f64::consts::FRAC_PI_2);
    }
    (a, b)
}

/// The surface with its frame's rounding residue cleared: a component of
/// the origin or an axis that is only a few ulps of the frame's own size
/// (the `6e-17` a constructor's `cos(pi/2)` leaves) becomes zero. Two
/// pipes built at right angles then meet as the Steinmetz pair they were
/// meant to be, not as skew cylinders 1e-16 apart whose section is
/// singular in all but name. The change is far below any tolerance.
pub(crate) fn cleaned(surface: &Surface) -> Surface {
    let clean_vec = |v: Vec3, scale: Scalar| {
        let floor = 8.0 * Scalar::EPSILON * scale;
        let c = |x: Scalar| if x.abs() <= floor { 0.0 } else { x };
        Vec3::new(c(v.x), c(v.y), c(v.z))
    };
    let frame = |f: axiolid_core::Frame3, size: Scalar| {
        let scale = f.origin.abs().max_element().max(size);
        axiolid_core::Frame3 {
            origin: clean_vec(f.origin, scale),
            x: clean_vec(f.x, f.x.abs().max_element()),
            y: clean_vec(f.y, f.y.abs().max_element()),
            z: clean_vec(f.z, f.z.abs().max_element()),
        }
    };
    match surface {
        Surface::Plane(p) => Surface::Plane(axiolid_surface::Plane {
            frame: frame(p.frame, 1.0),
        }),
        Surface::Cylinder(c) => Surface::Cylinder(axiolid_surface::Cylinder {
            frame: frame(c.frame, c.radius),
            ..*c
        }),
        Surface::EllipticalCylinder(c) => {
            Surface::EllipticalCylinder(axiolid_surface::EllipticalCylinder {
                frame: frame(c.frame, c.semi_axis_x.max(c.semi_axis_y)),
                ..*c
            })
        }
        Surface::Cone(c) => Surface::Cone(axiolid_surface::Cone {
            frame: frame(c.frame, c.radius.abs()),
            ..*c
        }),
        Surface::Sphere(s) => Surface::Sphere(axiolid_surface::Sphere {
            frame: frame(s.frame, s.radius),
            ..*s
        }),
        Surface::Torus(t) => Surface::Torus(axiolid_surface::Torus {
            frame: frame(t.frame, t.major_radius + t.minor_radius),
            ..*t
        }),
        other => other.clone(),
    }
}
