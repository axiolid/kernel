//! Exact predicates on the operands' own numbers (#236).
//!
//! A support decision ([`crate::report::support`]) asks one of these before
//! it reads anything within tolerance. Each evaluates a polynomial in the
//! surfaces' `f64` frame entries exactly, in dyadic arithmetic
//! (`axiolid-exact`), so "exactly coincident", "exactly parallel" and
//! "exactly perpendicular" mean what they say for the numbers given: an
//! opening placed with axis matrices whose entries are `0` and `+-1` meets
//! its wall in planes that are exactly coplanar or exactly perpendicular,
//! while the same opening under a general rotation does not.
//!
//! They are asked only for pairs whose `f64` measures are already within
//! tolerance or rounding of zero, so the big-integer cost is paid rarely.
//! A non-finite entry is never exact.

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_surface::{Cylinder, Plane, Surface};

/// An exact vector.
#[derive(Clone)]
struct V([Dyadic; 3]);

impl V {
    fn of(v: Vec3) -> Option<Self> {
        Some(Self([
            Dyadic::try_from_f64(v.x)?,
            Dyadic::try_from_f64(v.y)?,
            Dyadic::try_from_f64(v.z)?,
        ]))
    }

    fn point(p: Point3) -> Option<Self> {
        Self::of(p)
    }

    fn sub(&self, o: &Self) -> Self {
        Self([
            self.0[0].sub(&o.0[0]),
            self.0[1].sub(&o.0[1]),
            self.0[2].sub(&o.0[2]),
        ])
    }

    fn dot(&self, o: &Self) -> Dyadic {
        self.0[0]
            .mul(&o.0[0])
            .add(&self.0[1].mul(&o.0[1]))
            .add(&self.0[2].mul(&o.0[2]))
    }

    fn cross(&self, o: &Self) -> Self {
        let (a, b) = (&self.0, &o.0);
        Self([
            a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
            a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
            a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
        ])
    }

    fn is_zero(&self) -> bool {
        self.0.iter().all(zero)
    }
}

fn zero(x: &Dyadic) -> bool {
    *x == Dyadic::zero()
}

fn one(x: &Dyadic) -> bool {
    *x == Dyadic::from_f64(1.0)
}

/// A frame whose axes are exactly orthonormal: a cylinder on it is exactly
/// circular, its axis exactly `z`.
fn orthonormal(f: &Frame3) -> bool {
    let (Some(x), Some(y), Some(z)) = (V::of(f.x), V::of(f.y), V::of(f.z)) else {
        return false;
    };
    one(&x.dot(&x))
        && one(&y.dot(&y))
        && one(&z.dot(&z))
        && zero(&x.dot(&y))
        && zero(&y.dot(&z))
        && zero(&z.dot(&x))
}

/// A plane's exact normal, `x cross y` (the plane's points are
/// `origin + u x + v y`), and its origin.
fn plane(p: &Plane) -> Option<(V, V)> {
    let n = V::of(p.frame.x)?.cross(&V::of(p.frame.y)?);
    (!n.is_zero()).then_some(())?;
    Some((n, V::point(p.frame.origin)?))
}

/// Whether two supports are exactly one surface point set. Planes and
/// circular cylinders are decided by exact predicates; any other pair is
/// exactly one only when the two are the same numbers.
pub(crate) fn same_support(a: &Surface, b: &Surface) -> bool {
    match (a, b) {
        (Surface::Plane(p), Surface::Plane(q)) => (|| {
            let (np, op) = plane(p)?;
            let (nq, oq) = plane(q)?;
            Some(np.cross(&nq).is_zero() && zero(&oq.sub(&op).dot(&np)))
        })()
        .unwrap_or(false),
        (Surface::Cylinder(p), Surface::Cylinder(q)) => (|| {
            if p.radius != q.radius || !orthonormal(&p.frame) || !orthonormal(&q.frame) {
                return Some(false);
            }
            let (zp, zq) = (V::of(p.frame.z)?, V::of(q.frame.z)?);
            let offset = V::point(q.frame.origin)?.sub(&V::point(p.frame.origin)?);
            Some(zp.cross(&zq).is_zero() && offset.cross(&zp).is_zero())
        })()
        .unwrap_or(false),
        _ => a == b,
    }
}

/// A plane and an exactly circular cylinder, with the plane's exact normal
/// and the cylinder's exact axis and origin offset from the plane's.
fn plane_cylinder(p: &Plane, c: &Cylinder) -> Option<(V, V, V)> {
    if !orthonormal(&c.frame) {
        return None;
    }
    let (n, op) = plane(p)?;
    let axis = V::of(c.frame.z)?;
    let offset = V::point(c.frame.origin)?.sub(&op);
    Some((n, axis, offset))
}

/// Whether the plane is exactly perpendicular to the cylinder's axis (their
/// section is exactly a circle).
pub(crate) fn plane_perpendicular(p: &Plane, c: &Cylinder) -> bool {
    plane_cylinder(p, c).is_some_and(|(n, axis, _)| n.cross(&axis).is_zero())
}

/// Whether the plane is exactly parallel to the cylinder's axis (their
/// section is exactly rulings, or nothing).
pub(crate) fn plane_parallel(p: &Plane, c: &Cylinder) -> bool {
    plane_cylinder(p, c).is_some_and(|(n, axis, _)| zero(&n.dot(&axis)))
}

/// Whether the plane is exactly parallel to the cylinder's axis at exactly
/// the radius from it: `((o_c - o_p) . n)^2 = r^2 (n . n)`.
pub(crate) fn plane_touches(p: &Plane, c: &Cylinder) -> bool {
    let Some(r) = Dyadic::try_from_f64(c.radius) else {
        return false;
    };
    plane_cylinder(p, c).is_some_and(|(n, axis, offset)| {
        let d = offset.dot(&n);
        zero(&n.dot(&axis)) && d.mul(&d) == r.mul(&r).mul(&n.dot(&n))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(origin: Point3, x: Vec3, y: Vec3, z: Vec3) -> Frame3 {
        Frame3 { origin, x, y, z }
    }

    fn plane_at(origin: Point3, x: Vec3, y: Vec3) -> Plane {
        Plane {
            frame: frame(origin, x, y, x.cross(y)),
        }
    }

    #[test]
    fn axis_planes_are_exactly_coplanar_and_rotated_ones_are_not() {
        let a = Surface::Plane(plane_at(Point3::new(0.0, 0.125, 0.0), Vec3::X, Vec3::Z));
        // Same plane, other frame (origin elsewhere in it, axes swapped).
        let b = Surface::Plane(plane_at(Point3::new(3.0, 0.125, -2.0), Vec3::Z, -Vec3::X));
        assert!(same_support(&a, &b));
        // One ulp off.
        let c = Surface::Plane(plane_at(
            Point3::new(0.0, 0.125 + f64::EPSILON / 8.0, 0.0),
            Vec3::X,
            Vec3::Z,
        ));
        assert!(!same_support(&a, &c));
        // Tilted by the residue a constructor's `cos(pi / 2)` leaves: not
        // exactly the same plane.
        let d = Surface::Plane(plane_at(
            Point3::new(0.0, 0.125, 0.0),
            Vec3::new(1.0, (std::f64::consts::FRAC_PI_2).cos(), 0.0),
            Vec3::Z,
        ));
        assert!(!same_support(&a, &d));
    }

    #[test]
    fn exact_plane_cylinder_relations() {
        let c = Cylinder {
            frame: frame(Point3::new(1.0, 0.0, 0.0), Vec3::X, Vec3::Y, Vec3::Z),
            radius: 0.5,
        };
        let floor = plane_at(Point3::new(0.0, 0.0, 2.0), Vec3::X, Vec3::Y);
        assert!(plane_perpendicular(&floor, &c));
        assert!(!plane_parallel(&floor, &c));
        let jamb = plane_at(Point3::new(1.5, 0.0, 0.0), Vec3::Y, Vec3::Z);
        assert!(plane_parallel(&jamb, &c));
        assert!(plane_touches(&jamb, &c));
        let through = plane_at(Point3::new(1.25, 0.0, 0.0), Vec3::Y, Vec3::Z);
        assert!(plane_parallel(&through, &c) && !plane_touches(&through, &c));
    }
}
