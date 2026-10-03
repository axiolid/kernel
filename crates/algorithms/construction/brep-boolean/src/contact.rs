//! A plane read as touching a cylinder, sewn consistently (#243).
//!
//! A plane parallel to a cylinder's axis at the radius from it touches the
//! cylinder along a ruling. Placed separately, the two miss that by
//! rounding or by a fraction of the tolerance: the plane then cuts the
//! cylinder in two rulings `2 sqrt(2 r d)` apart (`d` how far it reaches
//! in), or misses it. `crate::section` reads such a pair as touching within
//! the tolerance (`PlaneTouchesCylinder`): the plane moved by `d` to touch
//! the cylinder along one ruling, the contact ruling.
//!
//! That reading has to hold everywhere the two surfaces meet, not only for
//! the pair of faces that carry them. Every other face pair that crosses
//! one of them near the ruling meets it there too: a cap disk of the
//! cylinder is cut by the plane in a chord `2 sqrt(2 r d)` long, a face
//! across the axis cuts the cylinder in a circle that crosses the plane
//! twice, and an edge in the plane passes the cylinder twice. Read pair by
//! pair, the exact roots of those crossings stay two points `2 sqrt(2 r d)`
//! apart, far more than the tolerance, while the plane/cylinder pair says
//! there is one ruling: the split faces then disagree and do not sew.
//!
//! So the reading is taken once per pair of surfaces, over the largest
//! window any of their face pairs needs, and kept here. Wherever a curve on
//! one of the two surfaces is cut by the other, the cut is where the curve
//! meets the contact ruling, which is what the moved plane gives:
//!
//! - a curve on the cylinder (a circle or ellipse cut from it by a third
//!   surface) meets the moved plane only on the ruling, where the ruling
//!   pierces the curve's own plane; that point lies on the cylinder
//!   exactly;
//! - a line in the plane meets the cylinder, for the moved plane, only on
//!   the ruling: where the line passes closest to it, at most `d` away.
//!
//! One point replaces each pair of roots, on the same ruling for every face
//! pair, so the perturbation is one and the same throughout: the result is
//! the exact boolean of operands whose plane moved by at most the recorded
//! `PlaneTouchesCylinder` distance `d`. A cut on a curve of the cylinder
//! lies on the ruling exactly; a cut on a line of the given plane lies at
//! most `d` from it (recorded as `TangentCrossing`), and the two are one
//! vertex within `d` (recorded as an incident or merged point), as the
//! moved plane makes them. Where the given and the moved plane disagree on
//! a side -- points of the cylinder beyond the given plane, points of the
//! plane inside the cylinder, both within `d` of the ruling -- no point
//! decides whether a piece lies in a face or a region in a solid
//! ([`Contacts::disputed`]).
//!
//! An exactly tangent pair (the exact predicate holds) needs none of this:
//! its crossings are exact double roots, decided without the tolerance.
//! A crossing this module cannot place on the ruling (a curve of the plane
//! that is no line, a traced curve on the cylinder) is refused by name
//! ([`BooleanError::UnsupportedContact`]): its exact roots would disagree
//! with the reading, and the result would not sew.

use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_curve::{Curve3, Line3};
use axiolid_evaluate::{curve::locate3, evaluate3};
use axiolid_surface::Surface;

use crate::predicate;
use crate::report::{self, ToleranceDecisionKind};
use crate::BooleanError;

/// One plane read as touching one cylinder within tolerance.
#[derive(Debug, Clone)]
pub(crate) struct Contact {
    pub(crate) plane: Surface,
    pub(crate) cylinder: Surface,
    /// The ruling of the cylinder the moved plane touches.
    pub(crate) ruling: Line3,
}

/// Every contact read within tolerance in one boolean.
#[derive(Debug, Clone, Default)]
pub(crate) struct Contacts {
    contacts: Vec<Contact>,
}

/// Whether two supports are one surface exactly: the same numbers, or (two
/// planes or circular cylinders) exactly one point set. Unrecorded: it
/// only names which surface a contact was read for.
fn is(a: &Surface, b: &Surface) -> bool {
    a == b || predicate::same_support(a, b)
}

impl Contacts {
    pub(crate) fn push(&mut self, contact: Contact) {
        self.contacts.push(contact);
    }

    /// The contact whose plane is `plane` and whose cylinder is `cylinder`.
    fn find(&self, plane: &Surface, cylinder: &Surface) -> Option<&Contact> {
        if !matches!(plane, Surface::Plane(_)) || !matches!(cylinder, Surface::Cylinder(_)) {
            return None;
        }
        self.contacts
            .iter()
            .find(|c| is(&c.plane, plane) && is(&c.cylinder, cylinder))
    }

    /// Where `curve`, lying on both `on` surfaces, is cut by `cutter`, when
    /// one of them and the cutter are a contact: parameters on `curve`
    /// where it meets the contact ruling (see the module docs). `None` when
    /// no contact applies, and the exact intersection decides.
    ///
    /// # Errors
    ///
    /// [`BooleanError::UnsupportedContact`] for a curve whose meeting with
    /// the ruling this stage cannot place (a curve of the plane that is no
    /// line, a curve of the cylinder that is no ruling or conic): its exact
    /// roots would disagree with the reading, so it is refused rather than
    /// sewn inconsistently. [`BooleanError::Evaluation`] for a curve that
    /// cannot be evaluated or inverted.
    pub(crate) fn crossing(
        &self,
        on: [&Surface; 2],
        cutter: &Surface,
        curve: &Curve3,
        tolerance: Tolerance,
    ) -> Result<Option<Vec<Scalar>>, BooleanError> {
        if self.contacts.is_empty() {
            return Ok(None);
        }
        // A curve on the cylinder, cut by the plane.
        if let Some(contact) = on.iter().find_map(|s| self.find(cutter, s)) {
            return on_cylinder(&contact.ruling, curve, tolerance);
        }
        // A curve in the plane, cut by the cylinder.
        if let Some(contact) = on.iter().find_map(|s| self.find(s, cutter)) {
            return in_plane(&contact.ruling, curve, tolerance);
        }
        Ok(None)
    }

    /// Whether `point`, on `surface`, lies where the given operands and the
    /// moved ones a contact stands for disagree: a point of the cylinder
    /// beyond the plane (it reaches through by at most the recorded
    /// distance), or a point of the plane inside the cylinder (the chord
    /// between the two exact rulings). Such a point cannot stand for its
    /// piece or region: the cuts around it were made for the moved plane.
    pub(crate) fn disputed(&self, surface: &Surface, point: Point3) -> bool {
        self.contacts.iter().any(|c| {
            let (Surface::Plane(plane), Surface::Cylinder(cylinder)) = (&c.plane, &c.cylinder)
            else {
                return false;
            };
            if is(&c.plane, surface) {
                let axis = cylinder.frame.z.normalize_or_zero();
                axis != axiolid_core::Vec3::ZERO
                    && (point - cylinder.frame.origin).cross(axis).length() < cylinder.radius
            } else if is(&c.cylinder, surface) {
                let n = plane.frame.z;
                let side = |p: Point3| (p - plane.frame.origin).dot(n);
                let (axis_side, point_side) = (side(cylinder.frame.origin), side(point));
                point_side != 0.0 && (point_side > 0.0) != (axis_side > 0.0)
            } else {
                false
            }
        })
    }
}

/// A curve on the cylinder meets the moved plane where it meets the
/// ruling: for a conic, where the ruling pierces its plane.
fn on_cylinder(
    ruling: &Line3,
    curve: &Curve3,
    tolerance: Tolerance,
) -> Result<Option<Vec<Scalar>>, BooleanError> {
    let frame = match curve {
        Curve3::Circle(c) => c.frame,
        Curve3::Ellipse(e) => e.frame,
        // Another ruling of the cylinder: parallel to the contact ruling,
        // so it meets the moved plane nowhere (the contact ruling itself
        // is no crossing).
        Curve3::Line(l) if parallel(l, ruling) => return Ok(Some(Vec::new())),
        _ => return Err(BooleanError::UnsupportedContact),
    };
    let n = frame.x.cross(frame.y).normalize_or_zero();
    let d = ruling.direction.normalize_or_zero();
    let along = d.dot(n);
    // The ruling in the curve's plane: a conic whose plane holds the axis
    // direction is no section of the cylinder.
    if along.abs() <= 1e-9 {
        return Err(BooleanError::UnsupportedContact);
    }
    let s = (frame.origin - ruling.origin).dot(n) / along;
    let point = ruling.origin + d * s;
    Ok(Some(
        on_curve(curve, point, tolerance)?.into_iter().collect(),
    ))
}

/// A line in the plane meets the cylinder, for the moved plane, where it
/// passes closest to the ruling.
fn in_plane(
    ruling: &Line3,
    curve: &Curve3,
    tolerance: Tolerance,
) -> Result<Option<Vec<Scalar>>, BooleanError> {
    let Curve3::Line(line) = curve else {
        return Err(BooleanError::UnsupportedContact);
    };
    if parallel(line, ruling) {
        // A line in the plane along the ruling meets the moved plane's
        // contact nowhere it crosses, or runs along it: no cut.
        return Ok(Some(Vec::new()));
    }
    let (d1, d2) = (line.direction, ruling.direction);
    let r = line.origin - ruling.origin;
    let (aa, bb, ab) = (d1.dot(d1), d2.dot(d2), d1.dot(d2));
    let denominator = aa * bb - ab * ab;
    let (c, f) = (d1.dot(r), d2.dot(r));
    let s = (ab * f - c * bb) / denominator;
    let t = (aa * f - ab * c) / denominator;
    let (p, q) = (line.origin + d1 * s, ruling.origin + d2 * t);
    if !(s.is_finite() && t.is_finite()) {
        return Err(BooleanError::Evaluation);
    }
    Ok(Some(
        if report::near(
            ToleranceDecisionKind::TangentCrossing,
            (p - q).length(),
            tolerance,
        ) {
            vec![s]
        } else {
            Vec::new()
        },
    ))
}

/// Whether two lines run parallel, to the slack of the two-line formula.
fn parallel(a: &Line3, b: &Line3) -> bool {
    let (d1, d2) = (a.direction, b.direction);
    let (aa, bb, ab) = (d1.dot(d1), d2.dot(d2), d1.dot(d2));
    aa * bb - ab * ab <= 1e-12 * aa * bb
}

/// The parameter of `point` on `curve`, when it lies on it within the
/// tolerance (recorded as a tangent crossing beyond rounding).
fn on_curve(
    curve: &Curve3,
    point: Point3,
    tolerance: Tolerance,
) -> Result<Option<Scalar>, BooleanError> {
    let Ok(t) = locate3(curve, point, report::floored(tolerance)) else {
        return Ok(None);
    };
    let on = evaluate3(curve, t).map_err(|_| BooleanError::Evaluation)?;
    Ok(report::near(
        ToleranceDecisionKind::TangentCrossing,
        (on - point).length(),
        tolerance,
    )
    .then_some(t))
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame3, Vec3};
    use axiolid_curve::Circle3;
    use axiolid_surface::{Cylinder, Plane};
    use core::f64::consts::FRAC_PI_2;

    const EPS: Scalar = 1e-6;

    fn frame(origin: Point3, x: Vec3, y: Vec3) -> Frame3 {
        Frame3 {
            origin,
            x,
            y,
            z: x.cross(y),
        }
    }

    /// The plane `y = 1 - d` and the unit cylinder about `z`, read as
    /// touching along the ruling `x = 0, y = 1`.
    fn contact(d: Scalar) -> (Contacts, Surface, Surface) {
        let plane = Surface::Plane(Plane {
            frame: frame(Point3::new(0.0, 1.0 - d, 0.0), Vec3::Z, Vec3::X),
        });
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame(Point3::ZERO, Vec3::X, Vec3::Y),
            radius: 1.0,
        });
        let mut contacts = Contacts::default();
        contacts.push(Contact {
            plane: plane.clone(),
            cylinder: cylinder.clone(),
            ruling: Line3 {
                origin: Point3::new(0.0, 1.0, 0.0),
                direction: Vec3::Z,
            },
        });
        (contacts, plane, cylinder)
    }

    #[test]
    fn crossings_of_either_surface_meet_the_contact_ruling() {
        let d = 0.5 * EPS;
        let (contacts, plane, cylinder) = contact(d);
        let tol = Tolerance::METRE;
        let cap = Surface::Plane(Plane {
            frame: frame(Point3::new(0.0, 0.0, 2.0), Vec3::X, Vec3::Y),
        });
        // The cap's circle, on the cylinder, cut by the plane: once, on the
        // ruling, where the exact roots are two, `2 sqrt(2 d)` apart.
        let circle = Curve3::Circle(Circle3 {
            frame: frame(Point3::new(0.0, 0.0, 2.0), Vec3::X, Vec3::Y),
            radius: 1.0,
        });
        let hits = contacts
            .crossing([&cylinder, &cap], &plane, &circle, tol)
            .unwrap()
            .expect("a contact applies");
        assert_eq!(hits.len(), 1);
        assert!((hits[0] - FRAC_PI_2).abs() < 1e-12, "{hits:?}");
        // The chord of the plane across the cap, cut by the cylinder: once,
        // closest to the ruling, `d` from it.
        let chord = Curve3::Line(Line3 {
            origin: Point3::new(-3.0, 1.0 - d, 2.0),
            direction: Vec3::X,
        });
        let hits = contacts
            .crossing([&plane, &cap], &cylinder, &chord, tol)
            .unwrap()
            .expect("a contact applies");
        assert_eq!(hits.len(), 1);
        assert!((hits[0] - 3.0).abs() < 1e-12, "{hits:?}");
        // A line of the plane along the ruling meets it nowhere it crosses.
        let along = Curve3::Line(Line3 {
            origin: Point3::new(0.5, 1.0 - d, 0.0),
            direction: Vec3::Z,
        });
        let hits = contacts.crossing([&plane, &cap], &cylinder, &along, tol);
        assert_eq!(hits, Ok(Some(Vec::new())));
        // Neither surface of the curve is in a contact with the cutter: the
        // exact intersection decides.
        assert_eq!(
            contacts.crossing([&cap, &cap], &cylinder, &chord, tol),
            Ok(None)
        );
        assert_eq!(
            Contacts::default().crossing([&cylinder, &cap], &plane, &circle, tol),
            Ok(None)
        );
        // A circle of the plane cut by the cylinder: its meeting with the
        // ruling is not placed here, so it is refused by name rather than
        // left to exact roots that disagree with the reading.
        let in_plane = Curve3::Circle(Circle3 {
            frame: frame(Point3::new(0.0, 1.0 - d, 0.0), Vec3::Z, Vec3::X),
            radius: 0.5,
        });
        assert_eq!(
            contacts.crossing([&plane, &cap], &cylinder, &in_plane, tol),
            Err(BooleanError::UnsupportedContact)
        );
    }

    #[test]
    fn only_points_the_given_and_the_moved_plane_disagree_on_are_disputed() {
        let d = 0.5 * EPS;
        let (contacts, plane, cylinder) = contact(d);
        // The plane's chord inside the cylinder, and the cylinder's sliver
        // beyond the plane.
        assert!(contacts.disputed(&plane, Point3::new(0.0, 1.0 - d, 5.0)));
        assert!(!contacts.disputed(&plane, Point3::new(0.1, 1.0 - d, 5.0)));
        assert!(contacts.disputed(&cylinder, Point3::new(0.0, 1.0, 5.0)));
        assert!(!contacts.disputed(&cylinder, Point3::new(1.0, 0.0, 5.0)));
        assert!(!contacts.disputed(&cylinder, Point3::new(0.0, -1.0, 5.0)));
        // A plane short of the cylinder disputes nothing.
        let (short, plane, cylinder) = contact(-d);
        assert!(!short.disputed(&plane, Point3::new(0.0, 1.0 + d, 5.0)));
        assert!(!short.disputed(&cylinder, Point3::new(0.0, 1.0, 5.0)));
        // Other surfaces are never disputed.
        let other = Surface::Plane(Plane {
            frame: frame(Point3::ZERO, Vec3::X, Vec3::Y),
        });
        assert!(!contacts.disputed(&other, Point3::new(0.0, 1.0, 0.0)));
    }
}
