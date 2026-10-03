//! Sound boxes around faces, for skipping work that cannot matter (#228).
//!
//! A face lies on its surface over its certified parameter box
//! (`FaceDomain::bounds`, whose boundary arcs are split at their turning
//! points, so their ends bound them). Each analytic surface maps that box
//! into a space box in closed form: a plane affinely (the four corners), a
//! cylinder, elliptical cylinder or cone as its axis segment widened by its
//! largest radius along each frame axis, a sphere or torus as its centre
//! widened likewise. Every box is then enlarged by the caller's linear
//! tolerance plus a rounding margin, so two faces that touch, or come within
//! tolerance, always have overlapping boxes: a pair is skipped only when no
//! decision about it could be taken. A B-spline face has no box and is never
//! skipped.

use axiolid_core::{Interval, Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::Curve3;
use axiolid_surface::Surface;

/// An axis-aligned box in space.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Aabb {
    lo: Point3,
    hi: Point3,
}

impl Aabb {
    /// Whether two boxes share a point.
    pub(crate) fn overlaps(&self, other: &Aabb) -> bool {
        self.lo.cmple(other.hi).all() && other.lo.cmple(self.hi).all()
    }

    /// Whether the ray `origin + t direction`, `t >= -slack`, can meet the
    /// box (the slab test, with the start widened by `slack` along the ray).
    pub(crate) fn met_by_ray(&self, origin: Point3, direction: Vec3, slack: Scalar) -> bool {
        let (mut near, mut far) = (-slack, Scalar::INFINITY);
        for axis in 0..3 {
            let (o, d) = (origin[axis], direction[axis]);
            let (lo, hi) = (self.lo[axis], self.hi[axis]);
            if d == 0.0 {
                if o < lo || o > hi {
                    return false;
                }
                continue;
            }
            let (a, b) = ((lo - o) / d, (hi - o) / d);
            near = near.max(a.min(b));
            far = far.min(a.max(b));
            if near > far {
                return false;
            }
        }
        true
    }

    /// The common part of two boxes, or `None` when they are apart.
    pub(crate) fn intersection(&self, other: &Aabb) -> Option<Aabb> {
        self.overlaps(other).then(|| Aabb {
            lo: self.lo.max(other.lo),
            hi: self.hi.min(other.hi),
        })
    }

    /// Parameters where a line, circle or ellipse crosses a face plane of
    /// the box; none for other curves. Between consecutive ones the curve
    /// is wholly inside the box or wholly outside it.
    pub(crate) fn crossings(&self, curve: &Curve3) -> Vec<Scalar> {
        let mut out = Vec::new();
        match curve {
            Curve3::Line(line) => {
                for axis in 0..3 {
                    let d = line.direction[axis];
                    if d != 0.0 {
                        for bound in [self.lo[axis], self.hi[axis]] {
                            out.push((bound - line.origin[axis]) / d);
                        }
                    }
                }
            }
            Curve3::Circle(_) | Curve3::Ellipse(_) => {
                let (centre, p, q) = conic_axes(curve);
                for axis in 0..3 {
                    // centre + p cos t + q sin t = bound, per axis.
                    let (a, b) = (p[axis], q[axis]);
                    let r = a.hypot(b);
                    if r == 0.0 {
                        continue;
                    }
                    let phase = b.atan2(a);
                    for bound in [self.lo[axis], self.hi[axis]] {
                        let k = (bound - centre[axis]) / r;
                        if k.abs() <= 1.0 {
                            let turn = k.acos();
                            out.extend(
                                [phase + turn, phase - turn]
                                    .map(|t| t.rem_euclid(core::f64::consts::TAU)),
                            );
                        }
                    }
                }
            }
            _ => {}
        }
        out.retain(|t| t.is_finite());
        out
    }

    /// The smallest box holding both.
    pub(crate) fn union(&self, other: &Aabb) -> Aabb {
        Aabb {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    /// Whether `point` lies in the box.
    pub(crate) fn contains(&self, point: Point3) -> bool {
        self.lo.cmple(point).all() && point.cmple(self.hi).all()
    }
}

/// A sound box around the face on `surface` whose parameters lie in
/// `[lo, hi]`, enlarged by the linear tolerance; `None` for a surface
/// without a closed-form bound (a B-spline) or non-finite input.
pub(crate) fn face_box(
    surface: &Surface,
    lo: Point2,
    hi: Point2,
    tolerance: Tolerance,
) -> Option<Aabb> {
    // Per frame axis, how far a radius `r` in the frame's x-y plane reaches.
    let radial = |x: Vec3, y: Vec3, a: Scalar, b: Scalar| x.abs() * a + y.abs() * b;
    let (centre_lo, centre_hi, reach) = match surface {
        Surface::Plane(p) => {
            let f = p.frame;
            let corners = [
                f.origin + f.x * lo.x + f.y * lo.y,
                f.origin + f.x * hi.x + f.y * lo.y,
                f.origin + f.x * lo.x + f.y * hi.y,
                f.origin + f.x * hi.x + f.y * hi.y,
            ];
            let mut a = corners[0];
            let mut b = corners[0];
            for c in corners {
                a = a.min(c);
                b = b.max(c);
            }
            (a, b, Vec3::ZERO)
        }
        Surface::Cylinder(c) => {
            let (a, b) = axis_segment(c.frame.origin, c.frame.z, lo.y, hi.y);
            (a, b, radial(c.frame.x, c.frame.y, c.radius, c.radius))
        }
        Surface::EllipticalCylinder(c) => {
            let (a, b) = axis_segment(c.frame.origin, c.frame.z, lo.y, hi.y);
            (
                a,
                b,
                radial(c.frame.x, c.frame.y, c.semi_axis_x, c.semi_axis_y),
            )
        }
        Surface::Cone(c) => {
            let (a, b) = axis_segment(c.frame.origin, c.frame.z, lo.y, hi.y);
            let slope = c.semi_angle.tan();
            let r = (c.radius + lo.y * slope)
                .abs()
                .max((c.radius + hi.y * slope).abs());
            (a, b, radial(c.frame.x, c.frame.y, r, r))
        }
        Surface::Sphere(s) => (
            s.frame.origin,
            s.frame.origin,
            radial(s.frame.x, s.frame.y, s.radius, s.radius) + s.frame.z.abs() * s.radius,
        ),
        Surface::Torus(t) => {
            let ring = t.major_radius + t.minor_radius;
            (
                t.frame.origin,
                t.frame.origin,
                radial(t.frame.x, t.frame.y, ring, ring) + t.frame.z.abs() * t.minor_radius,
            )
        }
        _ => return None,
    };
    padded(centre_lo - reach, centre_hi + reach, tolerance)
}

fn axis_segment(origin: Point3, axis: Vec3, v0: Scalar, v1: Scalar) -> (Point3, Point3) {
    let (p, q) = (origin + axis * v0, origin + axis * v1);
    (p.min(q), p.max(q))
}

/// A conic's centre and the vectors its cosine and sine terms run along.
fn conic_axes(curve: &Curve3) -> (Point3, Vec3, Vec3) {
    match curve {
        Curve3::Circle(c) => (c.frame.origin, c.frame.x * c.radius, c.frame.y * c.radius),
        Curve3::Ellipse(e) => (
            e.frame.origin,
            e.frame.x * e.semi_axis_x,
            e.frame.y * e.semi_axis_y,
        ),
        _ => (Point3::ZERO, Vec3::ZERO, Vec3::ZERO),
    }
}

/// A sound box around an edge, enlarged by the linear tolerance: a line's
/// span, or a whole circle or ellipse; `None` for other curves.
pub(crate) fn edge_box(curve: &Curve3, span: Interval, tolerance: Tolerance) -> Option<Aabb> {
    let (a, b) = match curve {
        Curve3::Line(line) => {
            let (p, q) = (
                line.origin + line.direction * span.start,
                line.origin + line.direction * span.end,
            );
            (p.min(q), p.max(q))
        }
        Curve3::Circle(_) | Curve3::Ellipse(_) => {
            let (centre, p, q) = conic_axes(curve);
            let reach = Vec3::new(p.x.hypot(q.x), p.y.hypot(q.y), p.z.hypot(q.z));
            (centre - reach, centre + reach)
        }
        _ => return None,
    };
    padded(a, b, tolerance)
}

fn padded(a: Point3, b: Point3, tolerance: Tolerance) -> Option<Aabb> {
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    // The tolerance, plus a margin for the rounding of the bound itself.
    let size = a.abs().max(b.abs()).max_element();
    let pad = tolerance.linear() + 8.0 * Scalar::EPSILON * (1.0 + size);
    Some(Aabb {
        lo: a - Vec3::splat(pad),
        hi: b + Vec3::splat(pad),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_brep::ExactBRep;
    use axiolid_construct::extrude::extrude_profile_exact;
    use axiolid_construct::revolve_exact::revolve_profile_exact;
    use axiolid_core::Transform3;
    use axiolid_evaluate::evaluate3;
    use axiolid_measure::FaceDomain;
    use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

    fn tol() -> Tolerance {
        Tolerance::METRE
    }

    fn placed(brep: ExactBRep) -> ExactBRep {
        brep.transformed(
            &(Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
                * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)),
        )
        .expect("rigid")
    }

    /// Every point of every edge lies in the box of each face it bounds and
    /// in its own box: the boxes are sound, so a pair they keep apart
    /// cannot touch.
    #[test]
    fn every_edge_point_lies_in_its_faces_and_its_own_box() {
        let circle = Profile::Circle(CircleProfile {
            radius: 0.4,
            thickness: None,
        });
        let wall = Profile::Rectangle(RectangleProfile {
            x: 2.0,
            y: 0.3,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        });
        let solids = [
            extrude_profile_exact(&wall, Vec3::Z, 3.0, tol()).unwrap(),
            extrude_profile_exact(&circle, Vec3::Z, 2.0, tol()).unwrap(),
            revolve_profile_exact(
                &circle,
                Point3::new(-1.0, 0.0, 0.0),
                Vec3::Y,
                core::f64::consts::TAU,
                tol(),
            )
            .unwrap(),
        ];
        for solid in solids.into_iter().map(placed) {
            let topology = solid.topology();
            for (index, face) in topology.faces().iter().enumerate() {
                let id = topology.face_id_at(index).unwrap();
                let surface = &solid.surfaces()[face.surface.unwrap().index()];
                let domain = FaceDomain::new(&solid, id, tol()).unwrap().unwrap();
                let (lo, hi) = domain.bounds();
                let Some(face_box) = face_box(surface, lo, hi, tol()) else {
                    continue;
                };
                for bound in &face.bounds {
                    for use_ in &topology.loops()[bound.loop_id.index()].edges {
                        let edge = &topology.edges()[use_.edge.index()];
                        let curve = &solid.curves3()[edge.curve.unwrap().index()];
                        let span = solid.edge_interval(use_.edge).unwrap();
                        let own = edge_box(curve, span, tol());
                        for k in 0..=32 {
                            let t = span.start + (span.end - span.start) * k as f64 / 32.0;
                            let p = evaluate3(curve, t).unwrap();
                            assert!(face_box.contains(p), "{p:?} outside {face_box:?}");
                            assert!(own.is_none_or(|b| b.contains(p)), "{p:?} outside {own:?}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn a_line_and_a_circle_are_cut_where_they_leave_a_box() {
        let unit = Aabb {
            lo: Point3::splat(-1.0),
            hi: Point3::splat(1.0),
        };
        let line = Curve3::Line(axiolid_curve::Line3 {
            origin: Point3::new(0.0, 0.5, 0.0),
            direction: Vec3::X,
        });
        let mut cuts = unit.crossings(&line);
        cuts.sort_by(f64::total_cmp);
        assert_eq!(cuts, vec![-1.0, 1.0]);
        let circle = Curve3::Circle(axiolid_curve::Circle3 {
            frame: axiolid_core::Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 2.0,
        });
        // Leaves the box across x = +-1 and y = +-1: eight crossings.
        let cuts = unit.crossings(&circle);
        assert_eq!(cuts.len(), 8);
        for t in cuts {
            let p = evaluate3(&circle, t).unwrap();
            assert!((p.x.abs() - 1.0).abs() < 1e-12 || (p.y.abs() - 1.0).abs() < 1e-12);
        }
    }
}
