//! Certified distance between the boundaries of two exact B-reps (#125, C18).
//!
//! # What is certified
//!
//! [`boundary_distance`] returns an interval `[lower, upper]` that contains
//! the true minimum distance between the two boundaries, and two points, one
//! on each boundary, exactly `upper` apart. Both ends are guarantees, not
//! estimates:
//!
//! - `upper` is the distance between two points that lie on the boundaries:
//!   points on edges, or surface points at parameters the face's domain is
//!   certified to contain (the crate-private `Domain` classifier).
//! - `lower` comes from bounding spheres that enclose every point of a patch
//!   of a face's parameter domain, from Lipschitz bounds on the exact
//!   surface, with a margin for rounding. Patches are dropped only when
//!   certified to lie outside the face.
//!
//! A rule that compares a distance with a limit asks [`boundary_clearance`],
//! which refines only until the interval clears the limit and otherwise says
//! [`Clearance::Indeterminate`] -- a value within rounding of the limit is
//! never reported as a pass or a fail.
//!
//! # Method
//!
//! Branch and bound over pairs of elements, one from each B-rep: face
//! patches (a rectangle of a face's parameters) and edge spans. The pair
//! with the smallest lower bound is refined by splitting its larger element.
//! Refinement stops when the smallest remaining lower bound is within the
//! requested accuracy of the best upper bound, or the step budget runs out;
//! either way the interval returned is sound.
//!
//! Bounds are second order where the family allows: the range of `d . x`
//! over a patch is exact for planes, cylinders, elliptical cylinders, cones,
//! spheres and tori and for line, circle and ellipse edges, and `d` is taken
//! along the line between centres and along each patch's normal. A face
//! patch whose normals cannot point at the other element is dropped from
//! that pair (`critical_possible`): the closest pair of two separated
//! boundaries is critical on each face it lies inside or lies on an edge,
//! and edges are elements of their own.
//!
//! # Convergence
//!
//! An isolated nearest pair (pole to pole, apex to sphere, a wall to a
//! block face) closes to `1e-9` in well under a second. Where the nearest
//! points form a whole line -- two parallel columns -- every slice along
//! that line is a near-minimal pair and refinement slows; ask such cases for
//! a looser accuracy, or use [`boundary_clearance`], which stops as soon as
//! the limit is cleared. The step budget is fixed; when it runs out the
//! interval returned is still sound, only wider.
//!
//! # Scope
//!
//! This is the distance between BOUNDARIES. Two solids that overlap measure
//! the distance between their surfaces where they cross (zero), but a solid
//! wholly inside another measures the gap between the two boundaries, not
//! zero. Containment is a separate classification.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use axiolid_brep::ExactBRep;
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::Curve3;
use axiolid_evaluate::evaluate3;
use axiolid_evaluate::surface::{evaluate, normal};
use axiolid_surface::Surface;

use crate::exact::ExactMeasureError;
use crate::exact_domain::Domain;

/// Pairs refined before a query stops and reports what it has.
const MAX_STEPS: usize = 400_000;

/// An interval certain to contain the distance between two boundaries.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceBounds {
    /// No two boundary points are closer than this.
    pub lower: Scalar,
    /// `point_a` and `point_b` are this far apart.
    pub upper: Scalar,
    /// A point on the first boundary.
    pub point_a: Point3,
    /// A point on the second boundary, `upper` from `point_a`.
    pub point_b: Point3,
}

/// How a certified distance compares with a limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Clearance {
    /// Certainly closer than the limit: `upper < limit`.
    Below,
    /// Certainly farther than the limit: `lower > limit`.
    Above,
    /// The interval still contains the limit.
    Indeterminate,
}

impl DistanceBounds {
    /// Compare with `limit` without rounding the interval to a verdict.
    #[must_use]
    pub fn against(&self, limit: Scalar) -> Clearance {
        if self.upper < limit {
            Clearance::Below
        } else if self.lower > limit {
            Clearance::Above
        } else {
            Clearance::Indeterminate
        }
    }
}

/// Distance between the boundaries of `a` and `b`, to within `accuracy`.
///
/// # Errors
///
/// A face whose domain cannot be bounded (an unbounded surface trimmed by a
/// pcurve family the domain classifier does not split), an evaluation
/// failure, or boundaries with no boundable edge or face.
pub fn boundary_distance(
    a: &ExactBRep,
    b: &ExactBRep,
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<DistanceBounds, ExactMeasureError> {
    let accuracy = accuracy.max(0.0);
    search(a, b, tolerance, &mut |lower, upper| {
        upper - lower <= accuracy
    })
}

/// Refine the distance only until it clears `limit`.
///
/// # Errors
///
/// As [`boundary_distance`].
pub fn boundary_clearance(
    a: &ExactBRep,
    b: &ExactBRep,
    limit: Scalar,
    tolerance: Tolerance,
) -> Result<(DistanceBounds, Clearance), ExactMeasureError> {
    let bounds = search(a, b, tolerance, &mut |lower, upper| {
        upper < limit || lower > limit
    })?;
    let clearance = bounds.against(limit);
    Ok((bounds, clearance))
}

/// What an element covers.
#[derive(Debug, Clone, Copy)]
enum Shape {
    /// A rectangle of face `face`'s parameters; `inside` once certified to
    /// lie wholly in the face.
    Face {
        face: usize,
        lo: Point2,
        hi: Point2,
        inside: bool,
    },
    /// A span of edge `edge`'s curve.
    Edge { edge: usize, t0: Scalar, t1: Scalar },
}

/// A piece of one boundary with a sphere that encloses it.
#[derive(Debug, Clone, Copy)]
struct Element {
    shape: Shape,
    centre: Point3,
    radius: Scalar,
    /// A point certainly on the boundary, when one is known.
    witness: Option<Point3>,
    /// The surface normal at the centre of a face patch.
    normal: Option<Vec3>,
    /// Half-angle of a cone about `normal` holding every normal of the
    /// patch, when the patch may be dropped from a pair whose directions
    /// its normals cannot meet (see [`critical_possible`]).
    spread: Option<Scalar>,
}

/// One B-rep, prepared for bounding.
struct Side<'a> {
    brep: &'a ExactBRep,
    domains: Vec<Option<Domain<'a>>>,
    /// Whether every edge on the face's boundary is an element, so a
    /// closest point on that boundary is found through the edges.
    edges_bounded: Vec<bool>,
    elements: Vec<Element>,
}

impl<'a> Side<'a> {
    fn new(brep: &'a ExactBRep, linear: Scalar) -> Result<Self, ExactMeasureError> {
        let topology = brep.topology();
        let mut side = Side {
            brep,
            domains: Vec::with_capacity(topology.faces().len()),
            edges_bounded: Vec::with_capacity(topology.faces().len()),
            elements: Vec::new(),
        };
        for face in topology.faces() {
            let mut bounded = true;
            for bound in &face.bounds {
                let wire = topology
                    .loops()
                    .get(bound.loop_id.index())
                    .ok_or(ExactMeasureError::DanglingReference)?;
                for use_ in &wire.edges {
                    let curve = topology.edges()[use_.edge.index()]
                        .curve
                        .and_then(|id| brep.curves3().get(id.index()));
                    bounded &= matches!(
                        curve,
                        Some(Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_))
                    );
                }
            }
            side.edges_bounded.push(bounded);
        }
        for face in topology.faces() {
            let surface = surface_of(brep, face.surface)?;
            side.domains.push(Domain::new(brep, face, surface, linear)?);
        }
        for (index, face) in topology.faces().iter().enumerate() {
            let surface = surface_of(brep, face.surface)?;
            let (lo, hi) = match &side.domains[index] {
                Some(domain) => (domain.min, domain.max),
                None => natural_range(surface).ok_or(ExactMeasureError::ParameterDomain(
                    "an unbounded face trimmed by a pcurve family the distance query cannot bound",
                ))?,
            };
            if let Some(element) = side.face_element(index, lo, hi, false)? {
                side.elements.push(element);
            }
        }
        for (index, edge) in topology.edges().iter().enumerate() {
            let Some(curve) = edge.curve.and_then(|id| brep.curves3().get(id.index())) else {
                continue;
            };
            let Some(span) = topology
                .edge_id_at(index)
                .and_then(|id| brep.edge_interval(id))
            else {
                continue;
            };
            if let Some(element) = edge_element(curve, index, span.start, span.end)? {
                side.elements.push(element);
            }
        }
        if side.elements.is_empty() {
            return Err(ExactMeasureError::Degenerate);
        }
        Ok(side)
    }

    /// A face patch, or `None` when it is certified to lie outside the face.
    fn face_element(
        &self,
        face: usize,
        lo: Point2,
        hi: Point2,
        inside: bool,
    ) -> Result<Option<Element>, ExactMeasureError> {
        let topology = self.brep.topology();
        let surface = surface_of(self.brep, topology.faces()[face].surface)?;
        let mid = (lo + hi) * 0.5;
        let mut inside = inside;
        if !inside {
            if let Some(domain) = &self.domains[face] {
                if !domain.touches(lo, hi)? {
                    match domain.contains(mid)? {
                        Some(false) => return Ok(None),
                        Some(true) => inside = true,
                        None => {}
                    }
                }
            }
        }
        let (centre, radius) = patch_sphere(surface, lo, hi)?;
        let witness = if inside {
            Some(point_on(surface, mid)?)
        } else {
            None
        };
        let normal = normal(surface, mid.x, mid.y).ok();
        let spread = if self.edges_bounded[face] {
            normal_spread(surface, lo, hi)
        } else {
            None
        };
        Ok(Some(Element {
            normal,
            spread,
            shape: Shape::Face {
                face,
                lo,
                hi,
                inside,
            },
            centre,
            radius,
            witness,
        }))
    }

    /// The element's two halves, less any certified outside the face: a
    /// face patch across its metrically longer side, an edge span in the
    /// middle.
    fn split(&self, element: &Element) -> Result<Vec<Element>, ExactMeasureError> {
        match element.shape {
            Shape::Face {
                face,
                lo,
                hi,
                inside,
            } => {
                let surface = surface_of(self.brep, self.brep.topology().faces()[face].surface)?;
                let (lu, lv) = lipschitz(surface, lo, hi)?;
                let mid = (lo + hi) * 0.5;
                let halves = if (hi.x - lo.x).abs() * lu >= (hi.y - lo.y).abs() * lv {
                    [
                        (lo, Point2::new(mid.x, hi.y)),
                        (Point2::new(mid.x, lo.y), hi),
                    ]
                } else {
                    [
                        (lo, Point2::new(hi.x, mid.y)),
                        (Point2::new(lo.x, mid.y), hi),
                    ]
                };
                let mut out = Vec::with_capacity(2);
                for (lo, hi) in halves {
                    if let Some(child) = self.face_element(face, lo, hi, inside)? {
                        out.push(child);
                    }
                }
                Ok(out)
            }
            Shape::Edge { edge, t0, t1 } => {
                let curve = self.brep.topology().edges()[edge]
                    .curve
                    .and_then(|id| self.brep.curves3().get(id.index()))
                    .ok_or(ExactMeasureError::DanglingReference)?;
                let tm = 0.5 * (t0 + t1);
                let mut out = Vec::with_capacity(2);
                for (a, b) in [(t0, tm), (tm, t1)] {
                    if let Some(child) = edge_element(curve, edge, a, b)? {
                        out.push(child);
                    }
                }
                Ok(out)
            }
        }
    }
}

fn surface_of(
    brep: &ExactBRep,
    id: Option<axiolid_brep::SurfaceId>,
) -> Result<&Surface, ExactMeasureError> {
    let id = id.ok_or(ExactMeasureError::MissingSurface)?;
    brep.surfaces()
        .get(id.index())
        .ok_or(ExactMeasureError::DanglingReference)
}

fn point_on(surface: &Surface, at: Point2) -> Result<Point3, ExactMeasureError> {
    evaluate(surface, at.x, at.y).map_err(|_| ExactMeasureError::Evaluation)
}

/// The whole parameter range of a closed surface, for a face whose trim
/// cannot be classified.
fn natural_range(surface: &Surface) -> Option<(Point2, Point2)> {
    use core::f64::consts::{FRAC_PI_2, TAU};
    match surface {
        Surface::Sphere(_) => Some((Point2::new(0.0, -FRAC_PI_2), Point2::new(TAU, FRAC_PI_2))),
        Surface::Torus(_) => Some((Point2::ZERO, Point2::new(TAU, TAU))),
        _ => None,
    }
}

/// Largest length a frame axis carries (1 for an orthonormal frame).
fn frame_scale(frame: &axiolid_core::Frame3) -> Scalar {
    frame.x.length().max(frame.y.length()).max(frame.z.length())
}

/// Bounds on `|S_u|` and `|S_v|` over the patch.
fn lipschitz(
    surface: &Surface,
    lo: Point2,
    hi: Point2,
) -> Result<(Scalar, Scalar), ExactMeasureError> {
    let bounds = match surface {
        Surface::Plane(p) => (p.frame.x.length(), p.frame.y.length()),
        Surface::Cylinder(c) => (c.radius.abs() * frame_scale(&c.frame), c.frame.z.length()),
        Surface::EllipticalCylinder(c) => (
            c.semi_axis_x.abs().max(c.semi_axis_y.abs()) * frame_scale(&c.frame),
            c.frame.z.length(),
        ),
        Surface::Cone(c) => {
            let slope = c.semi_angle.tan();
            let radius = (c.radius + lo.y * slope)
                .abs()
                .max((c.radius + hi.y * slope).abs());
            let scale = frame_scale(&c.frame);
            (radius * scale, scale * (1.0 + slope * slope).sqrt())
        }
        Surface::Sphere(s) => {
            // |S_u| = r cos v, which vanishes at the poles: bounding it by r
            // would split near-polar patches round the pole for nothing.
            let r = s.radius.abs() * frame_scale(&s.frame);
            (r * max_cos(lo.y, hi.y), r)
        }
        Surface::Torus(t) => {
            let scale = frame_scale(&t.frame);
            (
                (t.major_radius.abs() + t.minor_radius.abs() * max_cos(lo.y, hi.y)) * scale,
                t.minor_radius.abs() * scale,
            )
        }
        // A B-spline patch is bounded by its control net instead.
        Surface::BSpline(_) => (0.0, 0.0),
        _ => {
            return Err(ExactMeasureError::NonPlanarFace(crate::exact::family(
                surface,
            )))
        }
    };
    if bounds.0.is_finite() && bounds.1.is_finite() {
        Ok(bounds)
    } else {
        Err(ExactMeasureError::Evaluation)
    }
}

/// Largest `|cos v|` over `[a, b]`.
fn max_cos(a: Scalar, b: Scalar) -> Scalar {
    let (a, b) = (a.min(b), a.max(b));
    let k = (a / core::f64::consts::PI).ceil();
    if k * core::f64::consts::PI <= b {
        1.0
    } else {
        a.cos().abs().max(b.cos().abs())
    }
}

/// A sphere enclosing every surface point of the patch.
///
/// From the patch centre, any point is reached by a path along `u` then
/// along `v`, no longer than `du/2 |S_u|max + dv/2 |S_v|max`. A B-spline
/// patch lies in the convex hull of the surface's control net (positive
/// weights), which is not refined but is sound.
fn patch_sphere(
    surface: &Surface,
    lo: Point2,
    hi: Point2,
) -> Result<(Point3, Scalar), ExactMeasureError> {
    if let Surface::BSpline(spline) = surface {
        if let Some(weights) = &spline.weights {
            if weights.iter().flatten().any(|w| w.is_nan() || *w <= 0.0) {
                return Err(ExactMeasureError::NonPlanarFace(
                    "non-positive-weight B-spline",
                ));
            }
        }
        let points: Vec<Point3> = spline.control_points.iter().flatten().copied().collect();
        if points.is_empty() {
            return Err(ExactMeasureError::Degenerate);
        }
        let (mut min, mut max) = (points[0], points[0]);
        for p in &points {
            min = min.min(*p);
            max = max.max(*p);
        }
        let centre = (min + max) * 0.5;
        let radius = points
            .iter()
            .map(|p| (*p - centre).length())
            .fold(0.0, Scalar::max);
        return Ok((centre, pad(centre, radius)));
    }
    let centre = point_on(surface, (lo + hi) * 0.5)?;
    let (lu, lv) = lipschitz(surface, lo, hi)?;
    let radius = 0.5 * ((hi.x - lo.x).abs() * lu + (hi.y - lo.y).abs() * lv);
    Ok((centre, pad(centre, radius)))
}

/// Widen a bounding radius by the rounding its centre may carry.
fn pad(centre: Point3, radius: Scalar) -> Scalar {
    radius + 1e-12 * (centre.length() + radius) + Scalar::MIN_POSITIVE
}

/// An edge span with its enclosing sphere and its midpoint as a witness,
/// or `None` for a curve family with no derivative bound here.
fn edge_element(
    curve: &Curve3,
    edge: usize,
    t0: Scalar,
    t1: Scalar,
) -> Result<Option<Element>, ExactMeasureError> {
    let speed = match curve {
        Curve3::Line(line) => line.direction.length(),
        Curve3::Circle(circle) => circle.radius.abs() * frame_scale(&circle.frame),
        Curve3::Ellipse(ellipse) => {
            ellipse.semi_axis_x.abs().max(ellipse.semi_axis_y.abs()) * frame_scale(&ellipse.frame)
        }
        _ => return Ok(None),
    };
    let centre = evaluate3(curve, 0.5 * (t0 + t1)).map_err(|_| ExactMeasureError::Evaluation)?;
    let radius = pad(centre, 0.5 * (t1 - t0).abs() * speed);
    Ok(Some(Element {
        normal: None,
        spread: None,
        shape: Shape::Edge { edge, t0, t1 },
        centre,
        radius,
        witness: Some(centre),
    }))
}

/// Half-angle bound on how far the surface normal turns across the patch.
fn normal_spread(surface: &Surface, lo: Point2, hi: Point2) -> Option<Scalar> {
    let (du, dv) = ((hi.x - lo.x).abs(), (hi.y - lo.y).abs());
    let spread = match surface {
        Surface::Plane(_) => 0.0,
        // The normal turns with u at unit rate on a circular section.
        Surface::Cylinder(_) => 0.5 * du,
        Surface::Cone(c) => {
            // The apex is not smooth: a nearest point there is neither
            // critical nor on an edge, so a patch reaching it is never
            // dropped.
            let slope = c.semi_angle.tan();
            let apex = -c.radius / slope;
            if !apex.is_finite() || (apex >= lo.y.min(hi.y) - 1e-9 && apex <= lo.y.max(hi.y) + 1e-9)
            {
                return None;
            }
            0.5 * du
        }
        // On an ellipse it turns at most max/min times faster.
        Surface::EllipticalCylinder(c) => {
            let (a, b) = (c.semi_axis_x.abs(), c.semi_axis_y.abs());
            0.5 * du * a.max(b) / a.min(b)
        }
        Surface::Sphere(_) | Surface::Torus(_) => 0.5 * (du + dv),
        _ => return None,
    };
    spread.is_finite().then_some(spread + 1e-9)
}

/// Whether some pair of points of the two elements can be critical for
/// the distance on the face side: the segment joining them along the
/// surface normal there.
///
/// A closest pair of points of two separated boundaries is either critical
/// on each face it lies inside, or lies on an edge -- and every edge is an
/// element of its own. So a face patch whose normals cannot meet any
/// direction towards the other element holds no closest point that the
/// edges do not already hold, and the pair can be dropped. This is what
/// lets the bound close where the nearest points run along an edge: the
/// patches straddling that edge also hold points just outside the face,
/// closer than the true distance, that no bound on the patch can exclude.
fn critical_possible(face: &Element, other: &Element) -> bool {
    let (Some(normal), Some(spread)) = (face.normal, face.spread) else {
        return true;
    };
    let offset = other.centre - face.centre;
    let gap = offset.length();
    let reach = face.radius + other.radius;
    if gap.is_nan() || gap <= reach {
        return true;
    }
    let aperture = (reach / gap).asin();
    let angle = spread + aperture + 1e-9;
    if angle >= core::f64::consts::FRAC_PI_2 {
        return true;
    }
    (offset / gap).dot(normal).abs() >= angle.cos() * normal.length()
}

/// Range of `a cos t + b sin t` over `[t0, t1]`.
fn trig_range(a: Scalar, b: Scalar, t0: Scalar, t1: Scalar) -> (Scalar, Scalar) {
    let (t0, t1) = (t0.min(t1), t0.max(t1));
    let f = |t: Scalar| a * t.cos() + b * t.sin();
    let (mut lo, mut hi) = (f(t0).min(f(t1)), f(t0).max(f(t1)));
    let amplitude = a.hypot(b);
    let peak = b.atan2(a);
    let reaches = |angle: Scalar| {
        let k = ((t0 - angle) / core::f64::consts::TAU).ceil();
        angle + k * core::f64::consts::TAU <= t1
    };
    if reaches(peak) {
        hi = amplitude;
    }
    if reaches(peak + core::f64::consts::PI) {
        lo = -amplitude;
    }
    (lo, hi)
}

/// Range of `d . x` over an element's enclosing sphere.
fn sphere_range(element: &Element, d: Vec3) -> (Scalar, Scalar) {
    let c = element.centre.dot(d);
    (c - element.radius, c + element.radius)
}

impl Side<'_> {
    /// Exact range of `d . x` over the element where the family allows,
    /// else the enclosing sphere's.
    fn project(&self, element: &Element, d: Vec3) -> Result<(Scalar, Scalar), ExactMeasureError> {
        let sphere = sphere_range(element, d);
        let exact = match element.shape {
            Shape::Face { face, lo, hi, .. } => {
                let surface = surface_of(self.brep, self.brep.topology().faces()[face].surface)?;
                match surface {
                    Surface::Plane(p) => {
                        let base = p.frame.origin.dot(d);
                        let (x, y) = (p.frame.x.dot(d), p.frame.y.dot(d));
                        let values = [
                            base + x * lo.x + y * lo.y,
                            base + x * hi.x + y * lo.y,
                            base + x * lo.x + y * hi.y,
                            base + x * hi.x + y * hi.y,
                        ];
                        Some(
                            values
                                .iter()
                                .fold((Scalar::INFINITY, Scalar::NEG_INFINITY), |(a, b), v| {
                                    (a.min(*v), b.max(*v))
                                }),
                        )
                    }
                    Surface::Cylinder(c) => {
                        let (a, b) = trig_range(
                            c.radius * c.frame.x.dot(d),
                            c.radius * c.frame.y.dot(d),
                            lo.x,
                            hi.x,
                        );
                        let z = c.frame.z.dot(d);
                        let base = c.frame.origin.dot(d);
                        Some((
                            base + a + (z * lo.y).min(z * hi.y),
                            base + b + (z * lo.y).max(z * hi.y),
                        ))
                    }
                    Surface::EllipticalCylinder(c) => {
                        let (a, b) = trig_range(
                            c.semi_axis_x * c.frame.x.dot(d),
                            c.semi_axis_y * c.frame.y.dot(d),
                            lo.x,
                            hi.x,
                        );
                        let z = c.frame.z.dot(d);
                        let base = c.frame.origin.dot(d);
                        Some((
                            base + a + (z * lo.y).min(z * hi.y),
                            base + b + (z * lo.y).max(z * hi.y),
                        ))
                    }
                    Surface::Cone(c) => {
                        // Linear in v for fixed u: the extremes sit on the
                        // two v edges of the rectangle.
                        let slope = c.semi_angle.tan();
                        let base = c.frame.origin.dot(d);
                        let (x, y, z) = (c.frame.x.dot(d), c.frame.y.dot(d), c.frame.z.dot(d));
                        let mut range = (Scalar::INFINITY, Scalar::NEG_INFINITY);
                        for v in [lo.y, hi.y] {
                            let r = c.radius + v * slope;
                            let (a, b) = trig_range(r * x, r * y, lo.x, hi.x);
                            range = (range.0.min(base + a + z * v), range.1.max(base + b + z * v));
                        }
                        Some(range)
                    }
                    Surface::Sphere(sphere) => {
                        // d.S = d.c + r (cos v W(u) + sin v Z), W linear in
                        // the u-trig term: the extremes take W at its ends.
                        let r = sphere.radius;
                        let base = sphere.frame.origin.dot(d);
                        let (w_lo, w_hi) =
                            trig_range(sphere.frame.x.dot(d), sphere.frame.y.dot(d), lo.x, hi.x);
                        let z = sphere.frame.z.dot(d);
                        let mut range = (Scalar::INFINITY, Scalar::NEG_INFINITY);
                        for w in [w_lo, w_hi] {
                            let (a, b) = trig_range(r * w, r * z, lo.y, hi.y);
                            range = (range.0.min(base + a), range.1.max(base + b));
                        }
                        Some(range)
                    }
                    Surface::Torus(torus) => {
                        // d.S = d.c + R W(u) + r (cos v W(u) + sin v Z); for
                        // R > r the coefficient of W is positive, so again W
                        // takes its ends.
                        let (big, small) = (torus.major_radius, torus.minor_radius);
                        if big.is_nan() || big <= small.abs() {
                            return Ok(sphere_range(element, d));
                        }
                        let base = torus.frame.origin.dot(d);
                        let (w_lo, w_hi) =
                            trig_range(torus.frame.x.dot(d), torus.frame.y.dot(d), lo.x, hi.x);
                        let z = torus.frame.z.dot(d);
                        let mut range = (Scalar::INFINITY, Scalar::NEG_INFINITY);
                        for w in [w_lo, w_hi] {
                            let (a, b) = trig_range(small * w, small * z, lo.y, hi.y);
                            range = (
                                range.0.min(base + big * w + a),
                                range.1.max(base + big * w + b),
                            );
                        }
                        Some(range)
                    }
                    _ => None,
                }
            }
            Shape::Edge { edge, t0, t1 } => {
                let curve = self.brep.topology().edges()[edge]
                    .curve
                    .and_then(|id| self.brep.curves3().get(id.index()))
                    .ok_or(ExactMeasureError::DanglingReference)?;
                match curve {
                    Curve3::Line(line) => {
                        let (a, b) = (
                            (line.origin + line.direction * t0).dot(d),
                            (line.origin + line.direction * t1).dot(d),
                        );
                        Some((a.min(b), a.max(b)))
                    }
                    Curve3::Circle(c) => {
                        let (a, b) = trig_range(
                            c.radius * c.frame.x.dot(d),
                            c.radius * c.frame.y.dot(d),
                            t0,
                            t1,
                        );
                        let base = c.frame.origin.dot(d);
                        Some((base + a, base + b))
                    }
                    Curve3::Ellipse(e) => {
                        let (a, b) = trig_range(
                            e.semi_axis_x * e.frame.x.dot(d),
                            e.semi_axis_y * e.frame.y.dot(d),
                            t0,
                            t1,
                        );
                        let base = e.frame.origin.dot(d);
                        Some((base + a, base + b))
                    }
                    _ => None,
                }
            }
        };
        Ok(match exact {
            Some((lo, hi)) => {
                // Rounding in the projection itself.
                let pad =
                    1e-12 * (element.centre.length() + element.radius + lo.abs().max(hi.abs()));
                (lo.max(sphere.0) - pad, hi.min(sphere.1) + pad)
            }
            None => sphere,
        })
    }
}

/// The largest separation certified along the centre line or either
/// patch's normal, and never less than the spheres give.
fn lower_bound(
    side_a: &Side<'_>,
    a: &Element,
    side_b: &Side<'_>,
    b: &Element,
) -> Result<Scalar, ExactMeasureError> {
    let gap = (a.centre - b.centre).length();
    let rounding = 1e-12 * (a.centre.length() + b.centre.length() + gap);
    let mut best = (gap - a.radius - b.radius - rounding).max(0.0);
    let directions = [
        (gap > 0.0).then(|| (b.centre - a.centre) / gap),
        a.normal,
        b.normal,
    ];
    for d in directions.into_iter().flatten() {
        let (a_lo, a_hi) = side_a.project(a, d)?;
        let (b_lo, b_hi) = side_b.project(b, d)?;
        best = best.max(b_lo - a_hi).max(a_lo - b_hi);
    }
    Ok(best)
}

/// Total order on finite lower bounds for the heap.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Key(Scalar);

impl Eq for Key {}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

fn search(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
    done: &mut dyn FnMut(Scalar, Scalar) -> bool,
) -> Result<DistanceBounds, ExactMeasureError> {
    let linear = tolerance.linear().max(1e-12);
    let side_a = Side::new(a, linear)?;
    let side_b = Side::new(b, linear)?;
    let mut elements_a = side_a.elements.clone();
    let mut elements_b = side_b.elements.clone();

    let mut best: Option<(Scalar, Point3, Point3)> = None;
    let mut heap = BinaryHeap::new();
    let viable = |a: &Element, b: &Element| critical_possible(a, b) && critical_possible(b, a);
    for (i, ea) in elements_a.iter().enumerate() {
        for (j, eb) in elements_b.iter().enumerate() {
            if viable(ea, eb) {
                heap.push(Reverse((Key(lower_bound(&side_a, ea, &side_b, eb)?), i, j)));
            }
        }
    }

    let mut lower = 0.0;
    let mut steps = 0;
    while let Some(Reverse((Key(bound), i, j))) = heap.pop() {
        lower = bound;
        let (ea, eb) = (elements_a[i], elements_b[j]);
        if let (Some(wa), Some(wb)) = (ea.witness, eb.witness) {
            let d = (wa - wb).length();
            if best.is_none_or(|(current, _, _)| d < current) {
                best = Some((d, wa, wb));
            }
        }
        let upper = best.map_or(Scalar::INFINITY, |(d, _, _)| d);
        if upper.is_finite() && (bound >= upper || done(bound, upper)) {
            break;
        }
        steps += 1;
        if steps > MAX_STEPS {
            break;
        }
        // Refine the larger of the two; a pair that can shrink no further
        // holds the lower bound where it is.
        // Refine the larger of the two; a pair that can shrink no further
        // holds the lower bound where it is.
        let split_a = ea.radius >= eb.radius;
        let children = if split_a {
            side_a.split(&ea)?
        } else {
            side_b.split(&eb)?
        };
        let parent = if split_a { ea.radius } else { eb.radius };
        if children.iter().any(|child| child.radius >= parent) {
            heap.push(Reverse((Key(bound), i, j)));
            break;
        }
        for child in children {
            if split_a {
                elements_a.push(child);
                let index = elements_a.len() - 1;
                if viable(&child, &eb) {
                    heap.push(Reverse((
                        Key(lower_bound(&side_a, &child, &side_b, &eb)?),
                        index,
                        j,
                    )));
                }
            } else {
                elements_b.push(child);
                let index = elements_b.len() - 1;
                if viable(&ea, &child) {
                    heap.push(Reverse((
                        Key(lower_bound(&side_a, &ea, &side_b, &child)?),
                        i,
                        index,
                    )));
                }
            }
        }
    }
    // Every pair still queued bounds from below what is left unexplored.
    if let Some(Reverse((Key(bound), _, _))) = heap.peek() {
        lower = lower.min(*bound);
    }
    let (upper, point_a, point_b) = best.ok_or(ExactMeasureError::NotConverged)?;
    Ok(DistanceBounds {
        lower: lower.min(upper),
        upper,
        point_a,
        point_b,
    })
}

#[cfg(test)]
mod tests {
    //! The two claims the lower bound rests on, checked by dense sampling:
    //! every point of a patch lies in its sphere, and every normal lies in
    //! its cone.

    use super::{normal_spread, patch_sphere};
    use axiolid_core::{Frame3, Point2, Point3, Vec3};
    use axiolid_evaluate::surface::{evaluate, normal};
    use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface, Torus};

    fn frame() -> Frame3 {
        // Off the origin and turned, so no term vanishes by symmetry.
        let x = Vec3::new(0.6, 0.8, 0.0);
        let z = Vec3::new(0.0, 0.0, 1.0);
        Frame3 {
            origin: Point3::new(1.5, -2.0, 0.75),
            x,
            y: z.cross(x),
            z,
        }
    }

    fn families() -> Vec<(Surface, Point2, Point2)> {
        let f = frame();
        vec![
            (
                Surface::Plane(Plane { frame: f }),
                Point2::new(-1.0, 0.5),
                Point2::new(2.0, 1.25),
            ),
            (
                Surface::Cylinder(Cylinder {
                    frame: f,
                    radius: 2.5,
                }),
                Point2::new(0.3, -1.0),
                Point2::new(1.9, 2.0),
            ),
            (
                Surface::EllipticalCylinder(EllipticalCylinder {
                    frame: f,
                    semi_axis_x: 3.0,
                    semi_axis_y: 1.0,
                }),
                Point2::new(0.2, 0.0),
                Point2::new(1.4, 1.0),
            ),
            (
                Surface::Cone(Cone {
                    frame: f,
                    radius: 1.5,
                    semi_angle: 0.4,
                }),
                Point2::new(0.5, 0.2),
                Point2::new(2.5, 1.5),
            ),
            (
                Surface::Sphere(Sphere {
                    frame: f,
                    radius: 2.0,
                }),
                Point2::new(0.4, -0.3),
                Point2::new(2.0, 1.1),
            ),
            (
                Surface::Torus(Torus {
                    frame: f,
                    major_radius: 3.0,
                    minor_radius: 1.0,
                }),
                Point2::new(0.1, 0.5),
                Point2::new(1.7, 2.9),
            ),
        ]
    }

    fn samples(lo: Point2, hi: Point2) -> impl Iterator<Item = Point2> {
        const N: usize = 24;
        (0..=N).flat_map(move |i| {
            (0..=N).map(move |j| {
                Point2::new(
                    lo.x + (hi.x - lo.x) * i as f64 / N as f64,
                    lo.y + (hi.y - lo.y) * j as f64 / N as f64,
                )
            })
        })
    }

    #[test]
    fn every_patch_point_lies_in_its_sphere() {
        for (surface, lo, hi) in families() {
            let (centre, radius) = patch_sphere(&surface, lo, hi).expect("bounded");
            let mut reach: f64 = 0.0;
            for p in samples(lo, hi) {
                let point = evaluate(&surface, p.x, p.y).expect("point");
                reach = reach.max((point - centre).length());
            }
            assert!(
                reach <= radius,
                "{surface:?}: reach {reach} > radius {radius}"
            );
            // Not vacuous: the bound is within a small factor of the reach.
            assert!(
                radius <= 4.0 * reach + 1e-9,
                "{surface:?}: {radius} vs {reach}"
            );
        }
    }

    #[test]
    fn every_patch_normal_lies_in_its_cone() {
        let mid = |lo: Point2, hi: Point2| (lo + hi) * 0.5;
        for (surface, lo, hi) in families() {
            let Some(spread) = normal_spread(&surface, lo, hi) else {
                continue;
            };
            let c = mid(lo, hi);
            let axis = normal(&surface, c.x, c.y).expect("normal");
            let mut widest: f64 = 0.0;
            for p in samples(lo, hi) {
                let n = normal(&surface, p.x, p.y).expect("normal");
                widest = widest.max(n.dot(axis).clamp(-1.0, 1.0).acos());
            }
            assert!(
                widest <= spread,
                "{surface:?}: normals turn {widest} > {spread}"
            );
        }
    }
}
