//! Section edges: where two exact B-reps' faces cross (ADR 0075, step 2).
//!
//! For every pair of faces, one from each operand:
//!
//! 1. The support surfaces' exact intersection curves come from
//!    `exact_surface_intersection` (#119). Stage 1 takes lines, circles and
//!    ellipses; any other section is refused by name.
//! 2. Each curve is cut where it crosses a boundary edge of either face.
//!    Two curves on one surface meet only up to the rounding of their
//!    constructed doubles, so the crossing is found transversally instead:
//!    where the curve crosses the ADJACENT face's surface
//!    (`exact_curve_surface_intersection`), kept when the hit lies on the
//!    edge within tolerance and inside its span. Across a seam, where the
//!    adjacent face shares the surface, the cutter is the plane through the
//!    ruling and the surface normal.
//! 3. Each piece between consecutive cuts lies wholly inside or wholly
//!    outside each face, so its midpoint decides it, through the face's
//!    certified domain (`FaceDomain`). A piece inside both faces is a section
//!    edge.
//!
//! Faces that meet without crossing -- coincident or tangent supports, or a
//! section running along an existing edge -- are refused by name in this
//! stage rather than guessed.

use axiolid_brep::ExactBRep;
use axiolid_core::{Interval, Point2, Point3, Scalar, Tolerance};
use axiolid_curve::Curve3;
use axiolid_evaluate::surface::{invert, normal};
use axiolid_evaluate::{curve::invert3, evaluate3};
use axiolid_measure::FaceDomain;
use axiolid_nurbs::{
    exact_curve_surface_intersection, exact_surface_intersection, ExactCurveIntersection,
    ExactIntersectionRefusal,
};
use axiolid_surface::{Plane, Surface};
use axiolid_topology::FaceId;
use core::f64::consts::TAU;

use crate::BooleanError;

/// A piece of an intersection curve lying on a face of each operand.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionEdge {
    /// The face of the first operand the edge lies on.
    pub face_a: FaceId,
    /// The face of the second operand the edge lies on.
    pub face_b: FaceId,
    /// The exact curve.
    pub curve: Curve3,
    /// The curve's parameter span; for a circle or ellipse it may run past
    /// `2 pi`.
    pub span: Interval,
    /// The curve's point at `span.start`.
    pub start: Point3,
    /// The curve's point at `span.end`.
    pub end: Point3,
}

/// Every section edge between the faces of `a` and the faces of `b`.
///
/// # Errors
///
/// A face pair whose section is not a line, circle or ellipse, faces that
/// touch without crossing, a section running along an existing edge, or a
/// piece too close to a face boundary to classify.
pub fn section_edges(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
) -> Result<Vec<SectionEdge>, BooleanError> {
    let side_a = Side::new(a, tolerance)?;
    let side_b = Side::new(b, tolerance)?;
    let mut out = Vec::new();
    for fa in 0..side_a.faces.len() {
        for fb in 0..side_b.faces.len() {
            let (sa, sb) = (side_a.surface(fa)?, side_b.surface(fb)?);
            let curve = match exact_surface_intersection(sa, sb) {
                Ok(curve) => curve,
                Err(ExactIntersectionRefusal::Disjoint) => continue,
                Err(ExactIntersectionRefusal::NotRegularCurve) => {
                    return Err(BooleanError::NotTransverse {
                        face_a: side_a.faces[fa],
                        face_b: side_b.faces[fb],
                    })
                }
                Err(_) => return Err(BooleanError::UnsupportedSection),
            };
            for (branch, span) in curve.branches.iter().zip(&curve.spans) {
                if span.is_some()
                    || !matches!(
                        branch,
                        Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_)
                    )
                {
                    return Err(BooleanError::UnsupportedSection);
                }
                let mut cuts = side_a.cuts(fa, branch, tolerance)?;
                cuts.extend(side_b.cuts(fb, branch, tolerance)?);
                for span in pieces(branch, cuts) {
                    let mid = evaluate3(branch, 0.5 * (span.start + span.end))
                        .map_err(|_| BooleanError::Evaluation)?;
                    if side_a.inside(fa, mid, tolerance)? && side_b.inside(fb, mid, tolerance)? {
                        out.push(SectionEdge {
                            face_a: side_a.faces[fa],
                            face_b: side_b.faces[fb],
                            curve: branch.clone(),
                            span,
                            start: evaluate3(branch, span.start)
                                .map_err(|_| BooleanError::Evaluation)?,
                            end: evaluate3(branch, span.end)
                                .map_err(|_| BooleanError::Evaluation)?,
                        });
                    }
                }
            }
        }
    }
    Ok(out)
}

/// The spans between consecutive cuts: finite stretches of a line (its
/// unbounded ends leave every bounded face), the cyclic arcs of a conic.
fn pieces(curve: &Curve3, mut cuts: Vec<Scalar>) -> Vec<Interval> {
    cuts.sort_by(Scalar::total_cmp);
    cuts.dedup_by(|x, y| (*x - *y).abs() <= 1e-12 * (1.0 + x.abs()));
    match curve {
        Curve3::Line(_) => cuts
            .windows(2)
            .map(|pair| Interval::new(pair[0], pair[1]))
            .collect(),
        _ => {
            if cuts.is_empty() {
                return vec![Interval::new(0.0, TAU)];
            }
            let mut out: Vec<Interval> = cuts
                .windows(2)
                .map(|pair| Interval::new(pair[0], pair[1]))
                .collect();
            out.push(Interval::new(cuts[cuts.len() - 1], cuts[0] + TAU));
            out
        }
    }
}

/// One operand, with each face's certified domain built once.
struct Side<'a> {
    brep: &'a ExactBRep,
    faces: Vec<FaceId>,
    domains: Vec<FaceDomain<'a>>,
    /// For each edge, the faces whose loops use it.
    edge_faces: Vec<Vec<usize>>,
}

impl<'a> Side<'a> {
    fn new(brep: &'a ExactBRep, tolerance: Tolerance) -> Result<Self, BooleanError> {
        let topology = brep.topology();
        let mut faces = Vec::with_capacity(topology.faces().len());
        let mut domains = Vec::with_capacity(topology.faces().len());
        for index in 0..topology.faces().len() {
            let face = topology
                .face_id_at(index)
                .ok_or(BooleanError::DanglingReference)?;
            let domain = FaceDomain::new(brep, face, tolerance)
                .map_err(BooleanError::Measure)?
                .ok_or(BooleanError::UnsupportedTrim)?;
            faces.push(face);
            domains.push(domain);
        }
        let mut edge_faces = vec![Vec::new(); topology.edges().len()];
        for (index, face) in topology.faces().iter().enumerate() {
            for bound in &face.bounds {
                let wire = topology
                    .loops()
                    .get(bound.loop_id.index())
                    .ok_or(BooleanError::DanglingReference)?;
                for use_ in &wire.edges {
                    edge_faces[use_.edge.index()].push(index);
                }
            }
        }
        Ok(Self {
            brep,
            faces,
            domains,
            edge_faces,
        })
    }

    fn surface(&self, face: usize) -> Result<&'a Surface, BooleanError> {
        self.brep.topology().faces()[face]
            .surface
            .and_then(|id| self.brep.surfaces().get(id.index()))
            .ok_or(BooleanError::DanglingReference)
    }

    /// Whether `point`, on the face's support, lies inside the face.
    fn inside(
        &self,
        face: usize,
        point: Point3,
        tolerance: Tolerance,
    ) -> Result<bool, BooleanError> {
        let (u, v) =
            invert(self.surface(face)?, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
        self.domains[face]
            .contains(Point2::new(u, v))
            .map_err(BooleanError::Measure)?
            .ok_or(BooleanError::Undecided)
    }

    /// Parameters on `curve` where it crosses a boundary edge of the face.
    fn cuts(
        &self,
        face: usize,
        curve: &Curve3,
        tolerance: Tolerance,
    ) -> Result<Vec<Scalar>, BooleanError> {
        let topology = self.brep.topology();
        let mut out = Vec::new();
        let mut seen = Vec::new();
        for bound in &topology.faces()[face].bounds {
            let wire = topology
                .loops()
                .get(bound.loop_id.index())
                .ok_or(BooleanError::DanglingReference)?;
            for use_ in &wire.edges {
                if seen.contains(&use_.edge) {
                    continue;
                }
                seen.push(use_.edge);
                let edge = &topology.edges()[use_.edge.index()];
                let edge_curve = edge
                    .curve
                    .and_then(|id| self.brep.curves3().get(id.index()))
                    .ok_or(BooleanError::DanglingReference)?;
                let span = self
                    .brep
                    .edge_interval(use_.edge)
                    .ok_or(BooleanError::DanglingReference)?;
                let mut cutter =
                    self.cutter(face, use_.edge.index(), edge_curve, span, tolerance, false)?;
                let mut result = exact_curve_surface_intersection(curve, &cutter);
                if matches!(result, Ok(ExactCurveIntersection::Contained)) {
                    // The adjacent face may lie on the same surface under a
                    // different frame (a column's half-walls): cut across
                    // the seam instead.
                    if let Ok(seam) =
                        self.cutter(face, use_.edge.index(), edge_curve, span, tolerance, true)
                    {
                        cutter = seam;
                        result = exact_curve_surface_intersection(curve, &cutter);
                    }
                }
                match result {
                    Ok(ExactCurveIntersection::Points(hits)) => {
                        for hit in hits {
                            if on_edge(edge_curve, span, hit.point, tolerance)? {
                                out.push(hit.parameter.approx());
                            }
                        }
                    }
                    // The section lies in the adjacent face's surface: it
                    // runs along that face, or along the edge itself.
                    Ok(ExactCurveIntersection::Contained) => {
                        return Err(BooleanError::SectionAlongEdge)
                    }
                    Ok(_) => return Err(BooleanError::UnsupportedSection),
                    Err(_) => return Err(BooleanError::UnsupportedTrim),
                }
            }
        }
        Ok(out)
    }
}

impl Side<'_> {
    /// A surface the section crosses exactly where it crosses the edge: the
    /// adjacent face's support, or across a seam (the same support on both
    /// sides) the plane through the ruling and the surface normal.
    fn cutter(
        &self,
        face: usize,
        edge: usize,
        edge_curve: &Curve3,
        span: Interval,
        tolerance: Tolerance,
        across_seam: bool,
    ) -> Result<Surface, BooleanError> {
        let own = self.surface(face)?;
        let adjacent = self.edge_faces[edge]
            .iter()
            .copied()
            .find(|&other| other != face);
        if let (Some(other), false) = (adjacent, across_seam) {
            let theirs = self.surface(other)?;
            if !same_support(own, theirs, tolerance) {
                return Ok(theirs.clone());
            }
        }
        // A seam: the same support on both sides (or a free edge).
        let Curve3::Line(line) = edge_curve else {
            return Err(BooleanError::UnsupportedTrim);
        };
        let mid = evaluate3(edge_curve, 0.5 * (span.start + span.end))
            .map_err(|_| BooleanError::Evaluation)?;
        let (u, v) = invert(own, mid, tolerance).map_err(|_| BooleanError::Evaluation)?;
        let n = normal(own, u, v).map_err(|_| BooleanError::Evaluation)?;
        let across = line.direction.cross(n);
        let length = across.length();
        if length == 0.0 || !length.is_finite() {
            return Err(BooleanError::UnsupportedTrim);
        }
        let z = across / length;
        let x = line.direction.normalize();
        Ok(Surface::Plane(Plane {
            frame: axiolid_core::Frame3 {
                origin: mid,
                x,
                y: z.cross(x),
                z,
            },
        }))
    }
}

/// Whether two supports are the same surface, whatever their frames: a
/// column's half-walls lie on one cylinder parameterised from opposite
/// sides, and the edge between them is a seam, not a crease.
fn same_support(a: &Surface, b: &Surface, tolerance: Tolerance) -> bool {
    let eps = tolerance.linear().max(1e-9);
    let parallel = |x: axiolid_core::Vec3, y: axiolid_core::Vec3| {
        x.normalize().cross(y.normalize()).length() <= 1e-9
    };
    let on_axis = |o1: Point3, o2: Point3, z: axiolid_core::Vec3| {
        (o2 - o1).cross(z.normalize()).length() <= eps
    };
    match (a, b) {
        (Surface::Plane(p), Surface::Plane(q)) => {
            parallel(p.frame.z, q.frame.z)
                && (q.frame.origin - p.frame.origin)
                    .dot(p.frame.z.normalize())
                    .abs()
                    <= eps
        }
        (Surface::Cylinder(p), Surface::Cylinder(q)) => {
            (p.radius - q.radius).abs() <= eps
                && parallel(p.frame.z, q.frame.z)
                && on_axis(p.frame.origin, q.frame.origin, p.frame.z)
        }
        (Surface::Sphere(p), Surface::Sphere(q)) => {
            (p.radius - q.radius).abs() <= eps && (p.frame.origin - q.frame.origin).length() <= eps
        }
        _ => a == b,
    }
}

/// Whether `point` lies on the edge within tolerance and inside its span.
fn on_edge(
    curve: &Curve3,
    span: Interval,
    point: Point3,
    tolerance: Tolerance,
) -> Result<bool, BooleanError> {
    let Ok(t) = invert3(curve, point, tolerance) else {
        return Ok(false);
    };
    let on = evaluate3(curve, t).map_err(|_| BooleanError::Evaluation)?;
    if (on - point).length() > tolerance.linear().max(1e-9) {
        return Ok(false);
    }
    on_span(curve, span, point, tolerance)
}

/// Whether `point`, on `curve`, lies within the edge's span (a closed
/// conic's span may start anywhere and run past a full turn).
fn on_span(
    curve: &Curve3,
    span: Interval,
    point: Point3,
    tolerance: Tolerance,
) -> Result<bool, BooleanError> {
    let t = invert3(curve, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
    let (lo, hi) = (span.start.min(span.end), span.start.max(span.end));
    let slack = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
    let periodic = matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_));
    let candidates: &[Scalar] = if periodic {
        &[t - TAU, t, t + TAU, t + 2.0 * TAU]
    } else {
        &[t]
    };
    Ok(candidates
        .iter()
        .any(|c| *c >= lo - slack && *c <= hi + slack))
}
