//! Section edges: where two exact B-reps' faces cross (ADR 0075, step 2).
//!
//! For every pair of faces, one from each operand:
//!
//! 1. The support surfaces' exact intersection curves come from
//!    `exact_surface_intersection` (#119). Stage 1 takes lines, circles and
//!    ellipses; any other section is refused by name.
//! 2. Each curve is cut at every point where it crosses a boundary edge of
//!    either face, from `exact_curve_curve_intersection3`: exact hits on the
//!    curve, kept when they fall inside the edge's own span.
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
use axiolid_evaluate::surface::invert;
use axiolid_evaluate::{curve::invert3, evaluate3};
use axiolid_measure::FaceDomain;
use axiolid_nurbs::{
    exact_curve_curve_intersection3, exact_surface_intersection, ExactCurveIntersection,
    ExactIntersectionRefusal,
};
use axiolid_surface::Surface;
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
        Ok(Self {
            brep,
            faces,
            domains,
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
                match exact_curve_curve_intersection3(curve, edge_curve) {
                    Ok(ExactCurveIntersection::Points(hits)) => {
                        for hit in hits {
                            if on_span(edge_curve, span, hit.point, tolerance)? {
                                out.push(hit.parameter.approx());
                            }
                        }
                    }
                    // The section runs along the edge itself.
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
