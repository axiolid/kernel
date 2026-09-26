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
//! Faces that meet without crossing are handled, not refused:
//!
//! - **Tangent supports** touch in a point or along a curve without
//!   crossing, so they add no section: neither face changes sides there.
//! - **A section along an existing edge** (the curve lies in the adjacent
//!   face's surface too) is cut at the edge's ends and marked as lying
//!   along that face's boundary: it splits the other face, not this one.
//! - **Coincident supports** share a patch of surface. Each face's boundary
//!   edges are imprinted on the other where they run inside it, so the
//!   shared patch becomes a region of both faces; classification then sees
//!   it on the other solid's boundary (see `crate::boolean`).

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

use crate::support::same_support;
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
    /// Whether the edge runs along a boundary edge `face_a` already has, so
    /// it does not split `face_a`.
    pub along_a: bool,
    /// Whether the edge runs along a boundary edge `face_b` already has.
    pub along_b: bool,
}

/// Every section edge between the faces of `a` and the faces of `b`.
///
/// # Errors
///
/// A face pair whose section is not a line, circle or ellipse, or a piece
/// too close to a face boundary to classify.
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
            if same_support(sa, sb, tolerance) {
                imprint(&side_a, fa, &side_b, fb, true, tolerance, &mut out)?;
                imprint(&side_b, fb, &side_a, fa, false, tolerance, &mut out)?;
                continue;
            }
            let curve = match exact_surface_intersection(sa, sb) {
                Ok(curve) => curve,
                // Apart, or touching without crossing: no section.
                Err(
                    ExactIntersectionRefusal::Disjoint | ExactIntersectionRefusal::NotRegularCurve,
                ) => continue,
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
                let (mut cuts, along_a) = side_a.cuts(fa, branch, tolerance)?;
                let (more, along_b) = side_b.cuts(fb, branch, tolerance)?;
                cuts.extend(more);
                for span in pieces(branch, cuts) {
                    let mid = evaluate3(branch, 0.5 * (span.start + span.end))
                        .map_err(|_| BooleanError::Evaluation)?;
                    if touching(sa, sb, mid, tolerance)? {
                        continue;
                    }
                    let at_a = side_a.locate(fa, mid, &along_a, tolerance)?;
                    if at_a == Place::Outside {
                        continue;
                    }
                    let at_b = side_b.locate(fb, mid, &along_b, tolerance)?;
                    if at_b == Place::Outside {
                        continue;
                    }
                    out.push(SectionEdge {
                        face_a: side_a.faces[fa],
                        face_b: side_b.faces[fb],
                        curve: branch.clone(),
                        span,
                        start: evaluate3(branch, span.start)
                            .map_err(|_| BooleanError::Evaluation)?,
                        end: evaluate3(branch, span.end).map_err(|_| BooleanError::Evaluation)?,
                        along_a: at_a == Place::Boundary,
                        along_b: at_b == Place::Boundary,
                    });
                }
            }
        }
    }
    Ok(out)
}

/// Whether two surfaces through `point` share their tangent plane there: a
/// line where a cylinder rests on a plane, which it touches without
/// crossing. Analytic surfaces tangent along a whole curve lie on one side
/// of each other, so such a stretch splits neither face.
fn touching(
    a: &Surface,
    b: &Surface,
    point: Point3,
    tolerance: Tolerance,
) -> Result<bool, BooleanError> {
    let at = |s: &Surface| -> Result<axiolid_core::Vec3, BooleanError> {
        let (u, v) = invert(s, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
        Ok(normal(s, u, v)
            .map_err(|_| BooleanError::Evaluation)?
            .normalize())
    };
    Ok(at(a)?.cross(at(b)?).length() <= 1e-9)
}

/// The boundary edges of `from`'s face `face` that run inside `onto`'s face
/// `other`, on the same surface, as section edges (`first` when `from` is
/// the first operand).
#[allow(clippy::too_many_arguments)]
fn imprint(
    from: &Side<'_>,
    face: usize,
    onto: &Side<'_>,
    other: usize,
    first: bool,
    tolerance: Tolerance,
    out: &mut Vec<SectionEdge>,
) -> Result<(), BooleanError> {
    for (curve, span) in from.edges_of(face)? {
        let (cuts, along) = onto.cuts(other, &curve, tolerance)?;
        for piece in pieces_within(&curve, span, cuts) {
            let mid = evaluate3(&curve, 0.5 * (piece.start + piece.end))
                .map_err(|_| BooleanError::Evaluation)?;
            let at = onto.locate(other, mid, &along, tolerance)?;
            if at == Place::Outside {
                continue;
            }
            let (face_a, face_b, along_a, along_b) = if first {
                (
                    from.faces[face],
                    onto.faces[other],
                    true,
                    at == Place::Boundary,
                )
            } else {
                (
                    onto.faces[other],
                    from.faces[face],
                    at == Place::Boundary,
                    true,
                )
            };
            out.push(SectionEdge {
                face_a,
                face_b,
                curve: curve.clone(),
                span: piece,
                start: evaluate3(&curve, piece.start).map_err(|_| BooleanError::Evaluation)?,
                end: evaluate3(&curve, piece.end).map_err(|_| BooleanError::Evaluation)?,
                along_a,
                along_b,
            });
        }
    }
    Ok(())
}

/// Where a point on a face's support lies relative to the face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Place {
    Inside,
    Boundary,
    Outside,
}

/// The pieces of a bounded edge span between the cuts inside it.
fn pieces_within(curve: &Curve3, span: Interval, cuts: Vec<Scalar>) -> Vec<Interval> {
    let (lo, hi) = (span.start.min(span.end), span.start.max(span.end));
    let slack = 1e-12 * (1.0 + lo.abs().max(hi.abs()));
    let periodic = matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_));
    let mut inside = vec![lo, hi];
    for cut in cuts {
        let candidates: &[Scalar] = if periodic {
            &[cut - TAU, cut, cut + TAU, cut + 2.0 * TAU]
        } else {
            &[cut]
        };
        for &c in candidates {
            if c > lo + slack && c < hi - slack {
                inside.push(c);
            }
        }
    }
    inside.sort_by(Scalar::total_cmp);
    inside.dedup_by(|x, y| (*x - *y).abs() <= 1e-12 * (1.0 + x.abs()));
    inside
        .windows(2)
        .map(|pair| Interval::new(pair[0], pair[1]))
        .collect()
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

    /// Where `point`, on the face's support, lies: on one of the edges the
    /// curve through it runs along (`along`), or inside or outside the face.
    fn locate(
        &self,
        face: usize,
        point: Point3,
        along: &[(Curve3, Interval)],
        tolerance: Tolerance,
    ) -> Result<Place, BooleanError> {
        for (curve, span) in along {
            if on_edge(curve, *span, point, tolerance)? {
                return Ok(Place::Boundary);
            }
        }
        let (u, v) =
            invert(self.surface(face)?, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
        match self.domains[face]
            .contains(Point2::new(u, v))
            .map_err(BooleanError::Measure)?
        {
            Some(true) => Ok(Place::Inside),
            Some(false) => Ok(Place::Outside),
            None => Err(BooleanError::Undecided),
        }
    }

    /// The face's boundary edges, each once, with their curves and spans.
    fn edges_of(&self, face: usize) -> Result<Vec<(Curve3, Interval)>, BooleanError> {
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
                let curve = topology.edges()[use_.edge.index()]
                    .curve
                    .and_then(|id| self.brep.curves3().get(id.index()))
                    .ok_or(BooleanError::DanglingReference)?;
                let span = self
                    .brep
                    .edge_interval(use_.edge)
                    .ok_or(BooleanError::DanglingReference)?;
                out.push((curve.clone(), span));
            }
        }
        Ok(out)
    }

    /// Parameters on `curve` where it crosses a boundary edge of the face,
    /// and the edges it runs along (lying in the adjacent face's surface
    /// too), which cut it at their ends.
    #[allow(clippy::type_complexity)]
    fn cuts(
        &self,
        face: usize,
        curve: &Curve3,
        tolerance: Tolerance,
    ) -> Result<(Vec<Scalar>, Vec<(Curve3, Interval)>), BooleanError> {
        let topology = self.brep.topology();
        let mut out = Vec::new();
        let mut along = Vec::new();
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
                    // The curve lies in the adjacent face's surface as well
                    // as this one's: where it meets the edge it runs along
                    // it, so it is cut at the edge's ends.
                    Ok(ExactCurveIntersection::Contained) => {
                        for t in [span.start, span.end] {
                            let end =
                                evaluate3(edge_curve, t).map_err(|_| BooleanError::Evaluation)?;
                            if let Ok(s) = invert3(curve, end, tolerance) {
                                let on =
                                    evaluate3(curve, s).map_err(|_| BooleanError::Evaluation)?;
                                if (on - end).length() <= tolerance.linear().max(1e-9) {
                                    out.push(s);
                                }
                            }
                        }
                        along.push((edge_curve.clone(), span));
                    }
                    Ok(_) => return Err(BooleanError::UnsupportedSection),
                    Err(_) => return Err(BooleanError::UnsupportedTrim),
                }
            }
        }
        Ok((out, along))
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
