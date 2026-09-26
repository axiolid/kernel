//! Section edges: where two exact B-reps' faces cross (ADR 0075, step 2).
//!
//! For every pair of faces, one from each operand:
//!
//! 1. The support surfaces' exact intersection curves come from
//!    `exact_surface_intersection` (#119): lines and conics, ruled and
//!    torus sections (ADR 0076), traced sections (ADR 0077). A ruled or
//!    torus section piece is read as a traced section over its span, so
//!    every curve past this point is a line, a conic or an
//!    `ImplicitSection3`.
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
use axiolid_evaluate::surface::{locate, normal};
use axiolid_evaluate::{curve::locate3, evaluate3};
use axiolid_measure::FaceDomain;
use axiolid_nurbs::{
    exact_curve_curve_intersection3, exact_curve_surface_intersection, exact_surface_intersection,
    implicit_surface_intersection, ExactCurveIntersection, ExactIntersectionRefusal,
};
use axiolid_surface::{Plane, Surface};
use axiolid_topology::FaceId;
use core::f64::consts::TAU;

use crate::support::{cleaned, same_support, window};
use crate::BooleanError;

/// Where along a piece its inside/outside sample is taken: off-centre, so
/// the midpoint of a symmetric section (a meridian's pole) is never it.
const SAMPLE: Scalar = 0.414_213_562_373_095;

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
    /// The surface that meets `face_a`'s surface along this edge: `face_b`'s
    /// for a crossing, or for an edge imprinted from a coincident face the
    /// surface that bounds it there. Its equation on `face_a` is the edge's
    /// pcurve.
    pub other_a: Surface,
    /// The surface that meets `face_b`'s surface along this edge.
    pub other_b: Surface,
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
            // Lines and conics in closed form; every other section traced in
            // one face's parameter box (ADR 0077), which holds every part
            // of it that can matter.
            let closed_form = match exact_surface_intersection(&cleaned(sa), &cleaned(sb)) {
                Ok(curve) => {
                    let conic = curve.branches.iter().zip(&curve.spans).all(|(b, s)| {
                        matches!(
                            (b, s),
                            (Curve3::Line(_), _) | (Curve3::Circle(_) | Curve3::Ellipse(_), None)
                        )
                    });
                    conic.then_some(curve)
                }
                // Apart, or touching without crossing: no section.
                Err(
                    ExactIntersectionRefusal::Disjoint | ExactIntersectionRefusal::NotRegularCurve,
                ) => continue,
                Err(_) => None,
            };
            let branches: Vec<(Curve3, Option<Interval>)> = match closed_form {
                Some(curve) => curve.branches.into_iter().zip(curve.spans).collect(),
                None => {
                    // Only a B-spline can carry a section with a B-spline
                    // (it has no equation to read elsewhere); otherwise the
                    // compact and low-degree surfaces carry best.
                    let rank = |s: &Surface| match s {
                        Surface::BSpline(_) => 4,
                        Surface::Torus(_) => 3,
                        Surface::Sphere(_) => 2,
                        Surface::Plane(_) => 0,
                        _ => 1,
                    };
                    if matches!((sa, sb), (Surface::BSpline(_), Surface::BSpline(_))) {
                        return Err(BooleanError::UnsupportedSection);
                    }
                    let (carrier, other, side, face) = if rank(sa) >= rank(sb) {
                        (sa, sb, &side_a, fa)
                    } else {
                        (sb, sa, &side_b, fb)
                    };
                    let (lo, hi) = side.domains[face].bounds();
                    match implicit_surface_intersection(
                        carrier,
                        other,
                        Some(window(carrier, lo, hi)),
                    ) {
                        Ok(sections) => sections
                            .into_iter()
                            .map(|s| {
                                let end = s.curve.end();
                                (Curve3::ImplicitSection(s), Some(Interval::new(0.0, end)))
                            })
                            .collect(),
                        Err(ExactIntersectionRefusal::Disjoint) => continue,
                        Err(_) => return Err(BooleanError::UnsupportedSection),
                    }
                }
            };
            for (index, (branch, span)) in branches.iter().enumerate() {
                // A ray (a cone's ruling from its apex) is a line bounded
                // at its finite ends.
                let bounds = match (branch, span) {
                    (Curve3::Line(_), Some(span)) => Some(*span),
                    _ => None,
                };
                let branch = branch.clone();
                let branch = &branch;
                let (mut cuts, along_a) = side_a.cuts(fa, branch, sb, tolerance)?;
                let (more, along_b) = side_b.cuts(fb, branch, sa, tolerance)?;
                cuts.extend(more);
                // Branches of one section meet only where the surfaces touch
                // (the two ellipses of a Steinmetz pair): such a point splits
                // both, so the faces' graphs get a vertex there.
                for (other, _) in branches
                    .iter()
                    .skip(index + 1)
                    .chain(branches.iter().take(index))
                {
                    if let Ok(ExactCurveIntersection::Points(hits)) =
                        exact_curve_curve_intersection3(branch, other)
                    {
                        cuts.extend(hits.iter().map(|h| h.parameter.approx()));
                    }
                }
                if let Some(b) = bounds {
                    let (lo, hi) = (b.start.min(b.end), b.start.max(b.end));
                    cuts.retain(|&c| c >= lo && c <= hi);
                    cuts.extend([lo, hi].into_iter().filter(|x| x.is_finite()));
                }
                for (piece_curve, span) in pieces(branch, cuts)? {
                    // Off-centre, so a symmetric section's pole or seam
                    // crossing never becomes the sample.
                    let mid =
                        evaluate3(&piece_curve, span.start + SAMPLE * (span.end - span.start))
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
                        start: evaluate3(&piece_curve, span.start)
                            .map_err(|_| BooleanError::Evaluation)?,
                        end: evaluate3(&piece_curve, span.end)
                            .map_err(|_| BooleanError::Evaluation)?,
                        curve: piece_curve,
                        span,
                        along_a: at_a == Place::Boundary,
                        along_b: at_b == Place::Boundary,
                        other_a: sb.clone(),
                        other_b: sa.clone(),
                    });
                }
            }
        }
    }
    Ok(out)
}

/// The points where a surface's angle parameter has no value: a sphere's
/// poles, a cone's apex.
fn poles(surface: &Surface) -> Vec<Point3> {
    match surface {
        Surface::Sphere(s) => {
            let z = s.frame.z.normalize() * s.radius;
            vec![s.frame.origin + z, s.frame.origin - z]
        }
        Surface::Cone(c) => {
            let slope = c.semi_angle.tan();
            if slope == 0.0 {
                Vec::new()
            } else {
                vec![c.frame.origin - c.frame.z.normalize() * (c.radius / slope)]
            }
        }
        _ => Vec::new(),
    }
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
        let (u, v) = locate(s, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
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
    for (edge, curve, span) in from.edges_of(face)? {
        // The surface that bounds the imprinted edge in its own operand:
        // together with the shared surface it defines the edge's curve.
        let bounding = from.cutter(face, edge, &curve, span, tolerance, false)?;
        let (cuts, along) = onto.cuts(other, &curve, &bounding, tolerance)?;
        for piece in pieces_within(&curve, span, cuts) {
            let mid = evaluate3(&curve, piece.start + SAMPLE * (piece.end - piece.start))
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
                other_a: bounding.clone(),
                other_b: bounding.clone(),
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

/// The stretches between consecutive cuts, each with the curve it lies on:
/// finite stretches of a line (its unbounded ends leave every bounded
/// face), the cyclic arcs of a conic, and for a traced section the
/// stretches of its span. A traced loop's stretch across the loop's start
/// becomes a curve of its own.
fn pieces(curve: &Curve3, mut cuts: Vec<Scalar>) -> Result<Vec<(Curve3, Interval)>, BooleanError> {
    cuts.sort_by(Scalar::total_cmp);
    cuts.dedup_by(|x, y| (*x - *y).abs() <= 1e-12 * (1.0 + x.abs()));
    let plain = |spans: Vec<Interval>| spans.into_iter().map(|s| (curve.clone(), s)).collect();
    Ok(match curve {
        Curve3::Line(_) => plain(
            cuts.windows(2)
                .map(|pair| Interval::new(pair[0], pair[1]))
                .collect(),
        ),
        Curve3::ImplicitSection(section) => {
            let n = section.curve.end();
            let (pu, pv) = section.carrier.periodic();
            let closure = section.curve.closure(pu, pv);
            let slack = 1e-9 * (1.0 + n);
            let mut inner: Vec<Scalar> = cuts
                .into_iter()
                .filter(|&c| c > slack && c < n - slack)
                .collect();
            inner.dedup_by(|x, y| (*x - *y).abs() <= slack);
            let own = |a: Scalar, b: Scalar| -> Result<(Curve3, Interval), BooleanError> {
                // One cut on a loop: the whole loop, starting there.
                let sub = if (a - b).abs() <= slack {
                    closure.and_then(|c| section.curve.rotated(a, c))
                } else {
                    section.curve.sub(a, b, closure)
                };
                let sub = sub.ok_or(BooleanError::Evaluation)?;
                let end = sub.end();
                Ok((
                    Curve3::ImplicitSection(axiolid_curve::ImplicitSection3 {
                        carrier: section.carrier.clone(),
                        curve: sub,
                    }),
                    Interval::new(0.0, end),
                ))
            };
            if closure.is_some() {
                if inner.is_empty() {
                    return Ok(vec![(curve.clone(), Interval::new(0.0, n))]);
                }
                let mut out = Vec::new();
                for pair in inner.windows(2) {
                    out.push(own(pair[0], pair[1])?);
                }
                // Across the loop's own start.
                out.push(own(inner[inner.len() - 1], inner[0])?);
                out
            } else {
                let mut ends = vec![0.0];
                ends.extend(inner);
                ends.push(n);
                let mut out = Vec::new();
                for pair in ends.windows(2) {
                    out.push(own(pair[0], pair[1])?);
                }
                out
            }
        }
        _ => {
            if cuts.is_empty() {
                return Ok(vec![(curve.clone(), Interval::new(0.0, TAU))]);
            }
            let mut out: Vec<Interval> = cuts
                .windows(2)
                .map(|pair| Interval::new(pair[0], pair[1]))
                .collect();
            out.push(Interval::new(cuts[cuts.len() - 1], cuts[0] + TAU));
            plain(out)
        }
    })
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
            locate(self.surface(face)?, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
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
    fn edges_of(&self, face: usize) -> Result<Vec<(usize, Curve3, Interval)>, BooleanError> {
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
                out.push((use_.edge.index(), curve.clone(), span));
            }
        }
        Ok(out)
    }

    /// Parameters on `curve` where it crosses a boundary edge of the face,
    /// and the edges it runs along (lying in the adjacent face's surface
    /// too), which cut it at their ends.
    #[allow(clippy::type_complexity)]
    ///
    /// `meets` is the other surface `curve` lies on. Where the curve cannot
    /// be intersected with the adjacent face's surface (a B-spline), the
    /// edge is intersected with `meets` instead: the curve crosses the edge
    /// exactly where the edge crosses `meets`.
    fn cuts(
        &self,
        face: usize,
        curve: &Curve3,
        meets: &Surface,
        tolerance: Tolerance,
    ) -> Result<(Vec<Scalar>, Vec<(Curve3, Interval)>), BooleanError> {
        let topology = self.brep.topology();
        let mut out = Vec::new();
        let mut along = Vec::new();
        let mut seen = Vec::new();
        // A pole or apex the curve passes through cuts it: its pcurve jumps
        // round the angle there.
        for pole in poles(self.surface(face)?) {
            if let Ok(t) = locate3(curve, pole, tolerance) {
                let on = evaluate3(curve, t).map_err(|_| BooleanError::Evaluation)?;
                if (on - pole).length() <= tolerance.linear().max(1e-9) {
                    out.push(t);
                }
            }
        }
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
                            if let Ok(s) = locate3(curve, end, tolerance) {
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
                    Err(_) => {
                        let hits = match exact_curve_surface_intersection(edge_curve, meets) {
                            Ok(ExactCurveIntersection::Points(hits)) => hits,
                            _ => return Err(BooleanError::UnsupportedTrim),
                        };
                        for hit in hits {
                            if !on_span(edge_curve, span, hit.point, tolerance)? {
                                continue;
                            }
                            if let Ok(s) = locate3(curve, hit.point, tolerance) {
                                let on =
                                    evaluate3(curve, s).map_err(|_| BooleanError::Evaluation)?;
                                if (on - hit.point).length() <= tolerance.linear().max(1e-9) {
                                    out.push(s);
                                }
                            }
                        }
                    }
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
        // A seam: the same support on both sides (or a free edge). A seam
        // circle (a sphere's meridian or latitude, a torus's tube or ring
        // circle) is cut by the surface its normals sweep: coaxial with the
        // circle, a plane, cylinder or cone, which crosses the face
        // transversally along the whole circle.
        if let Curve3::Circle(circle) = edge_curve {
            return normal_sweep(own, circle, tolerance);
        }
        let Curve3::Line(line) = edge_curve else {
            return Err(BooleanError::UnsupportedTrim);
        };
        let mid = evaluate3(edge_curve, 0.5 * (span.start + span.end))
            .map_err(|_| BooleanError::Evaluation)?;
        let (u, v) = locate(own, mid, tolerance).map_err(|_| BooleanError::Evaluation)?;
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

/// The surface swept by `surface`'s normal lines along a circle on it,
/// when the surface is one of revolution about the circle's axis (or the
/// circle is a meridian, whose normals stay in its plane): a plane, a
/// cylinder or a cone coaxial with the circle.
fn normal_sweep(
    surface: &Surface,
    circle: &axiolid_curve::Circle3,
    tolerance: Tolerance,
) -> Result<Surface, BooleanError> {
    let f = circle.frame;
    let axis = f.x.cross(f.y).normalize();
    let x = f.x.normalize();
    let y = axis.cross(x);
    let p = f.origin + x * circle.radius;
    let (u, v) = locate(surface, p, tolerance).map_err(|_| BooleanError::Evaluation)?;
    let n = normal(surface, u, v)
        .map_err(|_| BooleanError::Evaluation)?
        .normalize();
    let frame = axiolid_core::Frame3 {
        origin: f.origin,
        x,
        y,
        z: axis,
    };
    let along = n.dot(axis);
    let radial = n.dot(x);
    let plane = || {
        Surface::Plane(Plane {
            frame: axiolid_core::Frame3 {
                origin: f.origin,
                x,
                y,
                z: axis,
            },
        })
    };
    // Normals within the circle's plane: the plane itself.
    if along.abs() <= 1e-12 {
        return Ok(plane());
    }
    // Normals along the axis: the cylinder through the circle.
    if radial.abs() <= 1e-12 {
        return Ok(Surface::Cylinder(axiolid_surface::Cylinder {
            frame,
            radius: circle.radius,
        }));
    }
    // Otherwise the cone of normal lines: at height `h` along the axis the
    // line is `radius + h * radial / along` from it.
    Ok(Surface::Cone(axiolid_surface::Cone {
        frame,
        radius: circle.radius,
        semi_angle: (radial / along).atan(),
    }))
}

/// Whether `point` lies on the edge within tolerance and inside its span.
fn on_edge(
    curve: &Curve3,
    span: Interval,
    point: Point3,
    tolerance: Tolerance,
) -> Result<bool, BooleanError> {
    let Ok(t) = locate3(curve, point, tolerance) else {
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
    let t = locate3(curve, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
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
