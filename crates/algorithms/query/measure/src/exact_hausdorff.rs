//! Certified Hausdorff distance between the boundaries of two exact B-reps
//! (#224).
//!
//! # What is certified
//!
//! [`one_sided_boundary_hausdorff`] returns an interval `[lower, upper]`
//! certain to contain `h(A, B) = max over p in dA of d(p, dB)`, how far the
//! boundary of `A` strays from the boundary of `B`, and
//! [`boundary_hausdorff_distance`] the two-sided `max(h(A, B), h(B, A))`
//! with both one-sided intervals: the contract of
//! [`crate::hausdorff_distance`] for meshes, on the exact boundaries instead
//! of a tessellation of them.
//! Both ends are guarantees:
//!
//! - `lower` is `d(p, dB)` bounded below for a point `p` certainly on `dA`
//!   (an edge point, or a surface point at parameters its face's domain is
//!   certified to contain), by the same branch and bound as
//!   [`crate::boundary_distance`] with one side shrunk to the point.
//! - `upper` bounds `d(p, dB)` for every point of `dA`. The faces of `A` are
//!   covered by patches of their parameter domains (a patch is dropped only
//!   when certified outside its face), and each patch's bound is the better
//!   of two:
//!   - a *matched* bound: when a face of `B` has the same surface family
//!     (for a cone the same radius and angle, for a B-spline the same knots,
//!     degrees and weights) and the very same trimming pcurves with the same
//!     intervals, both faces are one domain `D` mapped by `S_A` and `S_B`,
//!     so every `S_A(u, v)` of the patch has
//!     the boundary point `S_B(u, v)` of `B` within `|S_A - S_B|`. Written
//!     over the family's basis, `S_A - S_B = (o_A - o_B) + sum phi_k (a_k^A -
//!     a_k^B)`, bounded over the patch by its value at the centre plus each
//!     basis function's oscillation times its coefficient difference. A
//!     rigid move of a solid changes only frames, so for a pure translation
//!     the bound is the translation itself on every patch at once, and for
//!     an identical copy zero up to rounding. For a rotation it is the
//!     displacement of the rotation, which may exceed the distance: it then
//!     only caps what the Lipschitz bound must refine;
//!   - the same bound for a face of `B` that is a translate of `A`'s but
//!     trimmed in another chart (#227), as a builder that trims in world
//!     coordinates makes it: same family, axes and shape, and trims equal
//!     after the parameter shift the translation induces, to within a
//!     measured residue. `S_B` re-charted by the shift is again one domain
//!     with `S_A`, the residue is folded in, and the bound is `|t|`. A
//!     turned or resized face is never matched this way (the derivation and
//!     the gate are in the private `exact_hausdorff/translate.rs`);
//!   - a *Lipschitz* bound: `d(., dB)` is 1-Lipschitz, so over a patch inside
//!     a sphere of radius `r` about `c` it is at most `d(c, dB) + r`, with
//!     `d(c, dB)` bounded above by a point found on `dB`.
//!
//! The witnesses are `point_from`, the point realising `lower`, and
//! `point_to`, the nearest point found to it on the other boundary.
//!
//! # Method
//!
//! Branch and bound over `A`'s face patches, largest upper bound first: the
//! patch is split across its longer side, and each half is bounded and, if
//! its centre is certainly on the face, measured to raise `lower`. Edge
//! spans of `A` are split alongside, only to place witnesses (the faces'
//! closures already hold every edge point), so a farthest point on an edge
//! is approached along the edge. Point queries walk a tree of `B`'s
//! elements that is refined once and shared by every query.
//!
//! Before any split, each translation `t` a matched face carries seeds the
//! lower bound with `A`'s support point against `t` on its edges: for a
//! translate that point is exactly `|t|` from `B`'s boundary, whatever the
//! shape, so a farthest point isolated at a corner is not left to edge
//! bisection. For a body of several items ([`crate::exact_bodies`]) each
//! item's own support point is seeded: an item moved alone has its
//! farthest point there, wherever the union's support point lies.
//!
//! # Convergence and refusal
//!
//! Copies offset by a translation (and identical copies) close in a few
//! dozen splits whatever the accuracy, and usually in none: the matched
//! bound holds every patch at `|t|` at once, and the support point is `|t|`
//! from the other boundary when the solid's support against `t` is on an
//! edge (faces on planes, cylinders and cones). Elsewhere the Lipschitz bound closes at
//! first order: an accuracy `e` needs patches of radius about `e` wherever
//! the distance comes within `e` of the farthest, so the work grows with
//! that area over `e^2`. Refinement stops at a fixed budget, or at the
//! caller's ([`one_sided_boundary_hausdorff_with_budget`]); the interval
//! returned is then still sound, only wider than asked -- a caller decides
//! by its width, never by a guess. A face whose domain cannot be bounded is
//! refused as by [`crate::boundary_distance`].
//!
//! # Scope
//!
//! This is the distance between BOUNDARIES. A solid wholly inside another
//! is as far from it as its boundary is from the other's boundary, not
//! zero; containment is a separate classification.

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use axiolid_brep::ExactBRep;
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_surface::Surface;
use axiolid_topology::Face;

use crate::exact::ExactMeasureError;
use crate::exact_distance::{
    critical_toward, point_lower_bound, surface_of, Element, Metric, Shape, Side,
};
use crate::mesh_hausdorff::HausdorffBounds;

mod translate;

use axiolid_evaluate::surface::evaluate;
use translate::{support_point, translate, Shifted};

/// Splits of the measured boundary's patches and edge spans before a
/// one-sided query reports what it has.
pub(crate) const MAX_SPLITS: usize = 200_000;

/// Refinement steps one point query may take.
const POINT_STEPS: usize = 20_000;

/// Translations of matched faces whose support points seed the lower
/// bound; a translate gives one.
const MAX_DISPLACEMENTS: usize = 16;

/// The two-sided Hausdorff distance between two exact boundaries and both
/// one-sided ones.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryHausdorff {
    /// `max(h(A, B), h(B, A))`. The witness is that of the side with the
    /// larger lower bound; `point_from` lies on that side's boundary.
    pub distance: HausdorffBounds,
    /// `h(A, B)`: how far `A`'s boundary strays from `B`'s. Witness from
    /// `A` to `B`.
    pub forward: HausdorffBounds,
    /// `h(B, A)`: how far `B`'s boundary strays from `A`'s. Witness from
    /// `B` to `A`.
    pub backward: HausdorffBounds,
}

/// One-sided Hausdorff distance `h(from, to)` between two exact boundaries,
/// to within `accuracy` (#224).
///
/// Refines until `upper - lower <= accuracy` or the budget runs out; the
/// interval is sound either way. A negative or NaN accuracy is read as
/// zero, as by [`crate::boundary_distance`]: refine as far as the budget
/// allows.
///
/// # Errors
///
/// As [`crate::boundary_distance`]: a face whose domain cannot be bounded,
/// an evaluation failure, or a boundary with no boundable element; and
/// [`ExactMeasureError::NonPlanarFace`] (not converged) when no point of
/// `from` could be certified on its boundary.
pub fn one_sided_boundary_hausdorff(
    from: &ExactBRep,
    to: &ExactBRep,
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<HausdorffBounds, ExactMeasureError> {
    one_sided_boundary_hausdorff_with_budget(from, to, accuracy, tolerance, MAX_SPLITS)
}

/// [`one_sided_boundary_hausdorff`] with the caller's refinement budget: at
/// most `max_splits` splits of `from`'s face patches and edge spans (#227).
///
/// The interval is sound whatever the budget; a small one only leaves it
/// wider than `accuracy`. Translated and identical copies close in a few
/// dozen splits, so a caller comparing versions of a model can cap the work
/// per pair and read an interval left open as "not evaluated".
///
/// # Errors
///
/// As [`one_sided_boundary_hausdorff`].
pub fn one_sided_boundary_hausdorff_with_budget(
    from: &ExactBRep,
    to: &ExactBRep,
    accuracy: Scalar,
    tolerance: Tolerance,
    max_splits: usize,
) -> Result<HausdorffBounds, ExactMeasureError> {
    // One item: every edge.
    let every = 0..from.topology().edges().len();
    let edges = core::slice::from_ref(&every);
    Ok(witnessed(from, to, accuracy, tolerance, max_splits, edges)?.bounds)
}

/// A one-sided interval and the elements its witnesses lie on.
pub(crate) struct Witnessed {
    pub(crate) bounds: HausdorffBounds,
    /// The element of `from` holding `point_from`.
    pub(crate) from: Shape,
    /// The element of `to` holding `point_to`.
    pub(crate) to: Shape,
}

/// [`one_sided_boundary_hausdorff_with_budget`], naming the elements the
/// witnesses lie on. `items` partitions `from`'s edges by item: each item's
/// support point against a matched translation seeds the lower bound, so
/// an item moved alone is seeded even where another holds the union's.
pub(crate) fn witnessed(
    from: &ExactBRep,
    to: &ExactBRep,
    accuracy: Scalar,
    tolerance: Tolerance,
    max_splits: usize,
    items: &[core::ops::Range<usize>],
) -> Result<Witnessed, ExactMeasureError> {
    let accuracy = accuracy.max(0.0);
    let linear = tolerance.linear().max(1e-12);
    let source = Side::new(from, linear, Metric::Space)?;
    let mut target = Nearest::new(Side::new(to, linear, Metric::Space)?);
    search(&source, &mut target, accuracy, max_splits, items)
}

/// Two-sided Hausdorff distance between the boundaries of `a` and `b`, to
/// within `accuracy` (#224).
///
/// Both one-sided distances are refined to `accuracy`, so each of the three
/// intervals is at most that wide unless the budget ran out.
///
/// # Errors
///
/// As [`one_sided_boundary_hausdorff`].
pub fn boundary_hausdorff_distance(
    a: &ExactBRep,
    b: &ExactBRep,
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<BoundaryHausdorff, ExactMeasureError> {
    let forward = one_sided_boundary_hausdorff(a, b, accuracy, tolerance)?;
    let backward = one_sided_boundary_hausdorff(b, a, accuracy, tolerance)?;
    let mut distance = if backward.lower > forward.lower {
        backward
    } else {
        forward
    };
    distance.upper = forward.upper.max(backward.upper);
    Ok(BoundaryHausdorff {
        distance,
        forward,
        backward,
    })
}

/// Total order on finite bounds, ties broken first come, first served.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Key(Scalar, Reverse<usize>);

impl Eq for Key {}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0).then(self.1.cmp(&other.1))
    }
}

/// Point queries against one boundary, over a tree of its elements refined
/// on demand and kept for the next query.
struct Nearest<'a> {
    side: Side<'a>,
    nodes: Vec<Element>,
    children: Vec<Option<Vec<usize>>>,
    roots: Vec<usize>,
}

/// A certified interval on a point's distance to a boundary, and the
/// nearest boundary point found.
#[derive(Debug, Clone, Copy)]
struct PointDistance {
    lower: Scalar,
    upper: Scalar,
    nearest: Point3,
    /// The element `nearest` lies on.
    on: Shape,
}

impl<'a> Nearest<'a> {
    fn new(side: Side<'a>) -> Self {
        let nodes = side.elements.clone();
        let roots = (0..nodes.len()).collect();
        let children = vec![None; nodes.len()];
        Self {
            side,
            nodes,
            children,
            roots,
        }
    }

    fn children(&mut self, node: usize) -> Result<Vec<usize>, ExactMeasureError> {
        if let Some(known) = &self.children[node] {
            return Ok(known.clone());
        }
        let parent = self.side.size(&self.nodes[node]);
        let split = self.side.split(&self.nodes[node])?;
        // A pair that can shrink no further is a leaf.
        let ids = if split.iter().any(|child| self.side.size(child) >= parent) {
            Vec::new()
        } else {
            split
                .into_iter()
                .map(|child| {
                    self.nodes.push(child);
                    self.children.push(None);
                    self.nodes.len() - 1
                })
                .collect()
        };
        self.children[node] = Some(ids.clone());
        Ok(ids)
    }

    /// `d(p, boundary)` to within `accuracy` (never finer than the
    /// rounding the bounds already carry), or as far as the step budget
    /// allows; `None` when no boundary point was reached.
    fn query(
        &mut self,
        p: Point3,
        accuracy: Scalar,
    ) -> Result<Option<PointDistance>, ExactMeasureError> {
        let accuracy = accuracy.max(1e-12 * (1.0 + p.length()));
        let mut heap = BinaryHeap::new();
        let mut order = 0;
        for &root in &self.roots {
            let node = &self.nodes[root];
            if critical_toward(node, p) {
                let bound = point_lower_bound(p, &self.side, node)?;
                heap.push(Reverse((Key(bound, Reverse(order)), root)));
                order += 1;
            }
        }
        let mut best: Option<(Scalar, Point3, Shape)> = None;
        let mut lower = 0.0;
        let mut steps = 0;
        while let Some(Reverse((Key(bound, _), node))) = heap.pop() {
            lower = bound;
            if let Some(w) = self.nodes[node].witness {
                let d = (w - p).length();
                if best.is_none_or(|(current, ..)| d < current) {
                    best = Some((d, w, self.nodes[node].shape));
                }
            }
            let upper = best.map_or(Scalar::INFINITY, |(d, ..)| d);
            if upper.is_finite() && (bound >= upper || upper - bound <= accuracy) {
                break;
            }
            steps += 1;
            let children = self.children(node)?;
            if steps > POINT_STEPS || children.is_empty() {
                heap.push(Reverse((Key(bound, Reverse(order)), node)));
                break;
            }
            for child in children {
                let element = &self.nodes[child];
                if critical_toward(element, p) {
                    let bound = point_lower_bound(p, &self.side, element)?;
                    order += 1;
                    heap.push(Reverse((Key(bound, Reverse(order)), child)));
                }
            }
        }
        if let Some(Reverse((Key(bound, _), _))) = heap.peek() {
            lower = Scalar::min(lower, *bound);
        }
        Ok(best.map(|(upper, nearest, on)| PointDistance {
            lower: lower.min(upper),
            upper,
            nearest,
            on,
        }))
    }
}

/// A face of `B` matched to a face of `A`: on the same parameters, or on
/// `A`'s parameters once re-charted by a translate's shift.
#[derive(Debug, Clone)]
struct Match {
    face: usize,
    shifted: Option<Shifted>,
}

/// The faces of `B` matched to each face of `A` (see the module docs).
fn matches(a: &ExactBRep, b: &ExactBRep) -> Vec<Vec<Match>> {
    let fa = a.topology().faces();
    let fb = b.topology().faces();
    fa.iter()
        .map(|face_a| {
            fb.iter()
                .enumerate()
                .filter_map(|(face, face_b)| {
                    if same_face(a, face_a, b, face_b) {
                        Some(Match {
                            face,
                            shifted: None,
                        })
                    } else {
                        translate(a, face_a, b, face_b).map(|shifted| Match {
                            face,
                            shifted: Some(shifted),
                        })
                    }
                })
                .collect()
        })
        .collect()
}

/// Same surface family (see [`same_shape`]), and the very same trimming
/// pcurves with the same intervals in the same loops: one domain on both.
fn same_face(
    a: &ExactBRep,
    fa: &Face<axiolid_brep::SurfaceId>,
    b: &ExactBRep,
    fb: &Face<axiolid_brep::SurfaceId>,
) -> bool {
    let (Ok(sa), Ok(sb)) = (surface_of(a, fa.surface), surface_of(b, fb.surface)) else {
        return false;
    };
    if !same_shape(sa, sb) || fa.bounds.len() != fb.bounds.len() {
        return false;
    }
    let (ta, tb) = (a.topology(), b.topology());
    fa.bounds.iter().zip(&fb.bounds).all(|(ba, bb)| {
        if ba.orientation != bb.orientation || ba.outer != bb.outer {
            return false;
        }
        let (Some(la), Some(lb)) = (
            ta.loops().get(ba.loop_id.index()),
            tb.loops().get(bb.loop_id.index()),
        ) else {
            return false;
        };
        la.edges.len() == lb.edges.len()
            && la
                .edges
                .iter()
                .zip(&lb.edges)
                .enumerate()
                .all(|(i, (ua, ub))| {
                    let (Some(ca), Some(cb)) = (pcurve(a, ua.pcurve), pcurve(b, ub.pcurve)) else {
                        return false;
                    };
                    let ia = a.pcurve_interval(ba.loop_id, i);
                    let ib = b.pcurve_interval(bb.loop_id, i);
                    ia.is_some() && ia == ib && ca == cb
                })
    })
}

fn pcurve(brep: &ExactBRep, id: Option<axiolid_brep::Curve2Id>) -> Option<&axiolid_curve::Curve2> {
    id.and_then(|id| brep.curves2().get(id.index()))
}

/// Whether one parameter domain stands for a face on both surfaces, so
/// that `S_B(u, v)` is on `B`'s face wherever `S_A(u, v)` is on `A`'s: the
/// same family, whose chart (periods and poles) is then the same -- for a
/// cone the same radius and angle, which place its apex, and for a
/// B-spline the same degrees, knots and positive weights, so that the
/// difference is a convex combination of control point moves. Radii may
/// otherwise differ: they are folded into the coefficients.
fn same_shape(a: &Surface, b: &Surface) -> bool {
    match (a, b) {
        (Surface::Plane(_), Surface::Plane(_))
        | (Surface::Cylinder(_), Surface::Cylinder(_))
        | (Surface::EllipticalCylinder(_), Surface::EllipticalCylinder(_))
        | (Surface::Sphere(_), Surface::Sphere(_))
        | (Surface::Torus(_), Surface::Torus(_)) => true,
        (Surface::Cone(p), Surface::Cone(q)) => {
            p.radius == q.radius && p.semi_angle == q.semi_angle
        }
        (Surface::BSpline(p), Surface::BSpline(q)) => {
            let shape = |s: &axiolid_surface::BSplineSurface| {
                s.control_points.iter().map(Vec::len).collect::<Vec<_>>()
            };
            p.u_degree == q.u_degree
                && p.v_degree == q.v_degree
                && p.u_knots == q.u_knots
                && p.v_knots == q.v_knots
                && p.u_multiplicities == q.u_multiplicities
                && p.v_multiplicities == q.v_multiplicities
                && p.weights == q.weights
                && shape(p) == shape(q)
                && p.weights
                    .as_ref()
                    .is_none_or(|w| w.iter().flatten().all(|w| w.is_finite() && *w > 0.0))
        }
        _ => false,
    }
}

/// A surface as `o + sum phi_k(u, v) a_k`: its origin and coefficient
/// vectors, over the family's basis (see [`basis`]).
fn terms(surface: &Surface) -> Option<(Point3, Vec<Vec3>)> {
    Some(match surface {
        Surface::Plane(p) => (p.frame.origin, vec![p.frame.x, p.frame.y]),
        Surface::Cylinder(c) => (
            c.frame.origin,
            vec![c.frame.x * c.radius, c.frame.y * c.radius, c.frame.z],
        ),
        Surface::EllipticalCylinder(c) => (
            c.frame.origin,
            vec![
                c.frame.x * c.semi_axis_x,
                c.frame.y * c.semi_axis_y,
                c.frame.z,
            ],
        ),
        Surface::Cone(c) => {
            let slope = c.semi_angle.tan();
            (
                c.frame.origin,
                vec![
                    c.frame.x * c.radius,
                    c.frame.y * c.radius,
                    c.frame.x * slope,
                    c.frame.y * slope,
                    c.frame.z,
                ],
            )
        }
        Surface::Sphere(s) => (
            s.frame.origin,
            vec![
                s.frame.x * s.radius,
                s.frame.y * s.radius,
                s.frame.z * s.radius,
            ],
        ),
        Surface::Torus(t) => (
            t.frame.origin,
            vec![
                t.frame.x * t.major_radius,
                t.frame.y * t.major_radius,
                t.frame.x * t.minor_radius,
                t.frame.y * t.minor_radius,
                t.frame.z * t.minor_radius,
            ],
        ),
        _ => return None,
    })
}

/// The family's basis at `(u, v)`, each function's largest oscillation
/// from its value at the patch centre `m` over the patch, and each one's
/// largest magnitude there.
fn basis(
    surface: &Surface,
    lo: Point2,
    hi: Point2,
) -> Option<(Vec<Scalar>, Vec<Scalar>, Vec<Scalar>)> {
    let m = (lo + hi) * 0.5;
    let (hu, hv) = (0.5 * (hi.x - lo.x).abs(), 0.5 * (hi.y - lo.y).abs());
    // `cos` and `sin` are 1-Lipschitz and stay within [-1, 1].
    let (tu, tv) = (hu.min(2.0), hv.min(2.0));
    let reach_u = lo.x.abs().max(hi.x.abs());
    let reach_v = lo.y.abs().max(hi.y.abs());
    let (su, cu) = m.x.sin_cos();
    let (sv, cv) = m.y.sin_cos();
    Some(match surface {
        Surface::Plane(_) => (vec![m.x, m.y], vec![hu, hv], vec![reach_u, reach_v]),
        Surface::Cylinder(_) | Surface::EllipticalCylinder(_) => {
            (vec![cu, su, m.y], vec![tu, tu, hv], vec![1.0, 1.0, reach_v])
        }
        Surface::Cone(_) => {
            // |v cos u - m_v cos m_u| <= |v - m_v| + |m_v| |cos u - cos m_u|.
            let mixed = hv + m.y.abs() * tu;
            (
                vec![cu, su, m.y * cu, m.y * su, m.y],
                vec![tu, tu, mixed, mixed, hv],
                vec![1.0, 1.0, reach_v, reach_v, reach_v],
            )
        }
        Surface::Sphere(_) => (
            vec![cv * cu, cv * su, sv],
            vec![tv + tu, tv + tu, tv],
            vec![1.0, 1.0, 1.0],
        ),
        Surface::Torus(_) => (
            vec![cu, su, cv * cu, cv * su, sv],
            vec![tu, tu, tv + tu, tv + tu, tv],
            vec![1.0, 1.0, 1.0, 1.0, 1.0],
        ),
        _ => return None,
    })
}

/// A bound on `|S_A(u, v) - S_B(u, v)|` over the patch `[lo, hi]`, for two
/// surfaces matched by [`same_shape`].
fn matched_bound(a: &Surface, b: &Surface, lo: Point2, hi: Point2) -> Option<Scalar> {
    if let (Surface::BSpline(p), Surface::BSpline(q)) = (a, b) {
        // Equal knots and positive equal weights: `S_A - S_B` is a convex
        // combination of the control point differences, everywhere.
        let mut reach: Scalar = 0.0;
        let mut scale: Scalar = 0.0;
        for (row_a, row_b) in p.control_points.iter().zip(&q.control_points) {
            for (pa, pb) in row_a.iter().zip(row_b) {
                reach = reach.max((*pa - *pb).length());
                scale = scale.max(pa.length()).max(pb.length());
            }
        }
        return Some(reach + 1e-12 * scale + Scalar::MIN_POSITIVE);
    }
    let (oa, ka) = terms(a)?;
    let (ob, kb) = terms(b)?;
    let (phi, wobble, reach) = basis(a, lo, hi)?;
    let mut centre = oa - ob;
    let mut spread = 0.0;
    let mut scale = oa.length() + ob.length();
    for k in 0..phi.len() {
        let delta = ka[k] - kb[k];
        centre += delta * phi[k];
        spread += wobble[k] * delta.length();
        scale += reach[k] * (ka[k].length() + kb[k].length());
    }
    let bound = centre.length() + spread + 1e-12 * scale + Scalar::MIN_POSITIVE;
    bound.is_finite().then_some(bound)
}

/// Distinct translations matched faces carry: `S_B - S_A` where it is the
/// same at the centre and two corners of a matched face's parameter box,
/// at most [`MAX_DISPLACEMENTS`] of them. For a translate every face gives
/// the translation; a turned face gives none.
fn displacements(
    source: &Side<'_>,
    target: &ExactBRep,
    matched: &[Vec<Match>],
) -> Result<Vec<Vec3>, ExactMeasureError> {
    let mut found: Vec<Vec3> = Vec::new();
    for element in &source.elements {
        let Shape::Face { face, lo, hi, .. } = element.shape else {
            continue;
        };
        let surface = surface_of(source.brep, source.brep.topology().faces()[face].surface)?;
        for other in &matched[face] {
            let theirs = match &other.shifted {
                Some(shifted) => &shifted.surface,
                None => surface_of(target, target.topology().faces()[other.face].surface)?,
            };
            // The displacement at the centre and two corners of the box:
            // only a translation, the same at all three, is a seed.
            let moves: Vec<Vec3> = [(lo + hi) * 0.5, lo, hi]
                .iter()
                .filter_map(|q| {
                    let here = evaluate(surface, q.x, q.y).ok()?;
                    let there = evaluate(theirs, q.x, q.y).ok()?;
                    Some(there - here)
                })
                .collect();
            let [t, a, b] = moves[..] else {
                continue;
            };
            let rigid = (a - t).length().max((b - t).length()) <= 1e-9 * (1.0 + t.length());
            let fresh = rigid
                && t.is_finite()
                && t.length() > 0.0
                && found
                    .iter()
                    .all(|seen| (*seen - t).length() > 1e-9 * (1.0 + t.length()));
            if fresh && found.len() < MAX_DISPLACEMENTS {
                found.push(t);
            }
        }
    }
    Ok(found)
}

/// What the search has certified so far.
struct State {
    lower: Scalar,
    /// The witnesses, and the elements of each boundary they lie on.
    witness: Option<(Point3, Point3, Shape, Shape)>,
}

fn search(
    source: &Side<'_>,
    target: &mut Nearest<'_>,
    accuracy: Scalar,
    max_splits: usize,
    items: &[core::ops::Range<usize>],
) -> Result<Witnessed, ExactMeasureError> {
    let matched = matches(source.brep, target.side.brep);
    let point_accuracy = 0.25 * accuracy;
    let mut state = State {
        lower: 0.0,
        witness: None,
    };
    let mut faces: BinaryHeap<(Key, usize)> = BinaryHeap::new();
    let mut edges: BinaryHeap<(Key, usize)> = BinaryHeap::new();
    // Every element queued, by index; a queue holds its bound.
    let mut pieces: Vec<Element> = Vec::new();

    // A matched face moved by `t` puts the farthest point of a translate at
    // `A`'s support point against `t` (see `translate::support_point`).
    let seeds = displacements(source, target.side.brep, &matched)?
        .into_iter()
        .flat_map(|t| items.iter().map(move |edges| (t, edges.clone())));
    for (t, edges) in seeds {
        let Some((p, edge)) = support_point(source.brep, -t, edges) else {
            continue;
        };
        if let Some(found) = target.query(p, point_accuracy)? {
            if found.lower > state.lower || state.witness.is_none() {
                state.lower = state.lower.max(found.lower);
                let on = Shape::Edge {
                    edge,
                    t0: 0.0,
                    t1: 0.0,
                };
                state.witness = Some((p, found.nearest, on, found.on));
            }
        }
    }

    // Bound one element and queue it; measure its witness if it may raise
    // the lower bound.
    let mut admit = |element: Element,
                     state: &mut State,
                     faces: &mut BinaryHeap<(Key, usize)>,
                     edges: &mut BinaryHeap<(Key, usize)>,
                     pieces: &mut Vec<Element>|
     -> Result<(), ExactMeasureError> {
        let mut upper = Scalar::INFINITY;
        if let Shape::Face { face, lo, hi, .. } = element.shape {
            let surface = surface_of(source.brep, source.brep.topology().faces()[face].surface)?;
            for other in &matched[face] {
                let bound = match &other.shifted {
                    None => {
                        let target_surface = surface_of(
                            target.side.brep,
                            target.side.brep.topology().faces()[other.face].surface,
                        )?;
                        matched_bound(surface, target_surface, lo, hi)
                    }
                    Some(shifted) => matched_bound(surface, &shifted.surface, lo, hi)
                        .and_then(|bound| shifted.bound(bound, lo, hi)),
                };
                if let Some(bound) = bound {
                    upper = upper.min(bound);
                }
            }
        }
        // The witness raises the lower bound only if it can beat it.
        if let Some(w) = element.witness {
            if upper > state.lower {
                if let Some(found) = target.query(w, point_accuracy)? {
                    if found.lower > state.lower || state.witness.is_none() {
                        state.lower = state.lower.max(found.lower);
                        state.witness = Some((w, found.nearest, element.shape, found.on));
                    }
                    // `d(., dB)` is 1-Lipschitz about the witness.
                    upper = upper.min(found.upper + element.radius);
                }
            }
        }
        // Without a witness, about the enclosing sphere's centre, which
        // need not lie on the face.
        if element.witness.is_none() && upper > state.lower {
            if let Some(found) = target.query(element.centre, point_accuracy)? {
                upper = upper.min(found.upper + element.radius);
            }
        }
        if !upper.is_finite() {
            return Err(crate::exact::NOT_CONVERGED);
        }
        pieces.push(element);
        let key = Key(upper, Reverse(pieces.len()));
        match element.shape {
            Shape::Face { .. } => faces.push((key, pieces.len() - 1)),
            Shape::Edge { .. } => edges.push((key, pieces.len() - 1)),
        }
        Ok(())
    };

    for element in source.elements.clone() {
        admit(element, &mut state, &mut faces, &mut edges, &mut pieces)?;
    }

    let mut splits = 0;
    loop {
        // Every point of the boundary lies in a queued face patch, or in
        // one already bounded below `lower`.
        let face_top = faces.peek().map_or(Scalar::NEG_INFINITY, |(key, _)| key.0);
        let upper = face_top.max(state.lower);
        if upper - state.lower <= accuracy || splits >= max_splits {
            break;
        }
        splits += 1;
        // Alternate: an edge span that may hold a farther point than the
        // lower bound is split every other step, to place witnesses.
        let edge_worth = edges
            .peek()
            .is_some_and(|(key, _)| key.0 > state.lower + accuracy);
        let take_edge = edge_worth && (splits % 2 == 0 || faces.is_empty());
        let popped = if take_edge { edges.pop() } else { faces.pop() };
        let Some((key, index)) = popped else {
            break;
        };
        let element = pieces[index];
        let children = source.split(&element)?;
        let parent = source.size(&element);
        if children.iter().any(|child| source.size(child) >= parent) {
            // Can shrink no further: a face patch keeps its bound, and the
            // search stops there; an edge span is only a witness source.
            if !take_edge {
                faces.push((key, index));
                break;
            }
            continue;
        }
        for child in children {
            admit(child, &mut state, &mut faces, &mut edges, &mut pieces)?;
        }
    }

    let face_top = faces.peek().map_or(Scalar::NEG_INFINITY, |(key, _)| key.0);
    let (point_from, point_to, from, to) = state.witness.ok_or(crate::exact::NOT_CONVERGED)?;
    Ok(Witnessed {
        bounds: HausdorffBounds {
            lower: state.lower,
            upper: face_top.max(state.lower),
            point_from,
            point_to,
        },
        from,
        to,
    })
}

#[cfg(test)]
mod tests {
    //! The matched bound, checked by dense sampling: over a patch it holds
    //! `|S_A - S_B|` for every family, and it is not vacuous.

    use super::{matched_bound, same_shape};
    use axiolid_core::{Frame3, Point2, Point3, Vec3};
    use axiolid_curve::{BSplineSurface, KnotSpec};
    use axiolid_evaluate::surface::evaluate;
    use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface, Torus};

    fn frame(turn: f64, shift: Vec3) -> Frame3 {
        // A tilted frame, turned about its own z by `turn`, then shifted.
        let z = Vec3::new(0.0, -0.6, 0.8);
        let x0 = Vec3::new(0.6, 0.64, 0.48);
        let y0 = z.cross(x0);
        let (s, c) = turn.sin_cos();
        Frame3 {
            origin: Point3::new(1.5, -2.0, 0.75) + shift,
            x: x0 * c + y0 * s,
            y: y0 * c - x0 * s,
            z,
        }
    }

    fn families(f: Frame3) -> Vec<(Surface, Point2, Point2)> {
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

    /// The largest `|S_A - S_B|` over the samples, and the bound.
    fn reach_and_bound(a: &Surface, b: &Surface, lo: Point2, hi: Point2) -> (f64, f64) {
        let bound = matched_bound(a, b, lo, hi).expect("a matched pair");
        let mut reach: f64 = 0.0;
        for p in samples(lo, hi) {
            let pa = evaluate(a, p.x, p.y).expect("point");
            let pb = evaluate(b, p.x, p.y).expect("point");
            reach = reach.max((pa - pb).length());
        }
        (reach, bound)
    }

    #[test]
    fn a_translate_is_bounded_by_its_translation_exactly() {
        let shift = Vec3::new(0.3, -0.4, 1.2);
        let a = families(frame(0.0, Vec3::ZERO));
        let b = families(frame(0.0, shift));
        for ((sa, lo, hi), (sb, _, _)) in a.iter().zip(&b) {
            assert!(same_shape(sa, sb));
            let (reach, bound) = reach_and_bound(sa, sb, *lo, *hi);
            assert!(reach <= bound, "{sa:?}: {reach} > {bound}");
            assert!(
                (bound - shift.length()).abs() < 1e-9,
                "{sa:?}: {bound} is not the translation"
            );
        }
    }

    /// One frame axis moved alone, over patches tall in `v` or wide in
    /// `u` far from `v = 0`: each basis function's oscillation is needed
    /// on its own, which a rigid turn (moving two axes together) hides.
    #[test]
    fn each_axis_moved_alone_is_held_over_tall_and_wide_patches() {
        let base = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let tall = (Point2::new(0.3, 0.2), Point2::new(0.35, 1.4));
        let cases = |f: Frame3| {
            vec![
                (
                    Surface::Sphere(Sphere {
                        frame: f,
                        radius: 2.0,
                    }),
                    tall,
                ),
                (
                    Surface::Torus(Torus {
                        frame: f,
                        major_radius: 3.0,
                        minor_radius: 1.0,
                    }),
                    tall,
                ),
                (
                    Surface::Cone(Cone {
                        frame: f,
                        radius: 0.1,
                        semi_angle: 0.4,
                    }),
                    (Point2::new(0.0, 5.0), Point2::new(1.6, 5.2)),
                ),
            ]
        };
        for axis in 0..3 {
            let mut moved = base;
            let nudge = Vec3::new(0.3, 0.0, 0.0);
            match axis {
                0 => moved.x += nudge,
                1 => moved.y += nudge,
                _ => moved.z += nudge,
            }
            for ((sa, (lo, hi)), (sb, _)) in cases(base).iter().zip(cases(moved)) {
                let (reach, bound) = reach_and_bound(sa, &sb, *lo, *hi);
                assert!(reach <= bound, "{sa:?}, axis {axis}: {reach} > {bound}");
            }
        }
    }

    #[test]
    fn a_turned_and_moved_copy_is_held_and_not_vacuously() {
        let a = families(frame(0.0, Vec3::ZERO));
        let b = families(frame(0.05, Vec3::new(0.01, 0.0, -0.02)));
        for ((sa, lo, hi), (sb, _, _)) in a.iter().zip(&b) {
            let (reach, bound) = reach_and_bound(sa, sb, *lo, *hi);
            assert!(reach <= bound, "{sa:?}: {reach} > {bound}");
            assert!(bound <= 4.0 * reach + 1e-9, "{sa:?}: {bound} vs {reach}");
            // A quarter of the patch is held tighter.
            let mid = (*lo + *hi) * 0.5;
            let (small_reach, small) = reach_and_bound(sa, sb, *lo, mid);
            assert!(small_reach <= small, "{sa:?}: {small_reach} > {small}");
        }
    }

    #[test]
    fn a_different_family_cone_or_net_is_not_matched() {
        let f = frame(0.0, Vec3::ZERO);
        let thin = Surface::Cylinder(Cylinder {
            frame: f,
            radius: 1.0,
        });
        let thick = Surface::Cylinder(Cylinder {
            frame: f,
            radius: 1.5,
        });
        // Radii fold into the coefficients: still one family, and held.
        assert!(same_shape(&thin, &thick));
        let (lo, hi) = (Point2::new(0.2, -0.5), Point2::new(1.1, 0.5));
        let (reach, bound) = reach_and_bound(&thin, &thick, lo, hi);
        assert!(reach <= bound && bound <= 4.0 * reach, "{reach} vs {bound}");
        assert!(!same_shape(&thin, &Surface::Plane(Plane { frame: f })));
        let cone = |semi_angle| {
            Surface::Cone(Cone {
                frame: f,
                radius: 1.0,
                semi_angle,
            })
        };
        assert!(!same_shape(&cone(0.3), &cone(0.4)));
        // Another radius or angle moves a cone's apex, a pole of its chart.
        assert!(!same_shape(
            &cone(0.3),
            &Surface::Cone(Cone {
                frame: f,
                radius: 1.2,
                semi_angle: 0.3,
            })
        ));
    }

    #[test]
    fn a_spline_net_moved_point_by_point_is_held_by_its_largest_move() {
        let net = |lift: f64| BSplineSurface {
            u_degree: 1,
            v_degree: 1,
            control_points: vec![
                vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, lift)],
                vec![
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(1.0, 1.0, 0.5 * lift),
                ],
            ],
            u_knots: vec![0.0, 1.0],
            u_multiplicities: vec![2, 2],
            v_knots: vec![0.0, 1.0],
            v_multiplicities: vec![2, 2],
            weights: None,
            u_closed: false,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: None,
        };
        let (a, b) = (Surface::BSpline(net(0.0)), Surface::BSpline(net(0.2)));
        assert!(same_shape(&a, &b));
        let (lo, hi) = (Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
        let (reach, bound) = reach_and_bound(&a, &b, lo, hi);
        assert!(reach <= bound && bound <= 0.2 + 1e-9, "{reach} vs {bound}");
        // Other knots are another family member.
        let mut other = net(0.0);
        other.u_knots = vec![0.0, 2.0];
        assert!(!same_shape(&a, &Surface::BSpline(other)));
    }

    /// A one-face sheet on the plane through `origin`, its loop's three
    /// uses carrying `pcurve` over `interval`.
    fn sheet(
        origin: Point3,
        pcurve: axiolid_curve::Curve2,
        interval: axiolid_core::Interval,
    ) -> axiolid_brep::ExactBRep {
        use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};
        let mut builder = axiolid_brep::ExactBRepBuilder::default();
        let curve3 = builder.add_curve3(axiolid_curve::Curve3::Line(axiolid_curve::Line3 {
            origin,
            direction: Vec3::X,
        }));
        let curve2 = builder.add_curve2(pcurve);
        let surface = builder.add_surface(Surface::Plane(Plane {
            frame: Frame3 {
                origin,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        }));
        let topology = builder.topology_mut();
        let vertices: Vec<_> = [Point3::X, Point3::Y, Point3::Z]
            .into_iter()
            .map(|position| topology.add_vertex(Vertex { position }))
            .collect();
        let edges: Vec<_> = (0..3)
            .map(|i| {
                topology.add_edge(Edge {
                    start: vertices[i],
                    end: vertices[(i + 1) % 3],
                    curve: Some(curve3),
                })
            })
            .collect();
        let loop_id = topology.add_loop(Loop {
            edges: edges
                .iter()
                .map(|&edge| EdgeUse {
                    edge,
                    orientation: Orientation::Forward,
                    pcurve: Some(curve2),
                })
                .collect(),
        });
        topology.add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        });
        for edge in edges {
            builder.set_edge_interval(edge, axiolid_core::Interval::UNIT);
        }
        for use_index in 0..3 {
            builder.set_pcurve_interval(loop_id, use_index, interval);
        }
        builder.finish().expect("a valid sheet")
    }

    /// Faces match only on the very same trimming pcurves and intervals,
    /// wherever their planes are.
    #[test]
    fn faces_match_only_on_the_same_trim() {
        use axiolid_core::{Interval, Vec2};
        use axiolid_curve::{Curve2, Line2};
        let line = |dx: f64| {
            Curve2::Line(Line2 {
                origin: Point2::new(0.5, 0.0),
                direction: Vec2::new(dx, 1.0),
            })
        };
        let unit = Interval::UNIT;
        let a = sheet(Point3::ZERO, line(1.0), unit);
        let matched = |b: &axiolid_brep::ExactBRep| {
            super::same_face(&a, &a.topology().faces()[0], b, &b.topology().faces()[0])
        };
        assert!(matched(&sheet(
            Point3::new(3.0, 1.0, -2.0),
            line(1.0),
            unit
        )));
        assert!(!matched(&sheet(Point3::ZERO, line(2.0), unit)));
        assert!(!matched(&sheet(
            Point3::ZERO,
            line(1.0),
            Interval::new(0.0, 0.5)
        )));
        let found = super::matches(&a, &a);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].len(), 1);
        assert_eq!(found[0][0].face, 0);
        assert!(found[0][0].shifted.is_none());
    }
}
