//! Certified distance between the boundaries of two exact B-reps (#125, C18).
//!
//! # What is certified
//!
//! [`boundary_distance`] returns an interval `[lower, upper]` that contains
//! the true minimum distance between the two boundaries, and two computed
//! points, one on each boundary within its evaluation rounding, whose
//! boundary points `upper` covers (#439). Both ends are guarantees, not
//! estimates:
//!
//! - `upper` is at least the distance between two points that lie on the
//!   boundaries, which the two computed witnesses stand for: points on
//!   edges, or surface points at parameters the face's domain is certified
//!   to contain (the crate-private `Domain` classifier). It is their
//!   computed distance plus their rounding (see "Rounding").
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
//! Where a lower bound does not depend on a patch's size -- along a planar
//! patch's own normal, against a parallel planar patch or a line edge -- the
//! minimum is reached on a whole region and the bound is right from the
//! start; what never closes is the upper bound, since patch centres line up
//! only by chance. Such ties close on foot witnesses (`witness`): the foot,
//! on a planar face, of the other element's witness, and the closest points
//! of two line spans, each certified as any witness is (#318).
//!
//! # Touching and crossing boundaries
//!
//! Where the boundaries meet, the distance is zero, but witnesses at patch
//! centres and edge midpoints land on the contact only by chance, so the
//! refinement alone would close on it only when the patches round it were
//! within the accuracy -- often not within the step budget. A pair popped
//! with a lower bound of zero is therefore first asked for a point of each
//! element that the two share (#273): the closest points of two line edges,
//! or a point of a line edge and the point of a planar face at the same
//! place, each evaluated on its own boundary as any witness is, the face's
//! certified by the face's domain. Two such points no farther apart than
//! the rounding margin (128 ulps of the points' magnitudes: about `6e-14`
//! at unit scale) end the search with
//! `[0, d']`, `d'` their distance plus their rounding (#439), they being the witnesses: no lower bound can
//! rise above zero there, so this is the interval the search would only
//! converge to. Two boxes sharing a face, or crossing, stop within a few
//! steps. A gap wider than that margin is never taken for a touch, and boundaries that are apart measure exactly as
//! before. Curved faces and circle or ellipse edges show no touch this way
//! and are refined as before.
//!
//! # Rounding (#318)
//!
//! Every margin is derived from the magnitudes the computed numbers are
//! summed from, in units of `eps` (`f64::EPSILON`, twice the unit
//! roundoff), not from `1e-12` of the points' magnitudes, which left an
//! interval about `1e-5` wide at georeferenced coordinates (450000,
//! 5600000), where a point cannot be computed to better than a few `1e-9`.
//!
//! `M` is the summed magnitude of the terms a point of an element is
//! computed from, taken at the element's farthest parameters. It never
//! grows when an element is split, since a child's rectangle or span lies
//! in its parent's:
//!
//! | element | `M` |
//! |---|---|
//! | plane patch | `|o| + |u| Lu + |v| Lv` |
//! | cylinder patch | `|o| + 2 |r| s + |u| Lu + |v| Lv` |
//! | elliptical cylinder patch | `|o| + (|a| + |b|) s + ...` |
//! | cone patch | `|o| + 2 (|r| + |v| tan(angle)) s + ...` |
//! | sphere patch | `|o| + 3 |r| s + ...` |
//! | torus patch | `|o| + (2 R + 3 r) s + ...` |
//! | B-spline patch (sphere) | the largest control point |
//! | B-spline patch (witness) | `(p_u + p_v + 1) (w_max / w_min) max|P|` |
//! | line span | `|o| + |t| |d|` |
//! | circle span | `|o| + (2 + |t|) |r| s` |
//! | ellipse span | `|o| + (|a| + |b| + |t| max(|a|, |b|)) s` |
//!
//! where `o` is the frame origin, `s` the frame's scale, `Lu` and `Lv` the
//! Lipschitz bounds of the parameters. A side's scale `S` is the largest
//! `M` over its root elements, and every element of that side is padded by
//! it. The pad is therefore the same for every element, so it cannot
//! reorder the queue: a pair's lower bound is its geometric bound less a
//! constant of the whole search, and the rounding floor enters once.
//!
//! The margins, with `u = eps / 2` the unit roundoff and `gamma_n` the usual
//! `n u / (1 - n u)`:
//!
//! - `POINT` (16): a computed centre or witness against the exact point.
//!   A coordinate sums at most twelve rounded terms of magnitude at most
//!   `M` (`gamma_12 = 6 eps`), a vector's length takes a factor `sqrt 3`
//!   (10.4 `eps`), and the midpoint of the parameters and a cone's apex add
//!   the rest.
//! - `PROJECTION` (8): each end of a projected range `d . x`: a sum of
//!   rounded products of terms at most `M` against a unit direction.
//! - `RATIO` (8): relative: lengths, sums, differences, and directions that
//!   are unit only to rounding.
//! - `FLOOR` (128): the width a search settles at, in `eps (S_a + S_b)`.
//!   At zero size the computed bound may lie `2 (POINT + PROJECTION) eps
//!   (S_a + S_b) = 48 eps (S_a + S_b)` below the exact one, and the
//!   witnesses' own evaluation adds up to `32 eps (S_a + S_b)`: 80, with
//!   room to spare.
//!
//! A witness is a point evaluated at a parameter, so it is off its boundary
//! by `POINT eps W`, `W` the side's witness scale: the largest witness `M`
//! over its root elements, equal to `S` unless the side has a B-spline face.
//! The reported upper bound is
//! `|pa - pb| (1 + RATIO eps) + POINT eps (W_a + W_b)`, so it is at least
//! the distance between the boundary points the witnesses stand for, and
//! `Clearance::Below` is never claimed within rounding of a limit the
//! distance does not go below (#439). The search order and the witnesses
//! are those of the computed distance; only the stop test and the reported
//! value widen. A plan overlap shown by two planar faces still reports
//! `[0, 0]`. The floor adds `POINT eps ((W_a - S_a) + (W_b - S_b))`, nothing
//! without B-spline faces.
//!
//! A B-spline sphere's `M` is its control points, since its patch is
//! enclosed by them, but a witness runs `p_u + p_v` levels of de Boor in
//! homogeneous coordinates and divides by the weight. Each level adds at
//! most `4 eps` of the largest homogeneous coordinate `w_max |P|max`; the
//! division by a weight of at least `w_min` gives
//! `((8L + 0.5) rho + 0.5) eps |P|max` per coordinate (`L = p_u + p_v`,
//! `rho = w_max / w_min`). Times `sqrt 3` that is below
//! `16 (L + 1) rho eps |P|max`, which is `POINT eps M` for the witness `M`
//! of the table. A non-positive weight is refused.
//!
//! The `RATIO` guards of the pair and point bounds and a point query's own
//! pad `PROJECTION eps |p|` cover rounding that another margin covers too
//! (#447): `S` is at least the magnitude of every point of its side, so a
//! computed gap or projection is off by a few `eps (S_a + S_b)` at most,
//! inside the side pads. They stay, so each bound is sound by its own
//! arithmetic, but no test can observe their loss.
//!
//! The libm functions the evaluators call are assumed within one ulp.
//! The floor is about `3.2e-7` at (450000, 5600000) and about `6e-14` at
//! unit scale. It only stops refinement and changes neither bound; an
//! accuracy below it settles at it and is never reported as met, so the
//! width of the interval returned shows what was reached.
//!
//! # In plan
//!
//! [`plan_boundary_distance`] measures between the bodies' projections
//! onto the XY plane (#217). A solid's shadow is its boundary's, so the
//! same search applies with every gap, direction and enclosing disc
//! horizontal: elements carry a plan radius from their horizontal speeds,
//! lower bounds project onto horizontal directions only (a normal by its
//! horizontal part), and a face patch is dropped only where no normal of
//! it is horizontal and along the plan offset. [`plan_overlap`] shows two
//! shadows overlap over a patch of positive area from two planar faces,
//! not vertical, whose shadows share an open patch, and shows them apart
//! once the lower bound is positive. When the distance search shows
//! neither, a second search over patches of such faces alone, coarsest
//! pair first, looks for the overlap (#237): in the first, pairs of walls
//! and edges whose shadows merely cross take the budget, and a plan
//! distance of zero met between two boundary points ends it.
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

mod touch;
pub(crate) mod witness;

/// Pairs refined before a query stops and reports what it has.
const MAX_STEPS: usize = 400_000;

/// Pairs refined before [`plan_overlap`] reports it could not decide.
pub(crate) const OVERLAP_STEPS: usize = 20_000;

/// Rounding margins, in units of `eps * S`, where `S` is a side's rounding
/// scale: the largest summed magnitude `M` of the terms an element's points
/// are computed from (see "Rounding" in the module docs for the derivations).
///
/// A computed centre or witness against the exact point it stands for.
pub(crate) const POINT_EPSILONS: Scalar = 16.0;
/// Each end of a projected range `d . x`.
pub(crate) const PROJECTION_EPSILONS: Scalar = 8.0;
/// Relative: lengths, sums, differences and directions that are unit only
/// to rounding.
pub(crate) const RATIO_EPSILONS: Scalar = 8.0;
/// The width an interval settles at, in units of `eps * (S_a + S_b)`.
pub(crate) const FLOOR_EPSILONS: Scalar = 128.0;

/// An interval certain to contain the distance between two boundaries,
/// and two witnesses. `point_a` and `point_b` are computed points, each
/// within its evaluation rounding of a point of its boundary
/// (`POINT eps W`, see "Rounding" in the module docs), and `upper` is at
/// least the distance between those two boundary points (#439): the
/// computed `|point_a - point_b|` plus that rounding, never less. In plan,
/// an overlap shown by two planar faces reports `[0, 0]`.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceBounds {
    /// No two boundary points are closer than this.
    pub lower: Scalar,
    /// The boundary points `point_a` and `point_b` stand for are at most
    /// this far apart (in plan, for the plan queries).
    pub upper: Scalar,
    /// A computed point within rounding of the first boundary.
    pub point_a: Point3,
    /// A computed point within rounding of the second boundary.
    pub point_b: Point3,
}

/// How a certified distance compares with a limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Clearance {
    /// Certainly closer than the limit: `upper < limit`, `upper` covering
    /// the witnesses' rounding.
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
/// An accuracy below the rounding floor settles at the floor and is never
/// reported as met: the interval returned shows its width (see "Rounding"
/// in the module docs).
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
    Ok(search(
        a,
        b,
        tolerance,
        Metric::Space,
        &mut |lower, upper, floor| closed(lower, upper, accuracy, floor),
    )?
    .bounds)
}

/// Distance between the plan projections of `a` and `b` -- their shadows
/// on the XY plane -- to within `accuracy` (#217).
///
/// A solid's shadow is its boundary's: a vertical line through a point of
/// the solid meets the boundary. So this is the least horizontal distance
/// between two boundary points, found by the same certified search as
/// [`boundary_distance`] with horizontal gaps, directions and enclosing
/// discs. It is zero when the shadows overlap, and unlike the distance in
/// space it is zero for a body standing inside another's footprint.
/// `point_a` and `point_b` lie within rounding of the boundaries, their
/// projections at most `upper` apart. When two planar faces show the shadows overlap over a
/// patch of positive area, the interval is `[0, 0]` and the points project
/// into that patch (see [`plan_overlap`]).
///
/// # Errors
///
/// As [`boundary_distance`].
pub fn plan_boundary_distance(
    a: &ExactBRep,
    b: &ExactBRep,
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<DistanceBounds, ExactMeasureError> {
    let accuracy = accuracy.max(0.0);
    Ok(
        search(a, b, tolerance, Metric::Plan, &mut |lower, upper, floor| {
            closed(lower, upper, accuracy, floor)
        })?
        .bounds,
    )
}

/// [`plan_boundary_distance`] refined only until it clears `limit`.
///
/// # Errors
///
/// As [`boundary_distance`].
pub fn plan_boundary_clearance(
    a: &ExactBRep,
    b: &ExactBRep,
    limit: Scalar,
    tolerance: Tolerance,
) -> Result<(DistanceBounds, Clearance), ExactMeasureError> {
    let bounds = search(a, b, tolerance, Metric::Plan, &mut |lower, upper, floor| {
        upper < limit || lower > limit || closed(lower, upper, 0.0, floor)
    })?
    .bounds;
    let clearance = bounds.against(limit);
    Ok((bounds, clearance))
}

/// Whether the plan projections of two bodies overlap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlanOverlap {
    /// Over a patch of positive area: `at` and a disc around it lie in
    /// both shadows, shown by a planar face of each that is not vertical.
    Overlapping {
        /// A plan point inside both shadows.
        at: Point2,
    },
    /// Certainly apart: the plan distance is at least `gap`, positive.
    Disjoint {
        /// A certified lower bound on the plan distance.
        gap: Scalar,
    },
    /// Neither shown: the shadows may only touch, overlap where no planar
    /// face covers both, or lie closer than the search resolves.
    Undecided,
}

/// Whether the plan projections of `a` and `b` overlap over a patch of
/// positive area, are certainly apart, or neither could be shown (#217).
///
/// Overlap is shown by two planar faces, one of each body, not vertical,
/// whose shadows share an open patch: a slab under a column's base, two
/// slabs, a beam's soffit over a wall's top. Curved faces never show it;
/// their shadows are still measured by [`plan_boundary_distance`].
///
/// # Errors
///
/// As [`boundary_distance`].
pub fn plan_overlap(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
) -> Result<PlanOverlap, ExactMeasureError> {
    Ok(overlap_search(a, b, tolerance)?.0)
}

/// The elements of the two B-reps a result was shown on.
pub(crate) type Pair = (Shape, Shape);

/// [`plan_overlap`], with the elements an overlap was shown on.
///
/// First the distance search, which shows the shadows apart or meets an
/// overlap on its way; then, when it shows neither, a search over pairs of
/// patches of planar faces that are not vertical alone, coarsest first, so
/// that pairs of walls and edges whose shadows merely cross cannot use up
/// the budget, and a plan distance of zero shown by two boundary points
/// does not end it.
pub(crate) fn overlap_search(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
) -> Result<(PlanOverlap, Option<Pair>), ExactMeasureError> {
    // Overlap on curved faces alone is never shown, so the search is
    // bounded more tightly than a distance query's.
    let found = search_within(
        a,
        b,
        tolerance,
        Metric::Plan,
        OVERLAP_STEPS,
        &mut |lower, _, _| lower > 0.0,
    )?;
    if let Some(at) = found.overlap {
        return Ok((PlanOverlap::Overlapping { at }, Some(found.on)));
    }
    if found.bounds.lower > 0.0 {
        return Ok((
            PlanOverlap::Disjoint {
                gap: found.bounds.lower,
            },
            None,
        ));
    }
    Ok(match level_overlap(a, b, tolerance, OVERLAP_STEPS)? {
        Some((at, on)) => (PlanOverlap::Overlapping { at }, Some(on)),
        None => (PlanOverlap::Undecided, None),
    })
}

/// Search pairs of patches of planar faces that are not vertical, one of
/// each B-rep, whose shadows may meet, refining the coarsest pair first,
/// for two inside patches whose shadows share an open patch.
fn level_overlap(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
    max_steps: usize,
) -> Result<Option<(Point2, Pair)>, ExactMeasureError> {
    let linear = tolerance.linear().max(1e-12);
    let side_a = Side::new(a, linear, Metric::Plan)?;
    let side_b = Side::new(b, linear, Metric::Plan)?;
    let level = |side: &Side<'_>| -> Result<Vec<Element>, ExactMeasureError> {
        let mut out = Vec::new();
        for element in &side.elements {
            if let Shape::Face { face, .. } = element.shape {
                if level_frame(side, face)?.is_some() {
                    out.push(*element);
                }
            }
        }
        Ok(out)
    };
    let mut elements_a = level(&side_a)?;
    let mut elements_b = level(&side_b)?;
    // Coarsest first: the pair with the larger plan discs.
    let mut heap = BinaryHeap::new();
    let push = |heap: &mut BinaryHeap<(Key, usize, usize)>,
                ea: &Element,
                eb: &Element,
                i: usize,
                j: usize|
     -> Result<(), ExactMeasureError> {
        if lower_bound(&side_a, ea, &side_b, eb)? <= 0.0 {
            heap.push((Key(ea.plan_radius + eb.plan_radius), i, j));
        }
        Ok(())
    };
    for (i, ea) in elements_a.iter().enumerate() {
        for (j, eb) in elements_b.iter().enumerate() {
            push(&mut heap, ea, eb, i, j)?;
        }
    }
    let mut steps = 0;
    while let Some((_, i, j)) = heap.pop() {
        let (ea, eb) = (elements_a[i], elements_b[j]);
        if let Some((at, _, _)) = plan_patches_overlap(&side_a, &ea, &side_b, &eb)? {
            return Ok(Some((at, (ea.shape, eb.shape))));
        }
        steps += 1;
        if steps > max_steps {
            break;
        }
        let split_a = side_a.size(&ea) >= side_b.size(&eb);
        let (side, parent) = if split_a {
            (&side_a, ea)
        } else {
            (&side_b, eb)
        };
        let children = side.split(&parent)?;
        // A pair that can shrink no further shows nothing more.
        if children
            .iter()
            .any(|child| side.size(child) >= side.size(&parent))
        {
            continue;
        }
        for child in children {
            if split_a {
                elements_a.push(child);
                push(&mut heap, &child, &eb, elements_a.len() - 1, j)?;
            } else {
                elements_b.push(child);
                push(&mut heap, &ea, &child, i, elements_b.len() - 1)?;
            }
        }
    }
    Ok(None)
}

/// The frame of face `face` when it is planar and not vertical (nor within
/// rounding of it): a plane whose shadow has area.
fn level_frame(
    side: &Side<'_>,
    face: usize,
) -> Result<Option<axiolid_core::Frame3>, ExactMeasureError> {
    let surface = surface_of(side.brep, side.brep.topology().faces()[face].surface)?;
    let Surface::Plane(plane) = surface else {
        return Ok(None);
    };
    let f = plane.frame;
    let det = f.x.x * f.y.y - f.x.y * f.y.x;
    let scale = f.x.length() * f.y.length();
    // A vertical or near-vertical plane casts no area worth trusting.
    if det.abs().is_nan() || det.abs() <= 1e-9 * scale {
        return Ok(None);
    }
    Ok(Some(f))
}

/// Refine the distance only until it clears `limit`.
///
/// Like an accuracy below the rounding floor, a limit within the floor of
/// the distance settles at the floor and is never decided: the result is
/// [`Clearance::Indeterminate`].
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
    let bounds = search(
        a,
        b,
        tolerance,
        Metric::Space,
        &mut |lower, upper, floor| {
            upper < limit || lower > limit || closed(lower, upper, 0.0, floor)
        },
    )?
    .bounds;
    let clearance = bounds.against(limit);
    Ok((bounds, clearance))
}

/// Which distance a search measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Metric {
    /// In space.
    Space,
    /// Between projections onto the XY plane: every gap, direction and
    /// enclosing radius horizontal.
    Plan,
}

/// The horizontal part of a vector.
fn flat(v: Vec3) -> Vec3 {
    Vec3::new(v.x, v.y, 0.0)
}

/// Length of the horizontal part of a vector.
fn flat_length(v: Vec3) -> Scalar {
    v.x.hypot(v.y)
}

/// What an element covers.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Shape {
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
pub(crate) struct Element {
    pub(crate) shape: Shape,
    pub(crate) centre: Point3,
    pub(crate) radius: Scalar,
    /// A disc about the centre's projection holding the element's
    /// projection onto the XY plane; never above `radius`.
    pub(crate) plan_radius: Scalar,
    /// A point certainly on the boundary, when one is known.
    pub(crate) witness: Option<Point3>,
    /// The surface normal at the centre of a face patch.
    pub(crate) normal: Option<Vec3>,
    /// Half-angle of a cone about `normal` holding every normal of the
    /// patch, when the patch may be dropped from a pair whose directions
    /// its normals cannot meet (see [`critical_possible`]).
    pub(crate) spread: Option<Scalar>,
}

/// One B-rep, prepared for bounding.
pub(crate) struct Side<'a> {
    pub(crate) brep: &'a ExactBRep,
    metric: Metric,
    domains: Vec<Option<Domain<'a>>>,
    /// Whether every edge on the face's boundary is an element, so a
    /// closest point on that boundary is found through the edges.
    edges_bounded: Vec<bool>,
    pub(crate) elements: Vec<Element>,
    /// `S`: the largest `M` over the root elements, the rounding scale of a
    /// search over this side. Every element is padded by it, so the pads
    /// cannot reorder the queue (#318).
    pub(crate) scale: Scalar,
    /// `W`: the largest `M` of a point evaluated on this side, at least `S`.
    pub(crate) witness_scale: Scalar,
}

impl<'a> Side<'a> {
    pub(crate) fn new(
        brep: &'a ExactBRep,
        linear: Scalar,
        metric: Metric,
    ) -> Result<Self, ExactMeasureError> {
        let topology = brep.topology();
        let mut side = Side {
            brep,
            metric,
            domains: Vec::with_capacity(topology.faces().len()),
            edges_bounded: Vec::with_capacity(topology.faces().len()),
            elements: Vec::new(),
            scale: 0.0,
            witness_scale: 0.0,
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
        // First the scale, over every root element; then the elements, all
        // padded by it.
        let mut roots = Vec::with_capacity(topology.faces().len());
        for (index, face) in topology.faces().iter().enumerate() {
            let surface = surface_of(brep, face.surface)?;
            let (lo, hi) = match &side.domains[index] {
                Some(domain) => (domain.min, domain.max),
                None => natural_range(surface).ok_or(ExactMeasureError::NonPlanarFace(
                    "an unbounded face trimmed by a pcurve family the distance query cannot bound",
                ))?,
            };
            side.scale = side.scale.max(patch_scale(surface, lo, hi)?);
            side.witness_scale = side.witness_scale.max(witness_scale(surface, lo, hi)?);
            roots.push((index, lo, hi));
        }
        let mut spans = Vec::new();
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
            if let Some(scale) = span_scale(curve, span.start, span.end) {
                side.scale = side.scale.max(scale);
                side.witness_scale = side.witness_scale.max(scale);
            }
            spans.push((curve, index, span));
        }
        if !(side.scale.is_finite() && side.witness_scale.is_finite()) {
            return Err(crate::exact::EVALUATION);
        }
        for (index, lo, hi) in roots {
            if let Some(element) = side.face_element(index, lo, hi, false)? {
                side.elements.push(element);
            }
        }
        for (curve, index, span) in spans {
            if let Some(element) = edge_element(curve, index, span.start, span.end, side.scale)? {
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
        let (centre, radius) = patch_sphere(surface, lo, hi, self.scale)?;
        let plan_radius = patch_plan_radius(surface, lo, hi, centre, radius, self.scale)?;
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
            plan_radius,
            witness,
        }))
    }

    /// How far a point evaluated on this side may lie from the exact point.
    pub(crate) fn witness_error(&self) -> Scalar {
        POINT_EPSILONS * Scalar::EPSILON * self.witness_scale
    }

    /// What [`Side::witness_error`] adds to the floor beyond the
    /// `POINT eps S` it already holds: zero but for B-spline faces.
    pub(crate) fn witness_excess(&self) -> Scalar {
        POINT_EPSILONS * Scalar::EPSILON * (self.witness_scale - self.scale)
    }

    /// The radius a search by this side's metric refines by.
    pub(crate) fn size(&self, element: &Element) -> Scalar {
        match self.metric {
            Metric::Space => element.radius,
            Metric::Plan => element.plan_radius,
        }
    }

    /// The element's two halves, less any certified outside the face: a
    /// face patch across its metrically longer side, an edge span in the
    /// middle.
    pub(crate) fn split(&self, element: &Element) -> Result<Vec<Element>, ExactMeasureError> {
        match element.shape {
            Shape::Face {
                face,
                lo,
                hi,
                inside,
            } => {
                let surface = surface_of(self.brep, self.brep.topology().faces()[face].surface)?;
                let (lu, lv) = match self.metric {
                    Metric::Space => lipschitz(surface, lo, hi)?,
                    Metric::Plan => plan_lipschitz(surface, lo, hi)?,
                };
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
                    if let Some(child) = edge_element(curve, edge, a, b, self.scale)? {
                        out.push(child);
                    }
                }
                Ok(out)
            }
        }
    }
}

pub(crate) fn surface_of(
    brep: &ExactBRep,
    id: Option<axiolid_brep::SurfaceId>,
) -> Result<&Surface, ExactMeasureError> {
    let id = id.ok_or(ExactMeasureError::MissingSurface)?;
    brep.surfaces()
        .get(id.index())
        .ok_or(ExactMeasureError::DanglingReference)
}

fn point_on(surface: &Surface, at: Point2) -> Result<Point3, ExactMeasureError> {
    evaluate(surface, at.x, at.y).map_err(|_| crate::exact::EVALUATION)
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
        Err(crate::exact::EVALUATION)
    }
}

/// Bounds on `|P S_u|` and `|P S_v|` over the patch, `P` dropping the
/// vertical: the horizontal speeds, never above [`lipschitz`]'s.
fn plan_lipschitz(
    surface: &Surface,
    lo: Point2,
    hi: Point2,
) -> Result<(Scalar, Scalar), ExactMeasureError> {
    let full = lipschitz(surface, lo, hi)?;
    // |P (a X + b Y)| <= sqrt(a^2 |P X|^2 + b^2 |P Y|^2) when a^2 + b^2 <= 1.
    let turning = |frame: &axiolid_core::Frame3, a: Scalar, b: Scalar| {
        (a * flat_length(frame.x)).hypot(b * flat_length(frame.y))
    };
    let bounds = match surface {
        Surface::Plane(p) => (flat_length(p.frame.x), flat_length(p.frame.y)),
        Surface::Cylinder(c) => (
            turning(&c.frame, c.radius.abs(), c.radius.abs()),
            flat_length(c.frame.z),
        ),
        Surface::EllipticalCylinder(c) => (
            turning(&c.frame, c.semi_axis_x.abs(), c.semi_axis_y.abs()),
            flat_length(c.frame.z),
        ),
        Surface::Cone(c) => {
            let slope = c.semi_angle.tan();
            let radius = (c.radius + lo.y * slope)
                .abs()
                .max((c.radius + hi.y * slope).abs());
            (
                turning(&c.frame, radius, radius),
                turning(&c.frame, slope.abs(), slope.abs()) + flat_length(c.frame.z),
            )
        }
        _ => return Ok(full),
    };
    if bounds.0.is_finite() && bounds.1.is_finite() {
        Ok((bounds.0.min(full.0), bounds.1.min(full.1)))
    } else {
        Ok(full)
    }
}

/// A radius about the centre's projection holding the patch's projection:
/// as [`patch_sphere`], from the horizontal speeds, and never above the
/// sphere's own `radius`.
fn patch_plan_radius(
    surface: &Surface,
    lo: Point2,
    hi: Point2,
    centre: Point3,
    radius: Scalar,
    scale: Scalar,
) -> Result<Scalar, ExactMeasureError> {
    if let Surface::BSpline(spline) = surface {
        let reach = spline
            .control_points
            .iter()
            .flatten()
            .map(|p| flat_length(*p - centre))
            .fold(0.0, Scalar::max);
        return Ok(pad(reach, scale).min(radius));
    }
    let (lu, lv) = plan_lipschitz(surface, lo, hi)?;
    let plan = 0.5 * ((hi.x - lo.x).abs() * lu + (hi.y - lo.y).abs() * lv);
    Ok(pad(plan, scale).min(radius))
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
    scale: Scalar,
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
        return Ok((centre, pad(radius, scale)));
    }
    let centre = point_on(surface, (lo + hi) * 0.5)?;
    let (lu, lv) = lipschitz(surface, lo, hi)?;
    let radius = 0.5 * ((hi.x - lo.x).abs() * lu + (hi.y - lo.y).abs() * lv);
    Ok((centre, pad(radius, scale)))
}

/// Widen a bounding radius by the rounding its centre may carry: the radius
/// itself is a computed length, and the centre is a computed point whose
/// error scales with the side's `scale` (`S`), the same for every element,
/// not with the distance from the origin alone.
pub(crate) fn pad(radius: Scalar, scale: Scalar) -> Scalar {
    radius * (1.0 + RATIO_EPSILONS * Scalar::EPSILON)
        + POINT_EPSILONS * Scalar::EPSILON * scale
        + Scalar::MIN_POSITIVE
}

/// The width an interval settles at, whatever accuracy is asked: the
/// rounding of the lower bounds at zero size (about `32 eps S`, twice that
/// here) and of the upper, for a search of scale `scale` (`S_a + S_b`, or
/// `|p| + S` for a point) (#318).
pub(crate) fn floor(scale: Scalar, upper: Scalar) -> Scalar {
    FLOOR_EPSILONS * Scalar::EPSILON * scale + RATIO_EPSILONS * Scalar::EPSILON * upper
}

/// Whether `[lower, upper]` is as narrow as asked, or as the rounding
/// allows (#318). An accuracy below the floor is never reported as met: the
/// interval returned shows its width.
pub(crate) fn closed(lower: Scalar, upper: Scalar, accuracy: Scalar, floor: Scalar) -> bool {
    upper - lower <= accuracy.max(floor)
}

/// `M` of a patch: the magnitude of the terms `S(u, v)` is summed from. The
/// frame origin, the frame-scaled radii and the parameter terms
/// `|u| Lu + |v| Lv` at the patch's farthest parameters.
fn patch_scale(surface: &Surface, lo: Point2, hi: Point2) -> Result<Scalar, ExactMeasureError> {
    if let Surface::BSpline(spline) = surface {
        let scale = spline
            .control_points
            .iter()
            .flatten()
            .map(|p| p.length())
            .fold(0.0, Scalar::max);
        return if scale.is_finite() {
            Ok(scale)
        } else {
            Err(crate::exact::EVALUATION)
        };
    }
    let (lu, lv) = lipschitz(surface, lo, hi)?;
    let params = lo.x.abs().max(hi.x.abs()) * lu + lo.y.abs().max(hi.y.abs()) * lv;
    let scale = match surface {
        Surface::Plane(p) => p.frame.origin.length() + params,
        Surface::Cylinder(c) => {
            c.frame.origin.length() + 2.0 * c.radius.abs() * frame_scale(&c.frame) + params
        }
        Surface::EllipticalCylinder(c) => {
            c.frame.origin.length()
                + (c.semi_axis_x.abs() + c.semi_axis_y.abs()) * frame_scale(&c.frame)
                + params
        }
        Surface::Cone(c) => {
            let reach = lo.y.abs().max(hi.y.abs());
            c.frame.origin.length()
                + 2.0 * (c.radius.abs() + reach * c.semi_angle.tan().abs()) * frame_scale(&c.frame)
                + params
        }
        Surface::Sphere(s) => {
            s.frame.origin.length() + 3.0 * s.radius.abs() * frame_scale(&s.frame) + params
        }
        Surface::Torus(t) => {
            t.frame.origin.length()
                + (2.0 * t.major_radius.abs() + 3.0 * t.minor_radius.abs()) * frame_scale(&t.frame)
                + params
        }
        _ => {
            return Err(ExactMeasureError::NonPlanarFace(crate::exact::family(
                surface,
            )))
        }
    };
    if scale.is_finite() {
        Ok(scale)
    } else {
        Err(crate::exact::EVALUATION)
    }
}

/// `M` of the points evaluated on a patch as witnesses: [`patch_scale`]'s,
/// except for a B-spline (#439).
fn witness_scale(surface: &Surface, lo: Point2, hi: Point2) -> Result<Scalar, ExactMeasureError> {
    let Surface::BSpline(spline) = surface else {
        return patch_scale(surface, lo, hi);
    };
    let largest = spline
        .control_points
        .iter()
        .flatten()
        .map(|p| p.length())
        .fold(0.0, Scalar::max);
    let ratio = match &spline.weights {
        None => 1.0,
        Some(weights) => {
            let (mut lo, mut hi) = (Scalar::INFINITY, 0.0_f64);
            for w in weights.iter().flatten() {
                if w.is_nan() || *w <= 0.0 {
                    return Err(ExactMeasureError::NonPlanarFace(
                        "non-positive-weight B-spline",
                    ));
                }
                lo = lo.min(*w);
                hi = hi.max(*w);
            }
            hi / lo
        }
    };
    let levels = Scalar::from(spline.u_degree) + Scalar::from(spline.v_degree) + 1.0;
    let scale = levels * ratio * largest;
    if scale.is_finite() {
        Ok(scale)
    } else {
        Err(crate::exact::EVALUATION)
    }
}

/// The upper bound reported for witnesses computed `distance` apart.
pub(crate) fn witnessed_upper(distance: Scalar, error: Scalar) -> Scalar {
    distance * (1.0 + RATIO_EPSILONS * Scalar::EPSILON) + error
}

/// `M` of an edge span, or `None` for a curve family with no bound here.
fn span_scale(curve: &Curve3, t0: Scalar, t1: Scalar) -> Option<Scalar> {
    let t = t0.abs().max(t1.abs());
    match curve {
        Curve3::Line(line) => Some(line.origin.length() + t * line.direction.length()),
        Curve3::Circle(c) => {
            Some(c.frame.origin.length() + (2.0 + t) * c.radius.abs() * frame_scale(&c.frame))
        }
        Curve3::Ellipse(e) => {
            let (a, b) = (e.semi_axis_x.abs(), e.semi_axis_y.abs());
            Some(e.frame.origin.length() + (a + b + t * a.max(b)) * frame_scale(&e.frame))
        }
        _ => None,
    }
}

/// An edge span with its enclosing sphere and its midpoint as a witness,
/// or `None` for a curve family with no derivative bound here.
fn edge_element(
    curve: &Curve3,
    edge: usize,
    t0: Scalar,
    t1: Scalar,
    scale: Scalar,
) -> Result<Option<Element>, ExactMeasureError> {
    let speed = match curve {
        Curve3::Line(line) => line.direction.length(),
        Curve3::Circle(circle) => circle.radius.abs() * frame_scale(&circle.frame),
        Curve3::Ellipse(ellipse) => {
            ellipse.semi_axis_x.abs().max(ellipse.semi_axis_y.abs()) * frame_scale(&ellipse.frame)
        }
        _ => return Ok(None),
    };
    // The horizontal speed, as for a patch (see `plan_lipschitz`).
    let plan_speed = match curve {
        Curve3::Line(line) => flat_length(line.direction),
        Curve3::Circle(c) => {
            (c.radius.abs() * flat_length(c.frame.x)).hypot(c.radius.abs() * flat_length(c.frame.y))
        }
        Curve3::Ellipse(e) => (e.semi_axis_x.abs() * flat_length(e.frame.x))
            .hypot(e.semi_axis_y.abs() * flat_length(e.frame.y)),
        _ => speed,
    };
    let centre = evaluate3(curve, 0.5 * (t0 + t1)).map_err(|_| crate::exact::EVALUATION)?;
    let radius = pad(0.5 * (t1 - t0).abs() * speed, scale);
    let plan_radius = if plan_speed.is_finite() {
        pad(0.5 * (t1 - t0).abs() * plan_speed, scale).min(radius)
    } else {
        radius
    };
    Ok(Some(Element {
        normal: None,
        spread: None,
        shape: Shape::Edge { edge, t0, t1 },
        centre,
        radius,
        plan_radius,
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
///
/// Between plan projections the same holds with the offset horizontal: a
/// point inside a face is critical for the horizontal distance only where
/// its normal is horizontal and along the plan offset, as on a vertical
/// wall or at a silhouette.
fn critical_possible(face: &Element, other: &Element, metric: Metric) -> bool {
    let (offset, reach) = match metric {
        Metric::Space => (other.centre - face.centre, face.radius + other.radius),
        Metric::Plan => (
            flat(other.centre - face.centre),
            face.plan_radius + other.plan_radius,
        ),
    };
    critical_along(face, offset, reach)
}

/// [`critical_possible`] towards a single point `p`, in space: whether the
/// face patch can hold a point nearest to `p` that no edge holds.
pub(crate) fn critical_toward(face: &Element, p: Point3) -> bool {
    critical_along(face, p - face.centre, face.radius)
}

/// Whether some normal of the face patch can meet a direction from it to
/// a ball of `reach` about `face.centre + offset`.
fn critical_along(face: &Element, offset: Vec3, reach: Scalar) -> bool {
    let (Some(normal), Some(spread)) = (face.normal, face.spread) else {
        return true;
    };
    let gap = offset.length();
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
    pub(crate) fn project(
        &self,
        element: &Element,
        d: Vec3,
    ) -> Result<(Scalar, Scalar), ExactMeasureError> {
        // Rounding in the projection itself, on either range, the same for
        // every element of the side (#318).
        let pad = PROJECTION_EPSILONS * Scalar::EPSILON * self.scale;
        let sphere = match self.metric {
            Metric::Space => sphere_range(element, d),
            // A horizontal direction reads only the projection, which the
            // plan disc holds.
            Metric::Plan => {
                let c = element.centre.dot(d);
                (c - element.plan_radius, c + element.plan_radius)
            }
        };
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
                            return Ok((sphere.0 - pad, sphere.1 + pad));
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
        let (lo, hi) = exact.map_or(sphere, |(lo, hi)| (lo.max(sphere.0), hi.min(sphere.1)));
        Ok((lo - pad, hi + pad))
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
    let plan = side_a.metric == Metric::Plan;
    let (offset, reach) = if plan {
        (flat(b.centre - a.centre), a.plan_radius + b.plan_radius)
    } else {
        (b.centre - a.centre, a.radius + b.radius)
    };
    let gap = offset.length();
    let rounding = RATIO_EPSILONS * Scalar::EPSILON * (gap + reach);
    let mut best = (gap - reach - rounding).max(0.0);
    // Between projections only horizontal directions separate them: a
    // normal counts by its horizontal part.
    let direction = |n: Option<Vec3>| -> Option<Vec3> {
        let n = n?;
        if !plan {
            return Some(n);
        }
        let h = flat(n);
        let length = h.length();
        (length > 1e-6 * n.length()).then(|| h / length)
    };
    let directions = [
        (gap > 0.0).then(|| offset / gap),
        direction(a.normal),
        direction(b.normal),
    ];
    for d in directions.into_iter().flatten() {
        let (a_lo, a_hi) = side_a.project(a, d)?;
        let (b_lo, b_hi) = side_b.project(b, d)?;
        best = best.max(b_lo - a_hi).max(a_lo - b_hi);
    }
    Ok(best * (1.0 - RATIO_EPSILONS * Scalar::EPSILON))
}

/// A lower bound on the distance from the point `p` to every point of
/// `element`, in space: the enclosing sphere's, raised by the exact range
/// of `d . x` along the line to the centre and along the patch's normal.
pub(crate) fn point_lower_bound(
    p: Point3,
    side: &Side<'_>,
    element: &Element,
) -> Result<Scalar, ExactMeasureError> {
    let offset = element.centre - p;
    let gap = offset.length();
    let rounding = RATIO_EPSILONS * Scalar::EPSILON * (gap + element.radius);
    let mut best = (gap - element.radius - rounding).max(0.0);
    let directions = [(gap > 0.0).then(|| offset / gap), element.normal];
    for d in directions.into_iter().flatten() {
        let length = d.length();
        if length.is_nan() || length <= 0.0 {
            continue;
        }
        let d = d / length;
        let (lo, hi) = side.project(element, d)?;
        let x = p.dot(d);
        let pad = PROJECTION_EPSILONS * Scalar::EPSILON * p.length();
        best = best.max(lo - x - pad).max(x - hi - pad);
    }
    Ok(best * (1.0 - RATIO_EPSILONS * Scalar::EPSILON))
}

/// Total order on finite lower bounds for the heap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Key(pub(crate) Scalar);

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

/// What a search found: the bounds, the elements the witnesses lie on, and
/// a plan point shown inside both shadows when planar faces show them
/// overlapping.
pub(crate) struct Found {
    pub(crate) bounds: DistanceBounds,
    /// The element of each boundary holding `point_a` and `point_b`.
    pub(crate) on: (Shape, Shape),
    pub(crate) overlap: Option<Point2>,
    /// Pairs refined, which the tests count.
    #[cfg(test)]
    pub(crate) steps: usize,
}

pub(crate) fn search(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
    metric: Metric,
    done: &mut dyn FnMut(Scalar, Scalar, Scalar) -> bool,
) -> Result<Found, ExactMeasureError> {
    search_within(a, b, tolerance, metric, MAX_STEPS, done)
}

pub(crate) fn search_within(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
    metric: Metric,
    max_steps: usize,
    done: &mut dyn FnMut(Scalar, Scalar, Scalar) -> bool,
) -> Result<Found, ExactMeasureError> {
    search_with(a, b, tolerance, metric, max_steps, done, true)
}

/// [`search_within`], with the early end on a touch in space (#273) on or
/// off: off, it is the search before it, which the tests compare against.
fn search_with(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
    metric: Metric,
    max_steps: usize,
    done: &mut dyn FnMut(Scalar, Scalar, Scalar) -> bool,
    early_touch: bool,
) -> Result<Found, ExactMeasureError> {
    let linear = tolerance.linear().max(1e-12);
    let side_a = Side::new(a, linear, metric)?;
    let side_b = Side::new(b, linear, metric)?;
    let scale = side_a.scale + side_b.scale;
    let error = side_a.witness_error() + side_b.witness_error();
    let excess = side_a.witness_excess() + side_b.witness_excess();
    let apart = |p: Point3, q: Point3| match metric {
        Metric::Space => (p - q).length(),
        Metric::Plan => flat_length(p - q),
    };
    // A witness counts as soon as its pair is queued, not when it is popped
    // (#318): refining pairs whose bound cannot rise never runs ahead of the
    // witnesses in hand. Ties keep the first.
    let offer =
        |best: &mut Option<(Scalar, Point3, Point3, (Shape, Shape))>, a: &Element, b: &Element| {
            if let (Some(wa), Some(wb)) = (a.witness, b.witness) {
                let d = apart(wa, wb);
                if best.is_none_or(|(current, ..)| d < current) {
                    *best = Some((d, wa, wb, (a.shape, b.shape)));
                }
            }
        };
    let mut elements_a = side_a.elements.clone();
    let mut elements_b = side_b.elements.clone();

    let mut best: Option<(Scalar, Point3, Point3, (Shape, Shape))> = None;
    let mut heap = BinaryHeap::new();
    let viable = |a: &Element, b: &Element| {
        critical_possible(a, b, metric) && critical_possible(b, a, metric)
    };
    for (i, ea) in elements_a.iter().enumerate() {
        for (j, eb) in elements_b.iter().enumerate() {
            if viable(ea, eb) {
                heap.push(Reverse((Key(lower_bound(&side_a, ea, &side_b, eb)?), i, j)));
                offer(&mut best, ea, eb);
            }
        }
    }

    let mut lower = 0.0;
    let mut steps = 0;
    let mut overlap = None;
    while let Some(Reverse((Key(bound), i, j))) = heap.pop() {
        lower = bound;
        let (ea, eb) = (elements_a[i], elements_b[j]);
        if metric == Metric::Plan && bound == 0.0 {
            if let Some((at, wa, wb)) = plan_patches_overlap(&side_a, &ea, &side_b, &eb)? {
                overlap = Some(at);
                best = Some((0.0, wa, wb, (ea.shape, eb.shape)));
                break;
            }
        }
        // Boundaries shown to meet within rounding: no lower bound can rise
        // above zero, so the search would only close on this (#273).
        if early_touch && metric == Metric::Space && bound == 0.0 {
            if let Some((wa, wb)) = touch::touch(&side_a, &ea, &side_b, &eb)? {
                let d = apart(wa, wb);
                if best.is_none_or(|(current, ..)| d < current) {
                    best = Some((d, wa, wb, (ea.shape, eb.shape)));
                }
                break;
            }
        }
        // Where the bound is positive and still short of what is in hand, the
        // boundary point of each element nearest the other's witness may close
        // the gap at once (#318); a pair at zero is the touch check's.
        if bound > 0.0 && best.is_none_or(|(current, ..)| bound < current) {
            if let Some((wa, wb, on)) = witness::pair(&side_a, &ea, &side_b, &eb)? {
                let d = apart(wa, wb);
                if best.is_none_or(|(current, ..)| d < current) {
                    best = Some((d, wa, wb, on));
                }
            }
        }
        let upper = best.map_or(Scalar::INFINITY, |(d, ..)| witnessed_upper(d, error));
        if upper.is_finite() && (bound >= upper || done(bound, upper, floor(scale, upper) + excess))
        {
            break;
        }
        steps += 1;
        if steps > max_steps {
            break;
        }
        // Refine the larger of the two; a pair that can shrink no further
        // holds the lower bound where it is.
        let split_a = side_a.size(&ea) >= side_b.size(&eb);
        let children = if split_a {
            side_a.split(&ea)?
        } else {
            side_b.split(&eb)?
        };
        let (side, parent) = if split_a {
            (&side_a, side_a.size(&ea))
        } else {
            (&side_b, side_b.size(&eb))
        };
        if children.iter().any(|child| side.size(child) >= parent) {
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
                    offer(&mut best, &child, &eb);
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
                    offer(&mut best, &ea, &child);
                }
            }
        }
    }
    // Every pair still queued bounds from below what is left unexplored.
    if let Some(Reverse((Key(bound), _, _))) = heap.peek() {
        lower = lower.min(*bound);
    }
    let (distance, point_a, point_b, on) = best.ok_or(crate::exact::NOT_CONVERGED)?;
    let upper = if overlap.is_some() {
        lower = 0.0;
        0.0
    } else {
        witnessed_upper(distance, error)
    };
    Ok(Found {
        bounds: DistanceBounds {
            lower: lower.min(upper),
            upper,
            point_a,
            point_b,
        },
        on,
        overlap,
        #[cfg(test)]
        steps,
    })
}

/// When both elements are patches certified inside planar faces that are
/// not vertical, and their shadows -- parallelograms -- share an open
/// patch: a plan point in it with a margin, and the boundary points above
/// it on each face.
fn plan_patches_overlap(
    side_a: &Side<'_>,
    a: &Element,
    side_b: &Side<'_>,
    b: &Element,
) -> Result<Option<(Point2, Point3, Point3)>, ExactMeasureError> {
    // The shadow of an inside patch of a planar face: its corners in plan,
    // counter-clockwise, and the plane's plan map to invert.
    let shadow = |side: &Side<'_>,
                  e: &Element|
     -> Result<Option<([Point2; 4], axiolid_core::Frame3)>, ExactMeasureError> {
        let Shape::Face {
            face,
            lo,
            hi,
            inside: true,
        } = e.shape
        else {
            return Ok(None);
        };
        let Some(f) = level_frame(side, face)? else {
            return Ok(None);
        };
        let map = |u: Scalar, v: Scalar| {
            let p = f.origin + f.x * u + f.y * v;
            Point2::new(p.x, p.y)
        };
        let mut corners = [
            map(lo.x, lo.y),
            map(hi.x, lo.y),
            map(hi.x, hi.y),
            map(lo.x, hi.y),
        ];
        let det = f.x.x * f.y.y - f.x.y * f.y.x;
        if (hi.x - lo.x) * (hi.y - lo.y) * det < 0.0 {
            corners.reverse();
        }
        Ok(Some((corners, f)))
    };
    let (Some((pa, fa)), Some((pb, fb))) = (shadow(side_a, a)?, shadow(side_b, b)?) else {
        return Ok(None);
    };
    // Clip one parallelogram by the other; the clip's centroid, if it
    // clears every edge of both by a margin, sits in an open patch of both.
    let cross =
        |o: Point2, p: Point2, q: Point2| (p.x - o.x) * (q.y - o.y) - (p.y - o.y) * (q.x - o.x);
    let mut poly: Vec<Point2> = pa.to_vec();
    for k in 0..4 {
        let (e0, e1) = (pb[k], pb[(k + 1) % 4]);
        let input = core::mem::take(&mut poly);
        for m in 0..input.len() {
            let (p, q) = (input[m], input[(m + 1) % input.len()]);
            let (sp, sq) = (cross(e0, e1, p), cross(e0, e1, q));
            if sp >= 0.0 {
                poly.push(p);
            }
            if (sp >= 0.0) != (sq >= 0.0) {
                let t = sp / (sp - sq);
                poly.push(p + (q - p) * t);
            }
        }
        if poly.is_empty() {
            return Ok(None);
        }
    }
    let at = poly.iter().fold(Point2::ZERO, |s, p| s + *p) / poly.len() as Scalar;
    let scale = pa
        .iter()
        .chain(pb.iter())
        .fold(1.0f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));
    let margin = 1e-9 * scale;
    let clears = |ring: &[Point2; 4]| {
        (0..4).all(|k| {
            let (e0, e1) = (ring[k], ring[(k + 1) % 4]);
            let length = (e1 - e0).length();
            length > 0.0 && cross(e0, e1, at) / length > margin
        })
    };
    if !(clears(&pa) && clears(&pb)) {
        return Ok(None);
    }
    // The points above `at` on each plane, for witnesses.
    let lift = |f: axiolid_core::Frame3| {
        let det = f.x.x * f.y.y - f.x.y * f.y.x;
        let (dx, dy) = (at.x - f.origin.x, at.y - f.origin.y);
        let u = (dx * f.y.y - dy * f.y.x) / det;
        let v = (f.x.x * dy - f.x.y * dx) / det;
        f.origin + f.x * u + f.y * v
    };
    Ok(Some((at, lift(fa), lift(fb))))
}

#[cfg(test)]
pub(crate) mod fixtures {
    //! Boxes built by hand, shared by the tests of the distance search.

    use axiolid_brep::{ExactBRep, ExactBRepBuilder};
    use axiolid_core::{Frame3, Interval, Point2, Point3, Vec2, Vec3};
    use axiolid_curve::{Curve2, Curve3, Line2, Line3};
    use axiolid_surface::{Plane, Surface};
    use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};

    /// A box from corner `o` along the edge vectors `x`, `y`, `z`
    /// (orthogonal), with twelve shared line edges and six planar faces.
    pub(crate) fn cuboid(o: Point3, x: Vec3, y: Vec3, z: Vec3) -> ExactBRep {
        let mut builder = ExactBRepBuilder::default();
        let axes = [x, y, z];
        let corner = |i: usize| {
            o + x * (i & 1) as f64 + y * ((i >> 1) & 1) as f64 + z * ((i >> 2) & 1) as f64
        };
        let vertices: Vec<_> = (0..8)
            .map(|i| {
                builder.topology_mut().add_vertex(Vertex {
                    position: corner(i),
                })
            })
            .collect();
        // Edge from corner `i` along axis `k` (bit `k` of `i` clear).
        let mut edges = std::collections::HashMap::new();
        for i in 0..8usize {
            for (k, axis) in axes.iter().enumerate() {
                if i & (1 << k) != 0 {
                    continue;
                }
                let length = axis.length();
                let curve = builder.add_curve3(Curve3::Line(Line3 {
                    origin: corner(i),
                    direction: *axis / length,
                }));
                let edge = builder.topology_mut().add_edge(Edge {
                    start: vertices[i],
                    end: vertices[i | (1 << k)],
                    curve: Some(curve),
                });
                builder.set_edge_interval(edge, Interval::new(0.0, length));
                edges.insert((i, i | (1 << k)), edge);
            }
        }
        // Each face: base corner and its two axes; the third axis fixed.
        for (fixed, side) in [(2, 0), (2, 1), (1, 0), (1, 1), (0, 0), (0, 1)] {
            let (ka, kb) = match fixed {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            let base = side << fixed;
            let (la, lb) = (axes[ka].length(), axes[kb].length());
            let (ux, uy) = (axes[ka] / la, axes[kb] / lb);
            let surface = builder.add_surface(Surface::Plane(Plane {
                frame: Frame3 {
                    origin: corner(base),
                    x: ux,
                    y: uy,
                    z: ux.cross(uy),
                },
            }));
            // Counter-clockwise in (u, v): base, +a, +a+b, +b.
            let ring = [
                base,
                base | (1 << ka),
                base | (1 << ka) | (1 << kb),
                base | (1 << kb),
            ];
            let params = [
                Point2::new(0.0, 0.0),
                Point2::new(la, 0.0),
                Point2::new(la, lb),
                Point2::new(0.0, lb),
            ];
            let mut uses = Vec::new();
            let mut intervals = Vec::new();
            for m in 0..4 {
                let (from, to) = (ring[m], ring[(m + 1) % 4]);
                let (edge, orientation) = match edges.get(&(from, to)) {
                    Some(edge) => (*edge, Orientation::Forward),
                    None => (edges[&(to, from)], Orientation::Reversed),
                };
                let (p, q) = (params[m], params[(m + 1) % 4]);
                let length = (q - p).length();
                let pcurve = builder.add_curve2(Curve2::Line(Line2 {
                    origin: p,
                    direction: Vec2::new((q.x - p.x) / length, (q.y - p.y) / length),
                }));
                uses.push(EdgeUse {
                    edge,
                    orientation,
                    pcurve: Some(pcurve),
                });
                intervals.push(Interval::new(0.0, length));
            }
            let loop_id = builder.topology_mut().add_loop(Loop { edges: uses });
            for (index, interval) in intervals.into_iter().enumerate() {
                builder.set_pcurve_interval(loop_id, index, interval);
            }
            builder.topology_mut().add_face(Face {
                surface: Some(surface),
                bounds: vec![FaceBound {
                    loop_id,
                    orientation: Orientation::Forward,
                    outer: true,
                }],
                orientation: Orientation::Forward,
            });
        }
        builder.finish().expect("a valid box")
    }

    pub(crate) fn unit_box(o: Point3) -> ExactBRep {
        cuboid(o, Vec3::X, Vec3::Y, Vec3::Z)
    }
}

#[cfg(test)]
mod tests {
    //! The two claims the lower bound rests on, checked by dense sampling:
    //! every point of a patch lies in its sphere, and every normal lies in
    //! its cone.

    use super::{
        critical_possible, flat_length, floor, normal_spread, pad, patch_plan_radius, patch_scale,
        patch_sphere, plan_lipschitz, witness_scale, witnessed_upper, Element, Metric, Shape, Side,
        PROJECTION_EPSILONS,
    };
    use axiolid_core::{Frame3, Point2, Point3, Vec3};
    use axiolid_evaluate::surface::partials;
    use axiolid_evaluate::surface::{evaluate, normal};
    use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface, Torus};

    /// A georeferenced shift (#318).
    const FAR: Vec3 = Vec3::new(450_000.0, 5_600_000.0, 0.0);

    fn frame_at(shift: Vec3) -> Frame3 {
        // Off the origin and turned, so no term vanishes by symmetry.
        let x = Vec3::new(0.6, 0.8, 0.0);
        let z = Vec3::new(0.0, 0.0, 1.0);
        Frame3 {
            origin: Point3::new(1.5, -2.0, 0.75) + shift,
            x,
            y: z.cross(x),
            z,
        }
    }

    fn families() -> Vec<(Surface, Point2, Point2)> {
        families_at(Vec3::ZERO)
    }

    fn families_at(shift: Vec3) -> Vec<(Surface, Point2, Point2)> {
        let f = frame_at(shift);
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

    /// Walls side by side at georeferenced coordinates, turned off the
    /// axes so that the witnesses round across the gap (#439): the upper
    /// bound covers the distance, and a limit just below it is never
    /// `Below`.
    #[test]
    fn a_far_clearance_is_never_below_a_limit_the_distance_reaches() {
        use super::fixtures::cuboid;
        use super::{boundary_clearance, boundary_distance, Clearance};
        use axiolid_core::Tolerance;
        let (x, y, z) = (Vec3::new(0.6, 0.8, 0.0), Vec3::new(-0.8, 0.6, 0.0), Vec3::Z);
        for k in 0..40_i32 {
            let o = Point3::new(
                450_000.0 + f64::from(k) * 0.0137,
                5_600_000.0 + f64::from(k) * 0.0291,
                0.0,
            );
            let a = cuboid(o, x * 4.0, y, z * 2.0);
            let ob = o
                + y * 1.5
                + x * (0.5 + f64::from(k % 3) * 0.75)
                + z * (0.25 + f64::from(k % 4) * 0.125);
            let b = cuboid(ob, x, y, z * 0.5);
            // b's -y face lies over a's +y face, parallel: both planes carry
            // the unit axis of `x * 4`, and the far differences are exact.
            let ux = (x * 4.0) / (x * 4.0).length();
            let n = Vec3::new(ux.y, -ux.x, 0.0);
            let truth = (ob - (o + y)).dot(n / n.length()).abs();
            let d = boundary_distance(&a, &b, 1e-9, Tolerance::METRE).expect("bounded");
            assert!(d.lower <= truth && truth <= d.upper, "{k}: {truth} {d:?}");
            assert!((d.point_a - d.point_b).length() <= d.upper, "{k}: {d:?}");
            let verdict = |limit: f64| {
                boundary_clearance(&a, &b, limit, Tolerance::METRE)
                    .expect("bounded")
                    .1
            };
            assert_eq!(
                verdict(truth - 1e-15),
                Clearance::Indeterminate,
                "{k}: {truth}"
            );
            assert_eq!(verdict(truth + 1e-6), Clearance::Below, "{k}");
            assert_eq!(verdict(truth - 1e-6), Clearance::Above, "{k}");
        }
    }

    /// A B-spline's points are summed over `p_u + p_v` levels of de Boor
    /// and divided by the weight: its witnesses' `M` counts the levels and
    /// the weights' spread, not only the largest control point (#439).
    #[test]
    fn a_spline_witness_scale_counts_its_degrees_and_weights() {
        use axiolid_curve::{BSplineSurface, KnotSpec};
        let net = |weights: Option<Vec<Vec<f64>>>| {
            Surface::BSpline(BSplineSurface {
                u_degree: 1,
                v_degree: 1,
                control_points: vec![
                    vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
                    vec![Point3::new(1.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
                ],
                u_knots: vec![0.0, 1.0],
                u_multiplicities: vec![2, 2],
                v_knots: vec![0.0, 1.0],
                v_multiplicities: vec![2, 2],
                weights,
                u_closed: false,
                v_closed: false,
                knot_spec: KnotSpec::Unspecified,
                self_intersect: None,
            })
        };
        let (lo, hi) = (Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
        assert_eq!(patch_scale(&net(None), lo, hi).expect("bounded"), 5.0);
        assert_eq!(witness_scale(&net(None), lo, hi).expect("bounded"), 15.0);
        let rational = net(Some(vec![vec![1.0, 2.0], vec![4.0, 1.0]]));
        assert_eq!(witness_scale(&rational, lo, hi).expect("bounded"), 60.0);
        // Every other family's witnesses are its patches' points.
        for (surface, lo, hi) in families_at(FAR) {
            assert_eq!(
                witness_scale(&surface, lo, hi).expect("bounded"),
                patch_scale(&surface, lo, hi).expect("bounded")
            );
        }
    }

    /// A rational bicubic net moved far away by an exact shift: every point
    /// evaluated there lies within the witness error of the same point
    /// evaluated at the origin and shifted (#439).
    #[test]
    fn far_spline_points_lie_within_their_witness_error() {
        use axiolid_curve::{BSplineSurface, KnotSpec};
        let net = |shift: Vec3| {
            let control_points = (0..4)
                .map(|i| {
                    (0..4)
                        .map(|j| {
                            let (u, v) = (f64::from(i) * 0.37, f64::from(j) * 0.29);
                            Point3::new(u, v, 0.1 * (u * 3.1).sin() * (v * 2.3).cos()) + shift
                        })
                        .collect()
                })
                .collect();
            let weights = (0..4)
                .map(|i| {
                    (0..4)
                        .map(|j| 1.0 + 0.25 * f64::from((i * 3 + j) % 5))
                        .collect()
                })
                .collect();
            Surface::BSpline(BSplineSurface {
                u_degree: 3,
                v_degree: 3,
                control_points,
                u_knots: vec![0.0, 1.0],
                u_multiplicities: vec![4, 4],
                v_knots: vec![0.0, 1.0],
                v_multiplicities: vec![4, 4],
                weights: Some(weights),
                u_closed: false,
                v_closed: false,
                knot_spec: KnotSpec::Unspecified,
                self_intersect: None,
            })
        };
        let (near, far) = (net(Vec3::ZERO), net(FAR));
        let (lo, hi) = (Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
        let error = 16.0 * f64::EPSILON * witness_scale(&far, lo, hi).expect("bounded");
        let mut worst: f64 = 0.0;
        for p in samples(lo, hi) {
            let there = evaluate(&far, p.x, p.y).expect("point");
            let here = evaluate(&near, p.x, p.y).expect("point") + FAR;
            worst = worst.max((there - here).length());
        }
        // `here` carries the rounding of the shift itself, half an ulp a coordinate.
        assert!(worst <= error, "{worst} > {error}");
        eprintln!("worst {worst:e} error {error:e}");
    }

    #[test]
    fn rounding_margins_are_the_derived_constants() {
        assert_eq!(pad(0.0, 1.0), 16.0 * f64::EPSILON + f64::MIN_POSITIVE);
        assert_eq!(floor(3.0, 0.0), 128.0 * f64::EPSILON * 3.0);
        assert_eq!(
            witnessed_upper(0.5, 1e-9),
            0.5 * (1.0 + 8.0 * f64::EPSILON) + 1e-9
        );
        let origin = Point3::new(450_000.0, 5_600_000.0, 0.0);
        let plane = Surface::Plane(Plane {
            frame: Frame3 {
                origin,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        });
        let scale = patch_scale(&plane, Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
        assert_eq!(scale.expect("bounded"), origin.length() + 2.0);
    }

    /// `M(child) <= M(parent)`: the pad of a side, taken from its root
    /// elements, covers every element a split yields.
    #[test]
    fn child_scales_never_exceed_their_parent() {
        for shift in [Vec3::ZERO, FAR] {
            for (surface, lo, hi) in families_at(shift) {
                let whole = patch_scale(&surface, lo, hi).expect("bounded");
                let mid = (lo + hi) * 0.5;
                for (a, b) in [
                    (lo, mid),
                    (Point2::new(mid.x, lo.y), Point2::new(hi.x, mid.y)),
                    (Point2::new(lo.x, mid.y), Point2::new(mid.x, hi.y)),
                    (mid, hi),
                ] {
                    let part = patch_scale(&surface, a, b).expect("bounded");
                    assert!(part <= whole, "{surface:?}: {part} > {whole}");
                }
            }
        }
    }

    #[test]
    fn every_far_patch_point_lies_in_its_sphere() {
        for (surface, lo, hi) in families_at(FAR) {
            let scale = patch_scale(&surface, lo, hi).expect("bounded");
            let (centre, radius) = patch_sphere(&surface, lo, hi, scale).expect("bounded");
            let mut reach: f64 = 0.0;
            for p in samples(lo, hi) {
                let point = evaluate(&surface, p.x, p.y).expect("point");
                reach = reach.max((point - centre).length());
            }
            assert!(
                reach <= radius,
                "{surface:?}: reach {reach} > radius {radius}"
            );
            assert!(
                radius <= 4.0 * reach + 1e-9,
                "{surface:?}: {radius} vs {reach}"
            );
        }
    }

    #[test]
    fn a_far_projection_carries_its_pad() {
        use super::fixtures::unit_box;
        let b = unit_box(Point3::ZERO + FAR);
        let side = Side::new(&b, 1e-12, Metric::Space).expect("a side");
        let top = side
            .elements
            .iter()
            .find(|e| matches!(e.shape, Shape::Face { .. }) && e.centre.z == 1.0)
            .expect("the top face");
        let (lo, hi) = side.project(top, Vec3::Z).expect("projects");
        let pad = PROJECTION_EPSILONS * f64::EPSILON * side.scale;
        assert!(pad > 1e-9, "{pad}");
        assert!(hi - lo >= 2.0 * pad - 4.0 * f64::EPSILON, "[{lo}, {hi}]");
    }

    #[test]
    fn every_patch_point_lies_in_its_sphere() {
        for (surface, lo, hi) in families() {
            let scale = patch_scale(&surface, lo, hi).expect("bounded");
            let (centre, radius) = patch_sphere(&surface, lo, hi, scale).expect("bounded");
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

    /// Every point of a patch projects within its plan radius of the
    /// centre's projection -- for the frame of `families` and for one
    /// tilted off the vertical, where every axis has a horizontal part --
    /// and the plan radius is not vacuous.
    #[test]
    fn every_patch_projection_lies_in_its_plan_disc() {
        let tilt = |surface: Surface| -> Surface {
            let x = Vec3::new(0.6, 0.64, 0.48);
            let z = Vec3::new(0.0, -0.6, 0.8);
            let y = z.cross(x);
            let f = Frame3 {
                origin: Point3::new(1.5, -2.0, 0.75),
                x,
                y,
                z,
            };
            match surface {
                Surface::Plane(_) => Surface::Plane(Plane { frame: f }),
                Surface::Cylinder(c) => Surface::Cylinder(Cylinder { frame: f, ..c }),
                Surface::EllipticalCylinder(c) => {
                    Surface::EllipticalCylinder(EllipticalCylinder { frame: f, ..c })
                }
                Surface::Cone(c) => Surface::Cone(Cone { frame: f, ..c }),
                Surface::Sphere(c) => Surface::Sphere(Sphere { frame: f, ..c }),
                Surface::Torus(c) => Surface::Torus(Torus { frame: f, ..c }),
                other => other,
            }
        };
        let cases = families()
            .into_iter()
            .flat_map(|(s, lo, hi)| [(s.clone(), lo, hi), (tilt(s), lo, hi)]);
        for (surface, lo, hi) in cases {
            let scale = patch_scale(&surface, lo, hi).expect("bounded");
            let (centre, radius) = patch_sphere(&surface, lo, hi, scale).expect("bounded");
            let plan = patch_plan_radius(&surface, lo, hi, centre, radius, scale).expect("bounded");
            assert!(plan <= radius);
            let mut reach: f64 = 0.0;
            for p in samples(lo, hi) {
                let point = evaluate(&surface, p.x, p.y).expect("point");
                reach = reach.max((point.x - centre.x).hypot(point.y - centre.y));
            }
            assert!(
                reach <= plan,
                "{surface:?}: reach {reach} > plan radius {plan}"
            );
            assert!(plan <= 4.0 * reach + 1e-9, "{surface:?}: {plan} vs {reach}");
        }
    }

    /// The horizontal speeds bound every horizontal derivative, over a whole
    /// turn, on a frame tilted so the axes' shadows are not orthogonal:
    /// there the steepest horizontal turn exceeds either axis's shadow.
    #[test]
    fn plan_speeds_bound_every_horizontal_derivative() {
        let x = Vec3::new(0.6, 0.64, 0.48);
        let z = Vec3::new(0.0, -0.6, 0.8);
        let f = Frame3 {
            origin: Point3::new(1.5, -2.0, 0.75),
            x,
            y: z.cross(x),
            z,
        };
        let (lo, hi) = (
            Point2::new(0.0, 0.2),
            Point2::new(core::f64::consts::TAU, 1.5),
        );
        for surface in [
            Surface::Cylinder(Cylinder {
                frame: f,
                radius: 2.5,
            }),
            Surface::EllipticalCylinder(EllipticalCylinder {
                frame: f,
                semi_axis_x: 3.0,
                semi_axis_y: 1.0,
            }),
            Surface::Cone(Cone {
                frame: f,
                radius: 1.5,
                semi_angle: 0.4,
            }),
        ] {
            let (lu, lv) = plan_lipschitz(&surface, lo, hi).expect("bounded");
            let (mut su, mut sv) = (0.0f64, 0.0f64);
            for i in 0..=720 {
                for j in 0..=8 {
                    let u = lo.x + (hi.x - lo.x) * f64::from(i) / 720.0;
                    let v = lo.y + (hi.y - lo.y) * f64::from(j) / 8.0;
                    let (du, dv) = partials(&surface, u, v).expect("partials");
                    su = su.max(flat_length(du));
                    sv = sv.max(flat_length(dv));
                }
            }
            assert!(su <= lu * (1.0 + 1e-12), "{surface:?}: {su} > {lu}");
            assert!(sv <= lv * (1.0 + 1e-12), "{surface:?}: {sv} > {lv}");
        }
    }

    /// A wall patch is kept for an element across from it in plan, however
    /// far above: only the plan offset has to meet its normal.
    #[test]
    fn a_wall_patch_faces_an_element_high_above_in_plan() {
        let element = |centre: Point3, normal: Option<Vec3>, spread: Option<f64>| Element {
            shape: Shape::Edge {
                edge: 0,
                t0: 0.0,
                t1: 1.0,
            },
            centre,
            radius: 0.1,
            plan_radius: 0.1,
            witness: None,
            normal,
            spread,
        };
        let wall = element(Point3::ZERO, Some(Vec3::X), Some(0.01));
        let above = element(Point3::new(10.0, 0.0, 10.0), None, None);
        assert!(critical_possible(&wall, &above, Metric::Plan));
        // In space the offset turns 45 degrees off the normal: dropped.
        assert!(!critical_possible(&wall, &above, Metric::Space));
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
