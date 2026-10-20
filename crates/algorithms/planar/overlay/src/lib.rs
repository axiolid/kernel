#![forbid(unsafe_code)]
//! Validated, deterministic planar boolean overlay and offset.
//!
//! A straight ring's orientation and its `ZeroArea` check are exact signs
//! of its area taken as a fan of triangles from its first vertex (#274):
//! a shoelace over the coordinates as given rounds at the size of the
//! coordinates squared, which at georeferenced plan positions swamps the
//! area of a small ring. Areas reported as values are summed the same way,
//! so their rounding scales with the ring's extent, not its position.
//! `scripts/probe_ring_orientation_mutants.py` lists the faults the tests
//! in `tests/far_from_origin.rs` must catch.
//!
//! Whether two boundary edges meet, the `SelfIntersection` rule, is exact
//! too, and validation and settling share it (`contact`, #272): an end on
//! the other edge within the tolerance, or a proper crossing with each
//! segment's ends strictly on both sides of the other's line. So an
//! operation's settled output is a valid operand wherever it lies.
//!
//! A degenerate input ring refuses the operation by default. Asked
//! explicitly, [`overlay_with`] and [`union_soup_with`] leave out a ring
//! that is too short, has a repeated vertex, or has zero area under the
//! exact `ZeroArea` rule, and report each one with its area
//! (`degenerate`, #256); a self-intersecting ring is never dropped.
//! `scripts/probe_degenerate_drop_mutants.py` lists the faults the tests
//! in `tests/degenerate_drop.rs` must catch.
mod arc;
mod arc_overlay;
mod arrangement;
mod budget;
mod circle;
mod contact;
mod convex_parts;
mod cover;
mod degenerate;
mod exact_arc;
mod exact_overlay;
mod fill;
mod minkowski;
mod offset;
mod rectangle;
mod region;
mod segment_sweep;
mod settle;
mod visibility;

pub use arc::{
    arc_edge_radius, arc_ring_area, reverse_arc_ring, validate_arc_ring, ArcRing, ArcVertex,
};
pub use arc_overlay::{arc_overlay, ArcOverlayEvidence, ArcOverlayResult, ArcPolygon};
pub use arrangement::{
    ArcArrangement, ArrangementEdge, ArrangementRegion, EdgeSource, EdgeUse as ArrangementEdgeUse,
};
pub use budget::OverlayRefusal;
pub use circle::{
    minimum_enclosing_circle, CircleError, CircleEvidence, EnclosingCircle, MinimumCircle,
};
pub use degenerate::{
    DegenerateRings, DropReason, DroppedRing, OverlayOperand, OverlayOptions, ReportedOverlay,
    ReportedUnion, RingRole,
};
pub use minkowski::{BoundSide, MinkowskiError, MorphologyBound, MorphologyMode};
pub use offset::{
    offset_polygons, polygon_area, ring_area, stroke_polyline, total_area, CapStyle, JoinStyle,
    OffsetEvidence, OffsetResult,
};
pub use rectangle::{
    minimum_area_rectangle, MinimumRectangle, OrientedRectangle, RectangleError, RectangleEvidence,
};
pub use region::{Region, RegionEvidence};
pub use segment_sweep::{
    segment_intersections, ExactPoint2, Incidence, IntersectionPoint, SegmentIntersections,
    SegmentLocation, SegmentOverlap, SegmentSweepError, SweepEvidence,
};
pub use visibility::VisibilityError;

use std::borrow::Cow;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Polygon2, Tolerance};
use axiolid_exact::{certify, Arith, SignExpr};
use axiolid_guarantees::Sign;
use budget::Charge;
use contact::{on_segment, segments_intersect, side};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FillRule {
    EvenOdd,
    NonZero,
    Positive,
    Negative,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverlayOperation {
    Intersection,
    Union,
    Difference,
    Xor,
}
#[derive(Debug, Clone, PartialEq)]
pub struct Ring {
    pub points: Vec<Point2>,
}

/// A `Ring` is a closed boundary and so is [`Polygon2`]; converting between
/// them moves the points and nothing else.
///
/// The two exist separately because they are reached from different places:
/// `Polygon2` is a foundation value type usable without this crate, while
/// `Ring` is what the overlay consumes. Making them the same type
/// would drag the planar boolean vocabulary into `axiolid-core`.
impl From<Polygon2> for Ring {
    fn from(polygon: Polygon2) -> Self {
        Self {
            points: polygon.vertices,
        }
    }
}

impl From<Ring> for Polygon2 {
    fn from(ring: Ring) -> Self {
        Self::new(ring.points)
    }
}

impl From<&Ring> for Polygon2 {
    fn from(ring: &Ring) -> Self {
        Self::new(ring.points.clone())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    pub outer: Ring,
    pub holes: Vec<Ring>,
}

impl Polygon {
    /// The outer boundary as a [`Polygon2`], discarding any holes.
    ///
    /// Named `outline` rather than offered as a `From` impl because the
    /// conversion is lossy and the loss is silent: an annulus and a filled
    /// disc have the same outline, so a caller that reaches for this to
    /// compute area gets the wrong answer with no error. `Polygon2` models a
    /// simple polygon and cannot represent a hole, which is exactly why this
    /// has to be an explicit request rather than an implicit coercion.
    ///
    /// Use [`polygon_area`] when the holes matter.
    pub fn outline(&self) -> Polygon2 {
        Polygon2::from(&self.outer)
    }

    /// Whether the polygon has inner boundaries that [`outline`] would drop.
    ///
    /// [`outline`]: Polygon::outline
    pub fn has_holes(&self) -> bool {
        !self.holes.is_empty()
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct OverlayInput {
    pub frame: Frame2,
    pub polygons: Vec<Polygon>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OverlayError {
    InvalidFrame,
    NonFinitePoint,
    RingTooShort,
    RepeatedVertex,
    ZeroArea,
    /// Non-adjacent boundary segments meet or cross.
    SelfIntersection,
    HoleOutsideOuter,
    /// An offset distance or stroke width was not finite, or a width was not
    /// positive. A non-finite distance cannot produce a bounded region.
    InvalidOffsetDistance,
    /// A join or cap parameter was not a finite positive value.
    ///
    /// Separate from [`OverlayError::InvalidOffsetDistance`] because the fix is
    /// different: the caller passed a malformed style, not a malformed measure.
    InvalidOffsetStyle,
    /// An arc edge's implied radius was not above tolerance.
    ///
    /// Distinct from [`OverlayError::RepeatedVertex`]: the chord can be long
    /// enough while the bulge still implies a radius too small to be a real
    /// boundary. Such an edge is a point, not an arc.
    ZeroRadiusArc,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OverlayEvidence {
    pub subject_rings: usize,
    pub clip_rings: usize,
    pub output_polygons: usize,
    /// Number of inner boundary components across all result polygons.
    pub output_holes: usize,
}
#[derive(Debug, Clone, PartialEq)]
pub struct OverlayResult {
    pub polygons: Vec<Polygon>,
    pub evidence: OverlayEvidence,
}
/// A ring's signed area, positive counter-clockwise, as a measure.
///
/// Summed as triangles fanned from the first vertex, so each term is the
/// size of the ring's extent squared, not of its coordinates squared: a
/// shoelace over the coordinates as given loses a small ring's whole area
/// to rounding far from the origin (#274). Decisions use [`orientation`]
/// and [`zero_area`], which are exact; this is for magnitudes only.
fn signed(r: &Ring) -> f64 {
    let Some(&o) = r.points.first() else {
        return 0.0;
    };
    r.points
        .windows(2)
        .skip(1)
        .map(|w| (w[0] - o).perp_dot(w[1] - o))
        .sum::<f64>()
        * 0.5
}

/// The sign of `2 * area - offset` for the ring through `points`, exactly.
///
/// The area is the fan of triangles from the first vertex,
/// `sum((p[i] - p[0]) x (p[i + 1] - p[0]))`, which equals the shoelace in
/// exact arithmetic but keeps the interval filter tight far from the
/// origin (#274). Points must be finite; a refusal of the exact tier,
/// which finite points never cause, reads as zero.
struct TwiceArea<'a> {
    points: &'a [Point2],
    offset: f64,
}

impl SignExpr for TwiceArea<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let o = self.points.first()?;
        let local = |p: &Point2| (f(p.x).sub(&f(o.x)), f(p.y).sub(&f(o.y)));
        let mut sum = f(self.offset).neg();
        for w in self.points.windows(2).skip(1) {
            let ((ax, ay), (bx, by)) = (local(&w[0]), local(&w[1]));
            sum = sum.add(&ax.mul(&by).sub(&ay.mul(&bx)));
        }
        sum.sign()
    }
}

fn area_sign(points: &[Point2], offset: f64) -> Sign {
    certify(&TwiceArea { points, offset }).unwrap_or(Sign::Zero)
}

/// A ring's orientation, exactly: `Positive` counter-clockwise, `Negative`
/// clockwise, `Zero` when its signed area is exactly zero.
pub(crate) fn orientation(points: &[Point2]) -> Sign {
    area_sign(points, 0.0)
}

/// Whether a ring's area is at most `linear^2`, decided exactly against
/// that double (#274): the [`OverlayError::ZeroArea`] rule.
pub(crate) fn zero_area(points: &[Point2], linear: f64) -> bool {
    let bound = 2.0 * linear * linear;
    area_sign(points, bound) != Sign::Positive && area_sign(points, -bound) != Sign::Negative
}
fn self_intersects(r: &Ring, t: Tolerance) -> bool {
    first_touching_pair(&r.points, t.linear()).is_some()
}

/// Rings with at most this many edges compare every pair of them; longer
/// ones ask a box tree for the pairs that can meet.
const ALL_PAIRS_UP_TO: usize = 32;

/// The first pair of a ring's edges `(i, j)`, `i < j`, not neighbours,
/// that meet within `linear` ([`segments_intersect`]): the lowest `i`, then
/// the lowest `j`, as comparing every pair in that order finds it.
///
/// A long ring asks the crate's box tree (`exact_arc::boxes`, #341) for
/// the edges whose boxes, grown by `linear` and padded far beyond any
/// rounding of their own ([`exact_arc::Bounds::around`]), overlap edge
/// `i`'s. Two edges that meet within `linear` have boxes no farther apart
/// than `linear`, so every pair the exact test accepts is a candidate, and
/// the candidates come in ascending order: the pair found is the one the
/// all-pairs loop finds, at about `n log n` cost for a valid ring.
pub(crate) fn first_touching_pair(points: &[Point2], linear: f64) -> Option<(usize, usize)> {
    let n = points.len();
    let edge = |i: usize| (points[i], points[(i + 1) % n]);
    let meets = |i: usize, j: usize| {
        if j == i + 1 || (i == 0 && j + 1 == n) {
            return false;
        }
        let ((a, b), (c, d)) = (edge(i), edge(j));
        segments_intersect(a, b, c, d, linear)
    };
    if n <= ALL_PAIRS_UP_TO {
        return (0..n).find_map(|i| (i + 1..n).find(|&j| meets(i, j)).map(|j| (i, j)));
    }
    let boxes: Vec<exact_arc::Bounds> = (0..n)
        .map(|i| {
            let (a, b) = edge(i);
            exact_arc::Bounds::around(&[a, b], linear)
        })
        .collect();
    let tree = exact_arc::boxes::BoxTree::new(boxes.clone());
    let (mut candidates, mut stack) = (Vec::new(), Vec::new());
    for (i, own) in boxes.iter().enumerate() {
        tree.query_into(|b| b.overlaps(own), &mut candidates, &mut stack);
        if let Some(&j) = candidates.iter().find(|&&j| j > i && meets(i, j)) {
            return Some((i, j));
        }
    }
    None
}

/// Validate one ring: the contract every operand and every settled output
/// meets.
///
/// Refused, in this order: fewer than three points, a non-finite point,
/// an edge no longer than `t.linear()`, non-adjacent edges meeting within
/// it, and an area of at most `t.linear()^2` (`ZeroArea`). The area is
/// decided exactly, from the fan of triangles at the first vertex, so the
/// same ring decides the same at any distance from the origin (#274). A
/// ring whose vertices all lie exactly on one line is not refused: it
/// encloses nothing, and a boolean leaves it out (#219).
pub(crate) fn validate_ring(r: &Ring, t: Tolerance) -> Result<(), OverlayError> {
    #[cfg(test)]
    tests::VALIDATED.with(|count| count.set(count.get() + 1));
    ring_checks(r, t)
}

/// [`validate_ring`]'s checks, for debug assertions that must not count
/// as validating.
pub(crate) fn ring_checks(r: &Ring, t: Tolerance) -> Result<(), OverlayError> {
    if r.points.len() < 3 {
        return Err(OverlayError::RingTooShort);
    };
    if !r.points.iter().all(|p| p.is_finite()) {
        return Err(OverlayError::NonFinitePoint);
    };
    if r.points
        .iter()
        .zip(r.points.iter().cycle().skip(1))
        .take(r.points.len())
        .any(|(a, b)| (*a - *b).length() <= t.linear())
    {
        return Err(OverlayError::RepeatedVertex);
    };
    if self_intersects(r, t) {
        return Err(OverlayError::SelfIntersection);
    }
    // Both exact, so a ring decides the same wherever it lies (#274). A
    // ring whose vertices all lie on one line encloses nothing and a
    // boolean leaves it out (#219); it stays a valid operand, as the
    // shadow of a vertical face in a plan soup is.
    if zero_area(&r.points, t.linear()) && !exact_overlay::on_one_line(&r.points) {
        return Err(OverlayError::ZeroArea);
    };
    Ok(())
}
fn validate(input: &OverlayInput, t: Tolerance) -> Result<(), OverlayError> {
    validate_frame(input, t)?;
    for p in &input.polygons {
        validate_ring(&p.outer, t)?;
        for h in &p.holes {
            validate_ring(h, t)?;
            if hole_outside(&p.outer, h, t) {
                return Err(OverlayError::HoleOutsideOuter);
            }
        }
    }
    Ok(())
}
/// [`degenerate::admit`], or for an operand whose rings are `known` valid
/// and whose degenerate rings are refused, [`validate_known`] (#341).
fn admit_known<'a>(
    input: &'a OverlayInput,
    operand: OverlayOperand,
    tolerance: Tolerance,
    options: OverlayOptions,
    known: bool,
    dropped: &mut Vec<DroppedRing>,
) -> Result<Cow<'a, [Polygon]>, OverlayError> {
    if known && options.degenerate_rings == DegenerateRings::Refuse {
        validate_known(input, tolerance)?;
        return Ok(Cow::Borrowed(&input.polygons));
    }
    degenerate::admit(input, operand, tolerance, options, dropped)
}

/// [`validate`] for an input whose rings are known valid (#341): the
/// frame, and every hole inside its outer ring.
fn validate_known(input: &OverlayInput, t: Tolerance) -> Result<(), OverlayError> {
    validate_frame(input, t)?;
    for p in &input.polygons {
        if p.holes.iter().any(|h| hole_outside(&p.outer, h, t)) {
            return Err(OverlayError::HoleOutsideOuter);
        }
    }
    Ok(())
}
fn validate_frame(input: &OverlayInput, t: Tolerance) -> Result<(), OverlayError> {
    let f = input.frame;
    if !f.origin.is_finite()
        || !f.x.is_finite()
        || !f.y.is_finite()
        || (f.x.length() - 1.).abs() > t.linear()
        || (f.y.length() - 1.).abs() > t.linear()
        || f.x.dot(f.y).abs() > t.linear()
        || f.x.perp_dot(f.y) <= 0.
    {
        return Err(OverlayError::InvalidFrame);
    }
    Ok(())
}
/// Whether a hole leaves its outer ring: a vertex of it strictly outside.
/// Vertices on the outer boundary decide nothing -- a hole may touch its
/// outer ring at a vertex, which settled outputs do, and testing only the
/// first vertex called such a hole outside whenever that was the touching
/// one (axioval, #191).
fn hole_outside(outer: &Ring, hole: &Ring, t: Tolerance) -> bool {
    let n = outer.points.len();
    let on_boundary = |q: Point2| {
        (0..n).any(|i| {
            let (a, b) = (outer.points[i], outer.points[(i + 1) % n]);
            on_segment(a, b, q, t.linear())
        })
    };
    hole.points
        .iter()
        .any(|&q| !on_boundary(q) && !contains(outer, q))
}
fn contains(r: &Ring, p: Point2) -> bool {
    let mut inside = false;
    for (a, b) in r
        .points
        .iter()
        .zip(r.points.iter().cycle().skip(1))
        .take(r.points.len())
    {
        // The edge spans the ray's height and `p` lies strictly left of
        // where it crosses: left of an upward edge, right of a downward
        // one, decided exactly (#272). A point on the edge counts no
        // crossing, as the rounded comparison did.
        if (a.y > p.y) != (b.y > p.y) {
            let towards = if b.y > a.y {
                Sign::Positive
            } else {
                Sign::Negative
            };
            if side(*a, *b, p) == towards {
                inside = !inside
            }
        }
    }
    inside
}
pub(crate) fn canonical(mut r: Ring, want_positive: bool) -> Ring {
    if (orientation(&r.points) == Sign::Positive) != want_positive {
        r.points.reverse()
    };
    let k = r
        .points
        .iter()
        .enumerate()
        .min_by(|a, b| {
            a.1.x
                .total_cmp(&b.1.x)
                .then(a.1.y.total_cmp(&b.1.y))
                .then(a.0.cmp(&b.0))
        })
        .map(|x| x.0)
        .unwrap_or(0);
    r.points.rotate_left(k);
    r
}
/// Performs a neutral planar overlay. Output ordering is deterministic: polygons sort by outer-ring lexicographic start, rings are canonicalized CCW/CW.
///
/// Exact (#173): where boundaries cross and which side of each piece lies
/// in the result are exact decisions; the tolerance only validates the
/// operands and settles the output. An input vertex the operation does not
/// move comes back bit-identical, and a crossing of two edges is the double
/// nearest to the exact crossing point.
pub fn overlay(
    subject: &OverlayInput,
    clip: &OverlayInput,
    operation: OverlayOperation,
    fill: FillRule,
    tolerance: Tolerance,
) -> Result<OverlayResult, OverlayError> {
    overlay_with(
        subject,
        clip,
        operation,
        fill,
        tolerance,
        OverlayOptions::new(),
    )
    .map(|reported| reported.result)
}

/// [`overlay`] with explicit options; the default options make it
/// [`overlay`] exactly.
///
/// With [`OverlayOptions::dropping_degenerate_rings`] a ring that is too
/// short, has a repeated vertex, or has zero area under the exact
/// `ZeroArea` rule is left out instead of refusing the overlay, a hole
/// goes with its dropped outer ring, and each such ring is reported in
/// [`ReportedOverlay::dropped`] with its operand, polygon, ring, reason
/// and area (#256). Every other fault is refused as by [`overlay`]; a
/// self-intersecting ring is never dropped. The overlay of the kept rings
/// is exact as usual, and its evidence counts the kept rings. The rules
/// are in the docs of [`DropReason`] and of the `degenerate` module.
pub fn overlay_with(
    subject: &OverlayInput,
    clip: &OverlayInput,
    operation: OverlayOperation,
    fill: FillRule,
    tolerance: Tolerance,
    options: OverlayOptions,
) -> Result<ReportedOverlay, OverlayError> {
    overlay_charged(
        subject,
        clip,
        operation,
        fill,
        tolerance,
        options,
        &mut Charge::unlimited(),
    )
    .map_err(budget::unbudgeted)
}

/// [`overlay_with`] under the caller's execution options (#338): the
/// tolerance is [`ExecutionOptions::tolerance`], and with
/// [`ExecutionOptions::with_memory_budget`] the subdivision is refused
/// with [`OverlayRefusal::BudgetExceeded`] before it would hold more than
/// the budget. Without a budget it is [`overlay_with`] exactly, its
/// refusals wrapped in [`OverlayRefusal::Invalid`]. An admitted call gives
/// the same result as [`overlay_with`], bit for bit. No other option is
/// read: a cancellation token is not polled.
///
/// The charges and their bounds are in the docs of the `budget` module;
/// a budget covers the subdivision, not the operands or the output.
pub fn overlay_within(
    subject: &OverlayInput,
    clip: &OverlayInput,
    operation: OverlayOperation,
    fill: FillRule,
    options: OverlayOptions,
    execution: &ExecutionOptions,
) -> Result<ReportedOverlay, OverlayRefusal> {
    overlay_charged(
        subject,
        clip,
        operation,
        fill,
        execution.tolerance(),
        options,
        &mut Charge::within(execution.memory_budget_bytes()),
    )
}

/// [`overlay_with`], charged against `charge`.
fn overlay_charged(
    subject: &OverlayInput,
    clip: &OverlayInput,
    operation: OverlayOperation,
    fill: FillRule,
    tolerance: Tolerance,
    options: OverlayOptions,
    charge: &mut Charge,
) -> Result<ReportedOverlay, OverlayRefusal> {
    let known = [false, false];
    overlay_known(
        subject, clip, operation, fill, tolerance, options, charge, known,
    )
    .map(|(reported, _)| reported)
}

/// [`overlay_charged`] for operands some of which are known valid, and the
/// proof that the output is (#341).
///
/// An operand with `known` set, its rings proved valid under `tolerance`
/// by a [`settle::Validated`] the caller holds, has its frame and its
/// holes' places checked but no ring validated again; refused
/// degenerate rings are the only mode it applies to, the others validate
/// as usual. The proof returned is the settled output's.
#[allow(clippy::too_many_arguments)]
pub(crate) fn overlay_known(
    subject: &OverlayInput,
    clip: &OverlayInput,
    operation: OverlayOperation,
    fill: FillRule,
    tolerance: Tolerance,
    options: OverlayOptions,
    charge: &mut Charge,
    known: [bool; 2],
) -> Result<(ReportedOverlay, Option<settle::Validated>), OverlayRefusal> {
    let vertices = |input: &OverlayInput| -> usize {
        input
            .polygons
            .iter()
            .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
            .map(|ring| ring.points.len())
            .sum()
    };
    charge.add(vertices(subject) + vertices(clip), budget::PER_VERTEX)?;
    let mut dropped = Vec::new();
    let subject_polygons = admit_known(
        subject,
        OverlayOperand::Subject,
        tolerance,
        options,
        known[0],
        &mut dropped,
    )?;
    let clip_polygons = admit_known(
        clip,
        OverlayOperand::Clip,
        tolerance,
        options,
        known[1],
        &mut dropped,
    )?;
    if subject.frame != clip.frame {
        return Err(OverlayError::InvalidFrame.into());
    };
    let rings = exact_overlay::boolean(
        &subject_polygons,
        &clip_polygons,
        operation,
        fill,
        tolerance.linear(),
        charge,
    )?;
    let (polygons, validated) = settle::settle_validated(canonical_polygons(rings), tolerance);
    let evidence = OverlayEvidence {
        subject_rings: subject_polygons.iter().map(|p| 1 + p.holes.len()).sum(),
        clip_rings: clip_polygons.iter().map(|p| 1 + p.holes.len()).sum(),
        output_polygons: polygons.len(),
        output_holes: polygons.iter().map(|polygon| polygon.holes.len()).sum(),
    };
    Ok((
        ReportedOverlay {
            result: OverlayResult { polygons, evidence },
            dropped,
        },
        validated,
    ))
}

/// Union a set of individually-valid rings that may overlap each other.
///
/// [`overlay`] rejects a self-intersecting *operand*, which is correct for a
/// boolean between two shapes but wrong for a union of a triangle soup: a
/// projected mesh routinely overlaps itself, and mutual overlap is exactly
/// what a union is for. Each ring is still validated on its own, so malformed
/// geometry is refused rather than absorbed.
pub fn union_soup(rings: &[Ring], tolerance: Tolerance) -> Result<Vec<Polygon>, OverlayError> {
    for ring in rings {
        validate_ring(ring, tolerance)?;
    }
    union_valid(rings, tolerance, &mut Charge::unlimited()).map_err(budget::unbudgeted)
}

/// [`union_soup`] with explicit options; the default options make it
/// [`union_soup`] exactly.
///
/// With [`OverlayOptions::dropping_degenerate_rings`] a ring that is too
/// short, has a repeated vertex, or has zero area under the exact
/// `ZeroArea` rule is left out and reported in [`ReportedUnion::dropped`],
/// as [`overlay_with`] does (#256): operand
/// [`OverlayOperand::Subject`], polygon the ring's index, role
/// [`RingRole::Outer`]. Every other fault is refused as by
/// [`union_soup`].
pub fn union_soup_with(
    rings: &[Ring],
    tolerance: Tolerance,
    options: OverlayOptions,
) -> Result<ReportedUnion, OverlayError> {
    union_soup_charged(rings, tolerance, options, &mut Charge::unlimited())
        .map_err(budget::unbudgeted)
}

/// [`union_soup_with`] under the caller's execution options (#338), as
/// [`overlay_within`] is [`overlay_with`]: the tolerance is
/// [`ExecutionOptions::tolerance`], and a memory budget refuses the
/// subdivision with [`OverlayRefusal::BudgetExceeded`] before it would
/// hold more. Without a budget it is [`union_soup_with`] exactly.
pub fn union_soup_within(
    rings: &[Ring],
    options: OverlayOptions,
    execution: &ExecutionOptions,
) -> Result<ReportedUnion, OverlayRefusal> {
    union_soup_charged(
        rings,
        execution.tolerance(),
        options,
        &mut Charge::within(execution.memory_budget_bytes()),
    )
}

/// [`union_soup_with`], charged against `charge`.
fn union_soup_charged(
    rings: &[Ring],
    tolerance: Tolerance,
    options: OverlayOptions,
    charge: &mut Charge,
) -> Result<ReportedUnion, OverlayRefusal> {
    charge.add(
        rings.iter().map(|ring| ring.points.len()).sum(),
        budget::PER_VERTEX,
    )?;
    let mut dropped = Vec::new();
    let mut kept = Vec::with_capacity(rings.len());
    for (index, ring) in rings.iter().enumerate() {
        if options.degenerate_rings == DegenerateRings::Refuse {
            validate_ring(ring, tolerance)?;
            kept.push(ring.clone());
            continue;
        }
        match degenerate::checked(ring, tolerance)? {
            None => kept.push(ring.clone()),
            Some(reason) => {
                let place = degenerate::Place {
                    operand: OverlayOperand::Subject,
                    polygon: index,
                    ring: RingRole::Outer,
                };
                dropped.push(degenerate::report(ring, place, reason));
            }
        }
    }
    Ok(ReportedUnion {
        polygons: union_valid(&kept, tolerance, charge)?,
        dropped,
    })
}

/// The union of rings already validated one by one.
fn union_valid(
    rings: &[Ring],
    tolerance: Tolerance,
    charge: &mut Charge,
) -> Result<Vec<Polygon>, OverlayRefusal> {
    if rings.is_empty() {
        return Ok(Vec::new());
    }
    // All rings form one subject against an empty clip. The NonZero fill
    // then resolves the mutual overlaps in a single pass, which is both
    // correct and cheaper than folding pairwise.
    let subject: Vec<Polygon> = rings
        .iter()
        .map(|ring| Polygon {
            outer: ring.clone(),
            holes: Vec::new(),
        })
        .collect();
    let rings = exact_overlay::boolean(
        &subject,
        &[],
        OverlayOperation::Union,
        FillRule::NonZero,
        tolerance.linear(),
        charge,
    )?;
    // Settled like every other output, so the polygons are valid operands
    // (#191).
    Ok(settle::settle(canonical_polygons(rings), tolerance))
}

/// Canonical kernel polygons from a boolean's rings.
///
/// Shared by every operation so ring orientation, rotation and polygon
/// ordering cannot drift between them.
fn canonical_polygons(rings: Vec<(Ring, Vec<Ring>)>) -> Vec<Polygon> {
    let mut polygons: Vec<Polygon> = rings
        .into_iter()
        .map(|(outer, holes)| Polygon {
            outer: canonical(outer, true),
            holes: holes.into_iter().map(|h| canonical(h, false)).collect(),
        })
        .collect();
    polygons.sort_by(|a, b| {
        a.outer.points[0]
            .x
            .total_cmp(&b.outer.points[0].x)
            .then(a.outer.points[0].y.total_cmp(&b.outer.points[0].y))
    });
    polygons
}

#[cfg(test)]
mod tests {
    //! Ring validation's candidate search and the proof regions carry
    //! (#341).
    use super::*;
    use std::cell::Cell;

    thread_local! {
        /// Rings validated on this thread.
        pub(super) static VALIDATED: Cell<usize> = const { Cell::new(0) };
    }

    fn validated() -> usize {
        VALIDATED.with(Cell::get)
    }

    /// The search before #341: every pair, lowest `i` then lowest `j`.
    fn all_pairs(points: &[Point2], linear: f64) -> Option<(usize, usize)> {
        let n = points.len();
        for i in 0..n {
            for j in i + 1..n {
                if j == i + 1 || (i == 0 && j + 1 == n) {
                    continue;
                }
                let (a, b) = (points[i], points[(i + 1) % n]);
                let (c, d) = (points[j], points[(j + 1) % n]);
                if segments_intersect(a, b, c, d, linear) {
                    return Some((i, j));
                }
            }
        }
        None
    }

    /// A small deterministic generator, so failures reproduce.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// A star of `n` corners around the origin at seeded radii, with
    /// `faults` corners moved through the centre, onto another corner or
    /// within the tolerance of a far edge, so that some stars cross or
    /// touch themselves and others do not.
    fn star(rng: &mut Lcg, n: usize, faults: usize) -> Vec<Point2> {
        let mut points: Vec<Point2> = (0..n)
            .map(|i| {
                let angle = std::f64::consts::TAU * (i as f64 + 0.5 * rng.next()) / n as f64;
                let radius = 1.0 + 0.2 * rng.next();
                Point2::new(radius * angle.cos(), radius * angle.sin())
            })
            .collect();
        for _ in 0..faults {
            let k = (rng.next() * n as f64) as usize % n;
            points[k] = match (rng.next() * 3.0) as usize {
                0 => Point2::new(-0.9 * points[k].x, -0.9 * points[k].y),
                1 => points[(k + 2) % n],
                _ => {
                    let j = (k + n / 2) % n;
                    let (a, b) = (points[j], points[(j + 1) % n]);
                    Point2::new(0.5 * (a.x + b.x) + 1e-10, 0.5 * (a.y + b.y))
                }
            };
        }
        points
    }

    /// The box tree finds the pair the all-pairs loop finds, on rings long
    /// enough to use it, valid or crossing, touching or within the
    /// tolerance.
    #[test]
    fn the_first_touching_pair_is_the_all_pairs_one() {
        let mut rng = Lcg(341);
        let (mut found, mut clean) = (0, 0);
        for case in 0..400 {
            let n = [5, 33, 40, 64, 200][case % 5];
            let faults = case / 5 % 4;
            let points = star(&mut rng, n, faults);
            for linear in [1e-9, 1e-3, 0.05] {
                let pair = first_touching_pair(&points, linear);
                assert_eq!(pair, all_pairs(&points, linear), "case {case} at {linear}");
                match pair {
                    Some(_) => found += 1,
                    None => clean += 1,
                }
            }
        }
        assert!(found > 200 && clean > 200, "{found} found, {clean} clean");
    }

    /// A comb whose last corner comes within the tolerance of the first
    /// tooth's flank: a pair far apart in the ring, near in the plane.
    #[test]
    fn a_long_ring_touching_itself_far_along_is_found() {
        let mut points = Vec::new();
        for k in 0..2_000 {
            let x = k as f64;
            points.push(Point2::new(x, 0.0));
            points.push(Point2::new(x + 0.5, 1.0));
        }
        points.push(Point2::new(2_000.0, 0.0));
        points.push(Point2::new(2_000.0, -1.0));
        points.push(Point2::new(0.0, -1.0));
        let t = Tolerance::new(1e-9, 1e-9).expect("a tolerance");
        let ring = Ring {
            points: points.clone(),
        };
        assert_eq!(validate_ring(&ring, t), Ok(()));
        assert_eq!(first_touching_pair(&points, t.linear()), None);
        let n = points.len();
        points[n - 1] = Point2::new(0.25, 0.5 - 1e-10);
        let pair = first_touching_pair(&points, t.linear());
        assert!(pair.is_some());
        assert_eq!(pair, all_pairs(&points, t.linear()));
        assert_eq!(
            validate_ring(&Ring { points }, t),
            Err(OverlayError::SelfIntersection)
        );
    }

    fn square(x: f64, y: f64, side: f64) -> Polygon {
        Polygon {
            outer: Ring {
                points: vec![
                    Point2::new(x, y),
                    Point2::new(x + side, y),
                    Point2::new(x + side, y + side),
                    Point2::new(x, y + side),
                ],
            },
            holes: Vec::new(),
        }
    }

    /// Booleans on regions validate no ring a region carries a proof for
    /// under the call's tolerance, and the results carry one in turn; a
    /// call under another tolerance, or on a region without a proof,
    /// validates again and refuses what it refuses today.
    #[test]
    fn regions_are_validated_once() {
        let t = Tolerance::new(1e-9, 1e-9).expect("a tolerance");
        let a = Region::new(vec![square(0.0, 0.0, 2.0)], t).expect("valid");
        let b = Region::new(vec![square(1.0, 1.0, 2.0)], t).expect("valid");
        let before = validated();
        let union = a.union(&b, t).expect("the union");
        let cut = union.difference(&a, t).expect("the difference");
        let again = cut.intersection(&b, t).expect("the intersection");
        assert_eq!(validated(), before, "no ring validated again");
        assert!(!again.is_empty());
        let input = |region: &Region| OverlayInput {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Point2::X,
                y: Point2::Y,
            },
            polygons: region.polygons().to_vec(),
        };
        let raw = overlay(
            &input(&a),
            &input(&b),
            OverlayOperation::Union,
            FillRule::NonZero,
            t,
        )
        .expect("the union");
        assert_eq!(raw.polygons, union.polygons());
        assert!(validated() > before, "overlay validates every ring");

        // A ring with an edge of 1e-6 is valid at 1e-9, not at 1e-3.
        let mut short = square(5.0, 5.0, 1.0);
        short.outer.points.insert(1, Point2::new(5.5, 5.0));
        short.outer.points.insert(2, Point2::new(5.5 + 1e-6, 5.0));
        let fine = Region::new(vec![short], t).expect("valid at 1e-9");
        let coarse = Tolerance::new(1e-3, 1e-3).expect("a tolerance");
        assert_eq!(fine.union(&a, coarse), Err(OverlayError::RepeatedVertex));
        // A translated region carries no proof.
        let moved = a
            .translate(axiolid_core::Vec2::new(0.5, 0.5))
            .expect("moved");
        let before = validated();
        moved.union(&b, t).expect("the union");
        assert!(validated() > before, "a translate is validated again");

        // A disc morphology's settled result carries a proof too.
        let eroded = union.erode_inner(0.1, t).expect("the erosion");
        let before = validated();
        eroded.difference(&b, t).expect("the difference");
        assert_eq!(validated(), before, "an erosion is not validated again");

        // A hole outside its outer ring is valid ring by ring, so a region
        // holds it; a boolean refuses it, proof or not.
        let mut stray = square(10.0, 10.0, 1.0);
        let mut hole = square(20.0, 20.0, 0.5).outer;
        hole.points.reverse();
        stray.holes.push(hole);
        let stray = Region::new(vec![stray], t).expect("each ring valid");
        assert_eq!(stray.union(&a, t), Err(OverlayError::HoleOutsideOuter));
    }
}
