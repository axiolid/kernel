#![forbid(unsafe_code)]
//! Validated, deterministic planar boolean overlay and offset.
mod arc;
mod arc_overlay;
mod arrangement;
mod circle;
mod exact_arc;
mod exact_overlay;
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
pub use circle::{
    minimum_enclosing_circle, CircleError, CircleEvidence, EnclosingCircle, MinimumCircle,
};
pub use minkowski::{BoundSide, MinkowskiError, MorphologyBound};
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

use axiolid_core::{Frame2, Point2, Polygon2, Tolerance};

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
fn signed(r: &Ring) -> f64 {
    r.points
        .iter()
        .zip(r.points.iter().cycle().skip(1))
        .take(r.points.len())
        .map(|(a, b)| a.x * b.y - b.x * a.y)
        .sum::<f64>()
        * 0.5
}
fn cross(a: Point2, b: Point2, c: Point2) -> f64 {
    (b - a).perp_dot(c - a)
}

fn segments_intersect(a: Point2, b: Point2, c: Point2, d: Point2, epsilon: f64) -> bool {
    let ac = cross(a, b, c);
    let ad = cross(a, b, d);
    let ca = cross(c, d, a);
    let cb = cross(c, d, b);
    // Boundary contact is topology, not a repairable numerical nuisance. A
    // touching endpoint must lie ON the other segment, not merely on the
    // infinite line through it: two collinear edges of a U-shape or a comb
    // share a line without touching.
    (ac.abs() <= epsilon && within_extent(a, b, c, epsilon))
        || (ad.abs() <= epsilon && within_extent(a, b, d, epsilon))
        || (ca.abs() <= epsilon && within_extent(c, d, a, epsilon))
        || (cb.abs() <= epsilon && within_extent(c, d, b, epsilon))
        || ((ac > 0.0) != (ad > 0.0) && (ca > 0.0) != (cb > 0.0))
}

/// Whether `p` lies inside the axis-aligned box of segment `a`-`b`, widened
/// by `epsilon`. Combined with collinearity this puts `p` on the segment.
fn within_extent(a: Point2, b: Point2, p: Point2, epsilon: f64) -> bool {
    p.x >= a.x.min(b.x) - epsilon
        && p.x <= a.x.max(b.x) + epsilon
        && p.y >= a.y.min(b.y) - epsilon
        && p.y <= a.y.max(b.y) + epsilon
}

fn self_intersects(r: &Ring, t: Tolerance) -> bool {
    let n = r.points.len();
    for i in 0..n {
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j + 1 == n) {
                continue;
            }
            if segments_intersect(
                r.points[i],
                r.points[(i + 1) % n],
                r.points[j],
                r.points[(j + 1) % n],
                t.linear(),
            ) {
                return true;
            }
        }
    }
    false
}

pub(crate) fn validate_ring(r: &Ring, t: Tolerance) -> Result<(), OverlayError> {
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
    if signed(r).abs() <= t.linear().powi(2) {
        return Err(OverlayError::ZeroArea);
    };
    Ok(())
}
fn validate(input: &OverlayInput, t: Tolerance) -> Result<(), OverlayError> {
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
            cross(a, b, q).abs() <= t.linear() && within_extent(a, b, q, t.linear())
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
        if (a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x {
            inside = !inside
        }
    }
    inside
}
pub(crate) fn canonical(mut r: Ring, want_positive: bool) -> Ring {
    if (signed(&r) > 0.) != want_positive {
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
    validate(subject, tolerance)?;
    validate(clip, tolerance)?;
    if subject.frame != clip.frame {
        return Err(OverlayError::InvalidFrame);
    };
    let rings = exact_overlay::boolean(&subject.polygons, &clip.polygons, operation, fill)?;
    let polygons = settle::settle(canonical_polygons(rings), tolerance);
    let evidence = OverlayEvidence {
        subject_rings: subject.polygons.iter().map(|p| 1 + p.holes.len()).sum(),
        clip_rings: clip.polygons.iter().map(|p| 1 + p.holes.len()).sum(),
        output_polygons: polygons.len(),
        output_holes: polygons.iter().map(|polygon| polygon.holes.len()).sum(),
    };
    Ok(OverlayResult { polygons, evidence })
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
    let rings = exact_overlay::boolean(&subject, &[], OverlayOperation::Union, FillRule::NonZero)?;
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
