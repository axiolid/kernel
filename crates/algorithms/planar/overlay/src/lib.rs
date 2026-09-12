#![forbid(unsafe_code)]
//! Validated, deterministic planar boolean overlay and offset.
mod arc;
mod arc_overlay;
mod offset;
mod region;

pub use arc::{
    arc_edge_radius, arc_ring_area, reverse_arc_ring, validate_arc_ring, ArcRing, ArcVertex,
};
pub use arc_overlay::{arc_overlay, ArcOverlayEvidence, ArcOverlayResult, ArcPolygon};
pub use offset::{
    offset_polygons, polygon_area, ring_area, stroke_polyline, total_area, CapStyle, JoinStyle,
    OffsetEvidence, OffsetResult,
};
pub use region::{Region, RegionEvidence};

use axiolid_core::{Frame2, Point2, Tolerance};
use i_overlay::core::{fill_rule::FillRule as BackendFill, overlay_rule::OverlayRule};
use i_overlay::float::single::SingleFloatOverlay;

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
#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    pub outer: Ring,
    pub holes: Vec<Ring>,
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
    // Boundary contact is topology, not a repairable numerical nuisance.
    (ac.abs() <= epsilon || ad.abs() <= epsilon || ca.abs() <= epsilon || cb.abs() <= epsilon)
        || ((ac > 0.0) != (ad > 0.0) && (ca > 0.0) != (cb > 0.0))
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
            if !contains(&p.outer, h.points[0]) {
                return Err(OverlayError::HoleOutsideOuter);
            }
        }
    }
    Ok(())
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
fn backend(p: &[Polygon]) -> Vec<Vec<Vec<[f64; 2]>>> {
    p.iter()
        .map(|x| {
            core::iter::once(&x.outer)
                .chain(x.holes.iter())
                .map(|r| r.points.iter().map(|q| [q.x, q.y]).collect())
                .collect()
        })
        .collect()
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
    let rule = match operation {
        OverlayOperation::Intersection => OverlayRule::Intersect,
        OverlayOperation::Union => OverlayRule::Union,
        OverlayOperation::Difference => OverlayRule::Difference,
        OverlayOperation::Xor => OverlayRule::Xor,
    };
    let fill = match fill {
        FillRule::EvenOdd => BackendFill::EvenOdd,
        FillRule::NonZero => BackendFill::NonZero,
        FillRule::Positive => BackendFill::Positive,
        FillRule::Negative => BackendFill::Negative,
    };
    let shapes = backend(&subject.polygons).overlay(&backend(&clip.polygons), rule, fill);
    let polygons = shapes_to_polygons(shapes);
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
    // All rings go to the backend as one subject against an empty clip. The
    // NonZero fill then resolves the mutual overlaps in a single pass, which
    // is both correct and cheaper than folding pairwise.
    let subject: Vec<Vec<Vec<[f64; 2]>>> = rings
        .iter()
        .map(|ring| vec![ring.points.iter().map(|p| [p.x, p.y]).collect()])
        .collect();
    let empty: Vec<Vec<Vec<[f64; 2]>>> = Vec::new();
    let shapes = subject.overlay(&empty, OverlayRule::Union, BackendFill::NonZero);
    Ok(shapes_to_polygons(shapes))
}

/// Convert backend shapes to canonical kernel polygons.
///
/// Shared by every operation so ring orientation, rotation and polygon
/// ordering cannot drift between them.
fn shapes_to_polygons(shapes: Vec<Vec<Vec<[f64; 2]>>>) -> Vec<Polygon> {
    let mut polygons: Vec<Polygon> = shapes
        .into_iter()
        .filter_map(|shape| {
            let mut rings = shape.into_iter();
            let outer = rings.next()?;
            let outer = canonical(to_ring(outer), true);
            let holes = rings.map(|ring| canonical(to_ring(ring), false)).collect();
            Some(Polygon { outer, holes })
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

fn to_ring(points: Vec<[f64; 2]>) -> Ring {
    Ring {
        points: points
            .into_iter()
            .map(|p| Point2::new(p[0], p[1]))
            .collect(),
    }
}
