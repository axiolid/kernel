//! Planar offset: polygon inset/outset and polyline stroke offset (#42).
//!
//! # Boundary
//!
//! The kernel owns the geometric offset and its numeric contract. It does not
//! own *why* a distance was chosen. Clearance rules, code-mandated widths and
//! expansion direction are caller policy, so distance is always a parameter and
//! never a named constant here.
//!
//! # Sign convention
//!
//! A positive distance grows the region (outset), a negative distance shrinks
//! it (inset), and zero is the identity. This is fixed once here so consumers
//! do not each re-derive it from the backend's behaviour.
//!
//! # Collapse is a real answer
//!
//! Insetting further than a region's inradius removes it entirely. That returns
//! an empty result, never a degenerate ring: emitting a zero-area or
//! self-touching ring would hand the caller something that passes a ring-count
//! check while representing no region at all. `OffsetEvidence::collapsed`
//! reports it explicitly so a caller can distinguish "nothing left" from
//! "nothing given".
//!
//! # Validation is deliberately asymmetric
//!
//! Polygon input reuses the same validation as `overlay`, because offsetting a
//! self-intersecting or zero-area ring has no well-defined meaning. Polyline
//! input does NOT reject self-intersection: a stroke over a crossing path is
//! well defined — the crossing is resolved by the union of the swept region —
//! and rejecting it would refuse a case the operation genuinely handles.

use axiolid_core::{Point2, Tolerance};
use i_overlay::mesh::outline::offset::OutlineOffset;
use i_overlay::mesh::stroke::offset::StrokeOffset;
use i_overlay::mesh::style::{LineCap, LineJoin, OutlineStyle, StrokeStyle};

use crate::{canonical, validate_ring, OverlayError, Polygon, Ring};

/// How the outline turns a corner.
///
/// Named in kernel vocabulary rather than re-exported from the backend so the
/// choice stays a kernel contract if the backend is ever replaced.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum JoinStyle {
    /// Cut the corner off with a straight segment. Bounded by construction.
    Bevel,
    /// Extend the offset edges until they meet, limited by `angle_limit`
    /// radians: sharper corners than this fall back to a bevel.
    ///
    /// The limit is required rather than defaulted because an unlimited miter
    /// on a near-degenerate corner produces an arbitrarily distant spike.
    Miter { angle_limit: f64 },
    /// Approximate a circular arc, with `max_segment_ratio` bounding segment
    /// length over arc radius.
    Round { max_segment_ratio: f64 },
}

/// How an open stroke terminates.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum CapStyle {
    /// Stop flat at the endpoint.
    Butt,
    /// Extend flat by half the width past the endpoint.
    Square,
    /// Semicircular, with `max_segment_ratio` bounding the arc approximation.
    Round { max_segment_ratio: f64 },
}

/// What an offset actually did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct OffsetEvidence {
    /// Rings supplied by the caller.
    pub input_rings: usize,
    /// Polygons in the result.
    pub output_polygons: usize,
    /// Inner boundary components across all result polygons.
    pub output_holes: usize,
    /// The input was non-empty but the result is empty.
    ///
    /// Distinguishes a region inset out of existence from an empty input, which
    /// a bare empty vector cannot express.
    pub collapsed: bool,
}

/// Result of an offset.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct OffsetResult {
    /// Offset region, canonicalised exactly as `overlay` output is.
    pub polygons: Vec<Polygon>,
    /// What the operation did.
    pub evidence: OffsetEvidence,
}

impl JoinStyle {
    fn to_backend(self) -> Result<LineJoin<f64>, OverlayError> {
        match self {
            Self::Bevel => Ok(LineJoin::Bevel),
            Self::Miter { angle_limit } => {
                if !angle_limit.is_finite() || angle_limit <= 0.0 {
                    return Err(OverlayError::InvalidOffsetStyle);
                }
                Ok(LineJoin::Miter(angle_limit))
            }
            Self::Round { max_segment_ratio } => {
                if !max_segment_ratio.is_finite() || max_segment_ratio <= 0.0 {
                    return Err(OverlayError::InvalidOffsetStyle);
                }
                Ok(LineJoin::Round(max_segment_ratio))
            }
        }
    }
}

impl CapStyle {
    fn to_backend(self) -> Result<LineCap<[f64; 2], f64>, OverlayError> {
        match self {
            Self::Butt => Ok(LineCap::Butt),
            Self::Square => Ok(LineCap::Square),
            Self::Round { max_segment_ratio } => {
                if !max_segment_ratio.is_finite() || max_segment_ratio <= 0.0 {
                    return Err(OverlayError::InvalidOffsetStyle);
                }
                Ok(LineCap::Round(max_segment_ratio))
            }
        }
    }
}

/// Convert backend shapes into canonical kernel polygons.
///
/// Shares the winding and ordering convention with `overlay` so an offset
/// result and a boolean result are directly comparable.
fn to_kernel(shapes: Vec<Vec<Vec<[f64; 2]>>>) -> Vec<Polygon> {
    let mut polygons: Vec<Polygon> = shapes
        .into_iter()
        .filter_map(|shape| {
            let mut rings = shape.into_iter();
            let outer = rings.next()?;
            let outer = canonical(
                Ring {
                    points: outer.into_iter().map(|p| Point2::new(p[0], p[1])).collect(),
                },
                true,
            );
            // A backend ring with fewer than three points is not a region. It
            // is dropped rather than emitted, because a two-point "ring" would
            // satisfy a naive count check while bounding no area.
            if outer.points.len() < 3 {
                return None;
            }
            let holes = rings
                .filter_map(|ring| {
                    let ring = canonical(
                        Ring {
                            points: ring.into_iter().map(|p| Point2::new(p[0], p[1])).collect(),
                        },
                        false,
                    );
                    (ring.points.len() >= 3).then_some(ring)
                })
                .collect();
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

fn backend_shape(polygons: &[Polygon]) -> Vec<Vec<Vec<[f64; 2]>>> {
    polygons
        .iter()
        .map(|polygon| {
            core::iter::once(&polygon.outer)
                .chain(polygon.holes.iter())
                .map(|ring| ring.points.iter().map(|p| [p.x, p.y]).collect())
                .collect()
        })
        .collect()
}

/// Offset closed polygons by `distance`.
///
/// Positive grows, negative shrinks, zero is the identity. Holes are offset in
/// the opposite direction to the outer boundary automatically, which is what
/// makes an outset of a polygon-with-hole shrink the hole rather than grow it.
///
/// Returns [`OverlayError::ZeroArea`] and friends for malformed input via the
/// same ring validation `overlay` applies, so the two operations cannot
/// disagree about what a valid polygon is.
pub fn offset_polygons(
    polygons: &[Polygon],
    distance: f64,
    join: JoinStyle,
    tolerance: Tolerance,
) -> Result<OffsetResult, OverlayError> {
    if !distance.is_finite() {
        return Err(OverlayError::InvalidOffsetDistance);
    }
    for polygon in polygons {
        validate_ring(&polygon.outer, tolerance)?;
        for hole in &polygon.holes {
            validate_ring(hole, tolerance)?;
        }
    }
    let input_rings = polygons.iter().map(|p| 1 + p.holes.len()).sum();

    // Zero is the identity, and is handled without touching the backend so it
    // cannot pick up an incidental simplification pass.
    let result = if distance == 0.0 {
        to_kernel(backend_shape(polygons))
    } else {
        let style = OutlineStyle::new(distance).line_join(join.to_backend()?);
        crate::settle::settle(
            to_kernel(backend_shape(polygons).outline(&style)),
            tolerance,
        )
    };

    let evidence = OffsetEvidence {
        input_rings,
        output_polygons: result.len(),
        output_holes: result.iter().map(|p| p.holes.len()).sum(),
        collapsed: !polygons.is_empty() && result.is_empty(),
    };
    Ok(OffsetResult {
        polygons: result,
        evidence,
    })
}

/// Sweep an open polyline into a closed region of the given width.
///
/// `width` is the full stroke width, not a half-width: the region extends
/// `width / 2` either side of the path. Stating this explicitly matters because
/// both conventions are common and silently halving a clearance band is exactly
/// the kind of error this kernel refuses to make quietly.
///
/// Self-intersecting paths are accepted; the crossing is resolved by the union
/// of the swept region rather than rejected.
pub fn stroke_polyline(
    points: &[Point2],
    width: f64,
    join: JoinStyle,
    cap: CapStyle,
    closed: bool,
) -> Result<OffsetResult, OverlayError> {
    if points.len() < 2 {
        return Err(OverlayError::RingTooShort);
    }
    if !points.iter().all(|point| point.is_finite()) {
        return Err(OverlayError::NonFinitePoint);
    }
    if !width.is_finite() || width <= 0.0 {
        return Err(OverlayError::InvalidOffsetDistance);
    }

    let path: Vec<[f64; 2]> = points.iter().map(|p| [p.x, p.y]).collect();
    let style = StrokeStyle::new(width)
        .line_join(join.to_backend()?)
        .start_cap(cap.to_backend()?)
        .end_cap(cap.to_backend()?);

    let result = to_kernel(path.stroke(style, closed));
    let evidence = OffsetEvidence {
        input_rings: 1,
        output_polygons: result.len(),
        output_holes: result.iter().map(|p| p.holes.len()).sum(),
        collapsed: result.is_empty(),
    };
    Ok(OffsetResult {
        polygons: result,
        evidence,
    })
}

/// Absolute area of a ring.
///
/// Exposed because area monotonicity is the property callers most often want to
/// assert about an offset, and re-deriving the shoelace formula per consumer
/// invites sign-convention mistakes.
#[must_use]
pub fn ring_area(ring: &Ring) -> f64 {
    ring.points
        .iter()
        .zip(ring.points.iter().cycle().skip(1))
        .take(ring.points.len())
        .map(|(a, b)| a.x * b.y - b.x * a.y)
        .sum::<f64>()
        .abs()
        * 0.5
}

/// Net area of a polygon: outer boundary minus its holes.
#[must_use]
pub fn polygon_area(polygon: &Polygon) -> f64 {
    let holes: f64 = polygon.holes.iter().map(ring_area).sum();
    (ring_area(&polygon.outer) - holes).max(0.0)
}

/// Total net area across a set of polygons.
#[must_use]
pub fn total_area(polygons: &[Polygon]) -> f64 {
    polygons.iter().map(polygon_area).sum()
}
