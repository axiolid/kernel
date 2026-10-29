//! Opt-in cut regions: the closed regions where the solid's material lies on
//! both sides of the cutting plane (ADR 0092, #372).
//!
//! A request made with
//! [`SectionRequest::with_cut_regions`](crate::SectionRequest::with_cut_regions)
//! also reports the cut as [`CutRegion`]s. The cut is the closure of the plane
//! points that have the solid's material on both sides: its interior meets the
//! plane there. In-plane faces are not part of it; they stay
//! [`InPlaneRegion`](crate::InPlaneRegion)s with their side. A region edge
//! bounds the cut where the off-plane triangle next to it leaves the plane on
//! the side opposite the region's material.
//!
//! A consumer therefore needs no topology of its own: a cut just below the
//! plane is the cut regions together with the `Below` regions, and a cut just
//! above the plane is the cut regions together with the `Above` regions. Their
//! rings run through bit-identical vertices.
//!
//! # Invariants
//!
//! - Only a [`SectionRole::Solid`](crate::SectionRole::Solid) has a cut; a
//!   sheet request with cut regions is invalid input.
//! - An outer ring is counter-clockwise and a hole is clockwise, so the
//!   material lies on the left. Every ring has at least three points, no
//!   adjacent duplicates (wrap-around included), no repeated vertex, and
//!   starts at its lexicographically smallest point (`x`, then `y`, by
//!   `total_cmp`). Holes are sorted by point sequence, regions by `outer`.
//! - On the exact geometry, rings never cross and have no repeated vertex, and
//!   regions have disjoint interiors. The returned points are those exact
//!   points computed in binary64, so rings that touch may cross or overlap by
//!   rounding error in the returned coordinates.
//! - A ring keeps a collinear vertex only where something else does: a
//!   junction of a chain and a region edge, a vertex that a region ring keeps,
//!   or a chain vertex that chain simplification keeps (a vertex of an
//!   in-plane triangle, or a node whose neighbours are not exactly collinear).
//!   Any other vertex between two collinear neighbours is dropped.
//! - Every closed chain of a solid is a cut ring, bit for bit. Every open
//!   chain with both ends on a region is a contiguous run of exactly one cut
//!   ring, forward or reversed. Every cut ring vertex is bit-identical to a
//!   chain point or a region ring vertex.
//! - A provider that reports a [`Both`](crate::InPlaneSide::Both) region
//!   refuses a cut-region request with
//!   [`REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES`](crate::REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES).
//!
//! # Exactness
//!
//! Which edges bound the cut and how they link come from exact source
//! predicates. Outer or hole, nesting, contact and crossing are decided
//! exactly on source data too: the points are source vertices in the plane and
//! exact crossings of source edges with the plane through `origin`,
//! `origin + x` and `origin + y`, never the returned coordinates. A self-
//! intersecting mesh is not certified (ADR 0033); detected crossings are
//! refused by name. Shells that overlap in the input, even by the rounding of
//! the consumer's own transform, are still refused.

use axiolid_core::Point2;

/// Refusal detail: the cut boundary touches itself at a vertex.
pub const REFUSED_PINCHED_CUT_REGION: &str = "cut region boundary pinches at a vertex";
/// Refusal detail: rings of the cut cross or overlap each other.
pub const REFUSED_CROSSING_CUT_RINGS: &str = "cut rings cross or overlap";

/// Hard bounds on the cut-region work of one request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct CutRegionLimits {
    /// Maximum output cut regions.
    pub max_regions: usize,
    /// Maximum vertices across all output cut region rings.
    pub max_region_vertices: usize,
    /// Maximum candidate edge pairs of the crossing sweep plus ring pairs of
    /// the containment test.
    pub max_pairs: usize,
}

impl CutRegionLimits {
    /// Construct explicit cut-region limits.
    pub const fn new(max_regions: usize, max_region_vertices: usize, max_pairs: usize) -> Self {
        Self {
            max_regions,
            max_region_vertices,
            max_pairs,
        }
    }
}

/// A planar region of the plane where the solid has material on both sides.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct CutRegion {
    /// Counter-clockwise outer ring.
    pub outer: Vec<Point2>,
    /// Clockwise hole rings.
    pub holes: Vec<Vec<Point2>>,
}

impl CutRegion {
    /// Construct a cut region.
    pub fn new(outer: Vec<Point2>, holes: Vec<Vec<Point2>>) -> Self {
        Self { outer, holes }
    }
}
