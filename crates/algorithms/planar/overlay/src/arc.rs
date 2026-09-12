//! Arc-capable planar boundaries (ADR 0050, step 1).
//!
//! # Why a second ring type
//!
//! [`Ring`](crate::Ring) is a point list, so every edge between
//! consecutive points is a straight segment. A cylinder cross-section is a
//! circle, and approximating it by segments would destroy the exactness
//! that the exact boolean path exists to provide.
//!
//! [`ArcRing`] instead stores a bulge per edge, the DXF convention:
//! `bulge = tan(theta / 4)` where `theta` is the signed included angle of
//! the arc from one vertex to the next. `bulge == 0` is exactly a straight
//! segment, so a polygonal boundary is representable without a special
//! case, and a circle needs only two vertices.
//!
//! The sign carries direction: a positive bulge turns counter-clockwise
//! (bulging left of the chord), a negative one clockwise. That is what
//! makes a bite and a bump distinguishable on the same chord.
//!
//! # What this module does NOT do
//!
//! No boolean, no backend, no conversion to any third-party type. This is
//! the neutral contract and its validation only; ADR 0050 keeps the
//! adapter in a later step so that library types never reach this API.

use axiolid_core::{Point2, Tolerance};

use crate::OverlayError;

/// One boundary vertex and the bulge of the edge leaving it.
///
/// The bulge describes the edge from this vertex to the NEXT one, so a
/// ring of `n` vertices has exactly `n` edges and needs no separate edge
/// list. Storing it per departing vertex keeps insertion and reversal
/// local operations rather than re-indexing an edge array.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArcVertex {
    /// Vertex position.
    pub point: Point2,
    /// `tan(theta / 4)` for the edge leaving `point`; `0` means straight.
    pub bulge: f64,
}

impl ArcVertex {
    /// A vertex whose departing edge is a straight segment.
    pub const fn straight(point: Point2) -> Self {
        Self { point, bulge: 0.0 }
    }

    /// A vertex whose departing edge bulges by `bulge`.
    pub const fn bulged(point: Point2, bulge: f64) -> Self {
        Self { point, bulge }
    }

    /// True when the departing edge is straight.
    ///
    /// Exact comparison is deliberate: a bulge is either exactly zero, in
    /// which case the edge IS a segment, or it is an arc whose radius is
    /// whatever the value implies. Treating tiny bulges as straight would
    /// silently discard curvature the caller asked for.
    pub fn is_straight(&self) -> bool {
        self.bulge == 0.0
    }
}

/// A closed boundary whose edges may be arcs.
///
/// Always closed: the edge leaving the last vertex returns to the first.
/// An explicitly repeated closing vertex is therefore a zero-length edge
/// and is refused by [`validate_arc_ring`].
#[derive(Debug, Clone, PartialEq)]
pub struct ArcRing {
    /// Boundary vertices in order.
    pub vertices: Vec<ArcVertex>,
}

impl ArcRing {
    /// Build a ring from vertices.
    pub fn new(vertices: Vec<ArcVertex>) -> Self {
        Self { vertices }
    }

    /// Build a straight-edged ring, the polygonal case.
    pub fn from_points(points: &[Point2]) -> Self {
        Self {
            vertices: points.iter().copied().map(ArcVertex::straight).collect(),
        }
    }

    /// A full circle as two semicircular edges.
    ///
    /// Two vertices is the minimum honest encoding of a circle: one vertex
    /// would make the chord degenerate and the centre ambiguous, which is
    /// why [`validate_arc_ring`] requires at least two.
    pub fn circle(centre: Point2, radius: f64) -> Self {
        let left = Point2::new(centre.x - radius, centre.y);
        let right = Point2::new(centre.x + radius, centre.y);
        Self {
            vertices: vec![ArcVertex::bulged(left, 1.0), ArcVertex::bulged(right, 1.0)],
        }
    }

    /// Number of edges, which equals the number of vertices.
    pub fn edge_count(&self) -> usize {
        self.vertices.len()
    }

    /// True when every edge is straight.
    pub fn is_polygonal(&self) -> bool {
        self.vertices.iter().all(ArcVertex::is_straight)
    }
}

/// Signed area enclosed by an arc ring.
///
/// Positive is counter-clockwise. The value is the polygon area over the
/// vertices plus, for each arc edge, the signed area of its circular
/// segment:
///
/// ```text
/// theta = 4 * atan(bulge)                 signed included angle
/// R     = chord * (1 + bulge^2) / (4 |bulge|)
/// extra = R^2 * (theta - sin theta) / 2   signed by theta
/// ```
///
/// Keeping `theta` signed is what makes a reversed ring negate exactly: a
/// clockwise circle must return `-pi r^2`, not `+pi r^2`. An unsigned
/// segment term gets the disc right and the reversed disc wrong.
pub fn arc_ring_area(ring: &ArcRing) -> f64 {
    let count = ring.vertices.len();
    if count < 2 {
        return 0.0;
    }
    let mut area = 0.0;
    for index in 0..count {
        let from = ring.vertices[index];
        let to = ring.vertices[(index + 1) % count];
        area += from.point.x * to.point.y - to.point.x * from.point.y;
    }
    area *= 0.5;
    for index in 0..count {
        let from = ring.vertices[index];
        if from.is_straight() {
            continue;
        }
        let to = ring.vertices[(index + 1) % count];
        let chord = (to.point - from.point).length();
        if chord == 0.0 {
            continue;
        }
        let bulge = from.bulge;
        let theta = 4.0 * bulge.atan();
        let radius = chord * (1.0 + bulge * bulge) / (4.0 * bulge.abs());
        area += 0.5 * radius * radius * (theta - theta.sin());
    }
    area
}

/// Arc radius implied by an edge, or `None` for a straight edge.
pub fn arc_edge_radius(from: ArcVertex, to: ArcVertex) -> Option<f64> {
    if from.is_straight() {
        return None;
    }
    let chord = (to.point - from.point).length();
    if chord == 0.0 {
        return None;
    }
    let bulge = from.bulge;
    Some(chord * (1.0 + bulge * bulge) / (4.0 * bulge.abs()))
}

/// Validate an arc ring against the same contract as [`crate::Ring`],
/// plus the arc-specific degeneracies.
///
/// Refused, each with a reason the caller can act on:
///
/// - fewer than two vertices: a circle needs two, a polygon needs three,
///   and one vertex cannot describe either
/// - a polygonal ring with fewer than three vertices: two straight edges
///   enclose no area
/// - non-finite coordinates or a non-finite bulge
/// - a zero-length edge, which leaves the arc centre undefined
/// - an arc whose implied radius is below tolerance, the zero-radius case
///   ADR 0050 flagged: such an arc is a point, not a boundary
/// - a ring enclosing no measurable area
///
/// Self-intersection is NOT checked here. Arc/arc and arc/segment crossing
/// tests are genuinely part of the overlay algorithm, and a cheap
/// approximation would either reject valid input or pass invalid input.
/// The straight-edge contract checks it because there the test is exact;
/// claiming the same guarantee for arcs without the machinery would be
/// dishonest, so the gap is named instead.
pub fn validate_arc_ring(ring: &ArcRing, tolerance: Tolerance) -> Result<(), OverlayError> {
    let count = ring.vertices.len();
    if count < 2 {
        return Err(OverlayError::RingTooShort);
    }
    if ring.is_polygonal() && count < 3 {
        return Err(OverlayError::RingTooShort);
    }
    if !ring
        .vertices
        .iter()
        .all(|vertex| vertex.point.is_finite() && vertex.bulge.is_finite())
    {
        return Err(OverlayError::NonFinitePoint);
    }
    for index in 0..count {
        let from = ring.vertices[index];
        let to = ring.vertices[(index + 1) % count];
        if (to.point - from.point).length() <= tolerance.linear() {
            return Err(OverlayError::RepeatedVertex);
        }
        // A bulge implying a sub-tolerance radius is a point masquerading
        // as an arc: the chord is shorter than the tolerance relative to
        // the curvature, so no circle can be recovered from it.
        if let Some(radius) = arc_edge_radius(from, to) {
            if radius <= tolerance.linear() {
                return Err(OverlayError::ZeroRadiusArc);
            }
        }
    }
    if arc_ring_area(ring).abs() <= tolerance.linear().powi(2) {
        return Err(OverlayError::ZeroArea);
    }
    Ok(())
}

/// Reverse a ring's orientation, preserving its geometry.
///
/// Bulges are stored per departing edge, so reversing the vertex order
/// alone would attach each bulge to the wrong edge. The bulge list must
/// shift by one and negate: the edge that left vertex `i` toward `i + 1`
/// becomes the edge leaving `i + 1` toward `i`, curving the other way.
pub fn reverse_arc_ring(ring: &ArcRing) -> ArcRing {
    let count = ring.vertices.len();
    if count == 0 {
        return ArcRing::new(Vec::new());
    }
    let mut vertices = Vec::with_capacity(count);
    for index in (0..count).rev() {
        let point = ring.vertices[index].point;
        // The edge arriving at `index` came from its predecessor.
        let predecessor = (index + count - 1) % count;
        vertices.push(ArcVertex::bulged(point, -ring.vertices[predecessor].bulge));
    }
    ArcRing::new(vertices)
}
