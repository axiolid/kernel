//! Arc-aware planar boolean (ADR 0050, step 2).
//!
//! The polygon path keeps its integer-predicate backend. This path
//! handles boundaries that carry arcs, which that backend cannot
//! represent without tessellating them away.
//!
//! # Holes arrive separately
//!
//! The backend splits results into positive and negative loops. A
//! hole fully inside the subject comes back as a NEGATIVE loop, not
//! as a second positive one. Reading only positives silently drops
//! holes -- a wall would lose its opening -- so both are consumed
//! and negatives become holes.

use cavalier_contours::core::math::Vector2;
use cavalier_contours::polyline::{BooleanOp, PlineSource, PlineSourceMut, Polyline};

use axiolid_core::{Point2, Tolerance};

use crate::arc::{arc_ring_area, validate_arc_ring, ArcRing, ArcVertex};
use crate::{OverlayError, OverlayOperation};

/// An arc-aware region: one outer boundary and its holes.
#[derive(Debug, Clone, PartialEq)]
pub struct ArcPolygon {
    /// Outer boundary, counter-clockwise.
    pub outer: ArcRing,
    /// Inner boundaries, each clockwise.
    pub holes: Vec<ArcRing>,
}

/// Evidence that the arc path did what it claims.
///
/// `arc_edges` is the load-bearing number: a tessellating backend
/// would return zero here while still producing a plausible area.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArcOverlayEvidence {
    /// Outer boundaries in the result.
    pub regions: usize,
    /// Hole boundaries across all regions.
    pub holes: usize,
    /// Curved edges preserved as arcs, never tessellated.
    pub arc_edges: usize,
    /// Straight edges in the result.
    pub line_edges: usize,
}

/// The result of an arc-aware boolean.
#[derive(Debug, Clone, PartialEq)]
pub struct ArcOverlayResult {
    pub regions: Vec<ArcPolygon>,
    pub evidence: ArcOverlayEvidence,
}

/// Convert a neutral ring into a backend polyline.
fn to_backend(ring: &ArcRing) -> Polyline {
    let mut line: Polyline = Polyline::new_closed();
    for vertex in &ring.vertices {
        line.add(vertex.point.x, vertex.point.y, vertex.bulge);
    }
    line
}

/// Convert a backend polyline back into a neutral ring.
///
/// The bulge is carried across unchanged: both sides use the same
/// convention, where a vertex owns the bulge of the edge leaving it.
fn from_backend(line: &Polyline) -> ArcRing {
    let mut vertices = Vec::with_capacity(line.vertex_count());
    for index in 0..line.vertex_count() {
        let vertex = line.at(index);
        vertices.push(ArcVertex::bulged(
            Point2::new(vertex.x, vertex.y),
            vertex.bulge,
        ));
    }
    ArcRing { vertices }
}

/// Count arc and straight edges in a ring.
fn count_edges(ring: &ArcRing) -> (usize, usize) {
    let arcs = ring
        .vertices
        .iter()
        .filter(|vertex| vertex.bulge != 0.0)
        .count();
    (arcs, ring.vertices.len() - arcs)
}

/// Orient a ring so its signed area matches the wanted sign.
///
/// Outer boundaries are counter-clockwise and holes clockwise, so a
/// consumer can rely on winding without recomputing areas.
///
/// Applied on the way in AND on the way out, which is deliberately
/// redundant: mutation testing shows either one alone still yields the
/// right answer, because normalising the output repairs a mis-wound
/// input. Removing BOTH is caught. The input normalisation is kept so
/// the backend is never asked to interpret a clockwise operand, and the
/// output one so the result contract holds regardless of what the
/// backend chose to return.
fn oriented(ring: ArcRing, want_positive: bool) -> ArcRing {
    if (arc_ring_area(&ring) > 0.0) == want_positive {
        ring
    } else {
        crate::arc::reverse_arc_ring(&ring)
    }
}

/// Boolean of two arc-capable regions.
///
/// Both operands are validated by the same contract the polygon path
/// uses, so malformed input is refused by reason before any backend
/// work happens.
///
/// Returns regions with outer boundaries counter-clockwise and holes
/// clockwise. An empty result is not an error: an intersection of
/// disjoint shapes is legitimately empty.
pub fn arc_overlay(
    subject: &ArcRing,
    clip: &ArcRing,
    operation: OverlayOperation,
    tolerance: Tolerance,
) -> Result<ArcOverlayResult, OverlayError> {
    validate_arc_ring(subject, tolerance)?;
    validate_arc_ring(clip, tolerance)?;
    let operator = match operation {
        OverlayOperation::Intersection => BooleanOp::And,
        OverlayOperation::Union => BooleanOp::Or,
        OverlayOperation::Difference => BooleanOp::Not,
        OverlayOperation::Xor => BooleanOp::Xor,
    };

    // Operands are normalised to counter-clockwise first. The backend
    // reads subtraction from winding, so a clockwise operand would
    // invert the meaning of Not without reporting an error.
    let subject_line = to_backend(&oriented(subject.clone(), true));
    let clip_line = to_backend(&oriented(clip.clone(), true));
    let result = subject_line.boolean(&clip_line, operator);

    let mut regions: Vec<ArcPolygon> = result
        .pos_plines
        .iter()
        .map(|entry| ArcPolygon {
            outer: oriented(from_backend(&entry.pline), true),
            holes: Vec::new(),
        })
        .collect();

    // Negative loops are holes. Each is attached to the region that
    // contains it; with a single outer that is unambiguous, and with
    // several the containing one is found by point-in-ring.
    for entry in &result.neg_plines {
        let hole = oriented(from_backend(&entry.pline), false);
        let Some(probe) = hole.vertices.first().map(|vertex| vertex.point) else {
            continue;
        };
        // Containment uses the backend's arc-aware winding number: a
        // straight-edge point-in-polygon test would misjudge points
        // near a bulged edge, which is exactly where holes sit.
        let owner = regions.iter_mut().find(|region| {
            to_backend(&region.outer).winding_number(Vector2::new(probe.x, probe.y)) != 0
        });
        if let Some(region) = owner {
            region.holes.push(hole);
        } else if let Some(region) = regions.first_mut() {
            region.holes.push(hole);
        }
    }

    let mut arc_edges = 0;
    let mut line_edges = 0;
    let mut holes = 0;
    for region in &regions {
        let (arcs, lines) = count_edges(&region.outer);
        arc_edges += arcs;
        line_edges += lines;
        holes += region.holes.len();
        for hole in &region.holes {
            let (arcs, lines) = count_edges(hole);
            arc_edges += arcs;
            line_edges += lines;
        }
    }

    Ok(ArcOverlayResult {
        evidence: ArcOverlayEvidence {
            regions: regions.len(),
            holes,
            arc_edges,
            line_edges,
        },
        regions,
    })
}
