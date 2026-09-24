//! Arc-aware planar boolean (ADR 0050; exact since ADR 0070).
//!
//! The polygon path keeps its integer-predicate backend. This path handles
//! boundaries that carry arcs, which that backend cannot represent without
//! tessellating them away.
//!
//! # Exact, not tolerant
//!
//! Every topological decision is an exact sign over the given `f64` input
//! (see `exact_arc`): where boundaries cross, in which order, what lies
//! inside what, how the result links into rings, which ring is a hole of
//! which. None of them uses the tolerance, so a scene gives the same answer
//! in millimetres and in metres. Crossing points of two curves are in
//! general irrational; they are rounded to `f64` once, in the output.
//!
//! The tolerance still validates operands ([`validate_arc_ring`]), the
//! same contract the polygon path applies.

use axiolid_core::Tolerance;

use crate::arc::{arc_ring_area, validate_arc_ring, ArcRing};
use crate::exact_arc;
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
/// consumer can rely on winding without recomputing areas. Operands are
/// normalised on the way in (the exact core assumes the region lies left of
/// its boundary) and results on the way out.
fn oriented(ring: ArcRing, want_positive: bool) -> ArcRing {
    if (arc_ring_area(&ring) > 0.0) == want_positive {
        ring
    } else {
        crate::arc::reverse_arc_ring(&ring)
    }
}

/// Collapse edges no longer than `tolerance` after output rounding.
///
/// Exact topology keeps pieces of any length. A crossing that lies within
/// rounding of a vertex (a circle drawn through a corner, with coordinates
/// like 0.3 that binary cannot hold) leaves a piece far shorter than any
/// meaningful length, and once both ends are rounded to `f64` it can become
/// a repeated vertex. Such an edge is dropped: its end vertex goes, and the
/// following edge's bulge moves to its start, which lies within
/// `tolerance` of the dropped vertex. A ring left without a valid shape
/// (fewer vertices than a boundary needs, or no area) was a sliver below
/// the tolerance and is removed.
///
/// This is the only place the tolerance acts on a result, and it acts on
/// presentation only: which pieces exist and how they link were decided
/// exactly beforehand.
fn presented(mut ring: ArcRing, tolerance: Tolerance) -> Option<ArcRing> {
    loop {
        let count = ring.vertices.len();
        if count < 2 {
            return None;
        }
        let short = (0..count).find(|&i| {
            let (a, b) = (ring.vertices[i].point, ring.vertices[(i + 1) % count].point);
            (b - a).length() <= tolerance.linear()
        });
        let Some(i) = short else {
            break;
        };
        let j = (i + 1) % count;
        ring.vertices[i].bulge = ring.vertices[j].bulge;
        ring.vertices.remove(j);
    }
    match validate_arc_ring(&ring, tolerance) {
        Err(OverlayError::RingTooShort | OverlayError::ZeroArea) => None,
        _ => Some(ring),
    }
}

/// Boolean of two arc-capable regions.
///
/// Both operands are validated by the same contract the polygon path
/// uses, so malformed input is refused by reason before any work happens.
///
/// Returns regions with outer boundaries counter-clockwise and holes
/// clockwise. An empty result is not an error: an intersection of
/// disjoint shapes is legitimately empty.
///
/// # Errors
///
/// Any [`validate_arc_ring`] refusal, and
/// [`OverlayError::SelfIntersection`] when an operand's boundary crosses
/// itself, which the exact linking detects.
pub fn arc_overlay(
    subject: &ArcRing,
    clip: &ArcRing,
    operation: OverlayOperation,
    tolerance: Tolerance,
) -> Result<ArcOverlayResult, OverlayError> {
    validate_arc_ring(subject, tolerance)?;
    validate_arc_ring(clip, tolerance)?;
    let subject = oriented(subject.clone(), true);
    let clip = oriented(clip.clone(), true);

    let regions: Vec<ArcPolygon> = exact_arc::boolean(&subject, &clip, operation)?
        .into_iter()
        .filter_map(|(outer, holes)| {
            Some(ArcPolygon {
                outer: oriented(presented(outer, tolerance)?, true),
                holes: holes
                    .into_iter()
                    .filter_map(|h| presented(h, tolerance))
                    .map(|h| oriented(h, false))
                    .collect(),
            })
        })
        .collect();

    let mut arc_edges = 0;
    let mut line_edges = 0;
    let mut holes = 0;
    for region in &regions {
        holes += region.holes.len();
        for ring in std::iter::once(&region.outer).chain(&region.holes) {
            let (arcs, lines) = count_edges(ring);
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
