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

/// The backend's native position epsilon (`PlineBooleanOptions::new`).
///
/// Its other thresholds are fixed relative to this one: 1e-5 and 1e-3 in
/// debug self-validation, 1e-8 in its default fuzzy comparisons.
const BACKEND_POS_EPS: f64 = 1e-5;

/// Power-of-two factor that maps the caller's linear tolerance onto the
/// backend's native epsilon.
///
/// The backend hard-codes its thresholds in drawing units, so the same
/// scene gave different answers in millimetres and metres. Scaling the
/// geometry keeps every one of those thresholds at a fixed ratio to the
/// caller's tolerance, which setting a single option cannot: overriding
/// only `pos_equal_eps` trips the backend's own debug consistency checks.
/// A power of two makes the round trip bit-exact; the effective tolerance
/// is within a factor of sqrt(2) of the requested one.
///
/// The scale is capped so the scaled drawing stays within
/// [`MAX_BACKEND_EXTENT`]. Beyond that, f64 no longer resolves the
/// backend's finest (1e-8) comparisons and it fails its own consistency
/// checks. A tolerance finer than the input's extent allows is therefore
/// honoured only down to `extent * 1e-5 / MAX_BACKEND_EXTENT`, about 1e-11
/// of the extent. That is still far below anything the audit gate after
/// the boolean accepts.
///
/// `Tolerance::ZERO` asks for no scale-derived tolerance; the backend
/// then runs at its native thresholds (factor 1), as before.
fn backend_scale(tolerance: Tolerance, extent: f64) -> f64 {
    let linear = tolerance.linear();
    if linear <= 0.0 || !extent.is_finite() {
        return 1.0;
    }
    let mut exponent = (BACKEND_POS_EPS / linear).log2().round();
    if extent > 0.0 {
        exponent = exponent.min((MAX_BACKEND_EXTENT / extent).log2().floor());
    }
    2f64.powi(exponent.clamp(-1000.0, 1000.0) as i32)
}

/// Largest coordinate magnitude handed to the backend after scaling.
///
/// Its finest comparison is 1e-8 absolute; f64 resolves that only while
/// magnitudes stay below about 1e-8 / 2^-52 = 4.5e7. 1e6 leaves a margin.
const MAX_BACKEND_EXTENT: f64 = 1e6;

/// Largest absolute coordinate across both operands.
fn extent(rings: [&ArcRing; 2]) -> f64 {
    rings
        .iter()
        .flat_map(|ring| ring.vertices.iter())
        .map(|vertex| vertex.point.x.abs().max(vertex.point.y.abs()))
        .fold(0.0, f64::max)
}

/// Convert a neutral ring into a backend polyline, scaled by `scale`.
///
/// Bulges are dimensionless (`tan(theta/4)`), so only points scale.
fn to_backend(ring: &ArcRing, scale: f64) -> Result<Polyline, OverlayError> {
    let mut line: Polyline = Polyline::new_closed();
    for vertex in &ring.vertices {
        let (x, y) = (vertex.point.x * scale, vertex.point.y * scale);
        if !x.is_finite() || !y.is_finite() {
            return Err(OverlayError::NonFinitePoint);
        }
        line.add(x, y, vertex.bulge);
    }
    Ok(line)
}

/// Convert a backend polyline back into a neutral ring, undoing `scale`.
///
/// The bulge is carried across unchanged: both sides use the same
/// convention, where a vertex owns the bulge of the edge leaving it.
fn from_backend(line: &Polyline, scale: f64) -> ArcRing {
    let mut vertices = Vec::with_capacity(line.vertex_count());
    for index in 0..line.vertex_count() {
        let vertex = line.at(index);
        vertices.push(ArcVertex::bulged(
            Point2::new(vertex.x / scale, vertex.y / scale),
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
    let scale = backend_scale(tolerance, extent([subject, clip]));
    let subject_line = to_backend(&oriented(subject.clone(), true), scale)?;
    let clip_line = to_backend(&oriented(clip.clone(), true), scale)?;
    let result = subject_line.boolean(&clip_line, operator);

    let mut regions: Vec<ArcPolygon> = Vec::with_capacity(result.pos_plines.len());
    // Outers in backend coordinates, kept for hole containment below.
    let mut outer_lines = Vec::with_capacity(result.pos_plines.len());
    for entry in &result.pos_plines {
        regions.push(ArcPolygon {
            outer: oriented(from_backend(&entry.pline, scale), true),
            holes: Vec::new(),
        });
        outer_lines.push(&entry.pline);
    }

    // Negative loops are holes. Each is attached to the region that
    // contains it; with a single outer that is unambiguous, and with
    // several the containing one is found by point-in-ring.
    for entry in &result.neg_plines {
        if entry.pline.vertex_count() == 0 {
            continue;
        }
        let probe = entry.pline.at(0);
        let hole = oriented(from_backend(&entry.pline, scale), false);
        // Containment uses the backend's arc-aware winding number, in the
        // backend's own coordinates: a straight-edge point-in-polygon test
        // would misjudge points near a bulged edge, where holes sit.
        let owner = outer_lines
            .iter()
            .position(|line| line.winding_number(Vector2::new(probe.x, probe.y)) != 0);
        if let Some(index) = owner {
            regions[index].holes.push(hole);
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
