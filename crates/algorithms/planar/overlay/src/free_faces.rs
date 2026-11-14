//! Bounded free faces between obstacle polygons and open segments (#490).
//!
//! Given obstacle polygons (walls), barriers (open segments that also stop
//! the free space) and a closing width `gap`, [`free_faces`] returns the
//! bounded faces of the plane that no obstacle covers and no barrier cuts,
//! each with the input edges its boundary lies on, and [`FreeFaces::locate`]
//! says which face a point is in. Sections of a floor plan come as thin
//! wall pieces whose ends miss each other by float noise; a room must close
//! across such a joint, and the closure must be visible.
//!
//! # Gaps close by bridges, not by snapping
//!
//! Input is never moved. Where two boundary pieces that share no vertex come
//! closer than `gap` with nothing between them, a bridge, a straight segment
//! from a vertex to the nearest point of the other piece, joins them. A
//! bridge is a barrier the caller did not give: it is reported in
//! [`FreeFaces::bridges`] with the inputs it joins, and a ring edge on it
//! lists [`PieceSource::Bridge`]. The pipeline, on the exact arrangement of
//! `axiolid-exact-arrangement` (ADR 0094, ADR 0100):
//!
//! 1. Arrange the obstacle edges and barriers (`bridges`).
//! 2. For every vertex and every piece it is not on, whose closest point is
//!    strictly nearer than `gap` and whose open segment to it meets no piece,
//!    add a bridge. Every qualifying pair is bridged, so no choice is made
//!    and nothing can be ambiguous. Two edges that meet at a vertex are never
//!    bridged to each other there: a sharp corner stays sharp.
//! 3. Arrange inputs and bridges again, number the faces of the complement
//!    (`faces`) and read each face as obstacle or free by ring membership
//!    (`classify`, ADR 0104). A bridge inside an obstacle is discarded.
//! 4. A bounded free face that fits in a strip narrower than `gap` is a
//!    closed gap: reported in [`FreeFaces::closed_gaps`], not returned as a
//!    face. The strip rule is exact, on the convex hull of the face.
//! 5. Rings are traced, snap-rounded and certified as the visible faces are
//!    (`rings`, ADR 0100), so [`FreeFaceEvidence::snap_displacement`] bounds
//!    how far the binary64 ring moved from the exact one.
//!
//! A ring on one line is refused: a wall that vanished silently would make a
//! room leak. A piece with the same face on both sides, a dangling barrier or
//! bridge, is in no ring. `gap == 0` closes nothing and returns the exact
//! faces.
//!
//! The sources of a ring edge are lists, because inputs may coincide.
//! Everything is in `x`-then-`y` order (`axiolid_exact_arrangement`'s
//! `point_order`), and the output depends on the input order only through the
//! indices inside [`PieceSource`].
//!
//! `scripts/probe_free_faces_mutants.py` lists the faults the tests in
//! `tests/free_faces.rs` must catch.

mod bridges;
mod classify;
mod locate;
mod rings;

use axiolid_contracts::{DisplacementBound, ExecutionOptions, GeomError};
use axiolid_core::{Point2, Tolerance};
use axiolid_exact_arrangement::{arrange, faces, ExactPoints, PairBudget, PointId, Segment};

use crate::{exact_overlay, nesting, validate_ring, OverlayError, Polygon};

/// Bytes declared per segment of the arrangement: the pieces, vertices and
/// cut lists one input segment leads to (ADR 0088).
const SEGMENT_BYTES: usize = 256;

/// What blocks free space: obstacle polygons and open segments.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FreeSpaceInput {
    /// Walls and other solids. They may overlap or touch each other.
    pub obstacles: Vec<Polygon>,
    /// Open segments that stop free space like a wall of no thickness.
    pub barriers: Vec<[Point2; 2]>,
}

/// The input an edge of a free face lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum PieceSource {
    /// An edge of an obstacle. Ring 0 is the outer ring, `k >= 1` is
    /// `holes[k - 1]`; the edge leaves `points[edge]`.
    Obstacle {
        /// Index into [`FreeSpaceInput::obstacles`].
        polygon: usize,
        /// The ring of the polygon.
        ring: usize,
        /// The edge of the ring.
        edge: usize,
    },
    /// A barrier, by index into [`FreeSpaceInput::barriers`].
    Barrier {
        /// The index.
        barrier: usize,
    },
    /// A bridge the operation added, by index into [`FreeFaces::bridges`].
    Bridge {
        /// The index.
        bridge: usize,
    },
}

/// A closed boundary of a free face.
#[derive(Debug, Clone, PartialEq)]
pub struct FreeRing {
    /// The vertices, starting at the smallest by `x`, then `y`.
    pub points: Vec<Point2>,
    /// Per edge `points[i]` to `points[(i + 1) % n]`, the inputs it lies on:
    /// not empty, ascending, without repeats.
    pub sources: Vec<Vec<PieceSource>>,
}

/// A bounded free face.
#[derive(Debug, Clone, PartialEq)]
pub struct FreeFace {
    /// The outer ring, counter-clockwise.
    pub outer: FreeRing,
    /// The holes, clockwise, in the order of their points.
    pub holes: Vec<FreeRing>,
}

/// A segment added across a gap narrower than the closing width.
#[derive(Debug, Clone, PartialEq)]
pub struct Bridge {
    /// The ends, rounded once, the smaller by `x`, then `y`, first.
    pub ends: [Point2; 2],
    /// Per end, the inputs meeting or carrying it, ascending.
    pub joins: [Vec<PieceSource>; 2],
    /// The distance of the two ends as returned. The exact length is less
    /// than the closing width; this is `0.0` for a bridge shorter than the
    /// rounding.
    pub length: f64,
}

/// What the operation did besides the faces.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct FreeFaceEvidence {
    /// How far snap rounding moved a boundary point at most (ADR 0100).
    pub snap_displacement: DisplacementBound,
    /// Edges snap rounding routed through another vertex's cell.
    pub rerouted_edges: usize,
    /// Faces whose boundary rounded away.
    pub collapsed_faces: usize,
    /// Holes whose boundary rounded away.
    pub collapsed_holes: usize,
    /// Closed gaps whose boundary rounded away.
    pub collapsed_gaps: usize,
    /// Candidate pairs examined, over every step.
    pub candidate_pairs: usize,
}

/// Where a point lies relative to a [`FreeFaces`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FaceLocation {
    /// Inside the face of that index.
    Inside(usize),
    /// On a ring edge of a face. Ring 0 is the outer ring, `k >= 1` is
    /// `holes[k - 1]`; the edge's sources say what it lies on.
    OnBoundary {
        /// The face.
        face: usize,
        /// The ring of the face.
        ring: usize,
        /// The edge of the ring.
        edge: usize,
    },
    /// In or on the closed gap of that index.
    InClosedGap(usize),
    /// In or on the input obstacle of that index (the lowest).
    InObstacle(usize),
    /// In none of them: outside every bounded free face.
    Unbounded,
}

/// Why a barrier is refused.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum BarrierDefect {
    /// An end has a non-finite coordinate.
    NonFinitePoint,
    /// The ends are no further apart than the tolerance.
    TooShort,
}

/// Why no free faces were built.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum FreeFacesError {
    /// The closing width is not finite, or negative.
    InvalidGap,
    /// An obstacle ring is invalid, or lies on one line, or the rings of one
    /// obstacle are not apart and nested.
    InvalidObstacle {
        /// Index into [`FreeSpaceInput::obstacles`].
        polygon: usize,
        /// The first defect found.
        reason: OverlayError,
    },
    /// A barrier is invalid.
    InvalidBarrier {
        /// Index into [`FreeSpaceInput::barriers`].
        barrier: usize,
        /// The defect.
        reason: BarrierDefect,
    },
    /// A point to locate has a non-finite coordinate.
    NonFinitePoint,
    /// A resource limit was exceeded: `"memory"` or the name of a count.
    BudgetExceeded {
        /// The resource.
        resource: &'static str,
    },
    /// The execution options' token was cancelled.
    Cancelled,
    /// An exact step could not certify its result; the message names it.
    NotCertified(String),
}

impl core::fmt::Display for FreeFacesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidGap => f.write_str("the closing width is not finite, or is negative"),
            Self::InvalidObstacle { polygon, reason } => {
                write!(f, "obstacle {polygon} is invalid: {reason}")
            }
            Self::InvalidBarrier { barrier, reason } => {
                let why = match reason {
                    BarrierDefect::NonFinitePoint => "has a non-finite coordinate",
                    BarrierDefect::TooShort => "is not longer than the tolerance",
                };
                write!(f, "barrier {barrier} {why}")
            }
            Self::NonFinitePoint => f.write_str("a point has a non-finite coordinate"),
            Self::BudgetExceeded { resource } => write!(f, "the {resource} budget is exceeded"),
            Self::Cancelled => f.write_str("the operation was cancelled"),
            Self::NotCertified(message) => write!(f, "a result could not be certified: {message}"),
        }
    }
}

impl std::error::Error for FreeFacesError {}

impl From<GeomError> for FreeFacesError {
    fn from(error: GeomError) -> Self {
        match error {
            GeomError::BudgetExceeded { resource } => Self::BudgetExceeded { resource },
            GeomError::Cancelled => Self::Cancelled,
            other => Self::NotCertified(other.to_string()),
        }
    }
}

/// The bounded free faces of a [`FreeSpaceInput`] (module docs).
#[derive(Debug, Clone)]
pub struct FreeFaces {
    faces: Vec<FreeFace>,
    closed_gaps: Vec<FreeFace>,
    bridges: Vec<Bridge>,
    evidence: FreeFaceEvidence,
    obstacles: Vec<Polygon>,
}

impl FreeFaces {
    /// The faces, sorted by their outer rings' points.
    pub fn faces(&self) -> &[FreeFace] {
        &self.faces
    }

    /// The bounded free regions narrower than the closing width, in the same
    /// form and order as the faces.
    pub fn closed_gaps(&self) -> &[FreeFace] {
        &self.closed_gaps
    }

    /// The bridges, sorted by their ends.
    pub fn bridges(&self) -> &[Bridge] {
        &self.bridges
    }

    /// What the operation did besides the faces.
    pub fn evidence(&self) -> &FreeFaceEvidence {
        &self.evidence
    }
}

/// The bounded free faces of `input`, closing gaps narrower than `gap`.
///
/// `gap` is a distance in the units of the points: `0` closes nothing, there
/// is no upper limit, and a distance exactly `gap` is not closed.
/// `tolerance` validates the input, as everywhere in the crate; no decision
/// reads it.
///
/// # Errors
///
/// [`FreeFacesError`]; the first defect found in the order `gap`, obstacles
/// by index, barriers by index.
pub fn free_faces(
    input: &FreeSpaceInput,
    gap: f64,
    tolerance: Tolerance,
) -> Result<FreeFaces, FreeFacesError> {
    free_faces_within(input, gap, &ExecutionOptions::new(tolerance))
}

/// [`free_faces`] with execution options: the memory budget (ADR 0088) and
/// the cancellation token; the tolerance is `execution.tolerance()`.
pub fn free_faces_within(
    input: &FreeSpaceInput,
    gap: f64,
    execution: &ExecutionOptions,
) -> Result<FreeFaces, FreeFacesError> {
    validate(input, gap, execution.tolerance())?;
    let mut model = Model::assemble(input)?;
    let mut budget = PairBudget::new(usize::MAX, "free-face pairs");
    let inputs = model.segments.len();
    let first = arrange(
        &mut model.points,
        &model.segments,
        &mut budget,
        inputs * SEGMENT_BYTES,
        execution,
    )?;
    let found = bridges::find(
        &mut model.points,
        &first,
        &model.sources,
        gap,
        &mut budget,
        execution,
    )?;
    for (index, bridge) in found.iter().enumerate() {
        model.segments.push(Segment { ends: bridge.ends });
        model.sources.push(PieceSource::Bridge { bridge: index });
    }
    let arrangement = if found.is_empty() {
        first
    } else {
        arrange(
            &mut model.points,
            &model.segments,
            &mut budget,
            model.segments.len() * SEGMENT_BYTES,
            execution,
        )?
    };
    let numbered = faces(&model.points, &arrangement, &mut budget, execution)?;
    let classified = classify::classify(
        &model,
        &arrangement,
        &numbered,
        inputs,
        found.len(),
        gap,
        &mut budget,
        execution,
    )?;
    let built = rings::build(
        &model.points,
        &arrangement,
        &classified,
        &found,
        &mut budget,
        execution,
    )?;
    Ok(FreeFaces {
        faces: built.faces,
        closed_gaps: built.closed_gaps,
        bridges: built.bridges,
        evidence: FreeFaceEvidence {
            snap_displacement: built.snap_displacement,
            rerouted_edges: built.rerouted_edges,
            collapsed_faces: built.collapsed_faces,
            collapsed_holes: built.collapsed_holes,
            collapsed_gaps: built.collapsed_gaps,
            candidate_pairs: budget.used,
        },
        obstacles: input.obstacles.clone(),
    })
}

/// The first defect of `input`: gap, obstacles by index, barriers by index.
fn validate(input: &FreeSpaceInput, gap: f64, tolerance: Tolerance) -> Result<(), FreeFacesError> {
    if !gap.is_finite() || gap < 0.0 {
        return Err(FreeFacesError::InvalidGap);
    }
    for (index, polygon) in input.obstacles.iter().enumerate() {
        let invalid = |reason| FreeFacesError::InvalidObstacle {
            polygon: index,
            reason,
        };
        let rings = std::iter::once(&polygon.outer).chain(&polygon.holes);
        for ring in rings.clone() {
            validate_ring(ring, tolerance).map_err(invalid)?;
        }
        // A ring on one line is no defect of an operand, but a wall that
        // disappears would make a room leak.
        if rings
            .into_iter()
            .any(|ring| exact_overlay::on_one_line(&ring.points))
        {
            return Err(invalid(OverlayError::ZeroArea));
        }
        nesting::check(std::slice::from_ref(polygon)).map_err(invalid)?;
    }
    for (index, barrier) in input.barriers.iter().enumerate() {
        let invalid = |reason| FreeFacesError::InvalidBarrier {
            barrier: index,
            reason,
        };
        if !barrier[0].is_finite() || !barrier[1].is_finite() {
            return Err(invalid(BarrierDefect::NonFinitePoint));
        }
        if (barrier[1] - barrier[0]).length() <= tolerance.linear() {
            return Err(invalid(BarrierDefect::TooShort));
        }
    }
    Ok(())
}

/// The exact points and segments of the input, with what each segment is.
struct Model {
    points: ExactPoints,
    segments: Vec<Segment>,
    /// Per segment, the input it is: obstacle edges, barriers, then bridges.
    sources: Vec<PieceSource>,
}

impl Model {
    fn assemble(input: &FreeSpaceInput) -> Result<Self, FreeFacesError> {
        let mut points = ExactPoints::new();
        let mut segments = Vec::new();
        let mut sources = Vec::new();
        for (polygon, obstacle) in input.obstacles.iter().enumerate() {
            let rings = std::iter::once(&obstacle.outer).chain(&obstacle.holes);
            for (ring, found) in rings.enumerate() {
                let ids: Vec<PointId> = found
                    .points
                    .iter()
                    .map(|&point| points.push_flat(point))
                    .collect::<Result<_, _>>()?;
                for edge in 0..ids.len() {
                    segments.push(Segment {
                        ends: [ids[edge], ids[(edge + 1) % ids.len()]],
                    });
                    sources.push(PieceSource::Obstacle {
                        polygon,
                        ring,
                        edge,
                    });
                }
            }
        }
        for (barrier, ends) in input.barriers.iter().enumerate() {
            segments.push(Segment {
                ends: [points.push_flat(ends[0])?, points.push_flat(ends[1])?],
            });
            sources.push(PieceSource::Barrier { barrier });
        }
        Ok(Self {
            points,
            segments,
            sources,
        })
    }
}
