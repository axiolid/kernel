#![forbid(unsafe_code)]

//! General exact booleans between exact B-reps with analytic faces.
//!
//! The general-fuse pipeline of ADR 0075, built in stages:
//!
//! - [`section_edges`]: where two operands' faces cross, as exact curves
//!   trimmed exactly to both faces;
//! - [`split_face`]: one face cut along its section edges into regions, in
//!   the face's own parameters, with exact pcurves;
//! - [`boolean`]: the regions classified against the other solid, selected
//!   by operator and sewn into the result.
//!
//! Operands may touch as well as cross: faces sharing a patch of one
//! surface, sections running along existing edges, tangent contact, and
//! solids meeting along an edge or at a point are all handled, not
//! refused.

mod assemble;
mod classify;
mod section;
mod split;
mod support;

pub use section::{section_edges, SectionEdge};
pub use split::{split_face, Piece, PieceSource, Region};

use axiolid_brep::ExactBRep;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_evaluate::surface::{invert, normal};
use axiolid_surface::Surface;
use axiolid_topology::Orientation;

/// The exact boolean of two exact B-rep solids.
///
/// ADR 0075 stages 1 and 2: faces on planes, cylinders, elliptical
/// cylinders, cones, spheres and tori, meeting in any section
/// `exact_surface_intersection` builds (lines, conics, ruled, torus and
/// traced sections).
///
/// # Errors
///
/// Anything a step refuses (see [`BooleanError`]), or a result with nothing
/// left.
pub fn boolean(
    a: &ExactBRep,
    b: &ExactBRep,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> Result<ExactBRep, BooleanError> {
    let edges = section_edges(a, b, tolerance)?;
    let solid_a = classify::Solid::new(a, tolerance)?;
    let solid_b = classify::Solid::new(b, tolerance)?;
    let mut kept = Vec::new();
    let cuts: Vec<axiolid_core::Point3> = edges.iter().flat_map(|e| [e.start, e.end]).collect();
    for (operand, other, other_solid, first) in [(a, b, &solid_b, true), (b, a, &solid_a, false)] {
        let topology = operand.topology();
        for index in 0..topology.faces().len() {
            let face = topology
                .face_id_at(index)
                .ok_or(BooleanError::DanglingReference)?;
            let record = &topology.faces()[index];
            let surface: Surface = record
                .surface
                .and_then(|id| operand.surfaces().get(id.index()))
                .ok_or(BooleanError::DanglingReference)?
                .clone();
            let mine: Vec<SectionEdge> = edges
                .iter()
                .filter(|e| {
                    if first {
                        e.face_a == face
                    } else {
                        e.face_b == face
                    }
                })
                .cloned()
                .collect();
            // The other operand's faces on this face's surface: a region
            // may lie on one of them rather than inside or outside.
            let coincident: Vec<usize> = (0..other.topology().faces().len())
                .filter(|&f| {
                    other.topology().faces()[f]
                        .surface
                        .and_then(|s| other.surfaces().get(s.index()))
                        .is_some_and(|s| support::same_support(&surface, s, tolerance))
                })
                .collect();
            for region in split_face(operand, face, &mine, first, &cuts, tolerance)? {
                let sign = match record.orientation {
                    Orientation::Forward if !region.against => 1.0,
                    Orientation::Reversed if region.against => 1.0,
                    _ => -1.0,
                };
                let (keep, flip) = decide(
                    &region,
                    &surface,
                    sign,
                    other_solid,
                    &coincident,
                    operator,
                    first,
                    tolerance,
                )?;
                if keep {
                    kept.push(assemble::Kept {
                        surface: surface.clone(),
                        orientation: record.orientation,
                        region,
                        flip,
                    });
                }
            }
        }
    }
    assemble::assemble(&kept, tolerance)
}

/// Whether a region is kept, and whether it bounds the result from its
/// other side: classified at the first interior point that decides.
#[allow(clippy::too_many_arguments)]
fn decide(
    region: &Region,
    surface: &Surface,
    sign: f64,
    other_solid: &classify::Solid<'_>,
    coincident: &[usize],
    operator: BooleanOperator,
    first: bool,
    tolerance: Tolerance,
) -> Result<(bool, bool), BooleanError> {
    let mut last = BooleanError::Undecided;
    for point in classify::interior_points(region, surface)? {
        match classify_point(
            point,
            surface,
            sign,
            other_solid,
            coincident,
            operator,
            first,
            tolerance,
        ) {
            Err(BooleanError::Undecided) => last = BooleanError::Undecided,
            other => return other,
        }
    }
    Err(last)
}

#[allow(clippy::too_many_arguments)]
fn classify_point(
    point: axiolid_core::Point3,
    surface: &Surface,
    sign: f64,
    other_solid: &classify::Solid<'_>,
    coincident: &[usize],
    operator: BooleanOperator,
    first: bool,
    tolerance: Tolerance,
) -> Result<(bool, bool), BooleanError> {
    if let Some(theirs) = other_solid.on_face(point, coincident, tolerance)? {
        // On the other solid's boundary: kept once, from the first
        // operand, where the operator leaves a boundary.
        let (u, v) = invert(surface, point, tolerance).map_err(|_| BooleanError::Evaluation)?;
        let ours = normal(surface, u, v).map_err(|_| BooleanError::Evaluation)? * sign;
        let same = ours.dot(theirs) > 0.0;
        let keep = first
            && match operator {
                BooleanOperator::Union | BooleanOperator::Intersection => same,
                BooleanOperator::Difference => !same,
                _ => return Err(BooleanError::UnsupportedSection),
            };
        return Ok((keep, false));
    }
    let inside = other_solid.contains(point, tolerance)?;
    Ok(match (operator, first) {
        (BooleanOperator::Union, _) => (!inside, false),
        (BooleanOperator::Intersection, _) => (inside, false),
        (BooleanOperator::Difference, true) => (!inside, false),
        (BooleanOperator::Difference, false) => (inside, true),
        _ => return Err(BooleanError::UnsupportedSection),
    })
}

use axiolid_measure::ExactMeasureError;
use core::fmt;

/// Why a boolean step could not be carried out exactly.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BooleanError {
    /// Two faces' supports intersect in a curve this stage does not take
    /// (anything but a line, circle or ellipse).
    UnsupportedSection,
    /// A face is trimmed by an edge or pcurve family this stage cannot
    /// intersect or classify against.
    UnsupportedTrim,
    /// A point was too close to a face boundary to classify.
    Undecided,
    /// A curve or surface could not be evaluated or inverted.
    Evaluation,
    /// A handle referenced missing geometry.
    DanglingReference,
    /// A face domain could not be built.
    Measure(ExactMeasureError),
    /// A face whose surface, or a section on it, has no exact pcurve in
    /// this stage.
    UnsupportedSplit,
    /// The pieces of a split face do not close into loops.
    UnclosedSplit,
    /// Two pieces leave a vertex of a split face in the same direction and
    /// bend alike there, so no order between them can be read off.
    TangentSplit,
    /// The operation leaves nothing: an intersection of solids that only
    /// touch, or a difference that removes everything. An exact B-rep
    /// cannot be empty, so this is the answer, not a refusal.
    EmptyResult,
    /// A cavity of the result lies inside none of its solids.
    AmbiguousCavity,
    /// The kept faces did not sew into a valid exact B-rep.
    Assembly,
}

impl fmt::Display for BooleanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSection => {
                f.write_str("a face pair meets in a curve this stage does not build")
            }
            Self::UnsupportedTrim => {
                f.write_str("a face is trimmed by a curve family this stage cannot intersect")
            }
            Self::Undecided => f.write_str("a point lies too close to a face boundary to classify"),
            Self::Evaluation => f.write_str("a curve or surface could not be evaluated"),
            Self::DanglingReference => f.write_str("a handle references missing geometry"),
            Self::Measure(error) => write!(f, "a face domain could not be built: {error}"),
            Self::UnsupportedSplit => {
                f.write_str("a face or section has no exact pcurve in this stage")
            }
            Self::UnclosedSplit => f.write_str("the pieces of a split face do not close"),
            Self::TangentSplit => {
                f.write_str("two pieces leave a vertex of a split face in one direction")
            }
            Self::EmptyResult => f.write_str("the operation leaves nothing"),
            Self::AmbiguousCavity => {
                f.write_str("a cavity of the result lies inside none of its solids")
            }
            Self::Assembly => f.write_str("the kept faces did not sew into a valid exact B-rep"),
        }
    }
}

impl std::error::Error for BooleanError {}
