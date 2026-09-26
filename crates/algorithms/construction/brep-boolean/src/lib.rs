#![forbid(unsafe_code)]

//! General exact booleans between exact B-reps with analytic faces.
//!
//! The general-fuse pipeline of ADR 0075, built in stages:
//!
//! - [`section_edges`]: where two operands' faces cross, as exact curves
//!   trimmed exactly to both faces;
//! - [`split_face`]: one face cut along its section edges into regions, in
//!   the face's own parameters, with exact pcurves.
//!
//! Classification and assembly follow in later slices.

mod assemble;
mod classify;
mod section;
mod split;

pub use section::{section_edges, SectionEdge};
pub use split::{split_face, Piece, PieceSource, Region};

use axiolid_brep::ExactBRep;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_surface::Surface;

/// The exact boolean of two exact B-rep solids.
///
/// Stage 1 of ADR 0075: faces on planes and cylinders, meeting in lines,
/// circles and ellipses. Every other configuration is refused by name.
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
            let others = mine
                .iter()
                .map(|e| {
                    let id = if first { e.face_b } else { e.face_a };
                    other.topology().faces()[id.index()]
                        .surface
                        .and_then(|s| other.surfaces().get(s.index()))
                        .cloned()
                        .ok_or(BooleanError::DanglingReference)
                })
                .collect::<Result<Vec<_>, _>>()?;
            for region in split_face(operand, face, &mine, &others, tolerance)? {
                let point = classify::interior_point(&region, &surface)?;
                let inside = other_solid.contains(point, tolerance)?;
                let (keep, flip) = match (operator, first) {
                    (BooleanOperator::Union, _) => (!inside, false),
                    (BooleanOperator::Intersection, _) => (inside, false),
                    (BooleanOperator::Difference, true) => (!inside, false),
                    (BooleanOperator::Difference, false) => (inside, true),
                    _ => return Err(BooleanError::UnsupportedSection),
                };
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

use axiolid_measure::ExactMeasureError;
use axiolid_topology::FaceId;
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
    /// Two faces meet without crossing: coincident or tangent supports.
    NotTransverse {
        /// The face of the first operand.
        face_a: FaceId,
        /// The face of the second operand.
        face_b: FaceId,
    },
    /// A section curve runs along an existing boundary edge.
    SectionAlongEdge,
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
    /// Two pieces leave a vertex of a split face in the same direction.
    TangentSplit,
    /// The operation leaves nothing.
    EmptyResult,
    /// Cavities in a result of several solids: which holds each is not
    /// decided in this stage.
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
            Self::NotTransverse { face_a, face_b } => write!(
                f,
                "faces {} and {} meet without crossing (coincident or tangent)",
                face_a.index(),
                face_b.index()
            ),
            Self::SectionAlongEdge => f.write_str("a section runs along an existing edge"),
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
                f.write_str("cavities in a result of several solids are not assigned in this stage")
            }
            Self::Assembly => f.write_str("the kept faces did not sew into a valid exact B-rep"),
        }
    }
}

impl std::error::Error for BooleanError {}
