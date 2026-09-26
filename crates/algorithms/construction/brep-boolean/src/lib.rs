#![forbid(unsafe_code)]

//! General exact booleans between exact B-reps with analytic faces.
//!
//! The general-fuse pipeline of ADR 0075, built in stages. This first slice
//! provides the section edges: where two operands' faces cross, as exact
//! curves trimmed exactly to both faces. Face splitting, classification and
//! assembly follow in later slices.

mod section;

pub use section::{section_edges, SectionEdge};

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
        }
    }
}

impl std::error::Error for BooleanError {}
