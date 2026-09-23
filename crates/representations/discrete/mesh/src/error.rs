//! Structured mesh validation failures.

use core::fmt;

/// Cheap structural validation failure.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshValidationError {
    /// Triangle index buffer is not divisible by three.
    IncompleteTriangle { index_count: usize },
    /// Position index does not exist.
    PositionIndexOutOfRange { index: u32, position_count: usize },
    /// Non-indexed normals do not align with positions.
    NormalCount { expected: usize, actual: usize },
    /// Independently indexed normals do not align with corners.
    NormalIndexCount { expected: usize, actual: usize },
    /// Normal index does not exist.
    NormalIndexOutOfRange { index: u32, normal_count: usize },
    /// A channel does not carry one tuple per vertex.
    AttributeCount {
        /// Channel name.
        name: String,
        /// Vertices in the mesh.
        expected: usize,
        /// Vertices the channel covers.
        actual: usize,
    },
    /// A channel declares a zero tuple width, which covers no vertices.
    AttributeZeroWidth {
        /// Channel name.
        name: String,
    },
    /// Two channels share a name, so a lookup would be ambiguous.
    AttributeDuplicateName {
        /// The repeated name.
        name: String,
    },
    /// A corner-indexed channel's `values` is not whole tuples.
    AttributeRaggedValues {
        /// Channel name.
        name: String,
        /// Scalars found.
        values: usize,
        /// Declared tuple width.
        width: usize,
    },
    /// A corner-indexed channel does not have one entry per triangle corner.
    AttributeCornerCount {
        /// Channel name.
        name: String,
        /// Triangle corners in the mesh.
        expected: usize,
        /// Corner entries in the channel.
        actual: usize,
    },
    /// A corner index points past the channel's values.
    AttributeCornerIndexOutOfRange {
        /// Channel name.
        name: String,
        /// The offending index.
        index: u32,
        /// Tuples available.
        value_count: usize,
    },
    /// A triangle has some corners mapped and some unmapped.
    AttributePartiallyMapped {
        /// Channel name.
        name: String,
        /// The triangle.
        triangle: usize,
    },
}

impl fmt::Display for MeshValidationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IncompleteTriangle { index_count } => {
                write!(f, "index count {index_count} is not divisible by three")
            }
            Self::PositionIndexOutOfRange {
                index,
                position_count,
            } => write!(
                f,
                "position index {index} exceeds {position_count} positions"
            ),
            Self::NormalCount { expected, actual } => {
                write!(f, "expected {expected} normals, found {actual}")
            }
            Self::NormalIndexCount { expected, actual } => {
                write!(f, "expected {expected} normal indices, found {actual}")
            }
            Self::NormalIndexOutOfRange {
                index,
                normal_count,
            } => write!(f, "normal index {index} exceeds {normal_count} normals"),
            Self::AttributeCount {
                name,
                expected,
                actual,
            } => write!(
                f,
                "attribute channel {name} covers {actual} vertices, mesh has {expected}"
            ),
            Self::AttributeZeroWidth { name } => {
                write!(f, "attribute channel {name} declares a zero tuple width")
            }
            Self::AttributeDuplicateName { name } => {
                write!(f, "attribute channel name {name} is used more than once")
            }
            Self::AttributeRaggedValues {
                name,
                values,
                width,
            } => write!(
                f,
                "attribute channel {name} has {values} scalars, not a multiple of width {width}"
            ),
            Self::AttributeCornerCount {
                name,
                expected,
                actual,
            } => write!(
                f,
                "attribute channel {name} has {actual} corner indices, mesh has {expected} corners"
            ),
            Self::AttributeCornerIndexOutOfRange {
                name,
                index,
                value_count,
            } => write!(
                f,
                "attribute channel {name} corner index {index} exceeds {value_count} values"
            ),
            Self::AttributePartiallyMapped { name, triangle } => write!(
                f,
                "attribute channel {name} maps only some corners of triangle {triangle}"
            ),
        }
    }
}

impl std::error::Error for MeshValidationError {}
