#![forbid(unsafe_code)]

//! Scalar reference `MeshCompiler`.

mod bounded;
mod brep;
mod certify;
mod channels;
pub mod deviation;
mod directrix;
mod half_space_boundary;
mod pinch;
mod planar;
pub mod station;
mod weld;

use axiolid_contracts::BackendId;

/// This provider's identity.
pub const BACKEND_ID: BackendId = BackendId::new("scalar-compile");

mod compiler;
pub use compiler::ReferenceMeshCompiler;
pub use deviation::{DeviationBound, DeviationContribution, DeviationPath, DeviationReport};

mod exact;
/// What a boolean read within tolerance, as
/// [`ReferenceExactCompiler::compile_exact_with_report`] reports it (#236).
pub use axiolid_brep_boolean::{
    BooleanReport, ToleranceDecision, ToleranceDecisionKind, ROUNDING_FACTOR,
};
pub use directrix::ExactDirectrix;
pub use exact::{exact_directrix, ReferenceExactCompiler, SOLID_FAMILY_NAMES};
