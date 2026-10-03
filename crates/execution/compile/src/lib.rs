#![forbid(unsafe_code)]

//! Scalar reference `MeshCompiler`.

mod bounded;
mod brep;
mod certify;
mod channels;
pub mod deviation;
mod directrix;
mod pinch;
mod planar;

use axiolid_contracts::BackendId;

/// This provider's identity.
pub const BACKEND_ID: BackendId = BackendId::new("scalar-compile");

mod compiler;
pub use compiler::ReferenceMeshCompiler;
pub use deviation::{DeviationBound, DeviationContribution, DeviationPath, DeviationReport};

mod exact;
pub use directrix::ExactDirectrix;
pub use exact::{exact_directrix, ReferenceExactCompiler, SOLID_FAMILY_NAMES};
