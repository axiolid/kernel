#![forbid(unsafe_code)]

//! Scalar reference `MeshCompiler`.

mod brep;
mod channels;
mod directrix;

use axiolid_contracts::BackendId;

/// This provider's identity.
pub const BACKEND_ID: BackendId = BackendId::new("scalar-compile");

mod compiler;
pub use compiler::ReferenceMeshCompiler;

mod exact;
pub use exact::{ReferenceExactCompiler, SOLID_FAMILY_NAMES};
