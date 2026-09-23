#![forbid(unsafe_code)]
//! Explicit graph-to-triangle-mesh compilation contract.
//!
//! This operation never claims to preserve an exact B-rep result domain.

mod contract;
mod outcome;

pub use contract::MeshCompiler;
pub use outcome::CompileOutcome;

pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::GRAPH_TO_MESH;
