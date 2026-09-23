#![forbid(unsafe_code)]
//! Portable mesh-boolean operation contract, evidence, and conformance suite.

pub mod conformance;
mod contract;
mod evidence;

pub use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, ExecutionOptions,
    ExecutionTarget, GeomError, GeomResult, Operation, ScratchRequirement,
};
pub use contract::{symmetric_difference_via_composition, MeshBoolean};
pub use evidence::{merge_fates, BooleanEvidence, BooleanOutcome};

pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::MESH_BOOLEAN;
