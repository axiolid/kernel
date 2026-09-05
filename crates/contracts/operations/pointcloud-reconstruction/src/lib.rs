#![forbid(unsafe_code)]
//! Portable pointcloud-to-surface reconstruction contract, evidence, and
//! conformance suite.
//!
//! # Why a contract before an implementation
//!
//! Reconstruction has no single right answer: Poisson, ball pivoting, and
//! alpha shapes all produce different surfaces from identical input, and
//! each is correct for different data. Fixing the *seam* — request, result,
//! evidence, typed refusal — before any provider lands is what lets those
//! methods be swapped without every consumer rewriting itself.
//!
//! # The honesty requirement
//!
//! A reconstruction is an estimate. The contract therefore refuses to let a
//! provider present a guess as a measurement:
//!
//! - An absent surface is a typed [`ReconstructionRefusal`], never an empty
//!   mesh a caller might read as "the object is not there".
//! - Surface invented across gaps in the capture is counted in
//!   [`ReconstructionEvidence::interpolated_triangles`].
//! - The sample spacing the data actually resolves is reported, so a caller
//!   can tell measured detail from interpolated detail.
//! - `require_closed` is honoured or refused, never satisfied by
//!   fabricating unmeasured surface.

pub mod conformance;
mod contract;
mod evidence;

pub use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, Determinism, ExecutionOptions,
    ExecutionTarget, GeomError, GeomResult, Operation, ScratchRequirement,
};
pub use contract::{PointcloudReconstruction, Reconstruction, ReconstructionRefusal};
pub use evidence::{
    ReconstructionEvidence, ReconstructionOutcome, ReconstructionRequest, Resolution,
};

/// Capability this contract describes.
pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::POINTCLOUD_RECONSTRUCTION;
