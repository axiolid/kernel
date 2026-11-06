#![forbid(unsafe_code)]
//! Common provider-neutral execution and diagnostic contracts.
//!
//! Operation request/response schemas live in sibling operation packages;
//! provider selection and fallback policy live under `execution/`.

pub mod backend;
pub mod cancel;
pub mod capability;
pub mod capability_id;
pub mod displacement;
pub mod error;
pub mod execution;
pub mod integration;
pub mod plan;
pub mod profiles;

pub use axiolid_guarantees::{Certified, EscalationLadder, Precision, Sign};
pub use backend::Backend;
pub use cancel::{CancellationGranularity, CancellationToken};
pub use capability::{BackendDescriptor, BackendId, ExecutionTarget, Operation};
pub use capability_id::{capability_ids, CapabilityId};
pub use displacement::DisplacementBound;
pub use error::{GeomError, GeomResult};
pub use execution::{
    DataResidency, Determinism, DevicePreference, ExecutionOptions, OutputBound, Parallelism,
    Residency, ScratchRequirement,
};
pub use integration::{
    ApiVersion, BoundaryContract, CapabilityDescriptor, CapabilityRequirement, Exactness,
    IntegrationDescriptor, IntegrationProfile, Ownership, Representation, RequirementRefusal,
    ThreadSafety, INTEGRATION_API_VERSION, MINIMUM_RUST_VERSION,
};
pub use plan::{Plan, PlanStep};
pub use profiles::{ProfileContract, V04_PROFILE_CONTRACTS};
