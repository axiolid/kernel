#![forbid(unsafe_code)]
//! Portable curve-evaluation contract and conformance suite.
//!
//! Curve evaluation was the only kernel capability without a contract,
//! so a consumer that wanted to ask for a point at a distance had to
//! depend on a specific engine. This crate lets the capability be NAMED
//! without binding an implementation to it.
//!
//! See `docs/adr/0063-curve-evaluation-contract.md`.

pub mod conformance;
mod contract;
mod convention;
mod measure;

pub use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, Determinism, ExecutionTarget, GeomError, GeomResult,
    Operation,
};
pub use contract::CurveEvaluator;
pub use convention::DistanceConvention;
pub use measure::CurveMeasure;

/// Capability this contract names.
pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::CURVE_EVALUATE;
