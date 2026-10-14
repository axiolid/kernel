#![forbid(unsafe_code)]
//! Explicit graph-to-triangle-mesh compilation contract.
//!
//! This operation never claims to preserve an exact B-rep result domain.
//! Curve nodes, which have no area, compile to polylines through the
//! separate [`PolylineCompiler`] operation (#311).

mod contract;
mod outcome;
mod polyline;

pub use contract::MeshCompiler;
pub use outcome::{CompileOutcome, MeshClosure};
pub use polyline::{ChordDeviation, PolylineCompiler, PolylineOutcome};

pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::GRAPH_TO_MESH;
