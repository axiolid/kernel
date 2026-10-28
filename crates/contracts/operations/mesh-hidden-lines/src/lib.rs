#![forbid(unsafe_code)]
//! Portable orthographic hidden-line drawing operation contract and evidence.
//!
//! The request names triangle meshes, an orthographic view and the edge kinds
//! to draw. The outcome lists, per source edge, maximal pieces that are
//! visible or hidden in the view. Visibility is a geometric decision made
//! by the provider; the contract fixes what the decision means (ADR 0090):
//!
//! - a point is hidden exactly when a kept triangle point lies strictly in
//!   front of it at the same image point;
//! - triangles are closed sets, so touching an image boundary covers;
//! - equal depth never occludes;
//! - a face seen edge-on occludes through the depth range it spans.
//!
//! No tolerance enters a visibility decision; a provider refuses what it
//! cannot decide.

pub mod conformance;
mod contract;
mod error;
mod fixtures;

pub use contract::{
    EdgeKind, EdgeSelection, HiddenLineEvidence, HiddenLineLimits, HiddenLineMesh,
    HiddenLineOutcome, HiddenLineRequest, HiddenLineSegment, HiddenLineSource, MeshHiddenLines,
    MeshId, MeshRole, OrthographicView, SourceEdge, Visibility,
};
pub use error::HiddenLineError;

pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::MESH_HIDDEN_LINES;
