#![forbid(unsafe_code)]
//! Portable hidden-line drawing operation contract and evidence.
//!
//! The request names triangle meshes, a view (orthographic, or
//! perspective, ADR 0107) and the edge kinds
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
//! Visible face regions are a separate operation, [`MeshVisibleFaces`]
//! (ADR 0094), over the same view and meshes.
//!
//! No tolerance enters a visibility decision; a provider refuses what it
//! cannot decide.
//!
//! A refusal names the first refused mesh only ([`HiddenLineError`]). A
//! registry's audit lists every refused mesh ([`RefusedMesh`]) in one call, so
//! a caller can leave them out and ask again; no substitute drawing is made.

pub mod conformance;
mod contract;
mod error;
mod face_fixtures;
mod faces;
mod fixtures;

#[allow(deprecated)]
pub use contract::OrthographicView;
pub use contract::{
    EdgeKind, EdgeSelection, HiddenLineEvidence, HiddenLineLimits, HiddenLineMesh,
    HiddenLineOutcome, HiddenLineRequest, HiddenLineSegment, HiddenLineSource, MeshHiddenLines,
    MeshId, MeshRole, SourceEdge, View, ViewProjection, Visibility,
};
pub use error::{HiddenLineError, RefusedMesh};
pub use faces::{
    FaceGroup, FaceGroups, FaceSource, MeshVisibleFaces, VisibleFaceEvidence, VisibleFaceLimits,
    VisibleFaceOutcome, VisibleFaceRequest, VisibleRegion,
};

pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::MESH_HIDDEN_LINES;

/// Capability proven by implementing [`MeshVisibleFaces`] (ADR 0094).
pub const FACES_CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::MESH_VISIBLE_FACES;
