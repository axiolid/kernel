#![forbid(unsafe_code)]
//! Runtime provider registration and dispatch policy.
//!
//! Portable request/result schemas live in operation-contract packages. This
//! crate owns ordering, device matching, fallback, and budget admission.

#[cfg(any(
    feature = "mesh-boolean",
    feature = "mesh-section",
    feature = "mesh-hidden-lines",
    feature = "pointcloud-reconstruction"
))]
mod device;
#[cfg(any(feature = "mesh-section", feature = "mesh-hidden-lines"))]
mod frame;

#[cfg(feature = "mesh-boolean")]
mod boolean;
#[cfg(feature = "mesh-hidden-lines")]
mod hidden_lines;
#[cfg(feature = "pointcloud-reconstruction")]
mod reconstruction;
#[cfg(feature = "mesh-section")]
mod section;
#[cfg(feature = "mesh-section")]
mod section_request;

#[cfg(feature = "mesh-boolean")]
pub use boolean::MeshBooleanRegistry;
#[cfg(feature = "mesh-hidden-lines")]
pub use hidden_lines::MeshHiddenLinesRegistry;
#[cfg(feature = "pointcloud-reconstruction")]
pub use reconstruction::PointcloudReconstructionRegistry;
#[cfg(feature = "mesh-section")]
pub use section::MeshPlaneSectionRegistry;
