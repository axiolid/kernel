#![forbid(unsafe_code)]
//! Runtime provider registration and dispatch policy.
//!
//! Portable request/result schemas live in operation-contract packages. This
//! crate owns ordering, device matching, fallback, and budget admission.

#[cfg(any(
    feature = "mesh-boolean",
    feature = "mesh-section",
    feature = "pointcloud-reconstruction"
))]
mod device;

#[cfg(feature = "mesh-boolean")]
mod boolean;
#[cfg(feature = "pointcloud-reconstruction")]
mod reconstruction;
#[cfg(feature = "mesh-section")]
mod section;

#[cfg(feature = "mesh-boolean")]
pub use boolean::MeshBooleanRegistry;
#[cfg(feature = "pointcloud-reconstruction")]
pub use reconstruction::PointcloudReconstructionRegistry;
#[cfg(feature = "mesh-section")]
pub use section::MeshPlaneSectionRegistry;
