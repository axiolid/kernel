#![forbid(unsafe_code)]

//! Spatial acceleration contracts.
//!
//! BVH, octree, GPU broad phase, or a foreign spatial index can implement the
//! same callback API. Narrow-phase geometry remains outside the index.

pub mod bvh;
pub mod index;
pub mod points;

pub use bvh::{Bvh, CandidatePair, NearestCandidate, PairCandidates, SpatialQueryStats};
pub use index::{RayHit, SpatialIndex, SpatialItem};
pub use points::{PointHit, PointIndex, PointQueryError};
