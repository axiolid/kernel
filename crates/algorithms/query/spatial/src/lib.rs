#![forbid(unsafe_code)]

//! Spatial acceleration contracts.
//!
//! # What this crate provides
//!
//! Two indices, both deterministic and callback-based:
//!
//! - [`Bvh`] over bounded *objects* -- triangles, solids, anything with an
//!   AABB. Adapts to how geometry is distributed, so empty space costs
//!   nothing. This is what clash, ray casting, and healing use.
//! - [`PointIndex`] over *points*, on a uniform grid. Exact KNN and radius
//!   search for scattered samples, where every query has the same radius and
//!   cell arithmetic beats tree descent.
//!
//! An octree, k-d tree, GPU broad phase, or foreign index can implement the
//! same [`SpatialIndex`] callback API. None is provided here: the BVH covers
//! object queries and the grid covers point queries, and a third structure
//! should arrive with a measured workload that needs it, not before.
//!
//! Narrow-phase geometry remains outside the index: these answer *which
//! candidates*, never *what the intersection is*.
//!
//! Beside the indices, [`barycentric`] locates a point in a triangle,
//! tetrahedron or polygon by weights of its corners, for interpolating
//! values given there.

pub mod barycentric;
pub mod bvh;
pub mod index;
pub mod points;

pub use barycentric::{
    discrete_harmonic_coordinates2, mean_value_coordinates2, mean_value_coordinates3,
    tetrahedron_barycentric, triangle_barycentric2, triangle_barycentric3, wachspress_coordinates2,
    BarycentricError,
};
pub use bvh::{Bvh, CandidatePair, NearestCandidate, PairCandidates, SpatialQueryStats};
pub use index::{RayHit, SpatialIndex, SpatialItem};
pub use points::{PointHit, PointIndex, PointQueryError};
