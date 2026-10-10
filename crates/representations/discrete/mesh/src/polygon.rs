//! Polygonal face mesh that preserves n-gons and inner voids before triangulation.

use axiolid_core::Point3;

/// One polygonal face with an outer loop and zero or more inner loops.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct PolygonFace {
    /// Outer boundary position indices.
    pub outer: Vec<u32>,
    /// Inner boundary position indices.
    pub holes: Vec<Vec<u32>>,
}

/// Indexed polygon mesh in local coordinates.
#[derive(Debug, Clone, Default, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct PolygonMesh {
    /// Shared position list.
    pub positions: Vec<Point3>,
    /// Faces in source order.
    pub faces: Vec<PolygonFace>,
}
