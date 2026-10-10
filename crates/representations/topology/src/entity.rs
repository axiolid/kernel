//! B-rep entities. `G` is a caller-chosen curve/surface geometry handle.

use axiolid_core::Point3;

use crate::{EdgeId, FaceId, LoopId, ShellId, VertexId};

/// Topological orientation relative to the underlying geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum Orientation {
    /// Same parameter direction/normal.
    Forward,
    /// Reversed parameter direction/normal.
    Reversed,
}

/// Vertex with an explicit model-space position.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Vertex {
    /// Position.
    pub position: Point3,
}

/// Edge bounded by two vertices and optionally supported by exact curve data.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Edge<G> {
    /// Start vertex.
    pub start: VertexId,
    /// End vertex.
    pub end: VertexId,
    /// Exact support curve handle; absent for a straight topological edge whose
    /// endpoints are sufficient.
    pub curve: Option<G>,
}

/// One oriented use of an edge in a loop.
///
/// The `pcurve` is the edge's image in the parameter space of the face this
/// use belongs to. A 3D edge curve says where a boundary sits in model
/// space; it does not say where that boundary lies in a surface's `(u, v)`
/// domain, and inverting a surface to recover it is not generally solvable
/// in closed form. Trimming a curved face therefore needs it stated, which
/// is what exchange formats carry alongside the 3D edge curve.
///
/// It belongs to the USE rather than the edge: one edge bounds two faces
/// with different support surfaces, so it has a different parameter image
/// in each.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct EdgeUse<G> {
    /// Referenced edge.
    pub edge: EdgeId,
    /// Traversal direction.
    pub orientation: Orientation,
    /// Optional 2D curve handle in the owning face's surface parameters.
    pub pcurve: Option<G>,
}

/// Closed boundary wire.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Loop<G> {
    /// Consecutive oriented edges.
    pub edges: Vec<EdgeUse<G>>,
}

/// One oriented loop use on a face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct FaceBound {
    /// Referenced loop.
    pub loop_id: LoopId,
    /// Whether this bound has the same orientation as the face.
    pub orientation: Orientation,
    /// Whether this is the outer bound.
    pub outer: bool,
}

/// Face supported by an exact surface and bounded by loops.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Face<G> {
    /// Exact support surface handle. Planar polygonal faces may omit it.
    pub surface: Option<G>,
    /// Outer and inner bounds.
    pub bounds: Vec<FaceBound>,
    /// Orientation relative to the support surface normal.
    pub orientation: Orientation,
}

/// Connected collection of oriented faces.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Shell {
    /// Face handles.
    pub faces: Vec<(FaceId, Orientation)>,
    /// Whether the source asserts closure.
    pub closed: bool,
}

/// Solid with one outer shell and optional void shells.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Solid {
    /// Outer shell.
    pub outer: ShellId,
    /// Interior void shells.
    pub voids: Vec<ShellId>,
}
