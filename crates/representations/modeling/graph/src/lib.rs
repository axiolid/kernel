#![forbid(unsafe_code)]

//! Format-neutral geometry intermediate representation.
//!
//! Source adapters lower into an immutable [`GeometryGraph`]. Nodes preserve
//! exact curve, surface, topology, instancing, and construction intent; kernels
//! choose how to evaluate or tessellate them. Typed append-only handles replace
//! recursive `Box` trees and make mapped-item/CSG cycles impossible.
//!
//! Source identifiers stay outside the graph: an adapter keeps its own map
//! from source entities to [`NodeId`]s. The graph never lowers anything to a
//! mesh; compilation lives in the execution tier.

pub mod curve_relation;
pub mod graph;
pub mod id;
pub mod node;
pub mod solid_operation;
pub mod station;
pub mod surface_relation;
mod validation;
pub mod value;

pub use axiolid_core::BooleanOperator;
pub use curve_relation::{
    CurveRelation, CurveSegment, MasterRepresentation, SurfaceSides, Transition, TrimSelector,
    TrimmingPreference,
};
pub use graph::{GeometryGraph, GeometryGraphBuilder, GraphError};
pub use id::NodeId;
pub use node::{GeometryNode, Instance, OpenProfile, PointOnCurve, PointOnSurface};
pub use solid_operation::{Section, SolidOperation};
pub use station::{
    CurveStation, OrientedCurveStation, SectionAtStation, Station, StationFrame, StationOffsets,
    StationOrientation, StationedOpenSection, StationedSection, ORIENTATION_TOLERANCE,
};
pub use surface_relation::SurfaceRelation;
pub use value::BuiltInNode;
