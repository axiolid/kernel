#![forbid(unsafe_code)]

//! Analytic and spline curve/surface evaluation (ADR 0012, ADR 0036).
//!
//! This is the scalar evaluation oracle for parametric geometry: native-domain
//! evaluation, derivatives, jets, adaptive flattening, and elementary surface
//! inversion. It is deliberately separate from the `axiolid-reference`
//! umbrella so a parametric consumer (NURBS, CAD) acquires evaluation without
//! the umbrella's mesh, spatial, and measure dependencies.
//!
//! No intrinsics, no threading, no feature gates: it must stay obviously
//! correct in preference to being fast.
//!
//! `axiolid-reference` re-exports [`curve`] and [`surface`] unchanged, so
//! `axiolid_reference::curve::*` paths are part of this package's public
//! surface too: renaming or reshaping an item here breaks those callers.

pub mod arc_length;
pub mod arc_parameter;
pub mod banked;
pub mod bound;
pub mod chain;
pub mod curve;
pub mod elevation;
pub mod frenet;
pub mod intrinsic_relation;
mod nurbs;
pub mod polyline_length;
pub mod provider;
pub mod station;
pub mod surface;

pub use arc_length::{elevated_point, elevated_tangent, intrinsic_point, intrinsic_tangent};
pub use arc_parameter::{
    arc_length2, arc_length3, parameter_at_arc_length2, parameter_at_arc_length3,
    ARC_LENGTH_TOLERANCE,
};
pub use banked::{banked_derivative, banked_point, banked_section, banked_tangent, BankedSection};
pub use chain::{chain_point, chain_tangent};
pub use curve::{derivative2, derivative3, evaluate2, evaluate3, flatten2, ScalarCurve};
pub use elevation::{elevation_chord_bound, elevation_grade, elevation_height};
pub use frenet::{frenet_frame, frenet_point, frenet_tangent};
pub use intrinsic_relation::{join_intrinsic3, offset_intrinsic3, trim_intrinsic3};
pub use provider::ReferenceCurveEvaluator;
pub use station::{
    station_length2, station_length3, station_section2, station_section3, SectionFrame,
};
pub use surface::{partials, Patch, ScalarSurface};
