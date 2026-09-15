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

pub mod arc_length;
pub mod curve;
pub mod frenet;
pub mod intrinsic_relation;
mod nurbs;
pub mod surface;

pub use arc_length::{elevated_point, elevated_tangent, intrinsic_point, intrinsic_tangent};
pub use curve::{derivative2, derivative3, evaluate2, evaluate3, flatten2, ScalarCurve};
pub use frenet::{frenet_frame, frenet_point, frenet_tangent};
pub use intrinsic_relation::{join_intrinsic3, offset_intrinsic3, trim_intrinsic3};
pub use surface::{partials, Patch, ScalarSurface};
