#![forbid(unsafe_code)]

//! Certified intersections for linear geometry.
//!
//! This package exists so a line-query application can compile intersection
//! logic with a four-package internal closure — `axiolid-core`,
//! `axiolid-guarantees`, `axiolid-linear`, `axiolid-predicates` — and no
//! curves, surfaces, NURBS, meshes, B-rep, topology, providers, or execution
//! machinery (ADR 0036). The `linear-intersection-minimal` profile in
//! `docs/architecture/closure-profiles.toml` pins that closure, so any new
//! internal dependency fails `cargo xtask architecture closure check`.
//!
//! # Classification, not `Option<Point2>`
//!
//! Collapsing the answer to an optional point destroys facts a topology or
//! rule checker needs: crossing, endpoint contact, parallel-disjoint,
//! coincident, collinear-disjoint, and overlap are distinct outcomes and are
//! returned as distinct variants.
//!
//! # Certified topology
//!
//! Parallel-versus-coincident and on-segment decisions come from
//! `axiolid-predicates`, which escalates to exact arithmetic. The caller's
//! [`axiolid_core::Tolerance`] governs residual acceptance of the computed
//! coordinate, never the topological branch.

mod error;
mod line_line;
mod segment_segment;
mod validate;

/// Re-exported so a caller can name the required tolerance policy without a
/// separate dependency on `axiolid-core`.
pub use axiolid_core::Tolerance;
pub use error::{InputSide, LinearIntersectionError};
pub use line_line::{line_line2, LineLineIntersection2};
pub use segment_segment::{segment_segment2, SegmentSegmentIntersection2};
