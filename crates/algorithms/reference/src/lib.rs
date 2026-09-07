#![forbid(unsafe_code)]

//! Portable scalar reference implementation: the correctness oracle (ADR 0012).
//!
//! No intrinsics, no threading, no feature gates. Every optimized backend is
//! validated by differential test against this crate, so it must stay readable
//! and obviously correct in preference to being fast.
//!
//! This package is a **convenience umbrella** (ADR 0036). Certified predicates
//! now live in the focused `axiolid-predicates` package and are re-exported
//! here unchanged, so an existing `axiolid_reference::orient2d` caller is
//! unaffected while a narrow consumer can depend on the substrate directly
//! instead of acquiring this package's whole dependency graph.

pub mod assemble;
pub mod boolean;
pub mod clash;
pub mod convex_hull;
pub mod intersection;
pub mod polygon;
pub mod primitive;
pub mod retriangulate;
pub mod section;
pub mod segment_triangle;
pub mod tessellate;
pub mod triangle_triangle;

/// Analytic and spline evaluation moved to the focused `axiolid-evaluate`
/// package (ADR 0036). Re-exported unchanged so existing
/// `axiolid_reference::curve::*` and `::surface::*` callers are unaffected.
pub use axiolid_evaluate::{curve, surface};

pub use axiolid_predicates::{
    arithmetic, expansion, orient3, orientation, scene, sphere, static_filter,
};

pub use assemble::exact_boolean;
pub use axiolid_evaluate::{
    derivative2, derivative3, evaluate2, evaluate3, flatten2, partials, Patch, ScalarCurve,
    ScalarSurface,
};
pub use axiolid_predicates::{
    incircle, incircle_filter, insphere, insphere_filter, orient2d, orient2d_filter, orient3d,
    orient3d_filter, two_diff, two_product, two_sum, StaticFilter,
};
pub use boolean::ScalarBoolean;
pub use convex_hull::{minimum_area_rectangle, strict_convex_hull, OrientedRectangle2};
pub use intersection::{
    assemble_polylines, intersection_segments, EdgeKey, IntersectionCurve, IntersectionSegment,
    NodeKey, Operand, Polyline,
};
pub use polygon::{ring_orientation, signed_area2, triangulate_simple};
pub use retriangulate::{retriangulate_face, FacePatch};
pub use section::ScalarSection;
pub use segment_triangle::{segment_triangle_relation, SegmentTriangleRelation};
pub use triangle_triangle::{triangle_triangle_relation, TriangleTriangleRelation};
