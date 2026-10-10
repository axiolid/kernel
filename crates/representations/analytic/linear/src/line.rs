//! Unbounded parametric lines.

use axiolid_core::{Point2, Point3, Vec2, Vec3};

/// Infinite parametric line `origin + t * direction`.
///
/// The direction is deliberately not normalised on construction: import
/// adapters preserve authored vectors, and normalising here would silently
/// change the parameterisation a caller reasons about.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Line<P, V> {
    /// Point at parameter zero.
    pub origin: P,
    /// Parameter direction; import adapters may preserve a non-unit vector.
    pub direction: V,
}

/// Two-dimensional line.
pub type Line2 = Line<Point2, Vec2>;
/// Three-dimensional line.
pub type Line3 = Line<Point3, Vec3>;
