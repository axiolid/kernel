// SPDX-License-Identifier: MPL-2.0
#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Convex collision queries via the separating axis theorem.
//!
//! # What this answers, and what it deliberately does not
//!
//! Two convex shapes are disjoint if and only if some axis exists on which
//! their projections do not overlap. This crate finds such an axis, or proves
//! none exists.
//!
//! When the shapes are apart it also reports **how far** apart, because that
//! is the question a rule check asks: not "do these collide" but "is there
//! enough clearance". [`axiolid_inspect`]'s mesh clearance answers the same
//! question for triangle soups; this answers it for convex shapes without
//! building an index.
//!
//! It does **not** report penetration depth. That is a deliberate refusal,
//! consistent with `axiolid-measure`: a contact area and an interpenetration
//! depth are different measurements, and collapsing them behind one number
//! is how a caller ends up using the wrong one. Callers needing penetration
//! depth for physics want EPA, which belongs in a physics engine rather than
//! a certified geometry kernel.
//!
//! # Why SAT rather than GJK
//!
//! For the shapes a model-checking rule actually has -- boxes, extruded
//! profiles, hulls with tens of vertices -- SAT is direct, has no iteration
//! count to tune, and its failure mode is a clean answer rather than a
//! tolerance-dependent one. GJK wins on high vertex counts and on a generic
//! support function; neither is the common case here, and published
//! commercial model-checking geometry APIs ship convex hull and SAT-style
//! tests without GJK/EPA at all.
//!
//! # Robustness
//!
//! Projections are compared with an explicit tolerance rather than exactly.
//! An axis derived from a cross product of two nearly parallel edges is
//! numerically meaningless, so such axes are skipped instead of being
//! allowed to report a spurious separation of ~1e-17.

use axiolid_core::{Box3, Point3, Scalar, Vec3};

mod sat;
mod shape;

pub use sat::{separation, Contact, SeparationResult};
pub use shape::ConvexShape;

/// Shortest distance between two convex shapes, or zero if they overlap.
///
/// Convenience over [`separation`] for callers that only need the number.
#[must_use]
pub fn distance(a: &ConvexShape, b: &ConvexShape, tolerance: Scalar) -> Scalar {
    match separation(a, b, tolerance) {
        SeparationResult::Apart { distance, .. } => distance,
        SeparationResult::Overlapping => 0.0,
    }
}

/// Whether two convex shapes share any point.
#[must_use]
pub fn intersects(a: &ConvexShape, b: &ConvexShape, tolerance: Scalar) -> bool {
    matches!(separation(a, b, tolerance), SeparationResult::Overlapping)
}

/// Whether two oriented boxes overlap.
///
/// The case `Box3` could not answer before: it carried an orientation but had
/// no intersection test, so a caller had to fall back to an axis-aligned
/// bound and lose the tightness the orientation was for.
#[must_use]
pub fn boxes_intersect(a: &Box3, b: &Box3, tolerance: Scalar) -> bool {
    intersects(
        &ConvexShape::from_points(&a.corners()),
        &ConvexShape::from_points(&b.corners()),
        tolerance,
    )
}

/// Whether a point lies inside or on a convex shape.
#[must_use]
pub fn contains_point(shape: &ConvexShape, point: Point3, tolerance: Scalar) -> bool {
    intersects(shape, &ConvexShape::from_points(&[point]), tolerance)
}

/// Axis-aligned extent of a shape along a direction.
pub(crate) fn project(points: &[Point3], axis: Vec3) -> (Scalar, Scalar) {
    let mut min = Scalar::INFINITY;
    let mut max = Scalar::NEG_INFINITY;
    for point in points {
        let value = point.dot(axis);
        min = min.min(value);
        max = max.max(value);
    }
    (min, max)
}
