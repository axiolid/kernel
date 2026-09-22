// SPDX-License-Identifier: MPL-2.0

//! The separating axis test.
//!
//! # The degenerate-axis trap
//!
//! Candidate axes include cross products of edge pairs. When two edges are
//! nearly parallel their cross product is near zero, and normalising it
//! amplifies rounding into a direction that points essentially anywhere.
//! Projecting onto such an axis can report a separation of ~1e-16 for two
//! shapes that genuinely touch, which turns a clearance check into a
//! coin flip.
//!
//! Those axes are skipped by squared length before normalisation. Skipping is
//! safe: a degenerate cross product carries no separating information that
//! the two edges' own face normals do not already carry.

use axiolid_core::{Point3, Scalar, Vec3};

use crate::{project, ConvexShape};

/// The axis along which two shapes are furthest apart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// Unit direction of least overlap, pointing from `a` towards `b`.
    pub axis: Vec3,
    /// Gap along `axis`.
    pub distance: Scalar,
}

/// Outcome of a separating-axis query.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum SeparationResult {
    /// A separating axis exists: the shapes are disjoint.
    Apart {
        /// Unit direction that separates them.
        axis: Vec3,
        /// Gap along that axis. The largest gap found, which is the
        /// shapes' true separation distance for convex operands.
        distance: Scalar,
    },
    /// No separating axis exists: the shapes share at least one point.
    ///
    /// No depth is reported. See the crate docs: penetration depth and
    /// contact extent are different measurements and this crate refuses to
    /// conflate them.
    Overlapping,
}

/// Find a separating axis between two convex shapes, or prove none exists.
///
/// `tolerance` is the gap below which two projections count as touching
/// rather than apart, and the length below which a candidate axis is treated
/// as degenerate and skipped.
#[must_use]
pub fn separation(a: &ConvexShape, b: &ConvexShape, tolerance: Scalar) -> SeparationResult {
    if a.is_empty() || b.is_empty() {
        // An empty shape occupies no space, so it cannot overlap anything.
        // Reporting a direction would be inventing one.
        return SeparationResult::Apart {
            axis: Vec3::new(1.0, 0.0, 0.0),
            distance: Scalar::INFINITY,
        };
    }

    let mut best_axis = Vec3::ZERO;
    let mut best_gap = Scalar::NEG_INFINITY;

    for axis in candidate_axes(a, b, tolerance) {
        let (a_min, a_max) = project(&a.points, axis);
        let (b_min, b_max) = project(&b.points, axis);
        // Gap is positive when the projections are disjoint. Take the larger
        // of the two orderings so the sign convention is "b beyond a".
        let gap = (b_min - a_max).max(a_min - b_max);
        if gap > best_gap {
            best_gap = gap;
            best_axis = if b_min - a_max >= a_min - b_max {
                axis
            } else {
                -axis
            };
        }
        // A single axis with a real gap proves disjointness; no need to test
        // the rest.
        if gap > tolerance {
            return SeparationResult::Apart {
                axis: best_axis,
                distance: gap,
            };
        }
    }

    if best_gap > tolerance {
        SeparationResult::Apart {
            axis: best_axis,
            distance: best_gap,
        }
    } else {
        SeparationResult::Overlapping
    }
}

/// Every axis worth testing: both shapes' face normals, plus edge-pair cross
/// products, each normalised and filtered for degeneracy.
fn candidate_axes(a: &ConvexShape, b: &ConvexShape, tolerance: Scalar) -> Vec<Vec3> {
    // Below this squared length a cross product is numerical noise rather
    // than a direction.
    let floor = (tolerance * tolerance).max(1e-24);
    let mut axes = Vec::new();

    let push = |v: Vec3, axes: &mut Vec<Vec3>| {
        let length_squared = v.length_squared();
        if length_squared > floor {
            axes.push(v / length_squared.sqrt());
        }
    };

    for &normal in a.normals.iter().chain(b.normals.iter()) {
        push(normal, &mut axes);
    }
    for &edge_a in &a.edges {
        for &edge_b in &b.edges {
            push(edge_a.cross(edge_b), &mut axes);
        }
    }
    // With no face normals and no usable edge pairs (two coincident points,
    // say) fall back to the line between the shapes so the query still has an
    // answer rather than silently reporting overlap.
    if axes.is_empty() {
        let direction = centroid(&b.points) - centroid(&a.points);
        push(direction, &mut axes);
    }
    axes
}

fn centroid(points: &[Point3]) -> Point3 {
    let mut sum = Vec3::ZERO;
    for point in points {
        sum += *point;
    }
    sum / points.len() as Scalar
}
