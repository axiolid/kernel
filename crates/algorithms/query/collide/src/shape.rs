// SPDX-License-Identifier: MPL-2.0

//! The convex shape a query runs against.

use axiolid_core::{Point3, Vec3};

/// A convex shape as its vertices and face normals.
///
/// Storing normals alongside points is what keeps SAT honest: the candidate
/// axes for two polyhedra are their face normals plus the cross products of
/// their edge pairs. Deriving normals from a point cloud every call would
/// repeat work the caller usually already has.
///
/// A shape with no normals still works -- a point cloud's hull is implied by
/// the edge-pair axes -- but supplying them makes the test exact for
/// polyhedra rather than conservative.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvexShape {
    pub(crate) points: Vec<Point3>,
    pub(crate) normals: Vec<Vec3>,
    pub(crate) edges: Vec<Vec3>,
}

impl ConvexShape {
    /// A shape from its vertices alone.
    ///
    /// Face normals are not inferred: computing a hull here would make a
    /// cheap constructor expensive and would silently discard the caller's
    /// own topology. Use [`ConvexShape::with_faces`] when face normals are
    /// available.
    #[must_use]
    pub fn from_points(points: &[Point3]) -> Self {
        let edges = derive_edges(points);
        Self {
            points: points.to_vec(),
            normals: Vec::new(),
            edges,
        }
    }

    /// A shape with known face normals and edge directions.
    #[must_use]
    pub fn with_faces(points: &[Point3], normals: &[Vec3]) -> Self {
        let edges = derive_edges(points);
        Self {
            points: points.to_vec(),
            normals: normals.to_vec(),
            edges,
        }
    }

    /// An axis-aligned box as a convex shape.
    #[must_use]
    pub fn from_aabb(min: Point3, max: Point3) -> Self {
        let corners = [
            Point3::new(min.x, min.y, min.z),
            Point3::new(max.x, min.y, min.z),
            Point3::new(max.x, max.y, min.z),
            Point3::new(min.x, max.y, min.z),
            Point3::new(min.x, min.y, max.z),
            Point3::new(max.x, min.y, max.z),
            Point3::new(max.x, max.y, max.z),
            Point3::new(min.x, max.y, max.z),
        ];
        Self::with_faces(
            &corners,
            &[
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ],
        )
    }

    /// The shape's vertices.
    #[must_use]
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// Whether the shape has no vertices.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

/// Edge directions between consecutive and diagonal vertex pairs.
///
/// Capped deliberately. The full O(n^2) edge set makes SAT O(n^4) against
/// another shape, which is the wrong trade past a few dozen vertices; for the
/// box and prism cases this crate targets, the cap is never reached.
fn derive_edges(points: &[Point3]) -> Vec<Vec3> {
    const MAX_POINTS_FOR_FULL_EDGES: usize = 16;
    let mut edges = Vec::new();
    if points.len() <= MAX_POINTS_FOR_FULL_EDGES {
        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let direction = points[j] - points[i];
                if direction.length_squared() > 0.0 {
                    edges.push(direction);
                }
            }
        }
    } else {
        for window in points.windows(2) {
            let direction = window[1] - window[0];
            if direction.length_squared() > 0.0 {
                edges.push(direction);
            }
        }
    }
    edges
}
