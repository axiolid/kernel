//! Surfaces that nearly share a plane.
//!
//! Two meshes whose faces lie in one plane in exact arithmetic come out of
//! placement a rounding residue apart: each vertex rounds on its own. A
//! boolean that decides on the raw values reads those faces as crossing at a
//! sliver of an angle, or leaves a crack. This query finds such a pair: a
//! triangle of each, within the search distance of each other, whose corners
//! all lie within that distance of the other's plane, and which are not
//! exactly coplanar. Exact contact (every corner on the plane, certified)
//! is a decided case and is not reported; neither is a crossing at an angle.
//!
//! Distances are measurements (ADR 0045); only "exactly on the plane" is an
//! `orient3d` sign.

use axiolid_core::Point3;
use axiolid_guarantees::Sign;
use axiolid_mesh::TriMesh;
use axiolid_predicates::orient3d;
use axiolid_spatial::SpatialItem;

use crate::clearance::{build_index, triangle_at, triangle_distance};

/// A pair of triangles that nearly share a plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NearCoplanar {
    /// Triangle index in the first mesh.
    pub first: usize,
    /// Triangle index in the second mesh.
    pub second: usize,
    /// The largest distance from a corner of either triangle to the other's
    /// plane.
    pub offset: f64,
}

/// The first pair of triangles of `a` and `b` that nearly share a plane
/// within `linear`, in ascending `(first, second)` order.
///
/// A pair qualifies when both triangles are non-degenerate, the distance
/// between them is at most `linear`, every corner of each lies within
/// `linear` of the other's plane, and at least one corner of the second is
/// not exactly on the first's plane. `None` when no pair qualifies, and when
/// `linear` is negative or not finite (zero is allowed: then only an offset
/// of zero distance would qualify, which exact coplanarity excludes).
///
/// # Panics
///
/// Never. A malformed index buffer yields `None` rather than a panic.
pub fn near_coplanar(a: &TriMesh, b: &TriMesh, linear: f64) -> Option<NearCoplanar> {
    if !linear.is_finite() || linear < 0.0 {
        return None;
    }
    let index = build_index(b, linear);
    let mut candidates = Vec::new();
    for first in 0..a.indices.len() / 3 {
        let Some(first_triangle) = triangle_at(a, first) else {
            continue;
        };
        let probe = crate::clearance::grow(crate::clearance::bounds_of(&first_triangle), linear);
        candidates.clear();
        index.query_aabb(&probe, &mut candidates);
        let mut seconds: Vec<usize> = candidates
            .iter()
            .filter_map(|&candidate| {
                index
                    .item(candidate)
                    .map(|item: &SpatialItem<u32>| item.key as usize)
            })
            .collect();
        seconds.sort_unstable();
        for second in seconds {
            let Some(second_triangle) = triangle_at(b, second) else {
                continue;
            };
            if let Some(offset) = nearly_share_a_plane(&first_triangle, &second_triangle, linear) {
                return Some(NearCoplanar {
                    first,
                    second,
                    offset,
                });
            }
        }
    }
    None
}

/// The pair's offset when it qualifies.
fn nearly_share_a_plane(a: &[Point3; 3], b: &[Point3; 3], linear: f64) -> Option<f64> {
    let to_a = corner_distances(a, b)?;
    let to_b = corner_distances(b, a)?;
    let offset = to_a.max(to_b);
    if offset > linear || triangle_distance(a, b) > linear {
        return None;
    }
    let exactly_on_plane = b
        .iter()
        .all(|corner| orient3d(a[0], a[1], a[2], *corner).sign() == Some(Sign::Zero));
    (!exactly_on_plane).then_some(offset)
}

/// The largest distance from a corner of `corners` to the plane of `plane`,
/// or `None` if `plane` is degenerate.
fn corner_distances(plane: &[Point3; 3], corners: &[Point3; 3]) -> Option<f64> {
    let normal = (plane[1] - plane[0]).cross(plane[2] - plane[0]);
    let length = normal.length();
    if !length.is_finite() || length <= 0.0 {
        return None;
    }
    let normal = normal / length;
    Some(
        corners
            .iter()
            .map(|corner| normal.dot(*corner - plane[0]).abs())
            .fold(0.0, f64::max),
    )
}
