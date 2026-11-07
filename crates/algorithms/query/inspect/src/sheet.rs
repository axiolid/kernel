//! Two faces of one mesh that face each other within a linear distance.
//!
//! A boolean leaves such a thin sheet where an operand's face lands a
//! rounding residue off the face it should have met: the result is a Solid
//! whose lid is thinner than its coordinates can carry (#400). The query is
//! the mesh's own, so it needs no second operand.

use axiolid_core::{Point3, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_spatial::SpatialItem;

use crate::clearance::{build_index, triangle_at, triangle_distance};
use crate::coplanar::corner_distances;

/// A pair of triangles of one mesh that face each other within a linear
/// distance.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThinSheet {
    /// The lower triangle index.
    pub first: usize,
    /// The higher triangle index.
    pub second: usize,
    /// The largest distance from a corner of either triangle to the other's
    /// plane.
    pub offset: f64,
}

/// The first pair of triangles of `mesh` that form a sheet thinner than
/// `linear`, in ascending `(first, second)` order with `first < second`.
///
/// A pair qualifies when both triangles are non-degenerate, their unit
/// normals point opposite ways, the distance between them is at most
/// `linear`, every corner of each lies within `linear` of the other's plane,
/// and, projected onto the first's plane, they overlap in a polygon whose
/// doubled area exceeds `linear` times its perimeter: contact along a single
/// edge is not a sheet. Exactly coplanar pairs count. `None` when no pair
/// qualifies, and when `linear` is negative or not finite.
///
/// # Panics
///
/// Never. A malformed index buffer yields `None` rather than a panic.
pub fn thin_sheet(mesh: &TriMesh, linear: f64) -> Option<ThinSheet> {
    if !linear.is_finite() || linear < 0.0 {
        return None;
    }
    let index = build_index(mesh, linear);
    let mut candidates = Vec::new();
    for first in 0..mesh.indices.len() / 3 {
        let Some(first_triangle) = triangle_at(mesh, first) else {
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
            .filter(|&second| second > first)
            .collect();
        seconds.sort_unstable();
        for second in seconds {
            let Some(second_triangle) = triangle_at(mesh, second) else {
                continue;
            };
            if let Some(offset) = facing_within(&first_triangle, &second_triangle, linear) {
                return Some(ThinSheet {
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
fn facing_within(a: &[Point3; 3], b: &[Point3; 3], linear: f64) -> Option<f64> {
    let to_a = corner_distances(a, b)?;
    let to_b = corner_distances(b, a)?;
    let offset = to_a.max(to_b);
    if offset > linear {
        return None;
    }
    let na = (a[1] - a[0]).cross(a[2] - a[0]).normalize();
    let nb = (b[1] - b[0]).cross(b[2] - b[0]).normalize();
    if na.dot(nb) >= 0.0 || triangle_distance(a, b) > linear {
        return None;
    }
    let (area2, perimeter) = overlap(a, b, na);
    (area2 > linear * perimeter).then_some(offset)
}

/// Twice the area and the perimeter of the polygon in which `b`, projected
/// onto the plane of `a` (unit normal `n`), overlaps `a`.
fn overlap(a: &[Point3; 3], b: &[Point3; 3], n: Vec3) -> (f64, f64) {
    let u = (a[1] - a[0]).normalize();
    let v = n.cross(u);
    let flat = |p: Point3| [(p - a[0]).dot(u), (p - a[0]).dot(v)];
    let clip = a.map(flat);
    let mut polygon: Vec<[f64; 2]> = b.iter().map(|&p| flat(p)).collect();
    // `b` faces away from `a`: its projection winds the other way.
    if signed_area2(&polygon) < 0.0 {
        polygon.reverse();
    }
    for i in 0..3 {
        let (p, q) = (clip[i], clip[(i + 1) % 3]);
        let side = |x: [f64; 2]| (q[0] - p[0]) * (x[1] - p[1]) - (q[1] - p[1]) * (x[0] - p[0]);
        let mut kept = Vec::with_capacity(polygon.len() + 1);
        for j in 0..polygon.len() {
            let (s, e) = (polygon[j], polygon[(j + 1) % polygon.len()]);
            let (ds, de) = (side(s), side(e));
            if ds >= 0.0 {
                kept.push(s);
            }
            if (ds > 0.0 && de < 0.0) || (ds < 0.0 && de > 0.0) {
                let t = ds / (ds - de);
                kept.push([s[0] + t * (e[0] - s[0]), s[1] + t * (e[1] - s[1])]);
            }
        }
        polygon = kept;
        if polygon.len() < 3 {
            return (0.0, 0.0);
        }
    }
    let perimeter = (0..polygon.len())
        .map(|i| {
            let (s, e) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            (e[0] - s[0]).hypot(e[1] - s[1])
        })
        .sum();
    (signed_area2(&polygon).abs(), perimeter)
}

fn signed_area2(polygon: &[[f64; 2]]) -> f64 {
    (0..polygon.len())
        .map(|i| {
            let (s, e) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            s[0] * e[1] - e[0] * s[1]
        })
        .sum()
}
