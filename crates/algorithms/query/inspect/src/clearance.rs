//! Closest approach between two meshes.
//!
//! Clearance is the clash-detection question: not "do these collide" but "if
//! not, by how much do they miss". A boolean answers the first and discards
//! the second.
//!
//! The BVH supplies candidate triangle pairs within the search length; the
//! exact distance for each surviving pair is computed directly. Distance is
//! a measurement, not a predicate, so f64 is the right currency here and
//! ADR 0045's reasoning applies: the value is reported, not used to decide
//! topology.

use axiolid_core::{Aabb, Point3, Vec3};
use axiolid_guarantees::Sign;
use axiolid_mesh::TriMesh;
use axiolid_predicates::orient3d;
use axiolid_spatial::{Bvh, SpatialItem};

/// Smallest distance between the surfaces of `a` and `b`.
///
/// Returns `None` when nothing is found within `search_length`, and
/// `Some(0.0)` when the surfaces touch or interpenetrate. A caller doing
/// clash detection reads those as "clear beyond the search radius" and
/// "clash" respectively, and any positive value as the remaining clearance.
///
/// # Panics
///
/// Never. A malformed index buffer yields `None` rather than a panic.
pub fn min_gap(a: &TriMesh, b: &TriMesh, search_length: f64) -> Option<f64> {
    if !search_length.is_finite() || search_length <= 0.0 {
        return None;
    }
    // Interpenetration is a clash even when no triangle pair crosses: a
    // corner of one solid can sit strictly inside the other while every
    // surface-to-surface distance stays positive. Surface distance alone
    // would report that as clearance, which is exactly backwards. Reuse the
    // exact containment query rather than inventing a second answer.
    if any_vertex_inside(a, b) || any_vertex_inside(b, a) {
        return Some(0.0);
    }
    surface_gap(a, b, search_length)
}

/// Smallest distance between the surfaces of `a` and `b`, when it is at most
/// `search_length`.
///
/// Unlike [`min_gap`], containment is not tested: a mesh in another's cavity,
/// or wholly inside it, is as far from it as the surfaces are.
/// `Some(0.0)` when the surfaces touch or cross. `None` when nothing is within
/// `search_length`, or `search_length` is negative or not finite (zero is
/// allowed: only touching or crossing surfaces are found).
///
/// # Panics
///
/// Never. A malformed index buffer yields `None` rather than a panic.
pub fn surface_gap(a: &TriMesh, b: &TriMesh, search_length: f64) -> Option<f64> {
    if !search_length.is_finite() || search_length < 0.0 {
        return None;
    }
    let b_index = build_index(b, search_length);
    let mut best: Option<f64> = None;
    let mut candidates = Vec::new();

    for a_triangle in triangles(a) {
        let probe = grow(bounds_of(&a_triangle), search_length);
        candidates.clear();
        b_index.query_aabb(&probe, &mut candidates);
        for &candidate in &candidates {
            let Some(item) = b_index.item(candidate) else {
                continue;
            };
            let Some(b_triangle) = triangle_at(b, item.key as usize) else {
                continue;
            };
            let distance = triangle_distance(&a_triangle, &b_triangle);
            if best.is_none_or(|current| distance < current) {
                best = Some(distance);
            }
            if distance == 0.0 {
                return Some(0.0);
            }
        }
    }
    best.filter(|distance| *distance <= search_length)
}

/// Distance between two triangles, zero when they touch or cross.
///
/// Every closest pair of convex sets is realised on a boundary feature, so
/// checking all edge pairs and all vertex-to-face distances is exhaustive
/// for triangles.
fn triangle_distance(a: &[Point3; 3], b: &[Point3; 3]) -> f64 {
    // An edge through the other's interior has no closest pair on the
    // boundary features below: both sides of the crossing are positive
    // distances away (#390).
    if edge_crosses(a, b) || edge_crosses(b, a) {
        return 0.0;
    }
    let mut best = f64::INFINITY;
    for i in 0..3 {
        for j in 0..3 {
            let d = segment_distance(a[i], a[(i + 1) % 3], b[j], b[(j + 1) % 3]);
            best = best.min(d);
        }
    }
    for vertex in a {
        best = best.min(point_triangle_distance(*vertex, b));
    }
    for vertex in b {
        best = best.min(point_triangle_distance(*vertex, a));
    }
    best
}

/// Whether an edge of `edges` crosses or touches the triangle `face`.
///
/// Certified `orient3d` signs, as in `containment.rs`: the edge's ends lie on
/// opposite sides of, or on, the face's plane, and the three tetrahedra the
/// edge spans with the face's edges agree in sign, zeros allowed. An edge in
/// the face's plane is left to the edge and vertex distances; a sign that
/// cannot be certified is read as no crossing, for the same reason.
fn edge_crosses(edges: &[Point3; 3], face: &[Point3; 3]) -> bool {
    (0..3).any(|i| segment_crosses(edges[i], edges[(i + 1) % 3], face))
}

fn segment_crosses(p: Point3, q: Point3, face: &[Point3; 3]) -> bool {
    let [a, b, c] = *face;
    let sign = |value: Option<Sign>| match value {
        Some(Sign::Positive) => Some(1),
        Some(Sign::Negative) => Some(-1),
        Some(Sign::Zero) => Some(0),
        _ => None,
    };
    let (Some(sp), Some(sq)) = (
        sign(orient3d(a, b, c, p).sign()),
        sign(orient3d(a, b, c, q).sign()),
    ) else {
        return false;
    };
    if sp * sq > 0 || (sp == 0 && sq == 0) {
        return false;
    }
    let mut seen = 0;
    for (from, to) in [(a, b), (b, c), (c, a)] {
        let Some(side) = sign(orient3d(p, q, from, to).sign()) else {
            return false;
        };
        if side != 0 {
            if seen != 0 && seen != side {
                return false;
            }
            seen = side;
        }
    }
    true
}

/// Distance between two segments, handling the parallel case.
fn segment_distance(p0: Point3, p1: Point3, q0: Point3, q1: Point3) -> f64 {
    let u = p1 - p0;
    let v = q1 - q0;
    let w = p0 - q0;
    let a = u.dot(u);
    let b = u.dot(v);
    let c = v.dot(v);
    let d = u.dot(w);
    let e = v.dot(w);
    let denominator = a * c - b * b;

    // Parallel segments leave the parameters underdetermined; clamping to the
    // endpoints is what makes this total rather than a special case.
    let (mut s, mut t) = if denominator.abs() <= f64::EPSILON * a.max(c).max(1.0) {
        (0.0, if c > 0.0 { e / c } else { 0.0 })
    } else {
        ((b * e - c * d) / denominator, (a * e - b * d) / denominator)
    };
    s = s.clamp(0.0, 1.0);
    t = t.clamp(0.0, 1.0);

    // Re-solve each parameter against the other's clamped value, so a
    // clamped endpoint still finds the true nearest point on its partner.
    if c > 0.0 {
        t = ((s * b - e) / c).clamp(0.0, 1.0);
    }
    if a > 0.0 {
        s = (-((t * b + d) / a)).clamp(0.0, 1.0);
    }
    ((p0 + u * s) - (q0 + v * t)).length()
}

/// Distance from a point to a triangle, including its interior.
fn point_triangle_distance(point: Point3, triangle: &[Point3; 3]) -> f64 {
    let [a, b, c] = *triangle;
    let normal = (b - a).cross(c - a);
    let area_squared = normal.length_squared();
    if area_squared > 0.0 {
        // Project onto the plane and test the barycentric signs. If the
        // projection lands inside, the perpendicular IS the distance.
        let distance = normal.dot(point - a) / area_squared.sqrt();
        let projected = point - normal * (normal.dot(point - a) / area_squared);
        let inside = [(a, b), (b, c), (c, a)]
            .iter()
            .all(|(from, to)| normal.dot((*to - *from).cross(projected - *from)) >= 0.0);
        if inside {
            return distance.abs();
        }
    }
    // Outside the projection, or a degenerate triangle: the nearest point is
    // on an edge.
    let mut best = f64::INFINITY;
    for i in 0..3 {
        best = best.min(point_segment_distance(
            point,
            triangle[i],
            triangle[(i + 1) % 3],
        ));
    }
    best
}

/// Distance from a point to a segment.
fn point_segment_distance(point: Point3, from: Point3, to: Point3) -> f64 {
    let direction = to - from;
    let length_squared = direction.length_squared();
    if length_squared == 0.0 {
        return (point - from).length();
    }
    let t = ((point - from).dot(direction) / length_squared).clamp(0.0, 1.0);
    (point - (from + direction * t)).length()
}

/// Triangles of a mesh as coordinate triples.
fn triangles(mesh: &TriMesh) -> impl Iterator<Item = [Point3; 3]> + '_ {
    (0..mesh.indices.len() / 3).filter_map(move |index| triangle_at(mesh, index))
}

/// One triangle by index, or `None` if the index buffer is malformed.
pub(crate) fn triangle_at(mesh: &TriMesh, index: usize) -> Option<[Point3; 3]> {
    let corners = mesh.indices.get(index * 3..index * 3 + 3)?;
    Some([
        *mesh.positions.get(corners[0] as usize)?,
        *mesh.positions.get(corners[1] as usize)?,
        *mesh.positions.get(corners[2] as usize)?,
    ])
}

/// Bounding box of a triangle.
fn bounds_of(triangle: &[Point3; 3]) -> Aabb {
    let mut bounds = Aabb::from_point(triangle[0]);
    bounds.extend(triangle[1]);
    bounds.extend(triangle[2]);
    bounds
}

/// Grow a box by a uniform margin, so a query catches near misses.
fn grow(bounds: Aabb, margin: f64) -> Aabb {
    let padding = Vec3::splat(margin);
    let mut grown = Aabb::from_point(bounds.min - padding);
    grown.extend(bounds.max + padding);
    grown
}

/// BVH over the triangles of a mesh, padded by the search length.
fn build_index(mesh: &TriMesh, search_length: f64) -> Bvh<u32> {
    let items = (0..mesh.indices.len() / 3).filter_map(|index| {
        let triangle = triangle_at(mesh, index)?;
        Some(SpatialItem::new(
            index as u32,
            grow(bounds_of(&triangle), search_length),
        ))
    });
    Bvh::build(items)
}

/// Whether any vertex of `probe` lies strictly inside `solid`.
fn any_vertex_inside(probe: &TriMesh, solid: &TriMesh) -> bool {
    probe
        .positions
        .iter()
        .any(|point| crate::containment::contains(solid, *point).unwrap_or(false))
}
