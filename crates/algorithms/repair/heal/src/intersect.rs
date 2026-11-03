//! Triangle-triangle self-intersection over a mesh (#73).
//!
//! # Why adjacency is the hard part
//!
//! In a closed mesh, every triangle shares an edge with three neighbours and
//! a vertex with many more. Those touch by construction. A naive predicate
//! that answers "do these two triangles share a point" reports every closed
//! mesh as broken, so adjacency is excluded structurally: pairs sharing any
//! vertex INDEX are skipped before any arithmetic runs.
//!
//! Index-based exclusion, not coordinate comparison. Two distinct vertices
//! holding equal coordinates are a `DuplicateVertex` defect, reported
//! separately; treating them as adjacent here would hide it.
//!
//! # Narrow phase
//!
//! Exact end to end: every decision is the sign of a certified `orient3d` or
//! `orient2d`, so near-parallel planes (a rotation rounds vertices slightly
//! off their common plane) are decided as the inputs stand, not by a float
//! direction (#399). `triangles_intersect` gives the method and the contact
//! rule. The coplanar branch's `orient2d` has no exact fallback once its
//! products overflow (coordinates near 1e154 and above); TODO(#404).
//!
//! # Broad phase
//!
//! Candidate pairs come from the shared `Bvh`. The index is an accelerator
//! only: `self_intersections_brute_force` computes the same answer by
//! checking every pair, and the two must agree on every input. If they ever
//! disagree, the index is wrong -- that is a bug, not a tuning parameter.

use axiolid_core::{Aabb, Point2, Point3};
use axiolid_guarantees::Sign;
use axiolid_mesh::TriangleMeshView;
use axiolid_predicates::{orient2d, orient3d};
use axiolid_spatial::{Bvh, SpatialIndex, SpatialItem};
use std::ops::ControlFlow;

/// One intersecting triangle pair, lower index first.
///
/// Ordered and deduplicated so a caller can compare two runs directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct IntersectingPair {
    /// Lower triangle index.
    pub first: u32,
    /// Higher triangle index.
    pub second: u32,
}

impl IntersectingPair {
    /// Normalise so `first < second`, giving one canonical form per pair.
    fn new(a: usize, b: usize) -> Self {
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        Self {
            first: lo as u32,
            second: hi as u32,
        }
    }
}

/// Whether two triangles share at least one vertex index.
///
/// Adjacency is decided on indices alone, so it costs no arithmetic and
/// cannot be perturbed by coordinates.
fn share_a_vertex(a: [u64; 3], b: [u64; 3]) -> bool {
    a.iter().any(|i| b.contains(i))
}

/// Sign of `orient3d`, zero when the four points are coplanar.
///
/// The certified predicate decides the sign exactly, so a point lying in a
/// plane is reported as coplanar rather than as an arbitrary side chosen by
/// rounding.
fn side(a: Point3, b: Point3, c: Point3, d: Point3) -> i32 {
    match orient3d(a, b, c, d).sign() {
        Some(Sign::Positive) => 1,
        Some(Sign::Negative) => -1,
        _ => 0,
    }
}

/// Whether two non-adjacent triangles intersect, decided exactly.
///
/// Plane-side rejection alone is not sufficient. Two triangles can each
/// straddle the other's plane and still miss entirely -- the cube's opposite
/// faces do exactly that -- so passing the rejection test means only "not
/// separated by these two planes", which is not the same as intersecting.
///
/// # Method
///
/// Every decision is the sign of a certified `orient3d` or `orient2d`; no
/// direction, distance or intersection point is ever computed in floating
/// point. A float direction `n_p x n_q` is cancellation noise for nearly
/// parallel planes, so none is used (#399).
///
/// 1. The six vertex-against-plane signs reject triangles wholly on one side
///    of the other's plane.
/// 2. Otherwise the planes are distinct and meet in a line `L`. Each
///    triangle meets `L` in an interval, so the triangles meet exactly when
///    an end of one interval lies in the other triangle, and that end is on
///    an edge. So the pair intersects if and only if some edge of one
///    triangle meets the closed other triangle (`edge_meets_triangle`).
///    An edge lying in the other's plane is skipped: it lies on `L`, and the
///    interval end it carries is a vertex, reached through the next edge.
/// 3. Exactly coplanar pairs (all six signs zero) are projected onto a
///    coordinate plane that keeps the triangle non-degenerate, which is
///    exact, and decided with `orient2d` (`coplanar_triangles_overlap`).
///
/// Contact: outside the coplanar case the intersection is closed, so a vertex
/// or edge of one triangle on the other counts. Coplanar triangles are
/// reported only when they share interior area; sharing an edge or a vertex
/// is touching, not overlapping. A triangle of zero area has no interior and
/// is never reported.
fn triangles_intersect(p: [Point3; 3], q: [Point3; 3]) -> bool {
    let Some(p_sides) = plane_sides(q, p) else {
        return coplanar_triangles_overlap(p, q);
    };
    if separated(p_sides) {
        return false;
    }
    let Some(q_sides) = plane_sides(p, q) else {
        return coplanar_triangles_overlap(p, q);
    };
    if separated(q_sides) {
        return false;
    }
    (0..3).any(|i| edge_meets_triangle(p[i], p[(i + 1) % 3], [p_sides[i], p_sides[(i + 1) % 3]], q))
        || (0..3).any(|i| {
            edge_meets_triangle(q[i], q[(i + 1) % 3], [q_sides[i], q_sides[(i + 1) % 3]], p)
        })
}

/// Signed side of each vertex of `t` against the plane of `plane`.
fn plane_sides(plane: [Point3; 3], t: [Point3; 3]) -> Option<[i32; 3]> {
    let sides = t.map(|v| side(plane[0], plane[1], plane[2], v));
    // All three in the plane: the pair is coplanar, which the caller decides
    // in two dimensions.
    sides.iter().any(|&s| s != 0).then_some(sides)
}

/// Whether every vertex lies strictly on one side.
fn separated(sides: [i32; 3]) -> bool {
    sides.iter().all(|&s| s > 0) || sides.iter().all(|&s| s < 0)
}

/// Whether the closed segment `a`-`b` meets the closed triangle `t`.
///
/// `ends` are the signs of `a` and `b` against the plane of `t`. The segment
/// reaches the plane unless both ends are strictly on one side. A segment in
/// the plane is skipped (see [`triangles_intersect`]). Otherwise it crosses
/// the plane in one point `X`, and `X` lies in the closed triangle when the
/// three `orient3d(a, b, t_i, t_j)` signs, which are the orientations of `X`
/// against the three edges up to one common factor, are not both positive
/// and negative (the Guigue-Devillers sign test).
fn edge_meets_triangle(a: Point3, b: Point3, ends: [i32; 2], t: [Point3; 3]) -> bool {
    if ends[0] * ends[1] > 0 || (ends[0] == 0 && ends[1] == 0) {
        return false;
    }
    let around = [
        side(a, b, t[0], t[1]),
        side(a, b, t[1], t[2]),
        side(a, b, t[2], t[0]),
    ];
    !(around.iter().any(|&s| s > 0) && around.iter().any(|&s| s < 0))
}

/// Triangle corner positions and vertex indices, if the triangle is usable.
fn triangle_of<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    index: usize,
) -> Option<([Point3; 3], [u64; 3])> {
    let corners = mesh.triangle(index);
    let positions = corners.map(|j| mesh.position(j as usize));
    if positions.iter().any(|p| !p.is_finite()) {
        return None;
    }
    Some((positions, corners))
}

/// Self-intersecting triangle pairs, found through the spatial index.
///
/// Returns pairs in sorted order. An empty result means no pair intersects.
///
/// There is no tolerance parameter: the narrow phase decides with certified
/// `orient3d`, and the broad phase uses exact triangle bounds. A tolerance
/// here would only be able to make the answer *wrong*, by admitting or
/// rejecting pairs the exact test already decides.
#[must_use]
pub fn self_intersections<M: TriangleMeshView + ?Sized>(mesh: &M) -> Vec<IntersectingPair> {
    let count = mesh.triangle_count();
    let mut items = Vec::with_capacity(count);
    for index in 0..count {
        if let Some((positions, _)) = triangle_of(mesh, index) {
            // Pad by the caller's tolerance so the broad phase never rejects
            // a pair the exact narrow phase would have accepted.
            let mut bounds = Aabb::from_point(positions[0]);
            bounds.extend(positions[1]);
            bounds.extend(positions[2]);
            items.push(SpatialItem::new(index as u32, bounds));
        }
    }
    let bvh = Bvh::build(items);

    let mut found = Vec::new();
    for index in 0..count {
        let Some((positions, corners)) = triangle_of(mesh, index) else {
            continue;
        };
        let mut query = Aabb::from_point(positions[0]);
        query.extend(positions[1]);
        query.extend(positions[2]);
        bvh.visit_aabb(&query, &mut |other: &u32| {
            let other = *other as usize;
            // Each unordered pair is decided once, by its lower index.
            if other <= index {
                return ControlFlow::Continue(());
            }
            if let Some((other_positions, other_corners)) = triangle_of(mesh, other) {
                if !share_a_vertex(corners, other_corners)
                    && triangles_intersect(positions, other_positions)
                {
                    found.push(IntersectingPair::new(index, other));
                }
            }
            ControlFlow::Continue(())
        });
    }
    found.sort_unstable();
    found.dedup();
    found
}

/// The same answer without the spatial index, by checking every pair.
///
/// This exists to be compared against [`self_intersections`]. The index is an
/// optimisation, and an optimisation that changes the answer is a bug, so the
/// reference is kept in production code rather than in a test where it could
/// drift out of sync with the accelerated path.
#[must_use]
pub fn self_intersections_brute_force<M: TriangleMeshView + ?Sized>(
    mesh: &M,
) -> Vec<IntersectingPair> {
    let count = mesh.triangle_count();
    let mut found = Vec::new();
    for index in 0..count {
        let Some((positions, corners)) = triangle_of(mesh, index) else {
            continue;
        };
        for other in (index + 1)..count {
            let Some((other_positions, other_corners)) = triangle_of(mesh, other) else {
                continue;
            };
            if !share_a_vertex(corners, other_corners)
                && triangles_intersect(positions, other_positions)
            {
                found.push(IntersectingPair::new(index, other));
            }
        }
    }
    found.sort_unstable();
    found.dedup();
    found
}

/// Whether two coplanar triangles share interior area, decided exactly.
///
/// Both are projected onto a coordinate plane in which the first triangle has
/// non-zero area, so neither projection degenerates (both lie in one plane).
/// Dropping a coordinate is exact, and every test below is `orient2d`: no
/// tolerance and no constructed point.
///
/// Two convex regions have disjoint interiors if and only if an edge line of
/// one has the other wholly on its closed outer side (separating axis). That
/// is decided for the edges of both triangles, so the shapes need not cross
/// properly: equal triangles, a vertex on an edge and containment all overlap,
/// while triangles sharing only an edge or a vertex are separated by that
/// edge and are NOT reported.
fn coplanar_triangles_overlap(p: [Point3; 3], q: [Point3; 3]) -> bool {
    for axis in 0..3 {
        let pp = project_triangle(p, axis);
        let p_turn = orient_sign(pp[0], pp[1], pp[2]);
        if p_turn == 0 {
            continue;
        }
        let qq = project_triangle(q, axis);
        let q_turn = orient_sign(qq[0], qq[1], qq[2]);
        if q_turn == 0 {
            return false;
        }
        return !has_separating_edge(pp, p_turn, qq) && !has_separating_edge(qq, q_turn, pp);
    }
    false
}

/// Drop `axis`, keeping the other two coordinates in cyclic order.
fn project_triangle(t: [Point3; 3], axis: usize) -> [Point2; 3] {
    t.map(|p| match axis {
        0 => Point2::new(p.y, p.z),
        1 => Point2::new(p.z, p.x),
        _ => Point2::new(p.x, p.y),
    })
}

/// Sign of `orient2d` (zero when the three points are collinear).
fn orient_sign(a: Point2, b: Point2, c: Point2) -> i32 {
    match orient2d(a, b, c).sign() {
        Some(Sign::Positive) => 1,
        Some(Sign::Negative) => -1,
        _ => 0,
    }
}

/// Whether some edge line of `t` (turning `turn`, +1 or -1) has every vertex
/// of `other` on its outer side or on it.
fn has_separating_edge(t: [Point2; 3], turn: i32, other: [Point2; 3]) -> bool {
    (0..3).any(|i| {
        other
            .iter()
            .all(|&v| orient_sign(t[i], t[(i + 1) % 3], v) * turn <= 0)
    })
}
