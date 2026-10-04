//! Flat regions of a mesh: edge-connected triangles in one plane (#235).
//!
//! The per-triangle bound of the parent module needs one triangle to cover
//! a whole piece, so along every edge between two coplanar triangles (a fan
//! across a wall face, the diagonal of a cylinder facet) the search splits
//! down to the size of the bound it wants. A flat region removes those
//! edges.
//!
//! A region `R` is grown from a seed triangle across manifold edges (shared
//! by exactly two triangles, traversed in opposite directions) to
//! neighbours whose corners lie within a tiny distance of the seed's plane
//! `H` and whose normals agree with it. `thickness` is the largest distance
//! of any corner from `H`, so every point of `R` lies within it of `H`.
//!
//! # The bound
//!
//! Project onto `H`. Every region triangle projects with positive area
//! (its normal agrees with `H`'s), and every region edge whose two
//! triangles are both in `R` has them on opposite sides. So every point of
//! the projection `U` of `R` off the projected boundary edges `B` (edges
//! with one triangle in `R`) has a neighbourhood in `U`: inside a
//! triangle, on an interior edge (two triangles, opposite sides), or at a
//! vertex all of whose edges are interior (its triangles close a fan of
//! positive turn around it). For a convex piece `P` that misses `B`,
//! `P ∩ U` is then closed and open in `P`, so `P` lies in `U` as soon as one
//! of its points does.
//!
//! For `x` in such a piece, its projection `x'` is the projection of some
//! `y` in `R`, and `|y - x'| <= thickness`, so `d(x, R) <= d(x, H) +
//! thickness`. `d(., H)` is convex, so over the piece it is at most its
//! corners' largest. The test is conservative under rounding: a piece
//! counts as missing `B` only when it is further than a margin from every
//! boundary edge, and its centre as inside `U` only by a barycentric test
//! whose slack is far below that margin.

use std::collections::HashMap;

use axiolid_core::{Point2, Point3, Scalar, Vec3};

/// A flat region and its projection.
pub(super) struct Region {
    origin: Point3,
    e1: Vec3,
    e2: Vec3,
    normal: Vec3,
    /// Every point of the region is within this of its plane.
    thickness: Scalar,
    /// Projected edges with one region triangle.
    boundary: Vec<[Point2; 2]>,
    /// A piece closer than this to a boundary edge is not trusted inside.
    margin: Scalar,
}

/// The flat regions of a set of triangles, and each triangle's place in
/// them.
pub(super) struct Flat {
    regions: Vec<Region>,
    /// Region of each triangle, `u32::MAX` for none.
    region_of: Vec<u32>,
    /// Each region triangle projected onto its region's plane.
    projected: Vec<[Point2; 3]>,
}

/// Corners further than this fraction of the mesh's extent from a seed's
/// plane do not join its region.
const COPLANAR: Scalar = 1e-9;
/// Relative slack of the barycentric inside test.
const INSIDE_SLACK: Scalar = 1e-12;

impl Flat {
    /// Grow regions over `triangles`, whose corner indices are `ids`.
    #[cfg(test)]
    pub(super) fn new(triangles: &[[Point3; 3]], ids: &[[u32; 3]]) -> Self {
        Self::with_flatness(triangles, ids, 0.0)
    }

    /// [`Self::new`], also joining corners up to `flatness` from a seed's
    /// plane (#252). The bound `d(x, R) <= d(x, H) + thickness` holds for
    /// any thickness, which is measured, so a looser join only adds what
    /// it lets in: the two triangles of a slightly twisted tube quad then
    /// cover a piece together instead of forcing the search down to the
    /// width of their diagonal.
    pub(super) fn with_flatness(
        triangles: &[[Point3; 3]],
        ids: &[[u32; 3]],
        flatness: Scalar,
    ) -> Self {
        let count = triangles.len();
        let (lo, hi) = triangles.iter().flatten().fold(
            (
                Vec3::splat(Scalar::INFINITY),
                Vec3::splat(Scalar::NEG_INFINITY),
            ),
            |(lo, hi), p| (lo.min(*p), hi.max(*p)),
        );
        let scale = (hi - lo)
            .length()
            .max(lo.abs().max_element())
            .max(hi.abs().max_element());
        let mut flat = Self {
            regions: Vec::new(),
            region_of: vec![u32::MAX; count],
            projected: vec![[Point2::ZERO; 3]; count],
        };
        if !(scale.is_finite() && scale > 0.0) {
            return flat;
        }
        let near = (COPLANAR * scale).max(if flatness.is_finite() { flatness } else { 0.0 });
        // Directed edge -> triangles using it; an undirected edge is
        // manifold when used once each way and by no one else.
        let mut directed: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (t, id) in ids.iter().enumerate() {
            for k in 0..3 {
                directed
                    .entry((id[k], id[(k + 1) % 3]))
                    .or_default()
                    .push(t);
            }
        }
        let across = |t: usize, k: usize| -> Option<usize> {
            let (a, b) = (ids[t][k], ids[t][(k + 1) % 3]);
            let forward = directed.get(&(a, b))?;
            let backward = directed.get(&(b, a))?;
            (forward.len() == 1 && backward.len() == 1).then(|| backward[0])
        };
        let normals: Vec<Option<Vec3>> = triangles
            .iter()
            .map(|[a, b, c]| {
                let n = (*b - *a).cross(*c - *a);
                let len = n.length();
                (len > 1e-14 * scale * scale && len.is_finite()).then(|| n / len)
            })
            .collect();
        for seed in 0..count {
            if flat.region_of[seed] != u32::MAX {
                continue;
            }
            let Some(normal) = normals[seed] else {
                continue;
            };
            let origin = triangles[seed][0];
            let id = flat.regions.len() as u32;
            let mut members = vec![seed];
            flat.region_of[seed] = id;
            let mut next = 0;
            while next < members.len() {
                let t = members[next];
                next += 1;
                for k in 0..3 {
                    let Some(u) = across(t, k) else { continue };
                    if flat.region_of[u] != u32::MAX {
                        continue;
                    }
                    let Some(nu) = normals[u] else { continue };
                    let level = triangles[u]
                        .iter()
                        .all(|p| (*p - origin).dot(normal).abs() <= near);
                    if nu.dot(normal) > 0.5 && level {
                        flat.region_of[u] = id;
                        members.push(u);
                    }
                }
            }
            if members.len() < 2 {
                // A lone triangle gains nothing over the per-triangle bound.
                flat.region_of[seed] = u32::MAX;
                continue;
            }
            let e1 = perpendicular(normal);
            let e2 = normal.cross(e1);
            let project = |p: Point3| {
                let d = p - origin;
                Point2::new(d.dot(e1), d.dot(e2))
            };
            let mut thickness: Scalar = 0.0;
            let mut boundary = Vec::new();
            for &t in &members {
                for p in triangles[t] {
                    thickness = thickness.max((p - origin).dot(normal).abs());
                }
                flat.projected[t] = triangles[t].map(project);
                for k in 0..3 {
                    let interior = across(t, k).is_some_and(|u| flat.region_of[u] == id);
                    if !interior {
                        boundary.push([flat.projected[t][k], flat.projected[t][(k + 1) % 3]]);
                    }
                }
            }
            // Rounding of the projections and distances, at the mesh's
            // scale.
            let rounding = 64.0 * Scalar::EPSILON * scale;
            flat.regions.push(Region {
                origin,
                e1,
                e2,
                normal,
                thickness: thickness + rounding,
                boundary,
                margin: 1e3 * rounding,
            });
        }
        flat
    }

    /// The best region bound over the convex piece with corners `corners`
    /// (in order around it), from the regions of the `candidates`
    /// triangles, if any beats `best`.
    pub(super) fn bound(
        &self,
        corners: &[Point3; 4],
        candidates: &[usize],
        best: Scalar,
    ) -> Scalar {
        let mut best = best;
        let mut seen: Vec<u32> = Vec::new();
        for &t in candidates {
            let r = self.region_of[t];
            if r == u32::MAX || seen.contains(&r) {
                continue;
            }
            seen.push(r);
            let region = &self.regions[r as usize];
            let off = corners
                .iter()
                .map(|c| (*c - region.origin).dot(region.normal).abs())
                .fold(0.0, Scalar::max);
            let bound = off + region.thickness;
            if bound >= best {
                continue;
            }
            let quad = corners.map(|c| {
                let d = c - region.origin;
                Point2::new(d.dot(region.e1), d.dot(region.e2))
            });
            let centre = (quad[0] + quad[1] + quad[2] + quad[3]) * 0.25;
            let inside = candidates
                .iter()
                .any(|&u| self.region_of[u] == r && in_triangle(centre, &self.projected[u]));
            if inside && !touches(&quad, &region.boundary, region.margin) {
                best = bound;
            }
        }
        best
    }
}

/// A unit vector perpendicular to unit `n`.
fn perpendicular(n: Vec3) -> Vec3 {
    let axis = if n.x.abs() <= n.y.abs() && n.x.abs() <= n.z.abs() {
        Vec3::X
    } else if n.y.abs() <= n.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    (axis - n * n.dot(axis)).normalize()
}

fn in_triangle(p: Point2, [a, b, c]: &[Point2; 3]) -> bool {
    let area = (*b - *a).perp_dot(*c - *a);
    if area.is_nan() || area <= 0.0 {
        return false;
    }
    let slack = -INSIDE_SLACK * area;
    (*b - *a).perp_dot(p - *a) >= slack
        && (*c - *b).perp_dot(p - *b) >= slack
        && (*a - *c).perp_dot(p - *c) >= slack
}

/// Whether any boundary edge comes within `margin` of the convex `quad`.
fn touches(quad: &[Point2; 4], boundary: &[[Point2; 2]], margin: Scalar) -> bool {
    let lo = quad[0].min(quad[1]).min(quad[2]).min(quad[3]) - Point2::splat(margin);
    let hi = quad[0].max(quad[1]).max(quad[2]).max(quad[3]) + Point2::splat(margin);
    boundary.iter().any(|&[a, b]| {
        let (slo, shi) = (a.min(b), a.max(b));
        if shi.x < lo.x || shi.y < lo.y || slo.x > hi.x || slo.y > hi.y {
            return false;
        }
        segment_quad_distance(a, b, quad) <= margin
    })
}

/// Distance between segment `ab` and the convex quadrilateral `quad`: zero
/// when they meet, else the least distance between a vertex of one and the
/// other's boundary.
fn segment_quad_distance(a: Point2, b: Point2, quad: &[Point2; 4]) -> Scalar {
    if in_convex(a, quad) || in_convex(b, quad) {
        return 0.0;
    }
    let mut best = Scalar::INFINITY;
    for k in 0..4 {
        let (p, q) = (quad[k], quad[(k + 1) % 4]);
        if segments_meet(a, b, p, q) {
            return 0.0;
        }
        best = best
            .min(point_segment(a, p, q))
            .min(point_segment(b, p, q))
            .min(point_segment(p, a, b));
    }
    best
}

/// Inside a convex polygon (either winding), boundary included; a polygon
/// with no area contains nothing (its edges are tested separately).
fn in_convex(p: Point2, quad: &[Point2; 4]) -> bool {
    let mut positive = false;
    let mut negative = false;
    for k in 0..4 {
        let s = (quad[(k + 1) % 4] - quad[k]).perp_dot(p - quad[k]);
        positive |= s > 0.0;
        negative |= s < 0.0;
    }
    !(positive && negative) && (positive || negative)
}

fn point_segment(p: Point2, a: Point2, b: Point2) -> Scalar {
    let ab = b - a;
    let len = ab.length_squared();
    let t = if len > 0.0 {
        ((p - a).dot(ab) / len).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p - (a + ab * t)).length()
}

/// Whether closed segments `pq` and `rs` properly cross (touches are left
/// to the distance test).
fn segments_meet(p: Point2, q: Point2, r: Point2, s: Point2) -> bool {
    let orient = |a: Point2, b: Point2, c: Point2| (b - a).perp_dot(c - a);
    let (d1, d2) = (orient(r, s, p), orient(r, s, q));
    let (d3, d4) = (orient(p, q, r), orient(p, q, s));
    ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: Scalar, y: Scalar) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    /// Unit squares at the given lower-left corners, two triangles each,
    /// on a shared vertex lattice so neighbouring squares share edges.
    fn squares(at: &[(i32, i32)]) -> (Vec<[Point3; 3]>, Vec<[u32; 3]>) {
        let id = |x: i32, y: i32| (x + 10) as u32 * 100 + (y + 10) as u32;
        let mut triangles = Vec::new();
        let mut ids = Vec::new();
        for &(x, y) in at {
            let (fx, fy) = (Scalar::from(x), Scalar::from(y));
            triangles.push([p(fx, fy), p(fx + 1.0, fy), p(fx + 1.0, fy + 1.0)]);
            ids.push([id(x, y), id(x + 1, y), id(x + 1, y + 1)]);
            triangles.push([p(fx, fy), p(fx + 1.0, fy + 1.0), p(fx, fy + 1.0)]);
            ids.push([id(x, y), id(x + 1, y + 1), id(x, y + 1)]);
        }
        (triangles, ids)
    }

    fn piece(lo: (Scalar, Scalar), hi: (Scalar, Scalar), z: Scalar) -> [Point3; 4] {
        [
            Point3::new(lo.0, lo.1, z),
            Point3::new(hi.0, lo.1, z),
            Point3::new(hi.0, hi.1, z),
            Point3::new(lo.0, hi.1, z),
        ]
    }

    #[test]
    fn a_piece_across_a_diagonal_is_covered_by_its_square() {
        let (triangles, ids) = squares(&[(0, 0)]);
        let flat = Flat::new(&triangles, &ids);
        let bound = flat.bound(
            &piece((0.3, 0.3), (0.7, 0.7), 0.01),
            &[0, 1],
            Scalar::INFINITY,
        );
        assert!((bound - 0.01).abs() <= 1e-12, "{bound}");
    }

    /// An L of three squares: a piece over its notch is partly off the
    /// region, so the region must not cover it, while a piece across the
    /// inner corner's interior edges is covered.
    #[test]
    fn a_piece_over_a_notch_is_not_covered() {
        let (triangles, ids) = squares(&[(0, 0), (1, 0), (0, 1)]);
        let flat = Flat::new(&triangles, &ids);
        let all: Vec<usize> = (0..triangles.len()).collect();
        let over = flat.bound(&piece((0.8, 0.8), (1.2, 1.2), 0.0), &all, Scalar::INFINITY);
        assert_eq!(over, Scalar::INFINITY);
        let across = flat.bound(&piece((0.8, 0.2), (1.2, 0.9), 0.0), &all, Scalar::INFINITY);
        assert!(across <= 1e-12, "{across}");
    }

    /// A neighbour whose winding disagrees across a shared edge is no
    /// manifold neighbour: it does not join the region.
    #[test]
    fn a_neighbour_wound_against_its_edge_does_not_join() {
        let (mut triangles, mut ids) = squares(&[(0, 0)]);
        triangles[1].swap(1, 2);
        ids[1].swap(1, 2);
        let flat = Flat::new(&triangles, &ids);
        let bound = flat.bound(
            &piece((0.3, 0.3), (0.7, 0.7), 0.0),
            &[0, 1],
            Scalar::INFINITY,
        );
        assert_eq!(bound, Scalar::INFINITY);
    }

    /// A neighbour that folds back over the shared edge (consistent
    /// indices, opposite normal) lies on the same side of it: the union has
    /// that edge as boundary, so the two must not form a region that covers
    /// a piece across it.
    #[test]
    fn a_folded_neighbour_does_not_join() {
        let triangles = vec![
            [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0)],
            [p(0.0, 0.0), p(1.0, 1.0), p(0.9, 0.2)],
        ];
        let ids = vec![[0, 1, 2], [0, 2, 3]];
        let flat = Flat::new(&triangles, &ids);
        let bound = flat.bound(
            &piece((0.4, 0.35), (0.6, 0.55), 0.0),
            &[0, 1],
            Scalar::INFINITY,
        );
        assert_eq!(bound, Scalar::INFINITY);
    }
}
