//! The four passes of one remeshing iteration, each a sweep of guarded
//! local edits over a [`HalfedgeMesh`].
//!
//! Every edit is checked geometrically *before* it is applied: a collapse,
//! flip or vertex move that would turn a triangle's normal by a right angle
//! or more, or leave it flatter than [`MIN_SHAPE`], is skipped and the mesh
//! is left as it was; so is a collapse or flip that lowers the smallest
//! angle of the triangles it rewrites below [`QUALITY_ANGLE`]. Connectivity
//! refusals ([`HalfedgeMesh::collapse_edge`]'s
//! link condition, [`HalfedgeMesh::flip_edge`] onto an existing edge) are
//! skips too. Skipping is not degrading: the mesh stays a valid, consistently
//! oriented approximation of the input, and the report says how far its
//! edge lengths are from the target.
//!
//! Sweeps visit elements in id order and every choice is a pure function of
//! the current mesh, so the output is a function of the input alone.

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::{EdgeId, FaceId, HalfedgeMesh, VertexId};

use super::reference::Reference;
use super::RemeshError;

/// Smallest height-to-longest-edge ratio an edit may leave on a triangle it
/// creates or moves.
pub(super) const MIN_SHAPE: Scalar = 1e-3;

/// Two protected edges meeting at a vertex count as one straight feature
/// line when the sine of their angle is at most this.
const COLLINEAR: Scalar = 1e-9;

/// Smallest angle, in radians, a collapse or flip may leave when it lowers
/// the smallest angle of the triangles it rewrites. Without it valence and
/// length alone would cut ears off a fine feature polygon.
const QUALITY_ANGLE: Scalar = core::f64::consts::PI / 12.0;

/// Optimal valence of an interior vertex of a regular triangulation.
const INTERIOR_VALENCE: i64 = 6;
/// Optimal valence of a boundary vertex.
const BOUNDARY_VALENCE: i64 = 4;

/// Edit counts accumulated over every iteration.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Counts {
    pub(super) splits: usize,
    pub(super) collapses: usize,
    pub(super) flips: usize,
    pub(super) moves: usize,
}

pub(super) struct Remesher<'a> {
    pub(super) mesh: HalfedgeMesh,
    /// Per edge id: boundary or sharp feature, never flipped, collapsed only
    /// along a straight feature line.
    protected: Vec<bool>,
    /// Per face id: the input patch the face lies on.
    patch: Vec<u32>,
    reference: &'a Reference,
    low: Scalar,
    high: Scalar,
    max_triangles: usize,
    pub(super) counts: Counts,
}

impl<'a> Remesher<'a> {
    pub(super) fn new(
        mesh: HalfedgeMesh,
        protected: Vec<bool>,
        patch: Vec<u32>,
        reference: &'a Reference,
        target: Scalar,
        max_triangles: usize,
    ) -> Self {
        Self {
            mesh,
            protected,
            patch,
            reference,
            low: target * 4.0 / 5.0,
            high: target * 4.0 / 3.0,
            max_triangles,
            counts: Counts::default(),
        }
    }

    pub(super) fn is_protected(&self, e: EdgeId) -> bool {
        self.protected.get(e.index()).copied().unwrap_or(false)
    }

    fn set_protected(&mut self, e: EdgeId, value: bool) {
        if self.protected.len() <= e.index() {
            self.protected.resize(e.index() + 1, false);
        }
        self.protected[e.index()] = value;
    }

    fn face_patch(&self, f: FaceId) -> u32 {
        self.patch[f.index()]
    }

    fn set_patch(&mut self, f: FaceId, value: u32) {
        if self.patch.len() <= f.index() {
            self.patch.resize(f.index() + 1, 0);
        }
        self.patch[f.index()] = value;
    }

    fn length(&self, e: EdgeId) -> Scalar {
        let [a, b] = self.mesh.edge_vertices(e);
        (self.mesh.position(a) - self.mesh.position(b)).length()
    }

    /// Protected edges at a vertex.
    fn protected_degree(&self, v: VertexId) -> usize {
        self.mesh
            .outgoing_halfedges(v)
            .filter(|&h| self.is_protected(self.mesh.edge(h)))
            .count()
    }

    /// A vertex no protected edge touches: it may slide over its patch.
    fn is_free(&self, v: VertexId) -> bool {
        !self.mesh.is_isolated(v) && self.protected_degree(v) == 0
    }

    /// Patch of a free vertex: every face around it lies on the same one.
    fn vertex_patch(&self, v: VertexId) -> Option<u32> {
        self.mesh.vertex_faces(v).next().map(|f| self.face_patch(f))
    }

    // ---- split --------------------------------------------------------------

    /// Split every edge longer than `4/3` of the target at its midpoint,
    /// until none is. A free edge's midpoint is projected onto its patch; a
    /// protected edge's midpoint stays on the straight segment it divides.
    pub(super) fn split_long_edges(&mut self) -> Result<(), RemeshError> {
        loop {
            let long: Vec<EdgeId> = self
                .mesh
                .edges()
                .filter(|&e| self.length(e) > self.high)
                .collect();
            if long.is_empty() {
                return Ok(());
            }
            for e in long {
                self.split(e)?;
            }
        }
    }

    fn split(&mut self, e: EdgeId) -> Result<(), RemeshError> {
        if self.mesh.face_count() + 2 > self.max_triangles {
            return Err(RemeshError::BudgetExceeded {
                produced: self.mesh.face_count() + 2,
                limit: self.max_triangles,
            });
        }
        let h = self.mesh.edge_halfedge(e, 0);
        let o = self.mesh.opposite(h);
        let [a, b] = self.mesh.edge_vertices(e);
        let midpoint = self.mesh.position(a).midpoint(self.mesh.position(b));
        let protected = self.is_protected(e);
        let side_h = self.mesh.face(h).map(|f| self.face_patch(f));
        let side_o = self.mesh.face(o).map(|f| self.face_patch(f));

        let m = self
            .mesh
            .split_edge(e, midpoint)
            .map_err(RemeshError::Internal)?;
        self.counts.splits += 1;
        // `h` now runs a -> m and a new edge m -> b continues it.
        let g = self.mesh.find_halfedge(m, b).ok_or(RemeshError::Internal(
            axiolid_mesh::HalfedgeEditError::RemovedElement {
                element: "edge",
                index: e.get(),
            },
        ))?;
        self.set_protected(self.mesh.edge(g), protected);
        for (halfedge, side) in [
            (h, side_h),
            (g, side_h),
            (self.mesh.opposite(g), side_o),
            (o, side_o),
        ] {
            if let (Some(f), Some(p)) = (self.mesh.face(halfedge), side) {
                self.set_patch(f, p);
            }
        }
        // Diagonals the split added are new, unprotected edges.
        let spokes: Vec<EdgeId> = self
            .mesh
            .outgoing_halfedges(m)
            .map(|s| self.mesh.edge(s))
            .filter(|&s| s != e && s != self.mesh.edge(g))
            .collect();
        for s in spokes {
            self.set_protected(s, false);
        }

        if !protected {
            if let Some(p) = side_h.or(side_o) {
                let projected = self.reference.closest_on_patch(p, midpoint);
                if self.move_keeps_shape(m, projected) {
                    self.mesh.set_position(m, projected);
                }
            }
        }
        Ok(())
    }

    // ---- collapse -----------------------------------------------------------

    /// Collapse edges shorter than `4/5` of the target, where the collapse
    /// creates no edge longer than `4/3` of it and folds no triangle.
    pub(super) fn collapse_short_edges(&mut self) {
        loop {
            let short: Vec<EdgeId> = self
                .mesh
                .edges()
                .filter(|&e| self.length(e) < self.low)
                .collect();
            let mut changed = false;
            for e in short {
                if self.mesh.contains_edge(e) && self.length(e) < self.low && self.try_collapse(e) {
                    changed = true;
                }
            }
            if !changed {
                return;
            }
        }
    }

    fn try_collapse(&mut self, e: EdgeId) -> bool {
        let [a, b] = self.mesh.edge_vertices(e);
        let (pa, pb) = (self.mesh.position(a), self.mesh.position(b));
        let mut candidates: Vec<(VertexId, VertexId, Point3)> = Vec::with_capacity(2);
        if self.is_protected(e) {
            // Along a feature line only, and only removing a vertex the line
            // runs straight through, so the line keeps its exact shape.
            if self.is_straight_feature_vertex(a, e) {
                candidates.push((a, b, pb));
            }
            if self.is_straight_feature_vertex(b, e) {
                candidates.push((b, a, pa));
            }
        } else {
            match (self.is_free(a), self.is_free(b)) {
                (true, true) => {
                    let mid = pa.midpoint(pb);
                    let placed = self
                        .vertex_patch(a)
                        .map_or(mid, |p| self.reference.closest_on_patch(p, mid));
                    candidates.push((a, b, placed));
                    if placed != mid {
                        candidates.push((a, b, mid));
                    }
                }
                (true, false) => candidates.push((a, b, pb)),
                (false, true) => candidates.push((b, a, pa)),
                // Joining two feature vertices across a patch would pull the
                // patch onto its border.
                (false, false) => {}
            }
        }
        for (remove, keep, position) in candidates {
            if !self.collapse_is_acceptable(remove, keep, position) {
                continue;
            }
            let Some(h) = self.mesh.find_halfedge(remove, keep) else {
                continue;
            };
            let sides = self.collapse_sides(h);
            if self.mesh.collapse_edge(h).is_err() {
                continue;
            }
            self.mesh.set_position(keep, position);
            // Each removed triangle's two other sides merge into one edge,
            // which stays protected if either side was.
            for [x, y] in sides {
                let merged = self.is_protected(x) || self.is_protected(y);
                for side in [x, y] {
                    if self.mesh.contains_edge(side) {
                        self.set_protected(side, merged);
                    }
                }
            }
            self.counts.collapses += 1;
            return true;
        }
        false
    }

    /// The pairs of sides that merge when `h` collapses.
    fn collapse_sides(&self, h: axiolid_mesh::HalfedgeId) -> Vec<[EdgeId; 2]> {
        let mut sides = Vec::with_capacity(2);
        for side in [h, self.mesh.opposite(h)] {
            if self.mesh.face(side).is_some() {
                let next = self.mesh.next(side);
                let prev = self.mesh.prev(side);
                sides.push([self.mesh.edge(next), self.mesh.edge(prev)]);
            }
        }
        sides
    }

    /// Whether `v` lies on exactly two protected edges, one of them `e`,
    /// that continue each other in a straight line.
    fn is_straight_feature_vertex(&self, v: VertexId, e: EdgeId) -> bool {
        let protected: Vec<VertexId> = self
            .mesh
            .outgoing_halfedges(v)
            .filter(|&h| self.is_protected(self.mesh.edge(h)))
            .map(|h| self.mesh.target(h))
            .collect();
        let [x, y] = protected[..] else {
            return false;
        };
        let [ea, eb] = self.mesh.edge_vertices(e);
        if !(x == ea || x == eb || y == ea || y == eb) {
            return false;
        }
        let p = self.mesh.position(v);
        let incoming = p - self.mesh.position(x);
        let outgoing = self.mesh.position(y) - p;
        let scale = incoming.length() * outgoing.length();
        incoming.dot(outgoing) > 0.0 && incoming.cross(outgoing).length() <= COLLINEAR * scale
    }

    /// Whether merging `remove` into `keep` at `position` keeps every edge
    /// at most `4/3` of the target and every surviving triangle unfolded.
    fn collapse_is_acceptable(&self, remove: VertexId, keep: VertexId, position: Point3) -> bool {
        for v in [remove, keep] {
            for n in self.mesh.vertex_vertices(v) {
                if n != remove
                    && n != keep
                    && (self.mesh.position(n) - position).length() > self.high
                {
                    return false;
                }
            }
        }
        let moved = |u: VertexId| {
            if u == remove || u == keep {
                position
            } else {
                self.mesh.position(u)
            }
        };
        let mut worst_before = Scalar::INFINITY;
        let mut worst_after = Scalar::INFINITY;
        for v in [remove, keep] {
            for f in self.mesh.vertex_faces(v) {
                let corners: Vec<VertexId> = self.mesh.face_vertices(f).collect();
                let before = triangle(corners.iter().map(|&u| self.mesh.position(u)));
                worst_before = worst_before.min(min_angle(before));
                if corners.contains(&remove) && corners.contains(&keep) {
                    continue;
                }
                let after = triangle(corners.iter().map(|&u| moved(u)));
                if !replaces_without_folding(before, after) {
                    return false;
                }
                worst_after = worst_after.min(min_angle(after));
            }
        }
        keeps_quality(worst_before, worst_after)
    }

    // ---- flip ---------------------------------------------------------------

    /// Flip unprotected edges whose flip strictly lowers the total distance
    /// of the four corners' valences from optimal (6 inside, 4 on a
    /// boundary). That sum falls with every flip, so the sweep ends.
    pub(super) fn equalize_valences(&mut self) {
        loop {
            let edges: Vec<EdgeId> = self.mesh.edges().collect();
            let mut changed = false;
            for e in edges {
                if self.try_flip(e) {
                    changed = true;
                }
            }
            if !changed {
                return;
            }
        }
    }

    fn try_flip(&mut self, e: EdgeId) -> bool {
        if self.is_protected(e) || self.mesh.is_boundary_edge(e) {
            return false;
        }
        let h = self.mesh.edge_halfedge(e, 0);
        let o = self.mesh.opposite(h);
        let a = self.mesh.source(h);
        let b = self.mesh.target(h);
        let c = self.mesh.target(self.mesh.next(h));
        let d = self.mesh.target(self.mesh.next(o));
        let deviation = |v: VertexId, change: i64| {
            let optimal = if self.mesh.is_boundary_vertex(v) {
                BOUNDARY_VALENCE
            } else {
                INTERIOR_VALENCE
            };
            (self.mesh.degree(v) as i64 + change - optimal).abs()
        };
        let before = deviation(a, 0) + deviation(b, 0) + deviation(c, 0) + deviation(d, 0);
        let after = deviation(a, -1) + deviation(b, -1) + deviation(c, 1) + deviation(d, 1);
        if after >= before {
            return false;
        }
        if !self.flip_keeps_shape([a, b, c, d]) {
            return false;
        }
        if self.mesh.flip_edge(e).is_err() {
            return false;
        }
        self.counts.flips += 1;
        true
    }

    /// Whether flipping the edge `a -> b` of faces `(a, b, c)` and
    /// `(b, a, d)` to `(a, d, c)` and `(d, b, c)` folds or flattens neither
    /// new triangle and keeps quality.
    fn flip_keeps_shape(&self, corners: [VertexId; 4]) -> bool {
        let [pa, pb, pc, pd] = corners.map(|v| self.mesh.position(v));
        let old = [[pa, pb, pc], [pb, pa, pd]];
        let new = [[pa, pd, pc], [pd, pb, pc]];
        for created in new {
            for replaced in old {
                if !replaces_without_folding(replaced, created) {
                    return false;
                }
            }
        }
        let smallest = |pair: [[Point3; 3]; 2]| min_angle(pair[0]).min(min_angle(pair[1]));
        keeps_quality(smallest(old), smallest(new))
    }

    // ---- relax and project --------------------------------------------------

    /// Move each free vertex towards the area-weighted centroid of its
    /// triangles within its tangent plane, then onto its input patch. Every
    /// target is computed from the positions before the sweep; a move that
    /// would fold a triangle falls back to projecting the vertex where it
    /// stands, and failing that the vertex stays.
    pub(super) fn relax_and_project(&mut self) {
        let targets: Vec<(VertexId, Point3)> = self
            .mesh
            .vertices()
            .filter(|&v| self.is_free(v))
            .filter_map(|v| self.tangential_target(v).map(|t| (v, t)))
            .collect();
        for (v, target) in targets {
            let Some(patch) = self.vertex_patch(v) else {
                continue;
            };
            let current = self.mesh.position(v);
            let relaxed = self.reference.closest_on_patch(patch, target);
            let projected = self.reference.closest_on_patch(patch, current);
            for candidate in [relaxed, projected] {
                if candidate == current {
                    break;
                }
                if self.move_keeps_shape(v, candidate) {
                    self.mesh.set_position(v, candidate);
                    self.counts.moves += 1;
                    break;
                }
            }
        }
    }

    fn tangential_target(&self, v: VertexId) -> Option<Point3> {
        let p = self.mesh.position(v);
        let mut weighted = Vec3::ZERO;
        let mut total = 0.0;
        let mut normal = Vec3::ZERO;
        for f in self.mesh.vertex_faces(v) {
            let t = triangle(self.mesh.face_vertices(f).map(|u| self.mesh.position(u)));
            let n = (t[1] - t[0]).cross(t[2] - t[0]);
            let area = n.length();
            weighted += (t[0] + t[1] + t[2]) / 3.0 * area;
            total += area;
            normal += n;
        }
        let normal = normal.try_normalize()?;
        if total <= 0.0 {
            return None;
        }
        let shift = weighted / total - p;
        Some(p + shift - normal * shift.dot(normal))
    }

    /// Whether moving `v` to `position` leaves every triangle around it
    /// unfolded and not flat.
    fn move_keeps_shape(&self, v: VertexId, position: Point3) -> bool {
        self.mesh.vertex_faces(v).all(|f| {
            let corners: Vec<VertexId> = self.mesh.face_vertices(f).collect();
            let before = triangle(corners.iter().map(|&u| self.mesh.position(u)));
            let after = triangle(corners.iter().map(|&u| {
                if u == v {
                    position
                } else {
                    self.mesh.position(u)
                }
            }));
            replaces_without_folding(before, after)
        })
    }
}

fn triangle(mut corners: impl Iterator<Item = Point3>) -> [Point3; 3] {
    let mut next = || corners.next().unwrap_or(Point3::ZERO);
    [next(), next(), next()]
}

/// Whether an edit whose rewritten triangles had smallest angle `before`
/// and have `after` keeps quality: it does not lower the smallest angle, or
/// lowers it no further than [`QUALITY_ANGLE`].
fn keeps_quality(before: Scalar, after: Scalar) -> bool {
    after >= before.min(QUALITY_ANGLE)
}

/// Smallest interior angle of a triangle, in radians.
pub(super) fn min_angle([a, b, c]: [Point3; 3]) -> Scalar {
    let angle = |p: Point3, q: Point3, r: Point3| (q - p).angle_between(r - p);
    angle(a, b, c).min(angle(b, c, a)).min(angle(c, a, b))
}

/// Unnormalised normal: twice the area along the unit normal.
pub(super) fn normal([a, b, c]: [Point3; 3]) -> Vec3 {
    (b - a).cross(c - a)
}

/// Whether a triangle is no flatter than [`MIN_SHAPE`]: its height over its
/// longest edge, `|n| / longest^2`.
pub(super) fn is_well_shaped(t: [Point3; 3]) -> bool {
    let longest = [(t[1] - t[0]), (t[2] - t[1]), (t[0] - t[2])]
        .iter()
        .map(|e| e.length_squared())
        .fold(0.0, Scalar::max);
    longest > 0.0 && normal(t).length() > MIN_SHAPE * longest
}

/// Whether `created` may stand where `replaced` stood: well shaped and with
/// a normal less than a right angle from the old one.
fn replaces_without_folding(replaced: [Point3; 3], created: [Point3; 3]) -> bool {
    is_well_shaped(created) && normal(created).dot(normal(replaced)) > 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: Scalar, y: Scalar) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    /// A remesher over a planar mesh with nothing protected, one patch.
    fn with_remesher(
        positions: Vec<Point3>,
        faces: &[[u32; 3]],
        check: impl FnOnce(&mut Remesher<'_>),
    ) {
        let mesh = HalfedgeMesh::from_faces(positions.clone(), faces).expect("valid");
        let triangles = faces
            .iter()
            .map(|t| t.map(|i| positions[i as usize]))
            .collect();
        let reference = Reference::new(triangles, &vec![0; faces.len()]);
        let edges = mesh.edge_count();
        let mut remesher = Remesher::new(
            mesh,
            vec![false; edges],
            vec![0; faces.len()],
            &reference,
            1.0,
            1000,
        );
        check(&mut remesher);
    }

    #[test]
    fn an_edge_between_two_feature_lines_never_collapses() {
        // Rows y = -0.5, 0, 0.1, 0.6 of three columns; the rows y = 0 and
        // y = 0.1 are feature lines. The rung between them at x = 0.5 is
        // short, and collapsing it would move a whole stretch of one line
        // onto the other without folding anything.
        let rows = [-0.5, 0.0, 0.1, 0.6];
        let positions: Vec<Point3> = rows
            .iter()
            .flat_map(|&y| [0.0, 0.5, 1.0].map(|x| p(x, y)))
            .collect();
        let mut faces = Vec::new();
        for r in 0..3u32 {
            for c in 0..2u32 {
                let (a, b) = (r * 3 + c, r * 3 + c + 1);
                let (d, e) = (a + 3, b + 3);
                faces.extend([[a, b, e], [a, e, d]]);
            }
        }
        with_remesher(positions, &faces, |r| {
            for (from, to) in [(3, 4), (4, 5), (6, 7), (7, 8)] {
                let h = r
                    .mesh
                    .find_halfedge(VertexId::new(from), VertexId::new(to))
                    .expect("feature edge");
                r.set_protected(r.mesh.edge(h), true);
            }
            let rung = r
                .mesh
                .find_halfedge(VertexId::new(4), VertexId::new(7))
                .expect("rung");
            assert!(!r.try_collapse(r.mesh.edge(rung)));
            assert_eq!(r.mesh.vertex_count(), 12);
        });
    }

    #[test]
    fn a_flip_of_a_non_convex_quad_is_refused() {
        // Faces (a, b, c) and (b, a, d); the diagonal c-d misses a-b, so the
        // new triangle (a, d, c) would be inverted (yet well shaped).
        let positions = vec![p(0.0, 0.0), p(1.0, 0.0), p(-1.0, 1.0), p(-1.0, -0.5)];
        with_remesher(positions, &[[0, 1, 2], [1, 0, 3]], |r| {
            assert!(!r.flip_keeps_shape([0, 1, 2, 3].map(VertexId::new)));
        });
        // A convex square flips.
        let positions = vec![p(0.0, 0.0), p(1.0, 1.0), p(0.0, 1.0), p(1.0, 0.0)];
        with_remesher(positions, &[[0, 1, 2], [1, 0, 3]], |r| {
            assert!(r.flip_keeps_shape([0, 1, 2, 3].map(VertexId::new)));
        });
    }

    #[test]
    fn a_move_that_folds_or_flattens_a_triangle_is_refused() {
        // Centre vertex 0 of a square fan.
        let positions = vec![
            p(0.0, 0.0),
            p(-1.0, -1.0),
            p(1.0, -1.0),
            p(1.0, 1.0),
            p(-1.0, 1.0),
        ];
        let faces = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]];
        with_remesher(positions, &faces, |r| {
            let centre = VertexId::new(0);
            assert!(r.move_keeps_shape(centre, p(0.3, 0.2)));
            // Out of the fan's kernel: triangle (0, 1, 2) inverts.
            assert!(!r.move_keeps_shape(centre, p(0.0, -1.5)));
            // Onto (almost) the side 1-2: still upward, but flat.
            assert!(!r.move_keeps_shape(centre, p(0.0, -1.0 + 1e-6)));
        });
    }
}
