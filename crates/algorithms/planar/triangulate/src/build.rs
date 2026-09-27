// SPDX-License-Identifier: MPL-2.0

//! Incremental constrained Delaunay construction.

use axiolid_core::Point2;
use axiolid_guarantees::Sign;
use axiolid_predicates::orient2d;

use crate::mesh::{Triangulation, TriangulationError, NO_HALFEDGE};
use crate::recover::recover_constraints;
use crate::{decided, in_circumcircle, turns_left, Constraint};

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

/// Build a constrained Delaunay triangulation of `points`.
///
/// Every edge in `constraints` appears in the output as an edge of some
/// triangle. Away from constraints the result satisfies the Delaunay
/// empty-circumcircle property, decided exactly.
///
/// # Errors
///
/// Returns [`TriangulationError`] when the input has fewer than three points,
/// is entirely collinear, references a vertex that does not exist, or
/// contains constraints that cross one another.
pub fn triangulate(
    points: &[Point2],
    constraints: &[Constraint],
) -> Result<Triangulation, TriangulationError> {
    if points.len() < 3 {
        return Err(TriangulationError::TooFewPoints);
    }
    let count = u32::try_from(points.len()).unwrap_or(u32::MAX);
    for c in constraints {
        if c.a >= count || c.b >= count {
            let index = if c.a >= count { c.a } else { c.b };
            return Err(TriangulationError::ConstraintOutOfRange { index });
        }
    }
    if !has_non_collinear_triple(points) {
        return Err(TriangulationError::AllPointsCollinear);
    }

    let mut state = Builder::new(points);
    state.insert_all()?;
    let mut triangulation = state.finish(points, constraints);
    recover_constraints(&mut triangulation)?;
    Ok(triangulation)
}

/// Whether the input contains three points that are not collinear.
///
/// Checked with the exact predicate rather than by an area threshold: a
/// tolerance here would reject a legitimately thin but non-degenerate input.
fn has_non_collinear_triple(points: &[Point2]) -> bool {
    let a = points[0];
    // Find the first point distinct from `a`, then the first that is not on
    // the line through the two. Scanning is linear and runs once.
    let Some(b) = points.iter().copied().find(|p| *p != a) else {
        return false;
    };
    points
        .iter()
        .any(|&c| turns_left(a, b, c) || turns_left(b, a, c))
}

/// Incremental Delaunay construction over a super triangle.
struct Builder {
    points: Vec<Point2>,
    triangles: Vec<u32>,
    halfedges: Vec<u32>,
    real_count: usize,
}

impl Builder {
    fn new(points: &[Point2]) -> Self {
        let mut all = points.to_vec();
        // A super triangle that strictly contains every input point. Scaled
        // generously off the bounding box: a tight super triangle puts input
        // points on its edges, and a point exactly on an edge has no
        // containing triangle to split.
        let (min, max) = bounds(points);
        let dx = max.x - min.x;
        let dy = max.y - min.y;
        let span = if dx > dy { dx } else { dy };
        let span = if span > 0.0 { span } else { 1.0 };
        let cx = (min.x + max.x) * 0.5;
        let cy = (min.y + max.y) * 0.5;
        let far = span * 32.0;
        all.push(Point2::new(cx - far, cy - far));
        all.push(Point2::new(cx + far, cy - far));
        all.push(Point2::new(cx, cy + far));

        let n = points.len() as u32;
        Self {
            points: all,
            triangles: vec![n, n + 1, n + 2],
            halfedges: vec![NO_HALFEDGE, NO_HALFEDGE, NO_HALFEDGE],
            real_count: points.len(),
        }
    }

    fn insert_all(&mut self) -> Result<(), TriangulationError> {
        for v in 0..self.real_count as u32 {
            self.insert_point(v)?;
        }
        Ok(())
    }

    /// Insert one vertex, splitting its containing triangle and restoring the
    /// Delaunay property by flipping.
    fn insert_point(&mut self, v: u32) -> Result<(), TriangulationError> {
        // A vertex that cannot be located is a corrupted adjacency, not a
        // benign skip. Dropping it silently produced a triangulation missing
        // an input point -- the caller would get a valid-looking mesh with a
        // hole where their vertex should be.
        let Some(t) = self.locate(self.points[v as usize]) else {
            return Err(TriangulationError::VertexUnplaceable { index: v });
        };
        let (a, b, c) = (
            self.triangles[3 * t],
            self.triangles[3 * t + 1],
            self.triangles[3 * t + 2],
        );
        let p = self.points[v as usize];

        // On an edge: split the edge, and the triangles either side of it,
        // rather than the one triangle into three -- one of the three would
        // have no area, and flips over a triangle with no area never settle.
        for i in 0..3 {
            let (u, w) = (
                self.triangles[3 * t + i],
                self.triangles[3 * t + (i + 1) % 3],
            );
            if decided(orient2d(
                self.points[u as usize],
                self.points[w as usize],
                p,
            )) != Sign::Zero
            {
                continue;
            }
            let twin = self.halfedges[3 * t + i];
            if twin == NO_HALFEDGE {
                return Err(TriangulationError::VertexUnplaceable { index: v });
            }
            let x = self.triangles[3 * t + (i + 2) % 3];
            let ot = twin as usize / 3;
            let oi = twin as usize % 3;
            let y = self.triangles[3 * ot + (oi + 2) % 3];
            // (u, w, x) and (w, u, y) become (u, v, x), (v, w, x),
            // (w, v, y), (v, u, y).
            let t1 = self.triangles.len() / 3;
            let t2 = t1 + 1;
            self.triangles[3 * t] = x;
            self.triangles[3 * t + 1] = u;
            self.triangles[3 * t + 2] = v;
            self.triangles[3 * ot] = y;
            self.triangles[3 * ot + 1] = w;
            self.triangles[3 * ot + 2] = v;
            self.triangles.extend_from_slice(&[w, x, v]);
            self.triangles.extend_from_slice(&[u, y, v]);
            self.halfedges.extend_from_slice(&[NO_HALFEDGE; 6]);
            self.rebuild_adjacency();
            for tri in [t, ot, t1, t2] {
                self.legalize(3 * tri);
            }
            return Ok(());
        }

        // Reuse slot `t` for the first sub-triangle and append the other two,
        // so existing neighbour indices into `t` stay valid where possible.
        let t1 = self.triangles.len() / 3;
        let t2 = t1 + 1;
        self.triangles[3 * t] = a;
        self.triangles[3 * t + 1] = b;
        self.triangles[3 * t + 2] = v;
        self.triangles.extend_from_slice(&[b, c, v]);
        self.triangles.extend_from_slice(&[c, a, v]);
        self.halfedges.extend_from_slice(&[NO_HALFEDGE; 6]);

        // Rebuild adjacency from the triangle list rather than patching the
        // six affected twin pointers by hand. Hand-patching was wrong here:
        // reusing slot `t` for one sub-triangle leaves the OLD neighbours of
        // `t` pointing at halfedges that now belong to a different triangle,
        // and the corruption only surfaces later as a walk that leaves the
        // mesh. The rebuild is O(n) per insertion, which is the price of a
        // provably consistent structure; `locate` stays correct, and the
        // refinement loop above is what dominates runtime anyway.
        self.rebuild_adjacency();

        self.legalize(3 * t);
        self.legalize(3 * t1);
        self.legalize(3 * t2);
        Ok(())
    }

    /// Recompute every twin pointer from the triangle list.
    fn rebuild_adjacency(&mut self) {
        use std::collections::HashMap;
        let count = self.triangles.len();
        self.halfedges.clear();
        self.halfedges.resize(count, NO_HALFEDGE);
        let mut seen: HashMap<(u32, u32), u32> = HashMap::with_capacity(count);
        for e in 0..count {
            let t = e / 3;
            let i = e % 3;
            let from = self.triangles[3 * t + i];
            let to = self.triangles[3 * t + (i + 1) % 3];
            if let Some(&twin) = seen.get(&(to, from)) {
                self.halfedges[e] = twin;
                self.halfedges[twin as usize] = e as u32;
            } else {
                seen.insert((from, to), e as u32);
            }
        }
    }

    /// Find a triangle containing `p` by walking from the last triangle.
    ///
    /// A straight walk is O(sqrt n) on well-distributed input and needs no
    /// auxiliary structure. Returns `None` only if the walk leaves the mesh,
    /// which the super triangle makes impossible for in-range points.
    fn locate(&self, p: Point2) -> Option<usize> {
        let mut t = self.triangles.len() / 3 - 1;
        for _ in 0..self.triangles.len() {
            let (a, b, c) = (
                self.points[self.triangles[3 * t] as usize],
                self.points[self.triangles[3 * t + 1] as usize],
                self.points[self.triangles[3 * t + 2] as usize],
            );
            // Step across the first edge that `p` lies strictly outside of.
            // For a counter-clockwise triangle, "outside edge (u, w)" means
            // p is strictly to the RIGHT of u->w. Testing `turns_left(w, u, p)`
            // instead also fires when p is exactly ON the edge, which sends
            // the walk across a boundary it should have stopped at.
            let mut moved = false;
            for (i, (u, w)) in [(a, b), (b, c), (c, a)].into_iter().enumerate() {
                let right_of = decided(orient2d(u, w, p)) == Sign::Negative;
                if right_of {
                    let twin = self.halfedges[3 * t + i];
                    if twin == NO_HALFEDGE {
                        return None;
                    }
                    t = twin as usize / 3;
                    moved = true;
                    break;
                }
            }
            if !moved {
                return Some(t);
            }
        }
        None
    }

    /// Restore the Delaunay property across `edge`, whose triangle's third
    /// corner is the vertex just inserted, and propagate to the edges
    /// across from it.
    fn legalize(&mut self, edge: usize) {
        let mut stack = vec![edge];
        // Lawson's flips after one insertion are few; the bound only keeps
        // a corrupted mesh from spinning. Anything left is repaired after
        // constraint recovery.
        let mut budget = 4 * self.triangles.len() + 64;
        while let Some(e) = stack.pop() {
            if budget == 0 {
                break;
            }
            budget -= 1;
            let twin = self.halfedges[e];
            if twin == NO_HALFEDGE {
                continue;
            }
            let t = e / 3;
            let ti = e % 3;
            let o = twin as usize;
            let ot = o / 3;
            let oi = o % 3;

            let p0 = self.triangles[3 * t + ti];
            let p1 = self.triangles[3 * t + (ti + 1) % 3];
            let apex = self.triangles[3 * t + (ti + 2) % 3];
            let other = self.triangles[3 * ot + (oi + 2) % 3];

            if !in_circumcircle(
                self.points[p0 as usize],
                self.points[p1 as usize],
                self.points[apex as usize],
                self.points[other as usize],
            ) {
                continue;
            }
            // Only a strictly convex quadrilateral flips without inverting.
            let q = |k: u32| self.points[k as usize];
            if !(turns_left(q(apex), q(p0), q(other)) && turns_left(q(other), q(p1), q(apex))) {
                continue;
            }

            // Flip the shared edge to the other diagonal: (p0, other, apex)
            // and (p1, apex, other).
            self.triangles[3 * t + (ti + 1) % 3] = other;
            self.triangles[3 * ot + (oi + 1) % 3] = apex;

            // Same reasoning as insertion: recompute adjacency rather than
            // re-point twins by hand.
            self.rebuild_adjacency();

            // The two edges now across from the inserted vertex: p0 -> other
            // and other -> p1. (Pushing the new diagonal instead, as this
            // did, checked an edge that is Delaunay by construction and left
            // these unchecked; thin quads then flipped back and forth.)
            stack.push(3 * t + ti);
            stack.push(3 * ot + (oi + 2) % 3);
        }
    }

    /// Drop the super triangle and everything touching it.
    fn finish(self, original: &[Point2], constraints: &[Constraint]) -> Triangulation {
        let limit = self.real_count as u32;
        let mut triangles = Vec::with_capacity(self.triangles.len());
        for t in self.triangles.chunks_exact(3) {
            if t.iter().all(|&v| v < limit) {
                triangles.extend_from_slice(t);
            }
        }
        let mut sorted = constraints.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        let mut out = Triangulation {
            points: original.to_vec(),
            triangles,
            halfedges: Vec::new(),
            constraints: sorted,
        };
        out.rebuild_halfedges();
        out
    }
}

/// Axis-aligned bounds of a point set.
fn bounds(points: &[Point2]) -> (Point2, Point2) {
    let mut min = points[0];
    let mut max = points[0];
    for p in &points[1..] {
        min.x = min.x.min(p.x);
        min.y = min.y.min(p.y);
        max.x = max.x.max(p.x);
        max.y = max.y.max(p.y);
    }
    (min, max)
}
