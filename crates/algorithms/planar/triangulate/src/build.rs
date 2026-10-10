// SPDX-License-Identifier: MPL-2.0

//! Incremental constrained Delaunay construction.
//!
//! Points are inserted one at a time into a super triangle: located by a
//! walk, split one into three (or, on an edge, two into four), then
//! legalised by Lawson flips. Each split and flip re-links only the twin
//! pointers of the triangles it rewrites and of their outside neighbours,
//! so a step costs O(1) beyond the walk; recomputing the whole twin array
//! after each one, as 0.3.2 did, made n points cost O(n^2) or worse
//! (#268). Debug builds compare the result with a full recomputation
//! after every step (see the `adjacency` module). The walk starts next
//! to a grid hint (see the `hint` module), so it crosses a few triangles
//! rather than O(sqrt n).

use axiolid_core::Point2;
use axiolid_guarantees::Sign;
use axiolid_predicates::orient2d;

use crate::adjacency::{debug_check, link, note_update};
use crate::hint::CellHints;
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
    /// A triangle at each vertex, kept current by every split and flip,
    /// so a walk can start next to a hinted vertex.
    corner: Vec<u32>,
    hints: CellHints,
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
        let mut corner = vec![NO_HALFEDGE; all.len()];
        corner[points.len()..].fill(0);
        Self {
            points: all,
            triangles: vec![n, n + 1, n + 2],
            halfedges: vec![NO_HALFEDGE, NO_HALFEDGE, NO_HALFEDGE],
            real_count: points.len(),
            corner,
            hints: CellHints::new(min, max, points.len()),
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
        self.hints.set(self.points[v as usize], v);
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
            if self.halfedges[3 * t + i] == NO_HALFEDGE {
                return Err(TriangulationError::VertexUnplaceable { index: v });
            }
            let [t, ot, t1, t2] = self.split_edge(t, i, v);
            for tri in [t, ot, t1, t2] {
                self.legalize(3 * tri);
            }
            return Ok(());
        }

        let (t1, t2) = self.split_triangle(t, (a, b, c), v);
        self.legalize(3 * t);
        self.legalize(3 * t1);
        self.legalize(3 * t2);
        Ok(())
    }

    /// Split triangle `t` = `(a, b, c)` into three around `v`, keeping the
    /// twin pointers current, and return the two appended triangles.
    ///
    /// Slot `t` is reused for `(a, b, v)`, and `(b, c, v)` and `(c, a, v)`
    /// are appended, so existing neighbour indices into `t` stay valid
    /// where they can. Where they cannot is the trap: the old neighbours across
    /// `b -> c` and `c -> a` still point at halfedges `3t + 1` and `3t + 2`,
    /// which now run `b -> v` and `v -> a`. Those two neighbours are
    /// re-pointed at the appended triangles here; missing it corrupts the
    /// mesh in a way that only shows several insertions later, as a walk
    /// that leaves it.
    fn split_triangle(&mut self, t: usize, (a, b, c): (u32, u32, u32), v: u32) -> (usize, usize) {
        let across_bc = self.halfedges[3 * t + 1];
        let across_ca = self.halfedges[3 * t + 2];
        let t1 = self.triangles.len() / 3;
        let t2 = t1 + 1;
        self.triangles[3 * t] = a;
        self.triangles[3 * t + 1] = b;
        self.triangles[3 * t + 2] = v;
        self.triangles.extend_from_slice(&[b, c, v]);
        self.triangles.extend_from_slice(&[c, a, v]);
        self.halfedges.extend_from_slice(&[NO_HALFEDGE; 6]);

        // `a -> b` keeps its twin in slot `3t`.
        let h = &mut self.halfedges;
        link(h, 3 * t1, across_bc);
        link(h, 3 * t2, across_ca);
        link(h, 3 * t + 1, (3 * t1 + 2) as u32); // b -> v | v -> b
        link(h, 3 * t1 + 1, (3 * t2 + 2) as u32); // c -> v | v -> c
        link(h, 3 * t2 + 1, (3 * t + 2) as u32); // a -> v | v -> a
        self.corner[c as usize] = t1 as u32;
        self.corner[v as usize] = t as u32;
        note_update();
        debug_check(&self.triangles, &self.halfedges);
        (t1, t2)
    }

    /// Split edge `i` of triangle `t`, and the triangles either side of it,
    /// at `v`, keeping the twin pointers current. Returns the four
    /// triangles: the two reused slots, then the two appended.
    ///
    /// `(u, w, x)` and `(w, u, y)` become `(x, u, v)`, `(y, w, v)`,
    /// `(w, x, v)` and `(u, y, v)`. Every outside edge moves to a new slot,
    /// so all four outside neighbours are re-pointed.
    fn split_edge(&mut self, t: usize, i: usize, v: u32) -> [usize; 4] {
        let twin = self.halfedges[3 * t + i] as usize;
        let u = self.triangles[3 * t + i];
        let w = self.triangles[3 * t + (i + 1) % 3];
        let x = self.triangles[3 * t + (i + 2) % 3];
        let ot = twin / 3;
        let oi = twin % 3;
        let y = self.triangles[3 * ot + (oi + 2) % 3];
        let across_wx = self.halfedges[3 * t + (i + 1) % 3];
        let across_xu = self.halfedges[3 * t + (i + 2) % 3];
        let across_uy = self.halfedges[3 * ot + (oi + 1) % 3];
        let across_yw = self.halfedges[3 * ot + (oi + 2) % 3];

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

        let h = &mut self.halfedges;
        link(h, 3 * t, across_xu);
        link(h, 3 * ot, across_yw);
        link(h, 3 * t1, across_wx);
        link(h, 3 * t2, across_uy);
        link(h, 3 * t + 1, (3 * t2 + 2) as u32); // u -> v | v -> u
        link(h, 3 * t + 2, (3 * t1 + 1) as u32); // v -> x | x -> v
        link(h, 3 * ot + 1, (3 * t1 + 2) as u32); // w -> v | v -> w
        link(h, 3 * ot + 2, (3 * t2 + 1) as u32); // v -> y | y -> v
        for (vertex, tri) in [(u, t), (w, ot), (v, t)] {
            self.corner[vertex as usize] = tri as u32;
        }
        note_update();
        debug_check(&self.triangles, &self.halfedges);
        [t, ot, t1, t2]
    }

    /// Find a triangle containing `p`.
    ///
    /// The walk starts next to the last vertex inserted in `p`'s grid cell
    /// (see the `hint` module). A point strictly inside the triangle found
    /// has no other; a point on an edge or vertex, or a hinted walk that
    /// fails, is located again by the walk from the last triangle, the one
    /// the build always made, so the choice between neighbouring triangles,
    /// and with it the output, is unchanged.
    fn locate(&self, p: Point2) -> Option<usize> {
        let last = self.triangles.len() / 3 - 1;
        let start = self.hints.get(p).map_or(last, |v| {
            let t = self.corner[v as usize] as usize;
            // A stale entry would still be a valid start, only a slower
            // one; holding it exact keeps the bookkeeping honest.
            debug_assert!(
                self.triangles[3 * t..3 * t + 3].contains(&v),
                "vertex {v} is not a corner of its recorded triangle {t}"
            );
            t
        });
        if start != last {
            if let Some(t) = self.walk(start, p) {
                if self.strictly_inside(t, p) {
                    return Some(t);
                }
            }
        }
        self.walk(last, p)
    }

    /// Whether `p` lies strictly inside triangle `t`, off its edges.
    fn strictly_inside(&self, t: usize, p: Point2) -> bool {
        (0..3).all(|i| {
            let u = self.points[self.triangles[3 * t + i] as usize];
            let w = self.points[self.triangles[3 * t + (i + 1) % 3] as usize];
            decided(orient2d(u, w, p)) == Sign::Positive
        })
    }

    /// Walk from triangle `start` to one whose closure holds `p`.
    ///
    /// A straight walk needs no auxiliary structure. Returns `None` only if
    /// the walk leaves the mesh, which the super triangle makes impossible
    /// for in-range points.
    fn walk(&self, start: usize, p: Point2) -> Option<usize> {
        let mut t = start;
        for _ in 0..self.triangles.len() {
            #[cfg(test)]
            crate::hint::tests::STEPS.with(|c| c.set(c.get() + 1));
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
            // and (p1, apex, other). Slot `ti` of `t` now runs p0 -> other,
            // whose outside twin was across `ot`'s p0 -> other; slot `oi` of
            // `ot` runs p1 -> apex, whose twin was across `t`'s p1 -> apex.
            // The new diagonal pairs the two `+ 1` slots, and the remaining
            // edges keep both their slot and their twin.
            let across_p0_other = self.halfedges[3 * ot + (oi + 1) % 3];
            let across_p1_apex = self.halfedges[3 * t + (ti + 1) % 3];
            self.triangles[3 * t + (ti + 1) % 3] = other;
            self.triangles[3 * ot + (oi + 1) % 3] = apex;
            let h = &mut self.halfedges;
            link(h, 3 * t + ti, across_p0_other);
            link(h, 3 * ot + oi, across_p1_apex);
            link(h, 3 * t + (ti + 1) % 3, (3 * ot + (oi + 1) % 3) as u32);
            // p1 left `t` and p0 left `ot`.
            self.corner[p0 as usize] = t as u32;
            self.corner[p1 as usize] = ot as u32;
            note_update();
            debug_check(&self.triangles, &self.halfedges);

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
