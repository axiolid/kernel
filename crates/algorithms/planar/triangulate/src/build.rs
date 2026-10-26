// SPDX-License-Identifier: MPL-2.0

//! Incremental constrained Delaunay construction.
//!
//! Points are inserted one at a time: located by a walk, split one into
//! three (or, on an edge, two into four), then legalised by Lawson flips.
//! Each split and flip re-links only the twin pointers of the triangles it
//! rewrites and of their outside neighbours,
//! so a step costs O(1) beyond the walk; recomputing the whole twin array
//! after each one, as 0.3.2 did, made n points cost O(n^2) or worse
//! (#268). Debug builds compare the result with a full recomputation
//! after every step (see the `adjacency` module). The walk starts next
//! to a grid hint (see the `hint` module), so it crosses a few triangles
//! rather than O(sqrt n).
//!
//! # The vertex at infinity
//!
//! The build starts from one triangle of input points and closes the
//! convex hull with a single symbolic vertex at infinity, the standard
//! construction of Shewchuk's Triangle and CGAL's `Triangulation_2`:
//! every hull edge carries a *ghost* triangle with the infinite vertex on
//! its outer side, so the mesh is a closed surface and every halfedge has
//! a twin. The infinite vertex has no coordinates; the predicates decide
//! it symbolically:
//!
//! - a point lies in a ghost triangle when it lies strictly on the outer
//!   side of the ghost's hull edge: the ghost's "circumcircle" is that
//!   open half-plane;
//! - the infinite vertex lies in no finite triangle's circumcircle, so a
//!   hull edge never flips;
//! - an edge to the infinite vertex flips when the new point sees the hull
//!   edge beyond it strictly, which wraps the hull round a point inserted
//!   outside it; on that edge's line it does not, so collinear hull
//!   points stay hull vertices.
//!
//! Removing the ghosts leaves exactly the Delaunay triangulation of the
//! input, its convex hull covered. The finite super triangle this replaced
//! (32 times the input's span) lost thin triangles along the hull whose
//! circumcircles reached its corners, without a word (#299): 3,977
//! instead of 3,979 triangles for 2,000 random points.
//!
//! When the Delaunay triangulation is unique (no four input points
//! cocircular), the output is the same set of triangles as before for
//! every input the super triangle got right; their order in the list
//! differs, because the build starts from a different first triangle.
//! Where four or more points are cocircular the Delaunay triangulation is
//! not unique, and the two builds may choose differently.
//!
//! # Coincident points
//!
//! A point at exactly the position of one already inserted is merged into
//! it (#300): the walk ends in a triangle with that vertex as a corner,
//! and the point is recorded instead of split in. Splitting it in, as the
//! build did, left zero-area and overlapping triangles that used one
//! directed edge twice. The merge is reported through
//! [`Triangulation::merged_points`] (input index to kept vertex), so
//! nothing is lost. The kept vertex is the lowest input index at that
//! position: points are inserted in input order, and the first triangle
//! takes the first index at each of its positions. Constraints are carried
//! over to the kept vertices; one whose two endpoints coincide asks for an
//! edge from a vertex to itself and is refused as
//! [`TriangulationError::CoincidentConstraint`].

use axiolid_core::Point2;
use axiolid_guarantees::Sign;
use axiolid_predicates::orient2d;

use crate::adjacency::{debug_check, link, note_update};
use crate::hint::CellHints;
use crate::mesh::{MergedPoint, Triangulation, TriangulationError, NO_HALFEDGE};
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
/// is entirely collinear, references a vertex that does not exist, has a
/// constraint between two coincident points, or contains constraints that
/// cross one another.
///
/// Exactly coincident points are merged into the lowest input index at
/// their position and reported by [`Triangulation::merged_points`]; the
/// triangles and [`Triangulation::constraints`] use the kept indices.
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
    state.insert_all(first_triangle(points))?;
    let constraints = state.carry_over(constraints)?;
    let mut triangulation = state.finish(points, &constraints);
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

/// The first three points that span a triangle, counter-clockwise.
///
/// The first point, the first point distinct from it, and the first point
/// off the line through the two; the caller has checked that one exists.
/// For points in general position that is `(0, 1, 2)`, so the insertion
/// order is the input order, as it was with the super triangle.
fn first_triangle(points: &[Point2]) -> [u32; 3] {
    let a = 0;
    let b = points.iter().position(|p| *p != points[0]).unwrap_or(1);
    let (c, sign) = points
        .iter()
        .enumerate()
        .find_map(|(k, &q)| match decided(orient2d(points[a], points[b], q)) {
            Sign::Zero => None,
            sign => Some((k, sign)),
        })
        .unwrap_or((2, Sign::Positive));
    if sign == Sign::Positive {
        [a as u32, b as u32, c as u32]
    } else {
        [a as u32, c as u32, b as u32]
    }
}

/// Incremental Delaunay construction closed by a vertex at infinity.
struct Builder {
    /// The input points. The infinite vertex has no entry, so reading its
    /// position by mistake panics instead of computing with a stand-in.
    points: Vec<Point2>,
    triangles: Vec<u32>,
    halfedges: Vec<u32>,
    /// The infinite vertex: one past the last input point.
    infinite: u32,
    /// The vertex each input point is: itself, or the coincident point it
    /// was merged into.
    kept: Vec<u32>,
    /// A triangle at each vertex, kept current by every split and flip,
    /// so a walk can start next to a hinted vertex.
    corner: Vec<u32>,
    hints: CellHints,
}

impl Builder {
    fn new(points: &[Point2]) -> Self {
        let (min, max) = bounds(points);
        let n = points.len();
        Self {
            points: points.to_vec(),
            triangles: Vec::new(),
            halfedges: Vec::new(),
            infinite: n as u32,
            kept: (0..n as u32).collect(),
            corner: vec![NO_HALFEDGE; n + 1],
            hints: CellHints::new(min, max, n),
        }
    }

    /// Start from the triangle `first` and its three ghosts, then insert
    /// every other point in input order.
    fn insert_all(&mut self, first: [u32; 3]) -> Result<(), TriangulationError> {
        let [a, b, c] = first;
        let g = self.infinite;
        // Slot 0 is the finite triangle; slots 1, 2 and 3 are the ghosts
        // across its edges a -> b, b -> c and c -> a.
        self.triangles = vec![a, b, c, b, a, g, c, b, g, a, c, g];
        self.halfedges = vec![NO_HALFEDGE; 12];
        let h = &mut self.halfedges;
        link(h, 0, 3); // a -> b | b -> a
        link(h, 1, 6); // b -> c | c -> b
        link(h, 2, 9); // c -> a | a -> c
        link(h, 4, 11); // a -> g | g -> a
        link(h, 7, 5); // b -> g | g -> b
        link(h, 10, 8); // c -> g | g -> c
        for v in first {
            self.corner[v as usize] = 0;
            self.hints.set(self.points[v as usize], v);
        }
        self.corner[g as usize] = 1;
        note_update();
        debug_check(&self.triangles, &self.halfedges);
        for v in 0..self.infinite {
            if !first.contains(&v) {
                self.insert_point(v)?;
            }
        }
        Ok(())
    }

    /// The slot of the infinite vertex in triangle `t`, if it is a ghost.
    fn ghost_slot(&self, t: usize) -> Option<usize> {
        (0..3).find(|&i| self.triangles[3 * t + i] == self.infinite)
    }

    /// Insert one vertex, splitting its containing triangle and restoring the
    /// Delaunay property by flipping.
    fn insert_point(&mut self, v: u32) -> Result<(), TriangulationError> {
        // A vertex that cannot be located is a corrupted adjacency, not a
        // benign skip. Dropping it silently produced a triangulation missing
        // an input point -- the caller would get a valid-looking mesh with a
        // hole where their vertex should be.
        let p = self.points[v as usize];
        let Some(t) = self.locate(p) else {
            return Err(TriangulationError::VertexUnplaceable { index: v });
        };
        let (a, b, c) = (
            self.triangles[3 * t],
            self.triangles[3 * t + 1],
            self.triangles[3 * t + 2],
        );

        // At a vertex already inserted: merge, and report it (#300). A
        // point at a vertex lies in the closure of exactly the triangles
        // around it, so the triangle found has that vertex as a corner. A
        // ghost holds only points strictly outside the hull, and its
        // infinite corner has no position to compare.
        if self.ghost_slot(t).is_none() {
            if let Some(&at) = [a, b, c].iter().find(|&&u| self.points[u as usize] == p) {
                self.kept[v as usize] = at;
                return Ok(());
            }
        }
        self.hints.set(p, v);

        // Outside the hull: the ghost's hull edge is visible from `p`. The
        // split gives one finite triangle and two ghosts, and the flips
        // wrap the hull round `p`.
        if self.ghost_slot(t).is_some() {
            let (t1, t2) = self.split_triangle(t, (a, b, c), v);
            self.legalize(3 * t);
            self.legalize(3 * t1);
            self.legalize(3 * t2);
            return Ok(());
        }

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

    /// Whether `p` lies strictly inside finite triangle `t`, off its edges.
    ///
    /// Never for a ghost: a point outside the hull can see several hull
    /// edges, so the ghost a walk ends in need not be the only one.
    fn strictly_inside(&self, t: usize, p: Point2) -> bool {
        self.ghost_slot(t).is_none()
            && (0..3).all(|i| {
                let u = self.points[self.triangles[3 * t + i] as usize];
                let w = self.points[self.triangles[3 * t + (i + 1) % 3] as usize];
                decided(orient2d(u, w, p)) == Sign::Positive
            })
    }

    /// Walk from triangle `start` to a finite triangle whose closure holds
    /// `p`, or to a ghost whose hull edge `p` sees strictly.
    ///
    /// A straight walk needs no auxiliary structure. The ghosts close the
    /// mesh, so no halfedge lacks a twin; `None` means a corrupted mesh, or
    /// a walk that did not settle within one step per triangle.
    fn walk(&self, start: usize, p: Point2) -> Option<usize> {
        let mut t = start;
        for _ in 0..self.triangles.len() {
            #[cfg(test)]
            crate::hint::tests::STEPS.with(|c| c.set(c.get() + 1));
            if let Some(k) = self.ghost_slot(t) {
                match self.ghost_step(t, k, p) {
                    Some(next) => t = next,
                    None => return Some(t),
                }
                continue;
            }
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

    /// The next triangle of a walk standing in ghost `t`, whose infinite
    /// vertex is in slot `k`, or `None` when `p` lies in it.
    ///
    /// The ghost of hull edge `a -> b` holds the points strictly on the
    /// edge's outer side. Any other point steps across `a -> b` into the
    /// hull, and the walk goes on from there: one on the edge stops in
    /// the finite triangle beside it, and one on the edge's line beyond an
    /// end lies strictly outside another hull edge, which the walk then
    /// crosses into that edge's ghost. A walk only stands in a ghost it did
    /// not enter that way at its start (a hint or the last triangle), since
    /// it enters a ghost only across a hull edge `p` lies strictly outside.
    fn ghost_step(&self, t: usize, k: usize, p: Point2) -> Option<usize> {
        let edge = 3 * t + (k + 1) % 3;
        let a = self.points[self.triangles[edge] as usize];
        let b = self.points[self.triangles[3 * t + (k + 2) % 3] as usize];
        if decided(orient2d(a, b, p)) == Sign::Positive {
            None
        } else {
            Some(self.halfedges[edge] as usize / 3)
        }
    }

    /// Whether the edge `p0 -> p1` of the triangle `(p0, p1, apex)` must
    /// flip, with `other` the vertex across it.
    ///
    /// A finite quadrilateral flips when `other` lies strictly inside the
    /// circumcircle and the quadrilateral is strictly convex. With the
    /// infinite vertex: nothing lies beyond the hull, so a hull edge never
    /// flips (`other` infinite); an edge to the infinite vertex flips when
    /// `other` lies strictly on the outer side of the ghost's hull edge,
    /// which also makes the finite triangle the flip creates
    /// counter-clockwise. On that edge's line it does not flip, and the
    /// collinear vertex stays on the hull.
    fn must_flip(&self, p0: u32, p1: u32, apex: u32, other: u32) -> bool {
        let g = self.infinite;
        let q = |k: u32| self.points[k as usize];
        if other == g || apex == g {
            return false;
        }
        if p1 == g {
            // (p0, ∞, apex) is the ghost of hull edge apex -> p0.
            return turns_left(q(apex), q(p0), q(other));
        }
        if p0 == g {
            // (∞, p1, apex) is the ghost of hull edge p1 -> apex.
            return turns_left(q(p1), q(apex), q(other));
        }
        in_circumcircle(q(p0), q(p1), q(apex), q(other))
            // Only a strictly convex quadrilateral flips without inverting.
            && turns_left(q(apex), q(p0), q(other))
            && turns_left(q(other), q(p1), q(apex))
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

            if !self.must_flip(p0, p1, apex, other) {
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

    /// The constraints between kept vertices, or the first one whose two
    /// endpoints were merged into one vertex.
    fn carry_over(
        &self,
        constraints: &[Constraint],
    ) -> Result<Vec<Constraint>, TriangulationError> {
        constraints
            .iter()
            .map(|c| {
                let (a, b) = (self.kept[c.a as usize], self.kept[c.b as usize]);
                if a == b && c.a != c.b {
                    Err(TriangulationError::CoincidentConstraint { a: c.a, b: c.b })
                } else {
                    Ok(Constraint::new(a, b))
                }
            })
            .collect()
    }

    /// Drop the ghosts; the finite triangles cover the convex hull.
    fn finish(self, original: &[Point2], constraints: &[Constraint]) -> Triangulation {
        let g = self.infinite;
        let mut triangles = Vec::with_capacity(self.triangles.len());
        for t in self.triangles.as_chunks::<3>().0 {
            if !t.contains(&g) {
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
            merged: (0..self.kept.len() as u32)
                .filter(|&v| self.kept[v as usize] != v)
                .map(|v| MergedPoint {
                    input: v,
                    kept: self.kept[v as usize],
                })
                .collect(),
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
