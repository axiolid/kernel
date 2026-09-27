// SPDX-License-Identifier: MPL-2.0

//! Constraint recovery: force required edges into a Delaunay triangulation.
//!
//! After incremental insertion the triangulation is Delaunay but need not
//! contain the caller's constraint edges -- the Delaunay criterion and a
//! required edge can disagree. Recovery repairs that by repeatedly flipping
//! the edges that cross a missing constraint, the standard Anglada procedure.
//!
//! Flipping to recover an edge necessarily breaks the empty-circumcircle
//! property locally. That is the defining trade of a *constrained* Delaunay
//! triangulation and is why the crate documents the Delaunay guarantee as
//! holding "away from the constraints".

use axiolid_core::Point2;

use crate::mesh::{Triangulation, TriangulationError, NO_HALFEDGE};
use crate::{in_circumcircle, turns_left, Constraint};

impl Triangulation {
    /// Rebuild the halfedge twin array from the triangle list.
    ///
    /// Recovery rewrites triangles directly, so adjacency is recomputed
    /// rather than patched: an incrementally maintained twin array is easy to
    /// corrupt in a way that only shows up several flips later.
    pub(crate) fn rebuild_halfedges(&mut self) {
        use std::collections::HashMap;
        let count = self.triangles.len();
        self.halfedges = vec![NO_HALFEDGE; count];
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

    /// Whether some triangle already carries the edge `(a, b)`.
    fn has_edge(&self, a: u32, b: u32) -> bool {
        self.triangles
            .chunks_exact(3)
            .any(|t| (t[0] == a || t[1] == a || t[2] == a) && (t[0] == b || t[1] == b || t[2] == b))
    }
}

/// Force every constraint edge into the triangulation.
pub(crate) fn recover_constraints(tri: &mut Triangulation) -> Result<(), TriangulationError> {
    let constraints = tri.constraints.clone();
    for c in constraints {
        if tri.has_edge(c.a, c.b) {
            continue;
        }
        recover_one(tri, c)?;
    }
    tri.rebuild_halfedges();
    restore_delaunay(tri);
    Ok(())
}

/// Make every unconstrained edge locally Delaunay again.
///
/// Recovery's flips leave long, thin triangles whose circumcircles hold
/// points no constraint hides; without this pass the result was not
/// constrained Delaunay, though the crate promises it (#139: a skeleton
/// read from it ran across a corridor). Lawson's flips: an unconstrained
/// edge whose opposite vertex lies strictly inside the circumcircle is
/// flipped, until none is. Strictly inside, so cocircular ties never flip
/// and the loop ends; the bound only guards against a corrupted mesh.
fn restore_delaunay(tri: &mut Triangulation) {
    let mut budget = tri.triangles.len() * tri.triangles.len() + 64;
    loop {
        let mut changed = false;
        for t in 0..tri.triangle_count() {
            for i in 0..3 {
                let twin = tri.halfedges[3 * t + i];
                if twin == NO_HALFEDGE {
                    continue;
                }
                let (u, v) = (tri.triangles[3 * t + i], tri.triangles[3 * t + (i + 1) % 3]);
                if tri.is_constrained(u, v) {
                    continue;
                }
                let apex = tri.triangles[3 * t + (i + 2) % 3];
                let other = tri.triangles[3 * (twin as usize / 3) + (twin as usize % 3 + 2) % 3];
                let p = |k: u32| tri.points[k as usize];
                if in_circumcircle(p(u), p(v), p(apex), p(other)) && flip(tri, t, i) {
                    changed = true;
                    budget = budget.saturating_sub(1);
                }
            }
        }
        if !changed || budget == 0 {
            break;
        }
    }
}

/// Recover a single missing constraint by flipping the edges it crosses.
///
/// Anglada's queue: the crossing edges are taken in turn; one whose
/// quadrilateral is not convex yet goes to the back, and so does a new
/// diagonal that still crosses. Flipping whichever crossing edge came
/// first in a fresh scan, as this did, could flip a diagonal straight back
/// and cycle until the budget ran out (#190). There is always a flippable
/// edge while one crosses, so a pass over the whole queue with no flip
/// means another constraint blocks the way.
fn recover_one(tri: &mut Triangulation, c: Constraint) -> Result<(), TriangulationError> {
    use std::collections::VecDeque;
    let pa = tri.points[c.a as usize];
    let pb = tri.points[c.b as usize];
    let blocked = || TriangulationError::CrossingConstraints { a: c.a, b: c.b };
    let mut queue: VecDeque<(u32, u32)> = find_crossing(tri, pa, pb, c)
        .into_iter()
        .map(|(t, i)| (tri.triangles[3 * t + i], tri.triangles[3 * t + (i + 1) % 3]))
        .collect();
    if queue.is_empty() && !tri.has_edge(c.a, c.b) {
        return Err(blocked());
    }
    let mut idle = 0usize;
    let mut budget = tri.triangles.len() * tri.triangles.len() + 64;
    while let Some((u, v)) = queue.pop_front() {
        budget = budget.checked_sub(1).ok_or_else(blocked)?;
        let Some((t, i)) = halfedge(tri, u, v) else {
            continue;
        };
        if !flip(tri, t, i) {
            queue.push_back((u, v));
            idle += 1;
            if idle > queue.len() {
                return Err(blocked());
            }
            continue;
        }
        idle = 0;
        // The new diagonal runs apex -> other: slot 0 -> 1 of triangle `t`.
        let (x, y) = (tri.triangles[3 * t], tri.triangles[3 * t + 1]);
        let shares = x == c.a || x == c.b || y == c.a || y == c.b;
        if !shares
            && segments_properly_cross(pa, pb, tri.points[x as usize], tri.points[y as usize])
        {
            queue.push_back((x, y));
        }
    }
    if tri.has_edge(c.a, c.b) {
        Ok(())
    } else {
        Err(blocked())
    }
}

/// The triangle and slot of the edge between `u` and `v`, either way.
fn halfedge(tri: &Triangulation, u: u32, v: u32) -> Option<(usize, usize)> {
    for t in 0..tri.triangle_count() {
        for i in 0..3 {
            let (a, b) = (tri.triangles[3 * t + i], tri.triangles[3 * t + (i + 1) % 3]);
            if (a == u && b == v) || (a == v && b == u) {
                return Some((t, i));
            }
        }
    }
    None
}

/// Every triangle edge that properly crosses the segment `pa`-`pb`.
fn find_crossing(
    tri: &Triangulation,
    pa: Point2,
    pb: Point2,
    c: Constraint,
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for t in 0..tri.triangle_count() {
        for i in 0..3 {
            let u = tri.triangles[3 * t + i];
            let v = tri.triangles[3 * t + (i + 1) % 3];
            // An edge sharing an endpoint with the constraint cannot cross it.
            if u == c.a || u == c.b || v == c.a || v == c.b {
                continue;
            }
            // Never flip another constraint to satisfy this one.
            if tri.is_constrained(u, v) {
                continue;
            }
            let pu = tri.points[u as usize];
            let pv = tri.points[v as usize];
            if segments_properly_cross(pa, pb, pu, pv) {
                out.push((t, i));
            }
        }
    }
    out
}

/// Strict crossing test: shared endpoints and touching do not count.
fn segments_properly_cross(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let d1 = turns_left(a, b, c);
    let d2 = turns_left(a, b, d);
    let d3 = turns_left(c, d, a);
    let d4 = turns_left(c, d, b);
    d1 != d2 && d3 != d4
}

/// Flip the edge `i` of triangle `t` to the opposite diagonal.
///
/// Returns `false` when the edge is on the boundary or the quadrilateral is
/// not convex, in which case flipping would produce an inverted triangle.
fn flip(tri: &mut Triangulation, t: usize, i: usize) -> bool {
    let twin = tri.halfedges[3 * t + i];
    if twin == NO_HALFEDGE {
        return false;
    }
    let ot = twin as usize / 3;
    let oi = twin as usize % 3;

    let a = tri.triangles[3 * t + i];
    let b = tri.triangles[3 * t + (i + 1) % 3];
    let apex = tri.triangles[3 * t + (i + 2) % 3];
    let other = tri.triangles[3 * ot + (oi + 2) % 3];

    // Only a strictly convex quadrilateral can be flipped without inverting.
    let (pa, pb) = (tri.points[a as usize], tri.points[b as usize]);
    let (pap, pot) = (tri.points[apex as usize], tri.points[other as usize]);
    if !(turns_left(pap, pa, pot) && turns_left(pot, pb, pap)) {
        return false;
    }

    tri.triangles[3 * t] = apex;
    tri.triangles[3 * t + 1] = other;
    tri.triangles[3 * t + 2] = b;
    tri.triangles[3 * ot] = other;
    tri.triangles[3 * ot + 1] = apex;
    tri.triangles[3 * ot + 2] = a;
    tri.rebuild_halfedges();
    true
}
