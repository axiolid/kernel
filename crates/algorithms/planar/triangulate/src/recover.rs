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

use crate::adjacency::{debug_check, link, note_update, twins};
use crate::mesh::{Triangulation, TriangulationError, NO_HALFEDGE};
use crate::{collinear, in_circumcircle, turns_left, Constraint};

impl Triangulation {
    /// Rebuild the halfedge twin array from the triangle list.
    ///
    /// Run once, when the build hands over a triangulation with the super
    /// triangle removed. Recovery's flips then keep it current locally.
    pub(crate) fn rebuild_halfedges(&mut self) {
        self.halfedges = twins(&self.triangles);
    }

    /// Whether some triangle already carries the edge `(a, b)`.
    fn has_edge(&self, a: u32, b: u32) -> bool {
        self.triangles
            .as_chunks::<3>()
            .0
            .iter()
            .any(|t| (t[0] == a || t[1] == a || t[2] == a) && (t[0] == b || t[1] == b || t[2] == b))
    }
}

/// The constraints that were split, each with the vertices it runs through.
///
/// A constraint is a union of output edges, so one that passes exactly
/// through input vertices is recovered as the chain of pieces between them
/// (#360). Only constraints that were split are listed; any other is the
/// single edge between its endpoints.
#[derive(Debug, Clone, Default)]
pub(crate) struct Splits(Vec<(Constraint, Vec<u32>)>);

impl Splits {
    /// The vertices from `c.a` to `c.b` that `c` was split into, endpoints
    /// included, or `None` when `c` was not split.
    pub(crate) fn chain(&self, c: Constraint) -> Option<&[u32]> {
        self.0
            .binary_search_by_key(&c, |(key, _)| *key)
            .ok()
            .map(|i| self.0[i].1.as_slice())
    }
}

/// Force every constraint edge into the triangulation, splitting a
/// constraint at the kept vertices on it first.
///
/// `tri.constraints` becomes the sorted, deduplicated pieces before any is
/// recovered, so no recovery flips a piece of a later constraint.
pub(crate) fn recover_constraints(tri: &mut Triangulation) -> Result<Splits, TriangulationError> {
    let mut is_vertex = vec![true; tri.points.len()];
    for m in &tri.merged {
        is_vertex[m.input as usize] = false;
    }
    let mut splits = Vec::new();
    let mut pieces = Vec::new();
    for &c in &tri.constraints {
        if c.a == c.b || tri.has_edge(c.a, c.b) {
            pieces.push(c);
            continue;
        }
        let mut chain = vec![c.a];
        chain.extend(vertices_on(tri, &is_vertex, c.a, c.b));
        chain.push(c.b);
        pieces.extend(chain.windows(2).map(|w| Constraint::new(w[0], w[1])));
        if chain.len() > 2 {
            splits.push((c, chain));
        }
    }
    pieces.sort_unstable();
    pieces.dedup();
    tri.constraints = pieces;
    for c in tri.constraints.clone() {
        if tri.has_edge(c.a, c.b) {
            continue;
        }
        recover_one(tri, c)?;
    }
    restore_delaunay(tri);
    Ok(Splits(splits))
}

/// The kept vertices strictly inside the segment `a`-`b`, ordered from `a`
/// to `b`.
///
/// A merged duplicate is not a vertex (`is_vertex` is false): its kept twin
/// is the one on the segment. The test is the closed bounding box, then the
/// exact collinearity predicate, so a vertex exactly on the segment is
/// found and one off it by any amount is not.
///
/// TODO(#355): one pass over every vertex per missing constraint.
pub(crate) fn vertices_on(tri: &Triangulation, is_vertex: &[bool], a: u32, b: u32) -> Vec<u32> {
    let (pa, pb) = (tri.points[a as usize], tri.points[b as usize]);
    let (lo_x, hi_x) = (pa.x.min(pb.x), pa.x.max(pb.x));
    let (lo_y, hi_y) = (pa.y.min(pb.y), pa.y.max(pb.y));
    let mut hits: Vec<u32> = (0..tri.points.len() as u32)
        .filter(|&v| {
            let q = tri.points[v as usize];
            is_vertex[v as usize]
                && v != a
                && v != b
                && lo_x <= q.x
                && q.x <= hi_x
                && lo_y <= q.y
                && q.y <= hi_y
                && collinear(pa, pb, q)
        })
        .collect();
    // On one line the lexicographic order is the order along it.
    let key = |v: &u32| {
        let q = tri.points[*v as usize];
        (q.x, q.y)
    };
    hits.sort_by(|u, v| {
        key(u)
            .partial_cmp(&key(v))
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    if (pa.x, pa.y) > (pb.x, pb.y) {
        hits.reverse();
    }
    hits
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

    // Both triangles are rewritten from slot 0, so every outside edge
    // moves to a new slot: read all four outside twins before writing.
    let across_b_apex = tri.halfedges[3 * t + (i + 1) % 3];
    let across_apex_a = tri.halfedges[3 * t + (i + 2) % 3];
    let across_a_other = tri.halfedges[3 * ot + (oi + 1) % 3];
    let across_other_b = tri.halfedges[3 * ot + (oi + 2) % 3];

    tri.triangles[3 * t] = apex;
    tri.triangles[3 * t + 1] = other;
    tri.triangles[3 * t + 2] = b;
    tri.triangles[3 * ot] = other;
    tri.triangles[3 * ot + 1] = apex;
    tri.triangles[3 * ot + 2] = a;

    let h = &mut tri.halfedges;
    link(h, 3 * t, (3 * ot) as u32); // apex -> other | other -> apex
    link(h, 3 * t + 1, across_other_b);
    link(h, 3 * t + 2, across_b_apex);
    link(h, 3 * ot + 1, across_apex_a);
    link(h, 3 * ot + 2, across_a_other);
    note_update();
    debug_check(&tri.triangles, &tri.halfedges);
    true
}
