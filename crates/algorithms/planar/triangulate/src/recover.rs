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
use crate::{in_circumcircle, turns_left, Constraint};

/// "No triangle" in the vertex-to-triangle index.
const NO_TRIANGLE: u32 = u32::MAX;

/// The index or the twin array contradicts the triangles: a fan that does
/// not end, a missing corner, a walk that leaves the mesh. Reported as
/// [`TriangulationError::VertexUnplaceable`], the crate's corrupted-adjacency
/// error; valid input never raises it.
#[derive(Debug, Clone, Copy)]
struct Corrupt;

impl Triangulation {
    /// Rebuild the halfedge twin array from the triangle list.
    ///
    /// Run once, when the build hands over a triangulation with the super
    /// triangle removed. Recovery's flips then keep it current locally.
    pub(crate) fn rebuild_halfedges(&mut self) {
        self.halfedges = twins(&self.triangles);
    }
}

/// One triangle that has each vertex as a corner, [`NO_TRIANGLE`] for a
/// vertex no triangle uses.
///
/// Recovery reaches every edge through this index instead of scanning the
/// triangle array. [`flip`] keeps it current: a flip rewrites two triangles
/// and nothing else, and all four of their corners are updated.
fn incident_triangles(tri: &Triangulation) -> Vec<u32> {
    let mut incident = vec![NO_TRIANGLE; tri.points.len()];
    for (t, corners) in tri.triangles.as_chunks::<3>().0.iter().enumerate() {
        for &v in corners {
            incident[v as usize] = t as u32;
        }
    }
    incident
}

/// The halfedge after `h` in its triangle.
fn next_halfedge(h: usize) -> usize {
    3 * (h / 3) + (h % 3 + 1) % 3
}

/// The halfedge before `h` in its triangle.
fn prev_halfedge(h: usize) -> usize {
    3 * (h / 3) + (h % 3 + 2) % 3
}

/// Every halfedge that leaves `u`, in no particular order: the fan of
/// triangles around it, walked through the twins from one triangle in one
/// direction until it closes or reaches the boundary, then from the same
/// triangle in the other. A fan has at most one halfedge per triangle, so
/// a longer one is [`Corrupt`].
fn outgoing(tri: &Triangulation, incident: &[u32], u: u32) -> Result<Vec<usize>, Corrupt> {
    let t = incident[u as usize];
    if t == NO_TRIANGLE {
        return Ok(Vec::new());
    }
    let Some(slot) = (0..3).find(|&s| tri.triangles[3 * t as usize + s] == u) else {
        return Err(Corrupt);
    };
    let bound = tri.triangles.len() / 3;
    let start = 3 * t as usize + slot;
    let mut fan = vec![start];
    let mut h = start;
    let mut closed = false;
    loop {
        // The twin of the edge into `u` is the next edge out of it.
        let twin = tri.halfedges[prev_halfedge(h)];
        if twin == NO_HALFEDGE {
            break;
        }
        if twin as usize == start {
            closed = true;
            break;
        }
        fan.push(twin as usize);
        h = twin as usize;
        if fan.len() > bound {
            return Err(Corrupt);
        }
    }
    if !closed {
        h = start;
        loop {
            let twin = tri.halfedges[h];
            if twin == NO_HALFEDGE {
                break;
            }
            h = next_halfedge(twin as usize);
            fan.push(h);
            if fan.len() > bound {
                return Err(Corrupt);
            }
        }
    }
    Ok(fan)
}

/// The halfedge of the edge between `u` and `v`, either way: the lower
/// index of the two when both sides exist, which is the first one a scan of
/// the triangle array meets.
fn find_edge(
    tri: &Triangulation,
    incident: &[u32],
    u: u32,
    v: u32,
) -> Result<Option<usize>, Corrupt> {
    let mut best: Option<usize> = None;
    let mut keep = |h: usize| best = Some(best.map_or(h, |b| b.min(h)));
    for h in outgoing(tri, incident, u)? {
        if tri.triangles[next_halfedge(h)] == v {
            keep(h);
        }
        // A boundary edge into `u` has no twin among the outgoing ones.
        let into = prev_halfedge(h);
        if tri.triangles[into] == v {
            keep(into);
        }
    }
    Ok(best)
}

/// Whether some triangle already carries the edge `(a, b)`.
fn has_edge(tri: &Triangulation, incident: &[u32], a: u32, b: u32) -> Result<bool, Corrupt> {
    Ok(find_edge(tri, incident, a, b)?.is_some())
}

/// What a walk from `a` towards `b` met.
#[derive(Debug, Default)]
struct Walk {
    /// The vertices exactly on the segment, strictly inside it, in order.
    through: Vec<u32>,
    /// One halfedge of every edge the segment properly crosses, in order.
    crossed: Vec<usize>,
}

/// The sign of the exact orientation of `c` against `a -> b`: positive to
/// the left, zero on the line.
fn side(tri: &Triangulation, a: u32, b: u32, c: u32) -> i8 {
    let p = |k: u32| tri.points[k as usize];
    if turns_left(p(a), p(b), p(c)) {
        1
    } else if turns_left(p(b), p(a), p(c)) {
        -1
    } else {
        0
    }
}

/// Walk the segment `a`-`b` through the triangulation.
///
/// From `a`, the triangle of the fan whose wedge holds `b` is found with
/// the exact orientation predicate. Each edge the segment crosses is
/// recorded and the walk goes on into the triangle behind it, until it
/// reaches `b`. A vertex exactly on the segment is recorded and the walk
/// restarts from it. Every decision is an exact sign, so the result is the
/// set of edges that properly cross, with no tolerance.
///
/// [`Corrupt`] when the walk leaves the triangulation or does not end, which
/// a valid triangulation of a hull holding both ends does not do.
fn walk(tri: &Triangulation, incident: &[u32], a: u32, b: u32) -> Result<Walk, Corrupt> {
    let mut out = Walk::default();
    let mut budget = tri.triangles.len() + tri.points.len() + 8;
    let mut cur = a;
    while cur != b {
        // The wedge at `cur` that holds `b`: counter-clockwise (cur, p, q)
        // with `b` left of or on cur -> p and right of or on cur -> q.
        let wedge = outgoing(tri, incident, cur)?.into_iter().find_map(|h| {
            let (p, q) = (
                tri.triangles[next_halfedge(h)],
                tri.triangles[prev_halfedge(h)],
            );
            let (sp, sq) = (side(tri, cur, p, b), side(tri, cur, q, b));
            (sp >= 0 && sq <= 0).then_some((h, p, q, sp, sq))
        });
        let (h, p, q, sp, sq) = wedge.ok_or(Corrupt)?;
        // `b` on a ray from `cur`: that neighbour is on the segment (it
        // cannot lie beyond `b`, which is a vertex and so no edge's inside).
        let ahead = if sp == 0 {
            Some(p)
        } else if sq == 0 {
            Some(q)
        } else {
            None
        };
        if let Some(v) = ahead {
            if v != b {
                out.through.push(v);
            }
            cur = v;
            budget = budget.checked_sub(1).ok_or(Corrupt)?;
            continue;
        }
        // `b` is beyond the edge p -> q: cross it, and the ones after it.
        let mut e = next_halfedge(h);
        loop {
            budget = budget.checked_sub(1).ok_or(Corrupt)?;
            out.crossed.push(e);
            let twin = tri.halfedges[e];
            if twin == NO_HALFEDGE {
                return Err(Corrupt);
            }
            let twin = twin as usize;
            let apex = tri.triangles[prev_halfedge(twin)];
            if apex == b {
                cur = b;
                break;
            }
            match side(tri, a, b, apex) {
                0 => {
                    out.through.push(apex);
                    cur = apex;
                    break;
                }
                1 => {
                    // `apex` is left of the line: the segment leaves by
                    // the edge from the right end to `apex`.
                    e = next_halfedge(twin);
                }
                // Right of the line: by the edge from `apex` to the left end.
                _ => e = prev_halfedge(twin),
            }
        }
    }
    Ok(out)
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
    let mut incident = incident_triangles(tri);
    let mut splits = Vec::new();
    let mut pieces = Vec::new();
    for &c in &tri.constraints {
        if c.a == c.b || has_edge(tri, &incident, c.a, c.b).map_err(|e| unplaceable(c, e))? {
            pieces.push(c);
            continue;
        }
        let mut chain = vec![c.a];
        chain.extend(vertices_on(tri, &incident, c.a, c.b)?);
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
        if c.a == c.b || has_edge(tri, &incident, c.a, c.b).map_err(|e| unplaceable(c, e))? {
            continue;
        }
        recover_one(tri, &mut incident, c)?;
    }
    restore_delaunay(tri, &mut incident);
    #[cfg(test)]
    for (v, &t) in incident.iter().enumerate() {
        assert!(
            t == NO_TRIANGLE || tri.triangles[3 * t as usize..][..3].contains(&(v as u32)),
            "the index names a triangle without vertex {v}"
        );
    }
    Ok(Splits(splits))
}

/// The kept vertices strictly inside the segment `a`-`b`, ordered from `a`
/// to `b`.
///
/// A merged duplicate is in no triangle, so its kept twin is the one on the
/// segment. Found by walking the segment, with the exact orientation
/// predicate: a vertex exactly on the segment is found and one off it by
/// any amount is not.
pub(crate) fn vertices_on(
    tri: &Triangulation,
    incident: &[u32],
    a: u32,
    b: u32,
) -> Result<Vec<u32>, TriangulationError> {
    walk(tri, incident, a, b)
        .map(|w| w.through)
        .map_err(|e| unplaceable(Constraint::new(a, b), e))
}

/// The refusal for a [`Corrupt`] mesh found while recovering `c`.
fn unplaceable(c: Constraint, _: Corrupt) -> TriangulationError {
    TriangulationError::VertexUnplaceable { index: c.a }
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
fn restore_delaunay(tri: &mut Triangulation, incident: &mut [u32]) {
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
                if in_circumcircle(p(u), p(v), p(apex), p(other)) && flip(tri, incident, t, i) {
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
fn recover_one(
    tri: &mut Triangulation,
    incident: &mut [u32],
    c: Constraint,
) -> Result<(), TriangulationError> {
    use std::collections::VecDeque;
    let pa = tri.points[c.a as usize];
    let pb = tri.points[c.b as usize];
    let blocked = || TriangulationError::CrossingConstraints { a: c.a, b: c.b };
    let crossing = find_crossing(tri, incident, c)?;
    #[cfg(test)]
    assert_eq!(crossing, scan_crossing(tri, pa, pb, c));
    let mut queue: VecDeque<(u32, u32)> = crossing
        .into_iter()
        .map(|(t, i)| (tri.triangles[3 * t + i], tri.triangles[3 * t + (i + 1) % 3]))
        .collect();
    if queue.is_empty() && !has_edge(tri, incident, c.a, c.b).map_err(|e| unplaceable(c, e))? {
        return Err(blocked());
    }
    let mut idle = 0usize;
    let mut budget = tri.triangles.len() * tri.triangles.len() + 64;
    while let Some((u, v)) = queue.pop_front() {
        budget = budget.checked_sub(1).ok_or_else(blocked)?;
        let Some((t, i)) = halfedge(tri, incident, u, v).map_err(|e| unplaceable(c, e))? else {
            continue;
        };
        if !flip(tri, incident, t, i) {
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
    if has_edge(tri, incident, c.a, c.b).map_err(|e| unplaceable(c, e))? {
        Ok(())
    } else {
        Err(blocked())
    }
}

/// The triangle and slot of the edge between `u` and `v`, either way.
fn halfedge(
    tri: &Triangulation,
    incident: &[u32],
    u: u32,
    v: u32,
) -> Result<Option<(usize, usize)>, Corrupt> {
    let found = find_edge(tri, incident, u, v)?.map(|h| (h / 3, h % 3));
    #[cfg(test)]
    assert_eq!(found, scan_halfedge(tri, u, v));
    Ok(found)
}

/// The same edge by scanning every triangle: the oracle the index is held
/// to under test.
#[cfg(test)]
fn scan_halfedge(tri: &Triangulation, u: u32, v: u32) -> Option<(usize, usize)> {
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

/// Every triangle edge that properly crosses the segment of `c`, as a
/// triangle and slot, in the order of the triangle array: both halfedges of
/// an edge, so the queue holds it twice, as a scan of every triangle did.
///
/// Found by walking the segment. A failed walk is a corrupted mesh. A vertex
/// on the segment is refused as crossing constraints: the pieces are split
/// at those first, so none is, and a recovery that met one would otherwise
/// flip around it.
fn find_crossing(
    tri: &Triangulation,
    incident: &[u32],
    c: Constraint,
) -> Result<Vec<(usize, usize)>, TriangulationError> {
    let walked = walk(tri, incident, c.a, c.b).map_err(|e| unplaceable(c, e))?;
    if !walked.through.is_empty() {
        return Err(TriangulationError::CrossingConstraints { a: c.a, b: c.b });
    }
    let mut halfedges: Vec<usize> = Vec::with_capacity(2 * walked.crossed.len());
    for e in walked.crossed {
        halfedges.push(e);
        halfedges.push(tri.halfedges[e] as usize);
    }
    halfedges.sort_unstable();
    Ok(halfedges
        .into_iter()
        // Never flip another constraint to satisfy this one.
        .filter(|&h| !tri.is_constrained(tri.triangles[h], tri.triangles[next_halfedge(h)]))
        .map(|h| (h / 3, h % 3))
        .collect())
}

/// The same crossing edges by scanning every triangle: the oracle the walk
/// is held to under test.
#[cfg(test)]
fn scan_crossing(
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
            if u == c.a || u == c.b || v == c.a || v == c.b {
                continue;
            }
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
fn flip(tri: &mut Triangulation, incident: &mut [u32], t: usize, i: usize) -> bool {
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
    // The four corners are the only vertices of the two rewritten triangles.
    for v in [apex, other, b] {
        incident[v as usize] = t as u32;
    }
    incident[a as usize] = ot as u32;
    note_update();
    debug_check(&tri.triangles, &tri.halfedges);
    true
}

#[cfg(test)]
mod tests {
    use axiolid_core::Point2;

    use super::{incident_triangles, outgoing, prev_halfedge};
    use crate::{triangulate, triangulate_domain, Constraint, Ring, TriangulationError};

    /// A jittered `side` x `side` grid and zigzag breaklines over it, from
    /// a fixed seed (the shape of the `triangulate_constraints` benchmark).
    fn zigzag(side: usize, lines: usize) -> (Vec<Point2>, Vec<Constraint>) {
        let mut s: u64 = 0x3550_0000 + side as u64;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let mut points = Vec::new();
        for row in 0..side {
            for col in 0..side {
                points.push(Point2::new(
                    col as f64 + next() * 0.6,
                    row as f64 + next() * 0.6,
                ));
            }
        }
        let at = |col: usize, row: usize| (row * side + col) as u32;
        let mut constraints = Vec::new();
        for line in 0..lines {
            let row = |col: usize| 2 * line + col % 2;
            for col in 0..side - 1 {
                constraints.push(Constraint::new(
                    at(col, row(col)),
                    at(col + 1, row(col + 1)),
                ));
            }
        }
        (points, constraints)
    }

    fn edges(tri: &crate::Triangulation) -> std::collections::HashSet<(u32, u32)> {
        let mut out = std::collections::HashSet::new();
        for t in tri.triangles().as_chunks::<3>().0 {
            for i in 0..3 {
                out.insert((t[i].min(t[(i + 1) % 3]), t[i].max(t[(i + 1) % 3])));
            }
        }
        out
    }

    /// Recovery found its edges through the index and the walk; the result
    /// is the one the scans over every triangle gave, byte for byte (the
    /// hash was taken from the scanning implementation).
    #[test]
    fn zigzag_breaklines_are_recovered_as_the_scans_did() {
        let (points, constraints) = zigzag(40, 18);
        let tri = triangulate(&points, &constraints).expect("recoverable");
        let present = edges(&tri);
        assert!(constraints.iter().all(|c| present.contains(&(c.a, c.b))));
        let mut hash: u64 = 1469598103934665603;
        for &v in tri.triangles() {
            hash = (hash ^ u64::from(v)).wrapping_mul(1099511628211);
        }
        assert_eq!(hash, 18275924131169444898);
    }

    /// Long constraints through exactly collinear grid vertices are split
    /// at them: the walk finds each one, in order.
    #[test]
    fn constraints_through_collinear_vertices_are_chains() {
        let side = 9u32;
        let points: Vec<Point2> = (0..side * side)
            .map(|k| Point2::new(f64::from(k % side), f64::from(k / side)))
            .collect();
        let at = |x: u32, y: u32| y * side + x;
        let constraints = [
            Constraint::new(at(0, 0), at(8, 8)),
            Constraint::new(at(8, 0), at(0, 8)),
            Constraint::new(at(0, 4), at(8, 4)),
            Constraint::new(at(4, 0), at(4, 8)),
            Constraint::new(at(0, 2), at(6, 8)),
        ];
        let tri = triangulate(&points, &constraints).expect("recoverable");
        let present = edges(&tri);
        for k in 0..8 {
            assert!(present.contains(&(at(k, k), at(k + 1, k + 1))));
            assert!(present.contains(&(at(8 - k, k), at(8 - k - 1, k + 1))));
            assert!(present.contains(&(at(k, 4), at(k + 1, 4))));
            assert!(present.contains(&(at(4, k), at(4, k + 1))));
        }
    }

    /// Two required edges that cross away from any vertex are refused.
    #[test]
    fn crossing_constraints_are_refused() {
        let (points, _) = zigzag(12, 0);
        let at = |col: u32, row: u32| row * 12 + col;
        let constraints = [
            Constraint::new(at(0, 0), at(11, 7)),
            Constraint::new(at(0, 7), at(11, 0)),
        ];
        assert!(matches!(
            triangulate(&points, &constraints),
            Err(TriangulationError::CrossingConstraints { .. })
        ));
    }

    fn square_with_centre() -> Vec<Point2> {
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec()
    }

    /// A constraint from a vertex to itself is kept as it is given and
    /// needs no edge, as it did when every constraint was looked up by scan.
    #[test]
    fn a_constraint_from_a_vertex_to_itself_is_kept() {
        let points = square_with_centre();
        let tri = triangulate(&points, &[Constraint::new(4, 4)]).expect("a self-loop is no edge");
        assert_eq!(tri.constraints(), &[Constraint::new(4, 4)]);
        let tri = triangulate_domain(
            &points,
            &[Constraint::new(4, 4)],
            &[Ring::outer([0, 1, 2, 3])],
        )
        .expect("a self-loop in a domain");
        // The domain keeps only constraints that are edges of its triangles.
        assert!(!tri.constraints().contains(&Constraint::new(4, 4)));
    }

    /// A NaN or infinite coordinate reads as collinear in the predicates,
    /// so it is refused before any triangle exists.
    #[test]
    fn non_finite_points_are_refused() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for (x, y) in [(bad, 0.5), (0.5, bad)] {
                let mut points = square_with_centre();
                points[3] = Point2::new(x, y);
                let constraints = [Constraint::new(0, 3), Constraint::new(1, 2)];
                assert_eq!(
                    triangulate(&points, &constraints).unwrap_err(),
                    TriangulationError::NonFinitePoint { index: 3 }
                );
                assert_eq!(
                    triangulate_domain(&points, &[], &[Ring::outer([0, 1, 2, 4])]).unwrap_err(),
                    TriangulationError::NonFinitePoint { index: 3 }
                );
            }
        }
    }

    /// A fan that never closes or ends is a corrupted mesh and is refused,
    /// not walked for ever.
    #[test]
    fn a_fan_that_does_not_end_is_refused() {
        let points = square_with_centre();
        let mut tri = triangulate(&points, &[]).expect("triangulates");
        let incident = incident_triangles(&tri);
        assert!(outgoing(&tri, &incident, 4).is_ok());
        let t = incident[4] as usize;
        let slot = (0..3).find(|&k| tri.triangles[3 * t + k] == 4).unwrap();
        let start = 3 * t + slot;
        // The next halfedge out of vertex 4 is another triangle's own twin
        // cycle, which never returns to the start.
        let other = (start + 3) % tri.triangles.len();
        tri.halfedges[prev_halfedge(start)] = other as u32;
        tri.halfedges[prev_halfedge(other)] = other as u32;
        assert!(outgoing(&tri, &incident, 4).is_err());

        // The same on the hull, where the fan is walked the other way round
        // from the start once the first way ends at the boundary.
        let mut tri = triangulate(&points, &[]).expect("triangulates");
        let incident = incident_triangles(&tri);
        assert!(outgoing(&tri, &incident, 0).is_ok());
        let t = incident[0] as usize;
        let slot = (0..3).find(|&k| tri.triangles[3 * t + k] == 0).unwrap();
        let start = 3 * t + slot;
        // Make the first way end at once, then send the other way into a
        // cycle: `start` leads to `other`, whose next edge is twinned back.
        tri.halfedges[prev_halfedge(start)] = u32::MAX;
        let other = (start + 3) % tri.triangles.len();
        tri.halfedges[start] = other as u32;
        let after = super::next_halfedge(other);
        tri.halfedges[after] = other as u32;
        assert!(outgoing(&tri, &incident, 0).is_err());
    }
}
