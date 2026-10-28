// SPDX-License-Identifier: MPL-2.0

//! Triangulating a domain: the area inside an outline, minus its holes.
//!
//! [`triangulate`](crate::triangulate) covers the convex hull of its input.
//! A caller with a floor outline, a courtyard hole and an island in it
//! wants the triangles of that region only. [`triangulate_domain`] takes the
//! same points and constraints plus the rings that bound the region, and
//! keeps exactly the triangles inside.
//!
//! # The rule: winding with ring roles
//!
//! Each ring is tagged [`RingKind::Outer`] or [`RingKind::Hole`] and may be
//! given either way round. It is normalised (outer counter-clockwise, hole
//! clockwise), and then every piece of its boundary adds one to the winding
//! of the region on its left. The domain is the set of triangles with
//! winding exactly 1: inside one outer ring (or inside an island of a hole,
//! which is an outer ring again) and in no hole.
//!
//! The winding is computed by a flood fill over the triangles' shared
//! edges, seeded outside the hull at winding 0. Crossing a ring piece from
//! its left to its right changes the winding by minus the piece's net
//! contribution. Only ring pieces carry winding; breaklines never do, so a
//! breakline cannot change the domain. The fill is path independent: the
//! winding difference between two triangles is the signed number of ring
//! boundaries between them, whatever path joins them, so a revisit can only
//! disagree on corrupted input; it is asserted in debug builds and not
//! re-checked in release builds. A first visit with a winding outside
//! `{0, 1}` is refused as [`TriangulationError::MisnestedRing`], and so is
//! an empty result as [`TriangulationError::EmptyDomain`].
//!
//! Rejected alternatives:
//!
//! - Even-odd parity over unlabelled rings cannot tell a hole outside the
//!   outline from an island, and a hole listed twice cancels silently.
//! - Parity plus a per-edge role check still cancels a duplicate hole, and
//!   refuses valid input where three rings share an edge.
//!
//! # Rings
//!
//! A ring lists indices into the points in order and is closed implicitly;
//! repeating the first index at the end is refused as a repeated vertex.
//! Rings may share vertices and edges with each other, but a ring may not
//! visit one vertex twice. Its edges become constraints, so a ring edge
//! through input vertices is split at them like any constraint (#360), and
//! a breakline that crosses a ring edge between vertices is refused as
//! [`TriangulationError::CrossingConstraints`]. Any two crossing
//! breaklines are refused too, even outside the domain.
//!
//! The orientation of a ring is read exactly at its lexicographically
//! lowest vertex, which is always a convex corner. The fold check (both
//! neighbours of a vertex on one side of it on a common line) guarantees
//! that turn is non-zero.
//!
//! # Points outside the domain
//!
//! A point outside the domain, or inside a hole, is dropped from the
//! triangles and reported by [`Triangulation::outside_points`]; it is not an
//! error, because a caller's point cloud rarely respects the outline.
//! [`Triangulation::points`] still holds every point, so indices stay
//! stable. Inside the domain the triangulation is constrained Delaunay:
//! constraints block visibility, and so the dropped points do not change it.
//!
//! A constraint entirely outside the domain is removed from
//! [`Triangulation::constraints`], so every listed constraint is an edge of
//! an output triangle.

use axiolid_core::Point2;

use crate::build::build;
use crate::mesh::{Triangulation, TriangulationError, NO_HALFEDGE};
use crate::recover::Splits;
use crate::{collinear, turns_left, Constraint};

/// What a ring bounds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingKind {
    /// The interior of the ring is part of the domain. An island inside a
    /// hole is an outer ring again.
    Outer,
    /// The interior of the ring is excluded from the domain.
    Hole,
}

/// A closed boundary of the domain, as indices into the input points.
///
/// The ring is closed implicitly (the last vertex joins the first) and may
/// be listed in either orientation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Ring {
    /// The vertices in order, without repeating the first at the end.
    pub vertices: Vec<u32>,
    /// Whether the ring bounds the domain or a hole in it.
    pub kind: RingKind,
}

impl Ring {
    /// An outer ring through `vertices`.
    #[must_use]
    pub fn outer(vertices: impl Into<Vec<u32>>) -> Self {
        Self {
            vertices: vertices.into(),
            kind: RingKind::Outer,
        }
    }

    /// A hole bounded by `vertices`.
    #[must_use]
    pub fn hole(vertices: impl Into<Vec<u32>>) -> Self {
        Self {
            vertices: vertices.into(),
            kind: RingKind::Hole,
        }
    }
}

/// Build a constrained Delaunay triangulation of the domain bounded by
/// `rings`.
///
/// The result keeps exactly the triangles inside an outer ring and outside
/// every hole (see the module documentation for the rule). The ring edges
/// are constraints, as are `constraints`; the points that fall outside the
/// domain are listed by [`Triangulation::outside_points`].
///
/// # Errors
///
/// Everything [`triangulate`](crate::triangulate) refuses, and
/// [`TriangulationError::RingOutOfRange`],
/// [`TriangulationError::TooFewRingVertices`],
/// [`TriangulationError::RepeatedRingVertex`],
/// [`TriangulationError::FoldedRing`], [`TriangulationError::MisnestedRing`]
/// and [`TriangulationError::EmptyDomain`] (also for an empty ring list).
/// Checks run in this order, ring by ring and vertex by vertex.
pub fn triangulate_domain(
    points: &[Point2],
    constraints: &[Constraint],
    rings: &[Ring],
) -> Result<Triangulation, TriangulationError> {
    let built = build(points, constraints, Some(rings))?;
    let mut tri = built.tri;
    classify(&mut tri, &built.kept, &built.splits, rings)?;
    tri.domain = Some(rings.to_vec());
    Ok(tri)
}

/// Re-triangulate `points` the way `domain` says: the whole hull when it is
/// `None`, the domain's triangles otherwise. Refinement uses it so that a
/// Steiner point does not re-expand the mesh to the hull.
pub(crate) fn rebuild(
    points: &[Point2],
    constraints: &[Constraint],
    domain: Option<&[Ring]>,
) -> Result<Triangulation, TriangulationError> {
    match domain {
        None => crate::triangulate(points, constraints),
        Some(rings) => triangulate_domain(points, constraints, rings),
    }
}

/// Range and size of every ring, then the empty ring list.
pub(crate) fn check_ring_sizes(rings: &[Ring], count: u32) -> Result<(), TriangulationError> {
    for (r, ring) in rings.iter().enumerate() {
        let r = r as u32;
        if let Some(&index) = ring.vertices.iter().find(|&&v| v >= count) {
            return Err(TriangulationError::RingOutOfRange { ring: r, index });
        }
        if ring.vertices.len() < 3 {
            return Err(TriangulationError::TooFewRingVertices { ring: r });
        }
    }
    if rings.is_empty() {
        return Err(TriangulationError::EmptyDomain);
    }
    Ok(())
}

/// The ring's vertices as kept vertices.
fn kept_vertices(ring: &Ring, kept: &[u32]) -> Vec<u32> {
    ring.vertices.iter().map(|&v| kept[v as usize]).collect()
}

/// Lexicographic order of two positions; `-0.0` equals `0.0`.
fn lex(p: Point2, q: Point2) -> core::cmp::Ordering {
    (p.x, p.y)
        .partial_cmp(&(q.x, q.y))
        .unwrap_or(core::cmp::Ordering::Equal)
}

/// Check every ring's shape and return its edges as constraints between
/// kept vertices, ring by ring.
///
/// A repeated vertex (after merging coincident points) is refused first,
/// then a fold, per ring. The constraints are appended after the caller's.
pub(crate) fn check_rings(
    points: &[Point2],
    kept: &[u32],
    rings: &[Ring],
) -> Result<Vec<Constraint>, TriangulationError> {
    let at = |k: u32| points[k as usize];
    let mut stamp = vec![u32::MAX; points.len()];
    let mut edges = Vec::new();
    for (r, ring) in rings.iter().enumerate() {
        let r = r as u32;
        let seq = kept_vertices(ring, kept);
        for (i, &k) in seq.iter().enumerate() {
            if stamp[k as usize] == r {
                return Err(TriangulationError::RepeatedRingVertex {
                    ring: r,
                    index: ring.vertices[i],
                });
            }
            stamp[k as usize] = r;
        }
        let m = seq.len();
        for i in 0..m {
            let (prev, v, next) = (seq[(i + m - 1) % m], seq[i], seq[(i + 1) % m]);
            let folds = collinear(at(prev), at(v), at(next))
                && lex(at(prev), at(v)) == lex(at(next), at(v));
            if folds {
                return Err(TriangulationError::FoldedRing {
                    ring: r,
                    index: ring.vertices[i],
                });
            }
            edges.push(Constraint::new(v, next));
        }
    }
    Ok(edges)
}

/// Whether the ring `seq` runs counter-clockwise, read at its lowest vertex.
///
/// The lowest vertex by position is a convex corner of any simple ring, and
/// `check_rings` has already refused a fold there, so the turn is non-zero.
fn is_ccw(points: &[Point2], seq: &[u32]) -> bool {
    let m = seq.len();
    let at = |i: usize| points[seq[i % m] as usize];
    let lowest = (0..m)
        .min_by(|&i, &j| lex(at(i), at(j)))
        .unwrap_or_default();
    turns_left(at(lowest + m - 1), at(lowest), at(lowest + 1))
}

/// The net winding contribution of ring pieces, by piece `(min, max)`.
struct Pieces(Vec<((u32, u32), i32, u32)>);

impl Pieces {
    /// Net contribution along `u -> v` (crossing the piece from its left to
    /// its right lowers the winding by it) and the lowest ring carrying it.
    fn along(&self, u: u32, v: u32) -> (i32, u32) {
        let key = (u.min(v), u.max(v));
        match self.0.binary_search_by_key(&key, |e| e.0) {
            Ok(i) if u < v => (self.0[i].1, self.0[i].2),
            Ok(i) => (-self.0[i].1, self.0[i].2),
            Err(_) => (0, 0),
        }
    }
}

/// The ring pieces: each ring edge as the chain it was split into, with the
/// sign that normalises its ring.
fn ring_pieces(points: &[Point2], kept: &[u32], splits: &Splits, rings: &[Ring]) -> Pieces {
    let mut raw: Vec<((u32, u32), i32, u32)> = Vec::new();
    for (r, ring) in rings.iter().enumerate() {
        let seq = kept_vertices(ring, kept);
        let outer = ring.kind == RingKind::Outer;
        let sign = if is_ccw(points, &seq) == outer { 1 } else { -1 };
        let m = seq.len();
        for i in 0..m {
            let (u, v) = (seq[i], seq[(i + 1) % m]);
            let key = Constraint::new(u, v);
            let mut chain: Vec<u32> = splits
                .chain(key)
                .map_or_else(|| vec![key.a, key.b], <[u32]>::to_vec);
            if u != key.a {
                chain.reverse();
            }
            for w in chain.windows(2) {
                let (p, q) = (w[0], w[1]);
                if p < q {
                    raw.push(((p, q), sign, r as u32));
                } else {
                    raw.push(((q, p), -sign, r as u32));
                }
            }
        }
    }
    raw.sort_unstable();
    let mut net: Vec<((u32, u32), i32, u32)> = Vec::new();
    for (key, delta, ring) in raw {
        match net.last_mut() {
            Some(last) if last.0 == key => {
                last.1 += delta;
                last.2 = last.2.min(ring);
            }
            _ => net.push((key, delta, ring)),
        }
    }
    Pieces(net)
}

/// Keep the triangles with winding 1 and rebuild what depends on them.
fn classify(
    tri: &mut Triangulation,
    kept: &[u32],
    splits: &Splits,
    rings: &[Ring],
) -> Result<(), TriangulationError> {
    let pieces = ring_pieces(&tri.points, kept, splits, rings);
    let winding = flood(tri, &pieces)?;

    let mut triangles = Vec::new();
    for (t, corners) in tri.triangles.as_chunks::<3>().0.iter().enumerate() {
        if winding[t] == 1 {
            triangles.extend_from_slice(corners);
        }
    }
    if triangles.is_empty() {
        return Err(TriangulationError::EmptyDomain);
    }
    tri.triangles = triangles;
    tri.rebuild_halfedges();

    let mut edges: Vec<(u32, u32)> = tri
        .triangles
        .as_chunks::<3>()
        .0
        .iter()
        .flat_map(|c| (0..3).map(move |i| (c[i].min(c[(i + 1) % 3]), c[i].max(c[(i + 1) % 3]))))
        .collect();
    edges.sort_unstable();
    edges.dedup();
    tri.constraints
        .retain(|c| edges.binary_search(&(c.a, c.b)).is_ok());

    let mut used = vec![false; tri.points.len()];
    for &v in &tri.triangles {
        used[v as usize] = true;
    }
    let mut merged = vec![false; tri.points.len()];
    for m in &tri.merged {
        merged[m.input as usize] = true;
    }
    tri.outside = (0..tri.points.len() as u32)
        .filter(|&v| !used[v as usize] && !merged[v as usize])
        .collect();
    Ok(())
}

/// The winding of every triangle, by a first-in first-out flood fill.
fn flood(tri: &Triangulation, pieces: &Pieces) -> Result<Vec<i32>, TriangulationError> {
    use std::collections::VecDeque;

    let edge = |t: usize, i: usize| (tri.triangles[3 * t + i], tri.triangles[3 * t + (i + 1) % 3]);
    let seed = (0..tri.triangle_count())
        .find_map(|t| {
            (0..3)
                .find(|&i| tri.halfedges[3 * t + i] == NO_HALFEDGE)
                .map(|i| (t, i))
        })
        .ok_or(TriangulationError::EmptyDomain)?;

    // Entering the hull across its edge `u -> v` from outside (its right)
    // raises the winding by that edge's net contribution.
    let (u, v) = edge(seed.0, seed.1);
    let (net, ring) = pieces.along(u, v);
    check_winding(net, ring)?;
    let mut winding = vec![i32::MIN; tri.triangle_count()];
    winding[seed.0] = net;
    let mut queue = VecDeque::from([seed.0]);
    while let Some(t) = queue.pop_front() {
        for i in 0..3 {
            let twin = tri.halfedges[3 * t + i];
            let (u, v) = edge(t, i);
            let (net, ring) = pieces.along(u, v);
            if twin == NO_HALFEDGE {
                debug_assert_eq!(winding[t], net, "hull edge disagrees with the flood");
                continue;
            }
            let next = twin as usize / 3;
            let w = winding[t] - net;
            if winding[next] == i32::MIN {
                check_winding(w, ring)?;
                winding[next] = w;
                queue.push_back(next);
            } else {
                debug_assert_eq!(winding[next], w, "the winding is path dependent");
            }
        }
    }
    Ok(winding)
}

/// Refuse a region whose winding is neither outside (0) nor inside (1).
fn check_winding(winding: i32, ring: u32) -> Result<(), TriangulationError> {
    if (0..=1).contains(&winding) {
        Ok(())
    } else {
        Err(TriangulationError::MisnestedRing { ring, winding })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TriangulationError as E;

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    fn square() -> Vec<Point2> {
        vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)]
    }

    fn refused(points: &[Point2], rings: &[Ring]) -> E {
        triangulate_domain(points, &[], rings).unwrap_err()
    }

    #[test]
    fn a_ring_out_of_range_is_refused() {
        assert_eq!(
            refused(&square(), &[Ring::outer([0, 1, 9])]),
            E::RingOutOfRange { ring: 0, index: 9 }
        );
    }

    #[test]
    fn a_short_ring_is_refused() {
        assert_eq!(
            refused(&square(), &[Ring::outer([0, 1])]),
            E::TooFewRingVertices { ring: 0 }
        );
    }

    #[test]
    fn closing_a_ring_explicitly_is_a_repeat() {
        assert_eq!(
            refused(&square(), &[Ring::outer([0, 1, 2, 3, 0])]),
            E::RepeatedRingVertex { ring: 0, index: 0 }
        );
    }

    #[test]
    fn a_coincident_index_is_a_repeat() {
        let mut pts = square();
        pts.push(p(0.0, 0.0));
        assert_eq!(
            refused(&pts, &[Ring::outer([0, 1, 2, 3, 4])]),
            E::RepeatedRingVertex { ring: 0, index: 4 }
        );
    }

    #[test]
    fn three_collinear_vertices_fold() {
        let pts = vec![p(0.0, 0.0), p(1.0, 0.0), p(2.0, 0.0), p(0.0, 3.0)];
        assert_eq!(
            refused(&pts, &[Ring::outer([0, 1, 2])]),
            E::FoldedRing { ring: 0, index: 0 }
        );
    }

    #[test]
    fn a_spike_folds() {
        // A spike out to 4 and back along the same line through 5 folds
        // at 4.
        let pts = vec![
            p(0.0, 0.0),
            p(4.0, 0.0),
            p(4.0, 4.0),
            p(0.0, 4.0),
            p(6.0, 0.0),
            p(5.0, 0.0),
        ];
        assert_eq!(
            refused(&pts, &[Ring::outer([0, 1, 4, 5, 2, 3])]),
            E::FoldedRing { ring: 0, index: 4 }
        );
    }

    #[test]
    fn a_straight_vertex_is_accepted() {
        let mut pts = square();
        pts.push(p(2.0, 0.0));
        let tri = triangulate_domain(&pts, &[], &[Ring::outer([0, 4, 1, 2, 3])]).unwrap();
        assert_eq!(tri.triangle_count(), 3);
        assert!(tri.outside_points().is_empty());
    }

    #[test]
    fn orientation_is_read_at_a_reflex_first_vertex() {
        // A reflex first vertex: the orientation is read at the lowest
        // vertex (0, 0), a convex corner, whichever way the ring is listed.
        let pts = vec![
            p(2.0, 1.0),
            p(4.0, 0.0),
            p(4.0, 4.0),
            p(0.0, 4.0),
            p(0.0, 0.0),
        ];
        assert!(is_ccw(&pts, &[0, 1, 2, 3, 4]));
        assert!(!is_ccw(&pts, &[4, 3, 2, 1, 0]));
        let a = triangulate_domain(&pts, &[], &[Ring::outer([0, 1, 2, 3, 4])]).unwrap();
        let b = triangulate_domain(&pts, &[], &[Ring::outer([4, 3, 2, 1, 0])]).unwrap();
        assert_eq!(a.triangle_count(), 3);
        assert_eq!(b.triangle_count(), 3);
    }

    #[test]
    fn an_empty_ring_list_is_an_empty_domain() {
        assert_eq!(refused(&square(), &[]), E::EmptyDomain);
    }
}
