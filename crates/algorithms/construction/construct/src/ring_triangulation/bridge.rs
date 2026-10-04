//! Join every hole to the outer boundary, giving one weakly simple ring.
//!
//! Holes are joined in order of decreasing largest x, each from its
//! rightmost vertex `m`. That order is what guarantees a bridge exists: no
//! hole still waiting lies right of `m`, so Eberly's argument (the ray from
//! `m` to the right reaches the current boundary, and the nearest edge it
//! hits offers a mutually visible vertex) holds against the boundary built
//! so far. The bridge is not taken from that ray, whose hit point is not
//! representable: candidates on the current boundary are tried nearest
//! first and each is accepted only when exact predicates prove the open
//! segment touches no edge of any ring or earlier bridge, and that it
//! enters the polygon at its boundary end. A vertex lying exactly on the ray or
//! on the candidate segment therefore blocks that candidate instead of
//! being crossed.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;

use super::{locally_inside, segments_touch, Loop};

/// A ring of vertex indices as a doubly linked list. Bridge endpoints
/// appear as two nodes with the same vertex.
#[derive(Debug, Clone, Copy)]
pub(super) struct Node {
    pub(super) vertex: u32,
    pub(super) prev: usize,
    pub(super) next: usize,
}

/// The weakly simple ring: its nodes and one node still on it.
#[derive(Debug, Clone)]
pub(super) struct Polygon {
    pub(super) nodes: Vec<Node>,
    pub(super) start: usize,
    pub(super) len: usize,
}

impl Polygon {
    /// Append `vertices` as a closed ring of new nodes; returns the first.
    fn push_ring(&mut self, vertices: &[u32]) -> usize {
        let first = self.nodes.len();
        let n = vertices.len();
        for (k, &vertex) in vertices.iter().enumerate() {
            self.nodes.push(Node {
                vertex,
                prev: first + (k + n - 1) % n,
                next: first + (k + 1) % n,
            });
        }
        first
    }

    /// The node indices on the main ring, from `start`.
    fn walk(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.len);
        let mut node = self.start;
        for _ in 0..self.len {
            out.push(node);
            node = self.nodes[node].next;
        }
        out
    }

    /// Splice the ring through `m` into the main ring at `p`:
    /// `p -> m -> ... -> m' -> p' -> next(p)`.
    fn splice(&mut self, p: usize, m: usize, hole_len: usize) {
        let p_next = self.nodes[p].next;
        let m_prev = self.nodes[m].prev;
        let p2 = self.nodes.len();
        let m2 = p2 + 1;
        self.nodes.push(Node {
            vertex: self.nodes[p].vertex,
            prev: m2,
            next: p_next,
        });
        self.nodes.push(Node {
            vertex: self.nodes[m].vertex,
            prev: m_prev,
            next: p2,
        });
        self.nodes[p_next].prev = p2;
        self.nodes[m_prev].next = m2;
        self.nodes[p].next = m;
        self.nodes[m].prev = p;
        self.len += hole_len + 2;
    }
}

/// Bridge every hole of `loops` (outer first) into one ring.
pub(super) fn bridge_holes(points: &[Point2], loops: &[Loop]) -> GeomResult<Polygon> {
    let at = |v: u32| points[v as usize];
    let mut polygon = Polygon {
        nodes: Vec::with_capacity(points.len() + 2 * loops.len()),
        start: 0,
        len: loops[0].vertices.len(),
    };
    polygon.push_ring(&loops[0].vertices);

    // Each hole's rightmost vertex (largest x, then largest y, then first),
    // and the holes by decreasing x of it.
    let mut holes: Vec<(usize, usize)> = (1..loops.len())
        .map(|h| {
            let vertices = &loops[h].vertices;
            let mut best = 0;
            for (k, &v) in vertices.iter().enumerate() {
                let (p, q) = (at(v), at(vertices[best]));
                if p.x > q.x || (p.x == q.x && p.y > q.y) {
                    best = k;
                }
            }
            (h, best)
        })
        .collect();
    holes.sort_by(|&(h, k), &(g, j)| {
        let (p, q) = (at(loops[h].vertices[k]), at(loops[g].vertices[j]));
        q.x.total_cmp(&p.x)
            .then(q.y.total_cmp(&p.y))
            .then(h.cmp(&g))
    });

    let mut pending: Vec<bool> = vec![true; loops.len()];
    pending[0] = false;
    for &(h, k) in &holes {
        pending[h] = false;
        let hole = &loops[h].vertices;
        let n = hole.len();
        let m = at(hole[k]);

        // Every edge a bridge must not touch: the main ring (with its
        // earlier bridges), this hole, and the holes still waiting.
        let main = polygon.walk();
        let mut edges: Vec<(Point2, Point2)> = main
            .iter()
            .map(|&i| {
                let node = polygon.nodes[i];
                (at(node.vertex), at(polygon.nodes[node.next].vertex))
            })
            .collect();
        for (g, ring) in loops.iter().enumerate() {
            if g == h || pending[g] {
                let len = ring.vertices.len();
                edges.extend(
                    (0..len).map(|i| (at(ring.vertices[i]), at(ring.vertices[(i + 1) % len]))),
                );
            }
        }

        let mut candidates = main.clone();
        let distance = |i: usize| {
            let p = at(polygon.nodes[i].vertex);
            (p.x - m.x).powi(2) + (p.y - m.y).powi(2)
        };
        candidates.sort_by(|&i, &j| distance(i).total_cmp(&distance(j)).then(i.cmp(&j)));
        let found = candidates.into_iter().find(|&i| {
            let node = polygon.nodes[i];
            let p = at(node.vertex);
            let (p_prev, p_next) = (
                at(polygon.nodes[node.prev].vertex),
                at(polygon.nodes[node.next].vertex),
            );
            locally_inside(p_prev, p, p_next, m) && clear(m, p, &edges)
        });
        let Some(p) = found else {
            return Err(GeomError::Degenerate(format!(
                "profile hole {} has no vertex it can be bridged to",
                h - 1
            )));
        };
        let first = polygon.push_ring(hole);
        polygon.splice(p, first + k, n);
    }
    Ok(polygon)
}

/// Whether the segment `m p` meets no edge except at `m` or `p`.
///
/// Edges ending at `m` or `p` are not tested: one running along the
/// segment ends at a vertex on it, which the next edge from that vertex
/// touches. Nor is the hole's side tested locally: leaving `m` into the
/// hole, the segment would have to cross the hole's own edges to reach
/// `p`, which lies outside it. Only `p` can have several nodes (bridge
/// ends), so only there does the sector test pick one.
fn clear(m: Point2, p: Point2, edges: &[(Point2, Point2)]) -> bool {
    let (left, right) = (m.x.min(p.x), m.x.max(p.x));
    let (low, high) = (m.y.min(p.y), m.y.max(p.y));
    edges.iter().all(|&(a, b)| {
        if a.x.max(b.x) < left || a.x.min(b.x) > right || a.y.max(b.y) < low || a.y.min(b.y) > high
        {
            return true;
        }
        let at_end = |q: Point2| q == m || q == p;
        at_end(a) || at_end(b) || !segments_touch(m, p, a, b)
    })
}
