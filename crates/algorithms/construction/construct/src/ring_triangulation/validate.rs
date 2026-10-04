//! Refuse rings that do not bound a polygon with holes, and orient them.
//!
//! Every test is exact. Rings may not touch at all -- not even at a single
//! vertex -- because a pinch point extrudes to a non-manifold edge, and an
//! overlap has no single area to triangulate.

use core::ops::Range;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;

use super::{orient, segments_touch, Loop};

/// Validate `rings` (ranges into `points`, outer first) and return each as
/// a [`Loop`] with the polygon on its left.
pub(super) fn validate(points: &[Point2], rings: &[Range<usize>]) -> GeomResult<Vec<Loop>> {
    for (r, ring) in rings.iter().enumerate() {
        if ring.len() < 3 {
            return refuse(format!(
                "profile {} needs at least 3 vertices, got {}",
                name(r),
                ring.len()
            ));
        }
        if points[ring.clone()]
            .iter()
            .any(|p| !(p.x.is_finite() && p.y.is_finite()))
        {
            return refuse(format!("profile {} has a non-finite vertex", name(r)));
        }
        for k in 0..ring.len() {
            let next = (k + 1) % ring.len();
            if points[ring.start + k] == points[ring.start + next] {
                return refuse(format!("profile {} repeats vertex {k}", name(r)));
            }
        }
        for k in 0..ring.len() {
            let at = |i: usize| points[ring.start + i % ring.len()];
            if folds_back(at(k), at(k + 1), at(k + 2)) {
                return refuse(format!(
                    "profile {} folds back on itself at vertex {}",
                    name(r),
                    (k + 1) % ring.len()
                ));
            }
        }
    }
    check_edges(points, rings)?;

    let mut loops = Vec::with_capacity(rings.len());
    for (r, ring) in rings.iter().enumerate() {
        let mut vertices: Vec<u32> = ring.clone().map(|i| i as u32).collect();
        let counter_clockwise = ring_turns_left(&points[ring.clone()]);
        if counter_clockwise != (r == 0) {
            vertices.reverse();
        }
        loops.push(Loop { vertices });
    }

    let outer = &points[rings[0].clone()];
    for (h, hole) in rings.iter().enumerate().skip(1) {
        // No ring touches another, so one vertex places the whole hole.
        let probe = points[hole.start];
        if !strictly_inside(outer, probe) {
            return refuse(format!(
                "profile hole {} lies outside the outer ring",
                h - 1
            ));
        }
        for (g, other) in rings.iter().enumerate().skip(1) {
            if g != h && strictly_inside(&points[other.clone()], probe) {
                return refuse(format!("profile hole {} lies inside hole {}", h - 1, g - 1));
            }
        }
    }
    Ok(loops)
}

fn refuse<T>(message: String) -> GeomResult<T> {
    Err(GeomError::InvalidInput(message))
}

/// `outer ring` or `hole N`, for messages.
fn name(ring: usize) -> String {
    if ring == 0 {
        "outer ring".to_owned()
    } else {
        format!("hole {}", ring - 1)
    }
}

/// One ring edge: its ring, its position in the ring, and its endpoints.
#[derive(Debug, Clone, Copy)]
struct Edge {
    ring: usize,
    index: usize,
    a: Point2,
    b: Point2,
}

/// Refuse any two edges that share a point, except consecutive edges of one
/// ring meeting at their common vertex.
///
/// Edges are swept in order of their smallest x, so only pairs whose x
/// ranges overlap are compared.
fn check_edges(points: &[Point2], rings: &[Range<usize>]) -> GeomResult<()> {
    let mut edges = Vec::with_capacity(points.len());
    for (r, ring) in rings.iter().enumerate() {
        for k in 0..ring.len() {
            edges.push(Edge {
                ring: r,
                index: k,
                a: points[ring.start + k],
                b: points[ring.start + (k + 1) % ring.len()],
            });
        }
    }
    edges.sort_by(|e, f| e.a.x.min(e.b.x).total_cmp(&f.a.x.min(f.b.x)));
    for (i, e) in edges.iter().enumerate() {
        let right = e.a.x.max(e.b.x);
        let (low, high) = (e.a.y.min(e.b.y), e.a.y.max(e.b.y));
        for f in &edges[i + 1..] {
            if f.a.x.min(f.b.x) > right {
                break;
            }
            if f.a.y.max(f.b.y) < low || f.a.y.min(f.b.y) > high {
                continue;
            }
            check_pair(rings, e, f)?;
        }
    }
    Ok(())
}

fn check_pair(rings: &[Range<usize>], e: &Edge, f: &Edge) -> GeomResult<()> {
    if e.ring == f.ring {
        // Consecutive edges meet at their shared vertex, and folding back
        // there was refused already.
        let len = rings[e.ring].len();
        if (e.index + 1) % len == f.index || (f.index + 1) % len == e.index {
            return Ok(());
        }
        if segments_touch(e.a, e.b, f.a, f.b) {
            return refuse(format!("profile {} intersects itself", name(e.ring)));
        }
        return Ok(());
    }
    if !segments_touch(e.a, e.b, f.a, f.b) {
        return Ok(());
    }
    let (r, s) = (e.ring.min(f.ring), e.ring.max(f.ring));
    if r == 0 {
        refuse(format!(
            "profile hole {} touches or crosses the outer ring",
            s - 1
        ))
    } else {
        refuse(format!(
            "profile holes {} and {} overlap or touch",
            r - 1,
            s - 1
        ))
    }
}

/// Whether the path `a -> v -> b` turns back along itself: `b` collinear
/// with `a v` on the same side of `v` as `a`.
fn folds_back(a: Point2, v: Point2, b: Point2) -> bool {
    if orient(a, v, b) != 0 {
        return false;
    }
    let same = |p: f64, q: f64, o: f64| (p < o && q < o) || (p > o && q > o);
    same(a.x, b.x, v.x) || (a.x == v.x && same(a.y, b.y, v.y))
}

/// Whether a simple ring runs counter-clockwise, read exactly at its
/// lexicographically smallest vertex, which is strictly convex.
fn ring_turns_left(ring: &[Point2]) -> bool {
    let n = ring.len();
    let lowest = (0..n)
        .min_by(|&i, &j| {
            ring[i]
                .x
                .total_cmp(&ring[j].x)
                .then(ring[i].y.total_cmp(&ring[j].y))
        })
        .unwrap_or(0);
    orient(
        ring[(lowest + n - 1) % n],
        ring[lowest],
        ring[(lowest + 1) % n],
    ) > 0
}

/// Whether `p`, on no edge of the simple ring, lies inside it: the parity
/// of the edges crossing the horizontal ray to its right, each crossing
/// decided by `orient2d` with the half-open rule at vertices.
fn strictly_inside(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    for k in 0..ring.len() {
        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
        if (a.y > p.y) == (b.y > p.y) {
            continue;
        }
        // Upward edge: the crossing is right of `p` when `p` is left of it.
        let side = orient(a, b, p);
        if (b.y > a.y && side > 0) || (b.y < a.y && side < 0) {
            inside = !inside;
        }
    }
    inside
}
