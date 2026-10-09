//! Refuse rings that do not bound a polygon with holes, and orient them.
//!
//! Every test is exact. Under [`PinchPolicy::Refuse`] rings may not touch
//! at all -- not even at a single vertex -- because a pinch point extrudes
//! to a non-manifold edge, and an overlap has no single area to
//! triangulate. Under [`PinchPolicy::Accept`] two edges may meet at one
//! point that is a vertex of at least one of them (#262): the meeting is
//! returned as a [`Touch`] when it lies inside the other edge, and the
//! [`pinch`](super::pinch) path triangulates around it. Edges that cross,
//! or overlap along a stretch, are refused under both policies.

use core::ops::Range;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;

use super::{by_value, lexicographic, orient, segments_touch, within, Loop, PinchPolicy};

/// A vertex lying inside an edge of a ring, away from its ends: the edge
/// `index` of ring `ring` (from its vertex `index` to the next) passes
/// through point `vertex`. Only [`PinchPolicy::Accept`] returns these.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Touch {
    pub(super) ring: usize,
    pub(super) index: usize,
    pub(super) vertex: u32,
}

/// The checks of one ring on its own: enough vertices, all finite, none
/// repeated next to itself, no fold back.
pub(super) fn check_rings(points: &[Point2], rings: &[Range<usize>]) -> GeomResult<()> {
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
    Ok(())
}

/// Orient rings that touch nowhere, each as a [`Loop`] with the polygon on
/// its left, and check that every hole lies inside the outer ring and
/// outside every other hole.
pub(super) fn orient_and_place(points: &[Point2], rings: &[Range<usize>]) -> GeomResult<Vec<Loop>> {
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

pub(super) fn refuse<T>(message: String) -> GeomResult<T> {
    Err(GeomError::InvalidInput(message))
}

/// `outer ring` or `hole N`, for messages.
pub(super) fn name(ring: usize) -> String {
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
    /// The indices of `a` and `b` in the point list.
    ia: u32,
    ib: u32,
}

/// Refuse any two edges that share a point, except consecutive edges of one
/// ring meeting at their common vertex, and under [`PinchPolicy::Accept`]
/// edges meeting at one point that is a vertex of one of them; the latter
/// are returned when the point lies inside the other edge.
///
/// Edges are swept in order of their smallest x, so only pairs whose x
/// ranges overlap are compared.
pub(super) fn check_edges(
    points: &[Point2],
    rings: &[Range<usize>],
    policy: PinchPolicy,
) -> GeomResult<Vec<Touch>> {
    let mut edges = Vec::with_capacity(points.len());
    for (r, ring) in rings.iter().enumerate() {
        for k in 0..ring.len() {
            let (ia, ib) = (ring.start + k, ring.start + (k + 1) % ring.len());
            edges.push(Edge {
                ring: r,
                index: k,
                a: points[ia],
                b: points[ib],
                ia: ia as u32,
                ib: ib as u32,
            });
        }
    }
    let mut touches = Vec::new();
    edges.sort_by(|e, f| by_value(e.a.x.min(e.b.x), f.a.x.min(f.b.x)));
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
            match policy {
                PinchPolicy::Refuse => check_pair(rings, e, f)?,
                PinchPolicy::Accept => check_pair_pinched(rings, e, f, &mut touches)?,
            }
        }
    }
    Ok(touches)
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

/// [`check_pair`] under [`PinchPolicy::Accept`]: two edges may meet at a
/// single point that is an end of at least one of them. A shared end is a
/// pinch, found later from the coordinates; an end inside the other edge
/// is recorded in `touches`. Crossing and overlapping along a stretch stay
/// refused.
fn check_pair_pinched(
    rings: &[Range<usize>],
    e: &Edge,
    f: &Edge,
    touches: &mut Vec<Touch>,
) -> GeomResult<()> {
    if e.ring == f.ring {
        let len = rings[e.ring].len();
        if (e.index + 1) % len == f.index || (f.index + 1) % len == e.index {
            return Ok(());
        }
    }
    let touching = segments_touch(e.a, e.b, f.a, f.b);
    if !touching {
        return Ok(());
    }
    // A shared end: fine unless the other ends run back along one line
    // (an overlap) or coincide too (the same edge twice).
    for (p, far_e) in [(e.a, e.b), (e.b, e.a)] {
        for (q, far_f) in [(f.a, f.b), (f.b, f.a)] {
            if p == q {
                if far_e == far_f || folds_back(far_e, p, far_f) {
                    return overlap(e, f);
                }
                return Ok(());
            }
        }
    }
    // No shared end. Exactly one end inside the other edge is a touch;
    // anything else (a proper crossing, or both edges on one line, which
    // then overlap) is refused.
    let inside = |a: Point2, b: Point2, c: Point2| orient(a, b, c) == 0 && within(a, b, c);
    let hits = [
        (inside(e.a, e.b, f.a), e, f.ia),
        (inside(e.a, e.b, f.b), e, f.ib),
        (inside(f.a, f.b, e.a), f, e.ia),
        (inside(f.a, f.b, e.b), f, e.ib),
    ];
    let mut found = hits.iter().filter(|hit| hit.0);
    match (found.next(), found.next()) {
        (Some(&(_, edge, vertex)), None) => {
            touches.push(Touch {
                ring: edge.ring,
                index: edge.index,
                vertex,
            });
            Ok(())
        }
        (Some(_), Some(_)) => overlap(e, f),
        _ => crossing(e, f),
    }
}

fn overlap(e: &Edge, f: &Edge) -> GeomResult<()> {
    let (r, s) = (e.ring.min(f.ring), e.ring.max(f.ring));
    if r == s {
        refuse(format!("profile {} overlaps itself", name(r)))
    } else if r == 0 {
        refuse(format!("profile hole {} overlaps the outer ring", s - 1))
    } else {
        refuse(format!("profile holes {} and {} overlap", r - 1, s - 1))
    }
}

fn crossing(e: &Edge, f: &Edge) -> GeomResult<()> {
    let (r, s) = (e.ring.min(f.ring), e.ring.max(f.ring));
    if r == s {
        refuse(format!("profile {} intersects itself", name(r)))
    } else if r == 0 {
        refuse(format!("profile hole {} crosses the outer ring", s - 1))
    } else {
        refuse(format!("profile holes {} and {} cross", r - 1, s - 1))
    }
}

/// Whether the path `a -> v -> b` turns back along itself: `b` collinear
/// with `a v` on the same side of `v` as `a`.
pub(super) fn folds_back(a: Point2, v: Point2, b: Point2) -> bool {
    if orient(a, v, b) != 0 {
        return false;
    }
    let same = |p: f64, q: f64, o: f64| (p < o && q < o) || (p > o && q > o);
    same(a.x, b.x, v.x) || (a.x == v.x && same(a.y, b.y, v.y))
}

/// Whether a simple ring runs counter-clockwise, read exactly at its
/// lexicographically smallest vertex, which is strictly convex.
///
/// The order is by value ([`lexicographic`]). Under `total_cmp` a straight
/// corner at `(-0.0, 2)` came before the true lowest corner `(0, 0)`, its
/// zero turn read as clockwise, and the ring was reversed into one with no
/// ear (#269).
pub(super) fn ring_turns_left(ring: &[Point2]) -> bool {
    let n = ring.len();
    let lowest = (0..n)
        .min_by(|&i, &j| lexicographic(ring[i], ring[j]))
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
pub(super) fn strictly_inside(ring: &[Point2], p: Point2) -> bool {
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
