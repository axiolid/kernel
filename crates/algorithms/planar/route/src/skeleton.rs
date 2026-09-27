//! The skeleton of a region with holes: its corridors as a graph, with
//! path ends, junctions and the clearance at each (#139).
//!
//! # Approximate where it must be, certified where it can be
//!
//! The medial axis of a polygon with holes has parabolic arcs and
//! algebraic vertices. What a circulation check needs from it is its shape
//! as a graph -- where paths end, where they meet -- and how much room
//! there is along it. So the shape is approximated and the room is proven:
//!
//! - The boundary is sampled no coarser than `spacing` and triangulated
//!   with its edges as constraints. Each triangle inside the region gives
//!   a node at its circumcentre -- the centre of an empty circle touching
//!   three boundary samples, a Voronoi vertex of the samples -- joined to
//!   the nodes of the triangles it shares an unconstrained edge with. As
//!   the spacing shrinks these converge to the medial axis; how close they
//!   are is not proven, and nothing below depends on it. Nodes the
//!   rounding or a constrained triangle puts outside the region are
//!   dropped.
//! - Spurs -- branches into convex corners, which every corner has -- are
//!   pruned as in the lambda-medial axis: a node stays when the feet on its
//!   nearest walls lie at least `prune` times its clearance apart. Across
//!   a corridor they are two clearances apart; into a right-angled corner
//!   only the square root of two. So 1.5 keeps corridors, junctions and
//!   dead ends and drops the spurs into right-angled and sharper corners;
//!   leaf branches then left shorter than the clearance where they join
//!   go too. What remains ends where corridors and rooms end. The spacing
//!   should be under an eighth of the narrowest width.
//! - Every node is decided to lie in the region exactly, and its
//!   [`SkeletonNode::clearance`] -- the distance to the nearest wall -- is
//!   an interval proven to contain the true value (outward-rounded
//!   point-to-segment distances). A path end also names the wall ahead of
//!   it: the first wall the path, continued straight on, runs into,
//!   decided by exact crossing tests.

use axiolid_contracts::Sign;
use axiolid_core::Point2;
use axiolid_overlay::Polygon;
use axiolid_triangulate::{triangulate, Constraint};
use std::collections::HashMap;

use crate::{
    contains, crosses, dedup_points, ring_edges, side, validate_region, within, RouteError,
};

/// What a skeleton node is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum NodeKind {
    /// Where a path ends: one neighbour.
    End,
    /// Along a path: two neighbours.
    Path,
    /// Where paths meet: three or more.
    Junction,
    /// No neighbours: a region too small to hold a path.
    Isolated,
}

/// A wall edge: which polygon, which ring (0 the outer, `k` the `k`-th
/// hole) and which edge of it (from point `edge` to the next).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wall {
    /// Polygon index in the region.
    pub polygon: usize,
    /// Ring: 0 for the outer ring, `k` for hole `k - 1`.
    pub ring: usize,
    /// Edge index within the ring.
    pub edge: usize,
}

/// One node of the skeleton.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct SkeletonNode {
    /// Where it is; inside the region, decided exactly.
    pub point: Point2,
    /// Its kind, by its number of neighbours.
    pub kind: NodeKind,
    /// Proven bounds `(lower, upper)` on the distance to the nearest wall.
    pub clearance: (f64, f64),
    /// For an end: the wall the path runs into if continued straight on.
    pub ahead: Option<Wall>,
}

/// The skeleton: nodes and the edges joining them.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Skeleton {
    /// The nodes.
    pub nodes: Vec<SkeletonNode>,
    /// Edges as node index pairs, lower first, sorted.
    pub edges: Vec<(usize, usize)>,
    /// The boundary sample spacing used.
    pub spacing: f64,
}

impl Skeleton {
    /// Indices of the path ends.
    #[must_use]
    pub fn ends(&self) -> Vec<usize> {
        self.kinds(NodeKind::End)
    }

    /// Indices of the junctions.
    #[must_use]
    pub fn junctions(&self) -> Vec<usize> {
        self.kinds(NodeKind::Junction)
    }

    fn kinds(&self, kind: NodeKind) -> Vec<usize> {
        (0..self.nodes.len())
            .filter(|&i| self.nodes[i].kind == kind)
            .collect()
    }
}

/// Why no skeleton was built.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum SkeletonError {
    /// The region is malformed, or a predicate was undecidable.
    Route(RouteError),
    /// The spacing or prune factor is not positive and finite.
    InvalidParameter,
    /// The boundary could not be triangulated: rings cross.
    Triangulation,
}

impl From<RouteError> for SkeletonError {
    fn from(error: RouteError) -> Self {
        Self::Route(error)
    }
}

/// The skeleton of `region`, its boundary sampled at most `spacing` apart,
/// keeping the nodes whose nearest walls spread at least `prune` times
/// their clearance apart (1.5 drops the spurs into right-angled corners; 0
/// keeps every node).
///
/// # Errors
///
/// [`SkeletonError`] for a malformed region or parameters.
pub fn skeleton(region: &[Polygon], spacing: f64, prune: f64) -> Result<Skeleton, SkeletonError> {
    if !(spacing.is_finite() && spacing > 0.0 && prune.is_finite() && prune >= 0.0) {
        return Err(SkeletonError::InvalidParameter);
    }
    validate_region(region, &[])?;
    // Walls, and the boundary sampled along them.
    let mut walls: Vec<(Wall, Point2, Point2)> = Vec::new();
    for (pi, polygon) in region.iter().enumerate() {
        for (ri, ring) in std::iter::once(&polygon.outer)
            .chain(&polygon.holes)
            .enumerate()
        {
            for (ei, (a, b)) in ring_edges(ring).into_iter().enumerate() {
                if a != b {
                    let wall = Wall {
                        polygon: pi,
                        ring: ri,
                        edge: ei,
                    };
                    walls.push((wall, a, b));
                }
            }
        }
    }
    let mut points: Vec<Point2> = Vec::new();
    let mut pieces: Vec<(Point2, Point2)> = Vec::new();
    for &(_, a, b) in &walls {
        let count = ((b - a).length() / spacing).ceil().max(1.0) as usize;
        let mut prev = a;
        for k in 1..=count {
            let q = if k == count {
                b
            } else {
                a + (b - a) * (k as f64 / count as f64)
            };
            pieces.push((prev, q));
            points.push(prev);
            prev = q;
        }
    }
    dedup_points(&mut points);
    let mut index: HashMap<(u64, u64), u32> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        index.insert((p.x.to_bits(), p.y.to_bits()), i as u32);
    }
    let id = |p: Point2| index[&(p.x.to_bits(), p.y.to_bits())];
    let mut constraints: Vec<Constraint> = pieces
        .iter()
        .map(|&(p, q)| Constraint::new(id(p), id(q)))
        .collect();
    constraints.sort_unstable();
    constraints.dedup();
    let tri = triangulate(&points, &constraints).map_err(|_| SkeletonError::Triangulation)?;
    let at = tri.points();
    let fixed: std::collections::BTreeSet<Constraint> = tri.constraints().iter().copied().collect();
    // Triangles inside the region.
    let mut inside: Vec<[u32; 3]> = Vec::new();
    for t in tri.triangles().chunks_exact(3) {
        let (a, b, c) = (at[t[0] as usize], at[t[1] as usize], at[t[2] as usize]);
        let centroid = Point2::new((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0);
        if contains(region, centroid)? {
            inside.push([t[0], t[1], t[2]]);
        }
    }
    // Voronoi dual: a node per inside triangle at its circumcentre (the
    // centre of the empty circle through three boundary samples), joined
    // to the triangles it shares an unconstrained edge with.
    let mut nodes: Vec<Point2> = Vec::with_capacity(inside.len());
    let mut by_edge: HashMap<(u32, u32), usize> = HashMap::new();
    let mut links: Vec<(usize, usize)> = Vec::new();
    for (k, t) in inside.iter().enumerate() {
        let (a, b, c) = (at[t[0] as usize], at[t[1] as usize], at[t[2] as usize]);
        nodes.push(circumcentre(a, b, c));
        for e in 0..3 {
            let (u, v) = (t[e], t[(e + 1) % 3]);
            if fixed.contains(&Constraint::new(u, v)) {
                continue;
            }
            match by_edge.entry((u.min(v), u.max(v))) {
                std::collections::hash_map::Entry::Occupied(o) => links.push((*o.get(), k)),
                std::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert(k);
                }
            }
        }
    }
    let segments: Vec<(Point2, Point2)> = walls.iter().map(|&(_, a, b)| (a, b)).collect();
    let clearance: Vec<(f64, f64)> = nodes.iter().map(|&p| clearance_of(p, &segments)).collect();
    // Adjacency, then prune spurs.
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for &(a, b) in &links {
        if a != b && !adj[a].contains(&b) {
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    // Keep the nodes the walls around them spread wide enough about: the
    // feet on the walls nearest a node (within its clearance, and a slack
    // for the node being off the true axis) at least `prune` times its
    // clearance apart.
    let alive: Vec<bool> = (0..nodes.len())
        .map(|i| {
            let c = clearance[i].1;
            // A chordal node is off the true axis by a fraction of its
            // clearance (for a spacing under an eighth of the narrowest
            // width); a quarter allows for it.
            spread(nodes[i], c, 0.25 * c, &segments) >= prune * c
        })
        .collect();
    for (i, list) in adj.iter_mut().enumerate() {
        if !alive[i] {
            list.clear();
        } else {
            list.retain(|&j| alive[j]);
        }
    }
    // What the filter leaves of a spur is short: drop leaf branches
    // shorter than the clearance where they join a junction.
    let mut alive = alive;
    loop {
        let mut changed = false;
        for end in 0..nodes.len() {
            if !alive[end] || adj[end].len() != 1 {
                continue;
            }
            let mut path = vec![end];
            let mut length = 0.0;
            let (mut prev, mut here) = (end, adj[end][0]);
            while adj[here].len() == 2 {
                length += (nodes[here] - nodes[prev]).length();
                path.push(here);
                let next = if adj[here][0] == prev {
                    adj[here][1]
                } else {
                    adj[here][0]
                };
                prev = here;
                here = next;
            }
            length += (nodes[here] - nodes[prev]).length();
            if adj[here].len() >= 3 && length < clearance[here].1 {
                for &n in &path {
                    alive[n] = false;
                    for o in std::mem::take(&mut adj[n]) {
                        adj[o].retain(|&x| x != n);
                    }
                }
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    // Stray single nodes go, unless nothing else is left.
    if alive.iter().zip(&adj).any(|(&a, l)| a && !l.is_empty()) {
        for i in 0..nodes.len() {
            if adj[i].is_empty() {
                alive[i] = false;
            }
        }
    }
    // Renumber the survivors.
    let mut renumber = vec![usize::MAX; nodes.len()];
    let mut out_nodes = Vec::new();
    for i in 0..nodes.len() {
        if alive[i] && inside_exactly(region, nodes[i])? {
            renumber[i] = out_nodes.len();
            out_nodes.push(i);
        }
    }
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for (i, list) in adj.iter().enumerate() {
        for &j in list {
            let (a, b) = (renumber[i], renumber[j]);
            if i < j && a != usize::MAX && b != usize::MAX {
                edges.push((a.min(b), a.max(b)));
            }
        }
    }
    edges.sort_unstable();
    edges.dedup();
    let mut degree = vec![0usize; out_nodes.len()];
    for &(a, b) in &edges {
        degree[a] += 1;
        degree[b] += 1;
    }
    let mut result = Vec::with_capacity(out_nodes.len());
    for (k, &i) in out_nodes.iter().enumerate() {
        let kind = match degree[k] {
            0 => NodeKind::Isolated,
            1 => NodeKind::End,
            2 => NodeKind::Path,
            _ => NodeKind::Junction,
        };
        let ahead = if kind == NodeKind::End {
            // The path's direction into the end: back along it until the
            // nodes are half a clearance apart (cocircular samples put
            // neighbouring nodes on one point).
            let here = nodes[i];
            let (mut prev, mut at_node) = (usize::MAX, k);
            let mut from = None;
            for _ in 0..out_nodes.len() {
                let next = edges.iter().find_map(|&(a, b)| {
                    let o = if a == at_node {
                        b
                    } else if b == at_node {
                        a
                    } else {
                        return None;
                    };
                    (o != prev).then_some(o)
                });
                let Some(next) = next else { break };
                let p = nodes[out_nodes[next]];
                if (p - here).length() >= 0.5 * clearance[i].0 {
                    from = Some(p);
                    break;
                }
                prev = at_node;
                at_node = next;
            }
            match from {
                Some(from) => wall_ahead(from, here, &walls)?,
                None => None,
            }
        } else {
            None
        };
        result.push(SkeletonNode {
            point: nodes[i],
            kind,
            clearance: clearance[i],
            ahead,
        });
    }
    Ok(Skeleton {
        nodes: result,
        edges,
        spacing,
    })
}

/// The centre of the circle through three points (rounded).
fn circumcentre(a: Point2, b: Point2, c: Point2) -> Point2 {
    let (u, v) = (b - a, c - a);
    let d = 2.0 * u.perp_dot(v);
    let (uu, vv) = (u.dot(u), v.dot(v));
    a + Point2::new(v.y * uu - u.y * vv, u.x * vv - v.x * uu) / d
}

/// In some polygon of the region, closed, decided exactly.
fn inside_exactly(region: &[Polygon], p: Point2) -> Result<bool, RouteError> {
    for polygon in region {
        if crate::map::in_polygon(polygon, p)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The first wall the ray from `from` through `to`, beyond `to`, crosses,
/// by exact tests against a far point along it.
fn wall_ahead(
    from: Point2,
    to: Point2,
    walls: &[(Wall, Point2, Point2)],
) -> Result<Option<Wall>, RouteError> {
    let d = to - from;
    if d.length() == 0.0 {
        return Ok(None);
    }
    let reach = walls
        .iter()
        .map(|&(_, a, b)| (a - to).length().max((b - to).length()))
        .fold(0.0, f64::max);
    let far = to + d * (4.0 * reach / d.length() + 1.0);
    let mut best: Option<(f64, Wall)> = None;
    for &(wall, a, b) in walls {
        // Meets the segment from `to` to `far`, crossing or touching.
        let hit = crosses(to, far, a, b)?
            || (side(to, far, a)? == Sign::Zero && within(to, far, a))
            || (side(to, far, b)? == Sign::Zero && within(to, far, b));
        if !hit {
            continue;
        }
        // Distance along the ray to the wall's line (rounded; only the
        // order of nearby walls could be affected, and then either is a
        // wall the path runs into).
        let e = b - a;
        let den = d.perp_dot(e);
        let t = if den == 0.0 {
            (a - to).length()
        } else {
            (a - to).perp_dot(e) / den
        };
        if best.is_none_or(|(bt, _)| t < bt) {
            best = Some((t, wall));
        }
    }
    Ok(best.map(|(_, w)| w))
}

/// How far apart the feet of `p` on its nearest walls lie: the walls
/// within `clearance + slack` of it, each at its nearest point.
fn spread(p: Point2, clearance: f64, slack: f64, segments: &[(Point2, Point2)]) -> f64 {
    let feet: Vec<Point2> = segments
        .iter()
        .filter_map(|&(a, b)| {
            let e = b - a;
            let t = ((p - a).dot(e) / e.dot(e)).clamp(0.0, 1.0);
            let foot = a + e * t;
            ((p - foot).length() <= clearance + slack).then_some(foot)
        })
        .collect();
    let mut widest = 0.0f64;
    for (i, a) in feet.iter().enumerate() {
        for b in &feet[i + 1..] {
            widest = widest.max((*a - *b).length());
        }
    }
    widest
}

/// Proven bounds on the distance from `p` to the nearest of `segments`.
fn clearance_of(p: Point2, segments: &[(Point2, Point2)]) -> (f64, f64) {
    let mut low = f64::INFINITY;
    let mut high = f64::INFINITY;
    for &(a, b) in segments {
        let (lo, hi) = segment_distance(p, a, b);
        low = low.min(lo);
        high = high.min(hi);
    }
    (low, high)
}

/// Bounds on the distance from `p` to the segment `a b`: the squared
/// distance in outward-rounded intervals, square-rooted outward.
fn segment_distance(p: Point2, a: Point2, b: Point2) -> (f64, f64) {
    let iv = |x: f64| Iv { lo: x, hi: x };
    let (ex, ey) = (iv(b.x).sub(iv(a.x)), iv(b.y).sub(iv(a.y)));
    let (wx, wy) = (iv(p.x).sub(iv(a.x)), iv(p.y).sub(iv(a.y)));
    let dot = ex.mul(wx).add(ey.mul(wy));
    let len2 = ex.mul(ex).add(ey.mul(ey));
    let w2 = wx.mul(wx).add(wy.mul(wy));
    let (vx, vy) = (iv(p.x).sub(iv(b.x)), iv(p.y).sub(iv(b.y)));
    let v2 = vx.mul(vx).add(vy.mul(vy));
    // Which part of the segment is nearest, when decided; otherwise the
    // hull of the candidates.
    let perp = || {
        let c = ex.mul(wy).sub(ey.mul(wx));
        c.mul(c).div(len2)
    };
    let d2 = if dot.hi <= 0.0 {
        w2
    } else if dot.lo >= len2.hi {
        v2
    } else if dot.lo > 0.0 && dot.hi < len2.lo {
        perp()
    } else {
        let q = perp();
        Iv {
            lo: q.lo.min(w2.lo).min(v2.lo),
            hi: q.hi.max(w2.hi).max(v2.hi),
        }
    };
    let lo = d2.lo.max(0.0).sqrt().next_down().max(0.0);
    let hi = d2.hi.max(0.0).sqrt().next_up();
    (lo, hi)
}

/// An outward-rounded interval.
#[derive(Debug, Clone, Copy)]
struct Iv {
    lo: f64,
    hi: f64,
}

impl Iv {
    fn outward(lo: f64, hi: f64) -> Self {
        Self {
            lo: lo.next_down(),
            hi: hi.next_up(),
        }
    }

    fn add(self, o: Self) -> Self {
        Self::outward(self.lo + o.lo, self.hi + o.hi)
    }

    fn sub(self, o: Self) -> Self {
        Self::outward(self.lo - o.hi, self.hi - o.lo)
    }

    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self::outward(
            p.iter().copied().fold(f64::INFINITY, f64::min),
            p.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    }

    /// Quotient by a positive interval.
    fn div(self, o: Self) -> Self {
        if o.lo <= 0.0 {
            return Self {
                lo: 0.0,
                hi: f64::INFINITY,
            };
        }
        let q = [
            self.lo / o.lo,
            self.lo / o.hi,
            self.hi / o.lo,
            self.hi / o.hi,
        ];
        Self::outward(
            q.iter().copied().fold(f64::INFINITY, f64::min),
            q.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    }
}
