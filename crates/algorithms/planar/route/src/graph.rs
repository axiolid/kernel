//! The visibility graph, with each vertex split into its free sectors
//! (#189).
//!
//! # Why sectors
//!
//! Where obstacles meet at a point -- a barrier's foot on a wall, two holes
//! touching at a corner, a barrier bent at a vertex -- the free space
//! around that point falls into separate wedges. A route may come in
//! through one wedge and must leave through the same one: passing from one
//! to another would squeeze through a gap of zero width.
//!
//! So every graph vertex carries the obstacle rays leaving it, sorted by
//! angle exactly, and the wedges between consecutive rays -- its sectors --
//! each free or not by which side of each ray the region lies on. A route
//! is a path over (vertex, sector) states.
//!
//! An edge is travelled on one side: its left and right variants attach,
//! at each end, to the sector just beside it on that side, and where the
//! edge passes through another vertex, that side must be clear there -- no
//! obstacle ray strictly on it, and free space. An edge along a wall has
//! one usable side, an edge along a barrier two, kept apart.

use axiolid_contracts::Sign;
use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};

use crate::{ring_edges, side, visible, within, RouteError};

/// A ray from a vertex along an obstacle, with whether the region lies
/// just counter-clockwise of it and just clockwise, and whether a ring
/// edge (rather than only barriers) runs along it.
#[derive(Debug, Clone, Copy)]
struct Ray {
    to: Point2,
    ccw_free: bool,
    cw_free: bool,
    ring: bool,
}

impl Ray {
    /// Merge a ray in the same direction. Ring edges say where the region
    /// is, and two of them (polygons sharing an edge) put it on both
    /// sides; a barrier adds no region, and alone leaves both sides free.
    fn merge(&mut self, other: Ray) {
        match (self.ring, other.ring) {
            (true, true) => {
                self.ccw_free |= other.ccw_free;
                self.cw_free |= other.cw_free;
            }
            (false, true) => *self = other,
            _ => {}
        }
    }
}

/// The rays around one vertex, counter-clockwise from the positive x
/// direction, same-direction rays merged. Sector `i` runs from ray `i`
/// counter-clockwise to ray `i + 1`; with no rays there is one sector,
/// the whole disc.
#[derive(Debug, Clone)]
pub(crate) struct Star {
    at: Point2,
    rays: Vec<Ray>,
}

/// Where a direction lies among the rays.
enum Place {
    On(usize),
    In(usize),
}

/// Which side of travel an edge is taken on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Side {
    Left,
    Right,
}

impl Star {
    pub(crate) fn sectors(&self) -> usize {
        self.rays.len().max(1)
    }

    pub(crate) fn free(&self, sector: usize) -> bool {
        match self.rays.len() {
            0 => true,
            n => self.rays[sector].ccw_free && self.rays[(sector + 1) % n].cw_free,
        }
    }

    /// Half-plane of the direction to `p`: the upper half with the positive
    /// x-axis first, exactly.
    fn half(&self, p: Point2) -> u8 {
        u8::from(!(p.y > self.at.y || (p.y == self.at.y && p.x > self.at.x)))
    }

    /// Counter-clockwise order of the directions to `p` and `q`.
    fn order(&self, p: Point2, q: Point2) -> Result<core::cmp::Ordering, RouteError> {
        use core::cmp::Ordering;
        let (hp, hq) = (self.half(p), self.half(q));
        if hp != hq {
            return Ok(hp.cmp(&hq));
        }
        Ok(match side(self.at, p, q)? {
            Sign::Positive => Ordering::Less,
            Sign::Negative => Ordering::Greater,
            _ => Ordering::Equal,
        })
    }

    fn place(&self, toward: Point2) -> Result<Place, RouteError> {
        use core::cmp::Ordering;
        let n = self.rays.len();
        if n == 0 {
            return Ok(Place::In(0));
        }
        let mut before = None;
        for (i, ray) in self.rays.iter().enumerate() {
            match self.order(ray.to, toward)? {
                Ordering::Equal => return Ok(Place::On(i)),
                Ordering::Less => before = Some(i),
                Ordering::Greater => break,
            }
        }
        // Before the first ray means after the last, round the circle.
        Ok(Place::In(before.unwrap_or(n - 1)))
    }

    /// The sector just counter-clockwise of the direction to `toward`.
    pub(crate) fn ccw_of(&self, toward: Point2) -> Result<usize, RouteError> {
        Ok(match self.place(toward)? {
            Place::On(i) | Place::In(i) => i,
        })
    }

    /// The sector just clockwise of the direction to `toward`.
    pub(crate) fn cw_of(&self, toward: Point2) -> Result<usize, RouteError> {
        let n = self.sectors();
        Ok(match self.place(toward)? {
            Place::On(i) => (i + n - 1) % n,
            Place::In(i) => i,
        })
    }

    /// Whether a route along the line from `a` through this vertex to `b`
    /// may pass it on `side`: no obstacle ray strictly on that side, and
    /// free space there.
    fn passes(&self, a: Point2, b: Point2, on: Side) -> Result<bool, RouteError> {
        let blocking = match on {
            Side::Left => Sign::Positive,
            Side::Right => Sign::Negative,
        };
        for ray in &self.rays {
            if side(a, b, ray.to)? == blocking {
                return Ok(false);
            }
        }
        let sector = match on {
            Side::Left => self.ccw_of(b)?,
            Side::Right => self.cw_of(b)?,
        };
        Ok(self.free(sector))
    }
}

/// Whether the region lies to the left of a ring's edges in ring order.
fn region_left(ring: &Ring, hole: bool) -> Result<bool, RouteError> {
    let p = &ring.points;
    let n = p.len();
    let low = (0..n)
        .min_by(|&i, &j| p[i].x.total_cmp(&p[j].x).then(p[i].y.total_cmp(&p[j].y)))
        .unwrap_or(0);
    let ccw = match side(p[(low + n - 1) % n], p[low], p[(low + 1) % n])? {
        Sign::Positive => true,
        Sign::Negative => false,
        // Collinear at the extreme vertex: fall back on the area's sign.
        _ => {
            let twice: f64 = (0..n)
                .map(|i| p[i].x * p[(i + 1) % n].y - p[(i + 1) % n].x * p[i].y)
                .sum();
            twice > 0.0
        }
    };
    Ok(ccw != hole)
}

/// The star of every vertex: rays along each ring edge and barrier
/// segment through it.
pub(crate) fn stars(
    nodes: &[Point2],
    region: &[Polygon],
    barriers: &[Vec<Point2>],
) -> Result<Vec<Star>, RouteError> {
    // Every obstacle segment with the region's side: (from, to, left free,
    // right free) in the segment's own direction.
    let mut segments: Vec<(Point2, Point2, bool, bool, bool)> = Vec::new();
    for polygon in region {
        for (ring, hole) in
            core::iter::once((&polygon.outer, false)).chain(polygon.holes.iter().map(|h| (h, true)))
        {
            let left = region_left(ring, hole)?;
            for (p, q) in ring_edges(ring) {
                segments.push((p, q, left, !left, true));
            }
        }
    }
    for barrier in barriers {
        for pair in barrier.windows(2) {
            segments.push((pair[0], pair[1], true, true, false));
        }
    }
    let mut out = Vec::with_capacity(nodes.len());
    for &v in nodes {
        let mut star = Star {
            at: v,
            rays: Vec::new(),
        };
        let mut rays: Vec<Ray> = Vec::new();
        for &(p, q, left, right, ring) in &segments {
            if p == q {
                continue;
            }
            let on = v == p || v == q || (side(p, q, v)? == Sign::Zero && within(p, q, v));
            if !on {
                continue;
            }
            // Along the segment's direction the region's left is the ray's
            // counter-clockwise side; against it, the clockwise.
            if v != q {
                rays.push(Ray {
                    to: q,
                    ccw_free: left,
                    cw_free: right,
                    ring,
                });
            }
            if v != p {
                rays.push(Ray {
                    to: p,
                    ccw_free: right,
                    cw_free: left,
                    ring,
                });
            }
        }
        // Insertion sort by angle: the comparisons are fallible.
        for ray in rays {
            let mut at = star.rays.len();
            let mut merged = false;
            for (i, other) in star.rays.iter_mut().enumerate() {
                match star_order(v, other.to, ray.to)? {
                    core::cmp::Ordering::Equal => {
                        other.merge(ray);
                        merged = true;
                        break;
                    }
                    core::cmp::Ordering::Greater => {
                        at = i;
                        break;
                    }
                    core::cmp::Ordering::Less => {}
                }
            }
            if !merged {
                star.rays.insert(at, ray);
            }
        }
        out.push(star);
    }
    Ok(out)
}

fn star_order(at: Point2, p: Point2, q: Point2) -> Result<core::cmp::Ordering, RouteError> {
    Star {
        at,
        rays: Vec::new(),
    }
    .order(p, q)
}

/// The sides on which the segment from `a` to `b` may be travelled. Both
/// false when it leaves the region or crosses an obstacle; otherwise each
/// side must be clear at every vertex the segment passes through.
pub(crate) fn sides(
    a: Point2,
    b: Point2,
    region: &[Polygon],
    obstacles: &[(Point2, Point2)],
    nodes: &[Point2],
    stars: &[Star],
) -> Result<[bool; 2], RouteError> {
    if !visible(a, b, region, obstacles)? {
        return Ok([false, false]);
    }
    let mut ok = [true, true];
    for (v, star) in nodes.iter().zip(stars) {
        if *v == a || *v == b || star.rays.is_empty() {
            continue;
        }
        if side(a, b, *v)? != Sign::Zero || !within(a, b, *v) {
            continue;
        }
        for (k, on) in [Side::Left, Side::Right].into_iter().enumerate() {
            if ok[k] && !star.passes(a, b, on)? {
                ok[k] = false;
            }
        }
        if ok == [false, false] {
            break;
        }
    }
    Ok(ok)
}

/// The visibility graph over (vertex, sector) states.
#[derive(Debug, Clone)]
pub(crate) struct Graph {
    pub(crate) stars: Vec<Star>,
    /// First state of each vertex; its sectors follow.
    pub(crate) offset: Vec<usize>,
    pub(crate) adjacency: Vec<Vec<(usize, f64)>>,
}

impl Graph {
    pub(crate) fn build(
        nodes: &[Point2],
        region: &[Polygon],
        barriers: &[Vec<Point2>],
        obstacles: &[(Point2, Point2)],
    ) -> Result<Self, RouteError> {
        let stars = stars(nodes, region, barriers)?;
        let mut offset = Vec::with_capacity(nodes.len() + 1);
        let mut total = 0;
        for star in &stars {
            offset.push(total);
            total += star.sectors();
        }
        offset.push(total);
        let mut adjacency = vec![Vec::new(); total];
        for i in 0..nodes.len() {
            for j in i + 1..nodes.len() {
                let ok = sides(nodes[i], nodes[j], region, obstacles, nodes, &stars)?;
                let length = (nodes[i] - nodes[j]).length();
                let mut linked: Vec<(usize, usize)> = Vec::new();
                for (k, on) in [Side::Left, Side::Right].into_iter().enumerate() {
                    if !ok[k] {
                        continue;
                    }
                    let (si, sj) = match on {
                        Side::Left => (stars[i].ccw_of(nodes[j])?, stars[j].cw_of(nodes[i])?),
                        Side::Right => (stars[i].cw_of(nodes[j])?, stars[j].ccw_of(nodes[i])?),
                    };
                    if !stars[i].free(si) || !stars[j].free(sj) {
                        continue;
                    }
                    let (u, w) = (offset[i] + si, offset[j] + sj);
                    if !linked.contains(&(u, w)) {
                        linked.push((u, w));
                        adjacency[u].push((w, length));
                        adjacency[w].push((u, length));
                    }
                }
            }
        }
        Ok(Self {
            stars,
            offset,
            adjacency,
        })
    }

    /// The vertex of a state.
    pub(crate) fn node(&self, state: usize) -> usize {
        self.offset.partition_point(|&o| o <= state) - 1
    }

    /// The free states of a vertex.
    pub(crate) fn states(&self, node: usize) -> impl Iterator<Item = usize> + '_ {
        (0..self.stars[node].sectors())
            .filter(move |&s| self.stars[node].free(s))
            .map(move |s| self.offset[node] + s)
    }

    /// The state an edge from an outside point `from` arrives in at
    /// `node`, on `on` side of travel.
    pub(crate) fn arrival(
        &self,
        node: usize,
        from: Point2,
        on: Side,
    ) -> Result<Option<usize>, RouteError> {
        let star = &self.stars[node];
        let sector = match on {
            Side::Left => star.cw_of(from)?,
            Side::Right => star.ccw_of(from)?,
        };
        Ok(star.free(sector).then_some(self.offset[node] + sector))
    }
}

/// Dijkstra from several sources over states, returning distance and
/// predecessor. Ties break on the lowest state index.
pub(crate) fn dijkstra(
    adjacency: &[Vec<(usize, f64)>],
    sources: &[usize],
    stop: impl Fn(usize) -> bool,
) -> (Vec<f64>, Vec<usize>, Option<usize>) {
    let count = adjacency.len();
    let mut distance = vec![f64::INFINITY; count];
    let mut previous = vec![usize::MAX; count];
    let mut settled = vec![false; count];
    for &s in sources {
        distance[s] = 0.0;
    }
    for _ in 0..count {
        let mut current = None;
        for index in 0..count {
            if settled[index] || distance[index].is_infinite() {
                continue;
            }
            if current.is_none_or(|best: usize| distance[index] < distance[best]) {
                current = Some(index);
            }
        }
        let Some(current) = current else { break };
        if stop(current) {
            return (distance, previous, Some(current));
        }
        settled[current] = true;
        for &(neighbour, weight) in &adjacency[current] {
            let candidate = distance[current] + weight;
            if candidate < distance[neighbour] {
                distance[neighbour] = candidate;
                previous[neighbour] = current;
            }
        }
    }
    (distance, previous, None)
}
