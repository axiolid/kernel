//! Many-target shortest paths and the farthest point of a subregion (#186).
//!
//! # The map
//!
//! [`DistanceMap`] runs one Dijkstra from every target at once over the
//! visibility graph, so each graph vertex knows its distance to the nearest
//! target and the next vertex on the way. A query point then needs only the
//! vertices it sees: its distance is the least `|x - v| + d(v)` over them.
//! Visibility is decided exactly, as for [`crate::shortest_path`]; lengths
//! are sums of square roots in binary64.
//!
//! # The farthest point, bracketed
//!
//! The distance `D` to the nearest target is 1-Lipschitz along any segment
//! that stays in the free space: `D(y) <= D(a) + |y - a|`. So the free space
//! is cut into triangles that each lie inside it -- a constrained
//! triangulation whose constraints are every wall and barrier -- and on a
//! triangle `T` with a point `a` of known distance,
//!
//! ```text
//! max over T of D  <=  D(a) + (greatest distance from a to a corner of T).
//! ```
//!
//! Every distance evaluated at a point of the subregion is a lower bound on
//! the maximum. Triangles are split at their edge midpoints, the one with
//! the greatest upper bound first, and dropped once their upper bound falls
//! below the best lower one, until the two are within the tolerance or the
//! cell budget runs out. The result is an interval that contains the true
//! maximum: widened by the rounding of the lengths, and of the midpoints,
//! which leave gaps between cells no wider than a few ulps.
//!
//! Anchors -- the points whose distance bounds a cell -- are only ever
//! points proven inside the cell's original triangle, exactly: a rounded
//! midpoint on a wall may fall just outside the free space, and its
//! distance says nothing about the triangle. Nor may an anchor lie on a
//! barrier, unless strictly inside the triangle (which it then cannot be):
//! a point on a barrier is seen from both sides, so its distance is the
//! nearer side's, and bounds nothing on the farther.

use core::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};

use axiolid_contracts::Sign;
use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_triangulate::{triangulate, Constraint};

use crate::{
    contains, crosses, dedup_points, obstacle_segments, ring_edges, side, validate_region, visible,
    within, Route, RouteError, Unreachable, MAX_VERTICES,
};

/// Cells [`farthest_point`] refines at most.
pub const MAX_CELLS: usize = 20_000;

/// Why no distance map was built.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum MapError {
    /// The region, a barrier or a target is malformed, or the input is
    /// over budget ([`RouteError::TooManyVertices`], whose lower bound is
    /// zero here: there is no single pair of endpoints).
    Route(RouteError),
    /// No targets were given.
    NoTargets,
    /// A target lies outside the free-space region.
    TargetOutside {
        /// Its index in the targets given.
        index: usize,
    },
}

impl From<RouteError> for MapError {
    fn from(error: RouteError) -> Self {
        Self::Route(error)
    }
}

/// Shortest-path distances from every point of a region to the nearest of
/// several targets.
#[derive(Debug, Clone)]
pub struct DistanceMap {
    region: Vec<Polygon>,
    /// Barrier segments: zero-width, with free space on both sides.
    walls: Vec<(Point2, Point2)>,
    obstacles: Vec<(Point2, Point2)>,
    nodes: Vec<Point2>,
    /// Distance from each node to its nearest target; infinite if none.
    distance: Vec<f64>,
    /// The next node towards that target; `usize::MAX` at a target.
    next: Vec<usize>,
    /// Which target, as an index into the targets given.
    target: Vec<usize>,
    targets: usize,
}

/// The nearest target from a point, and the route there.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Reach {
    /// Index of the nearest target in the targets given. Ties go to the
    /// route through the lowest graph vertex, deterministically.
    pub target: usize,
    /// The route, from the query point to the target.
    pub route: Route,
}

/// A distance map from `targets` over `region`, avoiding `barriers`.
///
/// # Errors
///
/// [`MapError`] for malformed input, no targets, a target outside the
/// region, or more than [`MAX_VERTICES`] graph vertices.
pub fn distance_map(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    targets: &[Point2],
) -> Result<DistanceMap, MapError> {
    distance_map_within(region, barriers, targets, MAX_VERTICES)
}

/// [`distance_map`] with a caller-chosen vertex budget.
///
/// # Errors
///
/// As [`distance_map`], with `budget` for [`MAX_VERTICES`].
pub fn distance_map_within(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    targets: &[Point2],
    budget: usize,
) -> Result<DistanceMap, MapError> {
    validate_region(region, barriers)?;
    if targets.is_empty() {
        return Err(MapError::NoTargets);
    }
    if !targets.iter().all(|t| t.is_finite()) {
        return Err(RouteError::NonFinitePoint.into());
    }
    for (index, t) in targets.iter().enumerate() {
        if !contains(region, *t)? {
            return Err(MapError::TargetOutside { index });
        }
    }
    let mut nodes = targets.to_vec();
    for polygon in region {
        for ring in core::iter::once(&polygon.outer).chain(polygon.holes.iter()) {
            nodes.extend(ring.points.iter().copied());
        }
    }
    for barrier in barriers {
        nodes.extend(barrier.iter().copied());
    }
    dedup_points(&mut nodes);
    if nodes.len() > budget {
        return Err(RouteError::TooManyVertices {
            supplied: nodes.len(),
            budget,
            lower_bound: 0.0,
        }
        .into());
    }
    let obstacles = obstacle_segments(region, barriers);
    let mut adjacency = vec![Vec::new(); nodes.len()];
    for i in 0..nodes.len() {
        for j in i + 1..nodes.len() {
            if visible(nodes[i], nodes[j], region, &obstacles)? {
                let length = (nodes[i] - nodes[j]).length();
                adjacency[i].push((j, length));
                adjacency[j].push((i, length));
            }
        }
    }
    // Every target at distance zero, then one Dijkstra; ties go to the
    // lowest index, as in `shortest_path`.
    let count = nodes.len();
    let mut distance = vec![f64::INFINITY; count];
    let mut next = vec![usize::MAX; count];
    let mut target = vec![usize::MAX; count];
    for (i, node) in nodes.iter().enumerate() {
        if let Some(t) = targets.iter().position(|t| t == node) {
            distance[i] = 0.0;
            target[i] = t;
        }
    }
    let mut settled = vec![false; count];
    for _ in 0..count {
        let mut current = None;
        for i in 0..count {
            if settled[i] || distance[i].is_infinite() {
                continue;
            }
            if current.is_none_or(|best: usize| distance[i] < distance[best]) {
                current = Some(i);
            }
        }
        let Some(current) = current else { break };
        settled[current] = true;
        for &(neighbour, weight) in &adjacency[current] {
            let candidate = distance[current] + weight;
            if candidate < distance[neighbour] {
                distance[neighbour] = candidate;
                next[neighbour] = current;
                target[neighbour] = target[current];
            }
        }
    }
    Ok(DistanceMap {
        region: region.to_vec(),
        walls: barriers
            .iter()
            .flat_map(|b| b.windows(2).map(|w| (w[0], w[1])))
            .collect(),
        obstacles,
        nodes,
        distance,
        next,
        target,
        targets: targets.len(),
    })
}

impl DistanceMap {
    /// Vertices of the visibility graph.
    #[must_use]
    pub fn graph_vertices(&self) -> usize {
        self.nodes.len()
    }

    /// How many targets the map was built from.
    #[must_use]
    pub fn targets(&self) -> usize {
        self.targets
    }

    /// The nearest target from `point` and the route there.
    ///
    /// `Ok(Err(StartOutside))` for a point outside the region,
    /// `Ok(Err(DisconnectedComponents))` for one no target can reach.
    ///
    /// # Errors
    ///
    /// [`RouteError`] for a non-finite point or an undecidable predicate.
    pub fn nearest(&self, point: Point2) -> Result<Result<Reach, Unreachable>, RouteError> {
        if !point.is_finite() {
            return Err(RouteError::NonFinitePoint);
        }
        if !contains(&self.region, point)? {
            return Ok(Err(Unreachable::StartOutside));
        }
        let Some((length, first)) = self.via(point)? else {
            return Ok(Err(Unreachable::DisconnectedComponents));
        };
        let mut polyline = vec![point];
        let mut node = first;
        loop {
            if self.nodes[node] != point {
                polyline.push(self.nodes[node]);
            }
            if self.next[node] == usize::MAX {
                break;
            }
            node = self.next[node];
        }
        Ok(Ok(Reach {
            target: self.target[first],
            route: Route {
                polyline,
                length,
                graph_vertices: self.nodes.len(),
            },
        }))
    }

    /// The distance from `point`, known to be inside the region, and the
    /// first graph vertex on the way; `None` when no target is reachable.
    fn via(&self, point: Point2) -> Result<Option<(f64, usize)>, RouteError> {
        let mut best: Option<(f64, usize)> = None;
        for (i, node) in self.nodes.iter().enumerate() {
            if self.distance[i].is_infinite() {
                continue;
            }
            let leg = (*node - point).length();
            let length = leg + self.distance[i];
            if best.is_some_and(|(b, _)| length >= b) {
                continue;
            }
            if *node == point || visible(point, *node, &self.region, &self.obstacles)? {
                best = Some((length, i));
            }
        }
        Ok(best)
    }

    /// Distance to the nearest target, for a point inside the region.
    fn at(&self, point: Point2) -> Result<Option<f64>, RouteError> {
        Ok(self.via(point)?.map(|(length, _)| length))
    }

    /// The largest length of a path in the graph, in edges: bounds how many
    /// roundings a distance carries.
    fn hops(&self) -> usize {
        let mut most = 0;
        for start in 0..self.nodes.len() {
            let mut node = start;
            let mut hops = 0;
            while self.next[node] != usize::MAX && hops <= self.nodes.len() {
                node = self.next[node];
                hops += 1;
            }
            most = most.max(hops);
        }
        most
    }
}

/// A closed interval of lengths.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LengthInterval {
    /// No greater than the true value.
    pub lower: f64,
    /// No less than the true value.
    pub upper: f64,
}

/// The greatest distance from a subregion to the nearest target.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct Farthest {
    /// Contains the true maximum over the subregion's points in the free
    /// space.
    pub distance: LengthInterval,
    /// A point of the subregion, in the free space, at least
    /// `distance.lower` from every target; `None` only when the budget ran
    /// out before any point of the subregion was sampled.
    pub witness: Option<Point2>,
    /// Whether the interval is no wider than the tolerance asked for. When
    /// the cell budget runs out first it is still sound, only wider.
    pub converged: bool,
    /// Cells examined.
    pub cells: usize,
}

/// Why no bracket was produced.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum FarthestError {
    /// The subregion is malformed, or a predicate was undecidable.
    Route(RouteError),
    /// The tolerance is negative or not finite.
    InvalidTolerance,
    /// Two walls or barriers cross other than at a shared vertex. The free
    /// space is triangulated with them as constraints, and crossings would
    /// need rounded vertices; not supported yet.
    CrossingObstacles,
    /// The free space could not be triangulated, or a triangle of it is
    /// too thin to hold a point strictly inside.
    Triangulation,
    /// Part of the subregion lies in free space that no target reaches, so
    /// its farthest distance is infinite. The triangle lies in the free
    /// space and meets the subregion; no target reaches any of it.
    Unreachable {
        /// The triangle's corners.
        triangle: [Point2; 3],
    },
    /// The subregion does not meet the free space.
    Empty,
}

impl From<RouteError> for FarthestError {
    fn from(error: RouteError) -> Self {
        Self::Route(error)
    }
}

/// The greatest distance to the nearest target over the points of
/// `subregion` in the map's free space, bracketed to within `tolerance`.
///
/// # Errors
///
/// [`FarthestError`], among them [`FarthestError::Unreachable`] with a
/// triangle as evidence when part of the subregion reaches no target.
pub fn farthest_point(
    map: &DistanceMap,
    subregion: &Polygon,
    tolerance: f64,
) -> Result<Farthest, FarthestError> {
    farthest_point_within(map, subregion, tolerance, MAX_CELLS)
}

/// A cell: a triangle inside the free-space triangle `root`, with the best
/// anchor known for it and the upper bound that gives.
#[derive(Debug, Clone, Copy)]
struct Cell {
    corners: [Point2; 3],
    root: usize,
    anchor: (Point2, f64),
    upper: f64,
    depth: u32,
    order: usize,
}

impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Cell {}

impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Cell {
    /// Greatest upper bound first; the earlier cell on a tie.
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper
            .total_cmp(&other.upper)
            .then(other.order.cmp(&self.order))
    }
}

/// [`farthest_point`] with a caller-chosen cell budget.
///
/// # Errors
///
/// As [`farthest_point`].
pub fn farthest_point_within(
    map: &DistanceMap,
    subregion: &Polygon,
    tolerance: f64,
    max_cells: usize,
) -> Result<Farthest, FarthestError> {
    if !(tolerance.is_finite() && tolerance >= 0.0) {
        return Err(FarthestError::InvalidTolerance);
    }
    validate_region(core::slice::from_ref(subregion), &[])?;
    let roots = free_triangles(map)?;
    let scale = map
        .nodes
        .iter()
        .chain(subregion.outer.points.iter())
        .fold(0.0f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));
    // A distance is a sum of at most `hops + 1` rounded lengths, each a
    // rounded square root of rounded squares of rounded differences.
    let relative = (map.hops() as f64 + 8.0) * 2.0 * f64::EPSILON;
    // Rounded midpoints leave gaps of about an ulp per level between cells,
    // and the radius itself is rounded.
    let slack = |depth: u32, radius: f64| {
        f64::from(depth + 4) * 4.0 * f64::EPSILON * scale + 4.0 * f64::EPSILON * radius
    };
    let mut search = Search {
        map,
        subregion,
        roots: &roots,
        evaluated: HashMap::new(),
        lower: None,
    };
    let mut heap = BinaryHeap::new();
    let mut cells = 0usize;
    for (root, corners) in roots.iter().enumerate() {
        if search.outside(corners)? {
            continue;
        }
        cells += 1;
        // A root's centroid is strictly inside it, off every barrier,
        // unless the triangle is too thin to hold a representable point.
        let Some(cell) = search.cell(*corners, root, None, 0, cells, &slack)? else {
            return Err(FarthestError::Triangulation);
        };
        heap.push(cell);
    }
    if cells == 0 {
        return Err(FarthestError::Empty);
    }
    let best = |search: &Search| search.lower.map_or(f64::NEG_INFINITY, |(d, _)| d);
    while let Some(top) = heap.peek() {
        if top.upper <= best(&search) + tolerance || cells >= max_cells {
            break;
        }
        let cell = heap.pop().expect("peeked");
        let [a, b, c] = cell.corners;
        let mid = |p: Point2, q: Point2| Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y);
        let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
        for corners in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
            if search.outside(&corners)? {
                continue;
            }
            cells += 1;
            let child = search.cell(
                corners,
                cell.root,
                Some(cell.anchor),
                cell.depth + 1,
                cells,
                &slack,
            )?;
            if let Some(child) = child {
                if child.upper > best(&search) {
                    heap.push(child);
                }
            }
        }
    }
    let Some((lower, witness)) = search.lower else {
        if heap.is_empty() {
            return Err(FarthestError::Empty);
        }
        let upper = heap.peek().map_or(0.0, |c| c.upper);
        return Ok(Farthest {
            distance: LengthInterval {
                lower: 0.0,
                upper: upper * (1.0 + relative),
            },
            witness: None,
            converged: false,
            cells,
        });
    };
    let upper = heap.peek().map_or(lower, |c| c.upper.max(lower));
    Ok(Farthest {
        distance: LengthInterval {
            lower: lower * (1.0 - relative),
            upper: upper * (1.0 + relative),
        },
        witness: Some(witness),
        converged: upper - lower <= tolerance,
        cells,
    })
}

struct Search<'a> {
    map: &'a DistanceMap,
    subregion: &'a Polygon,
    roots: &'a [[Point2; 3]],
    evaluated: HashMap<(u64, u64), Option<f64>>,
    /// The greatest distance at a point of the subregion, and the point.
    lower: Option<(f64, Point2)>,
}

impl Search<'_> {
    /// Whether a triangle misses the subregion, decided exactly: no edge
    /// of the subregion meets it, and a corner lies outside.
    fn outside(&self, t: &[Point2; 3]) -> Result<bool, RouteError> {
        for ring in core::iter::once(&self.subregion.outer).chain(self.subregion.holes.iter()) {
            for (p, q) in ring_edges(ring) {
                if meets_triangle(p, q, t)? {
                    return Ok(false);
                }
            }
        }
        Ok(!in_polygon(self.subregion, t[0])?)
    }

    /// Whether `a` may anchor a cell of `root`: inside it, and not on a
    /// barrier.
    fn admissible(&self, root: usize, a: Point2) -> Result<bool, RouteError> {
        if !in_triangle(&self.roots[root], a)? {
            return Ok(false);
        }
        for &(p, q) in &self.map.walls {
            if side(p, q, a)? == Sign::Zero && within(p, q, a) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn distance(&mut self, p: Point2) -> Result<Option<f64>, RouteError> {
        let key = (p.x.to_bits(), p.y.to_bits());
        if let Some(d) = self.evaluated.get(&key) {
            return Ok(*d);
        }
        let d = self.map.at(p)?;
        self.evaluated.insert(key, d);
        Ok(d)
    }

    /// The cell's upper bound from its best anchor: its corners and
    /// centroid where proven inside the root triangle, and the parent's.
    /// Raises the lower bound from anchors in the subregion. `None` when
    /// the cell has no anchor (never for a root).
    fn cell(
        &mut self,
        corners: [Point2; 3],
        root: usize,
        inherited: Option<(Point2, f64)>,
        depth: u32,
        order: usize,
        slack: &dyn Fn(u32, f64) -> f64,
    ) -> Result<Option<Cell>, FarthestError> {
        let centroid = Point2::new(
            (corners[0].x + corners[1].x + corners[2].x) / 3.0,
            (corners[0].y + corners[1].y + corners[2].y) / 3.0,
        );
        let radius = |a: Point2| {
            corners
                .iter()
                .map(|c| (*c - a).length())
                .fold(0.0, f64::max)
        };
        let mut best: Option<(f64, (Point2, f64))> =
            inherited.map(|(a, d)| (d + radius(a), (a, d)));
        for a in [corners[0], corners[1], corners[2], centroid] {
            if !self.admissible(root, a)? {
                continue;
            }
            let Some(d) = self.distance(a)? else {
                return Err(FarthestError::Unreachable {
                    triangle: self.roots[root],
                });
            };
            if in_polygon(self.subregion, a)? && self.lower.is_none_or(|(l, _)| d > l) {
                self.lower = Some((d, a));
            }
            let bound = d + radius(a);
            if best.is_none_or(|(b, _)| bound < b) {
                best = Some((bound, (a, d)));
            }
        }
        Ok(best.map(|(bound, anchor)| Cell {
            corners,
            root,
            anchor,
            upper: bound + slack(depth, radius(anchor.0)),
            depth,
            order,
        }))
    }
}

/// The triangles of the free space: a constrained triangulation of every
/// wall and barrier, each segment first cut at the vertices lying on it,
/// keeping the triangles inside the region.
fn free_triangles(map: &DistanceMap) -> Result<Vec<[Point2; 3]>, FarthestError> {
    let mut points: Vec<Point2> = map.obstacles.iter().flat_map(|(p, q)| [*p, *q]).collect();
    dedup_points(&mut points);
    let mut pieces: Vec<(Point2, Point2)> = Vec::new();
    for &(p, q) in &map.obstacles {
        let d = q - p;
        let mut cuts = vec![p, q];
        for &v in &points {
            if v != p && v != q && side(p, q, v)? == Sign::Zero && within(p, q, v) {
                cuts.push(v);
            }
        }
        cuts.sort_by(|u, v| (*u - p).dot(d).total_cmp(&(*v - p).dot(d)));
        for pair in cuts.windows(2) {
            if pair[0] != pair[1] {
                pieces.push((pair[0], pair[1]));
            }
        }
    }
    for (i, &(p, q)) in pieces.iter().enumerate() {
        for &(r, s) in &pieces[i + 1..] {
            if crosses(p, q, r, s)? {
                return Err(FarthestError::CrossingObstacles);
            }
        }
    }
    let index = |v: Point2| {
        u32::try_from(points.iter().position(|p| *p == v).expect("an endpoint")).unwrap_or(u32::MAX)
    };
    let mut constraints: Vec<Constraint> = pieces
        .iter()
        .map(|(p, q)| Constraint::new(index(*p), index(*q)))
        .collect();
    constraints.sort_unstable();
    constraints.dedup();
    let triangulation =
        triangulate(&points, &constraints).map_err(|_| FarthestError::Triangulation)?;
    let at = triangulation.points();
    let mut out = Vec::new();
    for t in triangulation.triangles().chunks_exact(3) {
        let corners = [at[t[0] as usize], at[t[1] as usize], at[t[2] as usize]];
        let centroid = Point2::new(
            (corners[0].x + corners[1].x + corners[2].x) / 3.0,
            (corners[0].y + corners[1].y + corners[2].y) / 3.0,
        );
        if contains(&map.region, centroid)? {
            out.push(corners);
        }
    }
    Ok(out)
}

/// Whether `p` lies in the closed counter-clockwise triangle, exactly.
fn in_triangle(t: &[Point2; 3], p: Point2) -> Result<bool, RouteError> {
    for i in 0..3 {
        if side(t[i], t[(i + 1) % 3], p)? == Sign::Negative {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether the closed segment `pq` meets the closed triangle, exactly.
fn meets_triangle(p: Point2, q: Point2, t: &[Point2; 3]) -> Result<bool, RouteError> {
    if in_triangle(t, p)? || in_triangle(t, q)? {
        return Ok(true);
    }
    for i in 0..3 {
        if segments_meet(p, q, t[i], t[(i + 1) % 3])? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether two closed segments share a point, exactly.
fn segments_meet(p: Point2, q: Point2, r: Point2, s: Point2) -> Result<bool, RouteError> {
    let (d1, d2) = (side(p, q, r)?, side(p, q, s)?);
    let (d3, d4) = (side(r, s, p)?, side(r, s, q)?);
    if d1 != d2
        && d3 != d4
        && d1 != Sign::Zero
        && d2 != Sign::Zero
        && d3 != Sign::Zero
        && d4 != Sign::Zero
    {
        return Ok(true);
    }
    Ok((d1 == Sign::Zero && within(p, q, r))
        || (d2 == Sign::Zero && within(p, q, s))
        || (d3 == Sign::Zero && within(r, s, p))
        || (d4 == Sign::Zero && within(r, s, q)))
}

/// Whether `p` lies in the closed polygon (outer ring minus the open
/// holes), exactly.
fn in_polygon(polygon: &Polygon, p: Point2) -> Result<bool, RouteError> {
    if on_ring(&polygon.outer, p)? {
        return Ok(true);
    }
    if winding(&polygon.outer, p)? == 0 {
        return Ok(false);
    }
    for hole in &polygon.holes {
        if !on_ring(hole, p)? && winding(hole, p)? != 0 {
            return Ok(false);
        }
    }
    Ok(true)
}

fn on_ring(ring: &Ring, p: Point2) -> Result<bool, RouteError> {
    for (a, b) in ring_edges(ring) {
        if side(a, b, p)? == Sign::Zero && within(a, b, p) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Winding number of a ring around a point off it, from exact sides.
fn winding(ring: &Ring, p: Point2) -> Result<i32, RouteError> {
    let mut w = 0;
    for (a, b) in ring_edges(ring) {
        if a.y <= p.y && b.y > p.y && side(a, b, p)? == Sign::Positive {
            w += 1;
        } else if b.y <= p.y && a.y > p.y && side(a, b, p)? == Sign::Negative {
            w -= 1;
        }
    }
    Ok(w)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    fn square() -> Polygon {
        Polygon {
            outer: Ring {
                points: vec![p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0), p(0.0, 4.0)],
            },
            holes: Vec::new(),
        }
    }

    #[test]
    fn anchors_lie_in_their_triangle_and_off_every_barrier() {
        let barrier = vec![vec![p(1.0, 1.0), p(3.0, 3.0)]];
        let map = distance_map(&[square()], &barrier, &[p(0.5, 3.5)]).unwrap();
        let roots = [[p(0.0, 0.0), p(4.0, 0.0), p(4.0, 4.0)]];
        let subregion = square();
        let search = Search {
            map: &map,
            subregion: &subregion,
            roots: &roots,
            evaluated: HashMap::new(),
            lower: None,
        };
        assert!(search.admissible(0, p(3.0, 1.0)).unwrap());
        // On the barrier, though inside the triangle.
        assert!(!search.admissible(0, p(2.0, 2.0)).unwrap());
        // Off the barrier, but outside the triangle.
        assert!(!search.admissible(0, p(1.0, 3.0)).unwrap());
    }

    #[test]
    fn a_triangle_is_outside_only_when_no_subregion_edge_meets_it() {
        let map = distance_map(&[square()], &[], &[p(0.5, 0.5)]).unwrap();
        let subregion = Polygon {
            outer: Ring {
                points: vec![p(1.0, 1.0), p(3.0, 1.0), p(3.0, 3.0), p(1.0, 3.0)],
            },
            holes: Vec::new(),
        };
        let search = Search {
            map: &map,
            subregion: &subregion,
            roots: &[],
            evaluated: HashMap::new(),
            lower: None,
        };
        // First corner outside, but the triangle reaches in.
        assert!(!search
            .outside(&[p(0.0, 0.0), p(2.0, 0.0), p(2.0, 2.0)])
            .unwrap());
        assert!(search
            .outside(&[p(0.0, 0.0), p(0.9, 0.0), p(0.0, 0.9)])
            .unwrap());
        assert!(!search
            .outside(&[p(1.5, 1.5), p(2.0, 1.5), p(2.0, 2.0)])
            .unwrap());
    }
}
