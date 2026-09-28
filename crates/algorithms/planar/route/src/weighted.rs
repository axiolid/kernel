//! Weighted distance maps: travel costing a factor inside cost polygons
//! (#195).
//!
//! # The cost
//!
//! Each [`CostRegion`] is a polygon with a factor of at least 1; a metre
//! travelled inside it counts `factor` metres. Outside every polygon the
//! factor is 1, and where polygons overlap the greatest applies. The
//! weighted distance is the least cost of a walk in the free space to the
//! nearest target. On a cost polygon's edge the cheaper side applies: a
//! walk along the edge can always be moved off it, to that side, for as
//! little extra as wanted.
//!
//! Weighted shortest paths bend where they cross a cost edge (Snell's
//! law), so there is no exact visibility graph. The map brackets the
//! distance instead.
//!
//! # The upper bound
//!
//! Every cost edge carries points at most `spacing` apart. The visibility
//! graph over the region's vertices, the cost polygons' vertices, the
//! targets and those points is exact, as for [`crate::distance_map`], and
//! each edge costs its segment's weighted length, integrated piece by piece
//! across the cost polygons and rounded up. Every path in the graph is a
//! walk, so its cost bounds the distance from above; the query's route is
//! that walk. A walk crossing a cost edge between two of its points is only
//! longer by the square of the offset, so the upper bound is close.
//!
//! # The lower bound
//!
//! An optimal walk is straight between the points where it bends or
//! crosses a cost edge. Cut it there: each piece is a segment inside one
//! cell of constant factor, or along one line, from a vertex or a point of
//! a cost edge to another. Each cost edge is split into intervals, even
//! ones at most `spacing` long, graded down to 1/64 of that toward the
//! edge's vertices; a point of an edge lies in some interval. So every
//! piece runs between two nodes -- vertices, or intervals -- and costs at
//! least
//!
//! ```text
//! factor  x  least distance between the two nodes,
//! ```
//!
//! where the factor is the greatest of: that of every polygon holding the
//! whole hull of the two nodes, and, at an interval, that of the side the
//! piece leaves or enters it on. A piece along a line costs at least the
//! weighted length of the gap between its nodes. The graph over the nodes
//! with those weights has a path no dearer than the walk, so its distances
//! bound the walk from below. A hop is left out only when one obstacle --
//! or, across a cell, one cost edge -- certainly crosses every such piece
//! (a piece from an interval's end that is a vertex is the vertex's own).
//!
//! Relaxed like that, a path could enter an interval at one end and leave
//! at the other, and slide along an edge from interval to interval for
//! nothing. Three facts about optimal walks bar it; the graph's states
//! remember the line and side of each first hop to apply them:
//!
//! - no two pieces in a row along one line (they would be one piece);
//! - no touching a cost edge and turning back to the side it came from:
//!   the chord between the two pieces is shorter;
//! - no running along a cost edge between two touches on one side unless
//!   along is cheaper than that side (critical reflection): else the chord
//!   is no dearer.
//!
//! Cost edges along a wall get no intervals: no walk crosses them, and
//! none bends at a wall except at a vertex, the free side having one
//! factor.
//!
//! What is left is first order in the interval length where a walk
//! crosses a cost edge obliquely or wraps its corner; square crossings
//! lose nothing. Halving `spacing` about halves the gap.
//!
//! Every exact decision -- sides, crossings, containment -- is made with
//! certified predicates; lengths are rounded down for the lower bound and
//! up for the upper. The points cutting a cost edge at an angle are
//! interpolated, so they lie off its line by rounding: the sides of an
//! interval are taken of its edge's exact ends, an edge never blocks a
//! hop from its own intervals, and a hop along an edge is costed above as
//! the walk along the edge itself, which is that close to it.
//!
//! # At the walls
//!
//! A cost region clipped to the free region meets the walls only up to
//! rounding. A vertex left a hair inside would open a sliver along the
//! wall, costing 1; so a vertex that near a wall, on its free side, is
//! moved just beyond it, where the region costs nothing (see
//! [`weighted_distance_map`]). An edge left so on or beyond a wall is
//! wall-borne: it gets no intervals, and along the wall its free side's
//! factor counts.
//!
//! # Seeded targets
//!
//! Each target may start at its own cost (see
//! [`weighted_distance_map_seeded`]): both searches start there, so both
//! bounds count it.

use core::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap};

use axiolid_contracts::Sign;
use axiolid_core::Point2;
use axiolid_overlay::Polygon;

use crate::graph::{region_left, Graph, Side};
use crate::map::{
    free_triangles_in, in_polygon, in_triangle, meets_triangle, outside, segments_meet,
};
use crate::{
    contains, crosses, dedup_points, obstacle_segments, ring_edges, side, validate_region, within,
    Farthest, FarthestError, LengthInterval, MapError, Route, RouteError, Unreachable, MAX_CELLS,
};

/// Halvings of the end intervals of each cost edge toward its vertices.
const GRADING: u32 = 6;

/// Graph nodes a weighted map builds at most: region, barrier, target and
/// cost-polygon vertices, and the points along cost edges.
pub const MAX_WEIGHTED_NODES: usize = 2048;

/// A polygon inside which travel costs `factor` times its length.
#[derive(Debug, Clone, PartialEq)]
pub struct CostRegion {
    /// The polygon, closed.
    pub polygon: Polygon,
    /// At least 1.
    pub factor: f64,
}

impl CostRegion {
    /// A cost region.
    #[must_use]
    pub fn new(polygon: Polygon, factor: f64) -> Self {
        Self { polygon, factor }
    }
}

/// The nearest target by weighted distance, bracketed, and a walk there.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct WeightedReach {
    /// Index of the target the walk reaches, in the targets given.
    pub target: usize,
    /// Contains the weighted distance to the nearest target, counting the
    /// target's start weight (see [`weighted_distance_map_seeded`]).
    /// `upper` is the weighted cost of `route` plus that weight, rounded
    /// up.
    pub cost: LengthInterval,
    /// The walk, from the query point to the target. Its `length` is its
    /// plain Euclidean length.
    pub route: Route,
}

/// What a graph node is.
#[derive(Debug, Clone, Copy)]
enum Kind {
    /// A region, barrier, target or cost-polygon vertex: exact.
    Vertex,
    /// A point on a cost edge standing for the interval `a`-`b` of it, on
    /// line `line` and cost piece `piece` (no vertex lies strictly inside
    /// a piece), with the factors just left and right of `a`-`b`, and the
    /// factor of travel along it: the cheaper free side.
    Interval {
        a: Point2,
        b: Point2,
        /// The piece's own ends, exact: `a` and `b` are interpolated and
        /// may lie off its line by rounding, so sides are taken of this.
        p: Point2,
        q: Point2,
        line: u32,
        piece: u32,
        left: f64,
        right: f64,
        along: f64,
        a_vertex: bool,
        b_vertex: bool,
    },
}

const NO_LINE: u32 = u32::MAX;

/// Weighted distances to the nearest of several targets, bracketed.
#[derive(Debug, Clone)]
pub struct WeightedMap {
    region: Vec<Polygon>,
    walls: Vec<(Point2, Point2)>,
    obstacles: Vec<(Point2, Point2)>,
    weights: Weights,
    nodes: Vec<Point2>,
    kinds: Vec<Kind>,
    /// Lines of cost edges through each node.
    lines: Vec<Vec<u32>>,
    graph: Graph,
    /// Per state: upper-bound distance, next state, target.
    upper: Vec<f64>,
    next: Vec<usize>,
    target: Vec<usize>,
    /// Per state: lower-bound distances by the line of the hop that
    /// reached the state (`NO_LINE` for none).
    lower: Vec<Vec<(Tag, f64)>>,
    sites: Vec<Point2>,
    /// The distance each target starts at.
    seeds: Vec<f64>,
    /// The cost regions, as snapped (see [`weighted_distance_map`]).
    costs: Vec<CostRegion>,
    scale: f64,
}

/// A weighted distance map from `targets` over `region`, avoiding
/// `barriers`, with travel inside `costs` weighted by their factors and
/// points along cost edges at most `spacing` apart.
///
/// A cost region may touch the region's boundary up to rounding, as one
/// clipped to the free region does (#198): a cost vertex within 2^-24 of
/// the region's extent (and a few ulps) of a region edge, on its free
/// side, is moved just across it, and a cost edge may cross a region edge
/// that near one of either edge's ends. The map is the one for the
/// regions so moved: no sliver along a wall is left for a walk to slip
/// through at factor 1.
///
/// # Errors
///
/// [`MapError`] for malformed input, no targets, a target outside the
/// region, a factor below 1 or not finite ([`MapError::InvalidFactor`]), a
/// spacing not positive and finite ([`MapError::InvalidSpacing`]), a cost
/// edge that crosses an obstacle or another cost edge other than by such
/// a touch, or lies along a barrier ([`MapError::CostCrossing`]), or more
/// than [`MAX_WEIGHTED_NODES`] graph nodes.
pub fn weighted_distance_map(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    targets: &[Point2],
    costs: &[CostRegion],
    spacing: f64,
) -> Result<WeightedMap, MapError> {
    weighted_distance_map_within(
        region,
        barriers,
        targets,
        costs,
        spacing,
        MAX_WEIGHTED_NODES,
    )
}

/// [`weighted_distance_map`] with a caller-chosen node budget.
///
/// # Errors
///
/// As [`weighted_distance_map`], with `budget` for [`MAX_WEIGHTED_NODES`].
pub fn weighted_distance_map_within(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    targets: &[Point2],
    costs: &[CostRegion],
    spacing: f64,
    budget: usize,
) -> Result<WeightedMap, MapError> {
    let seeded: Vec<(Point2, f64)> = targets.iter().map(|t| (*t, 0.0)).collect();
    weighted_distance_map_seeded_within(region, barriers, &seeded, costs, spacing, budget)
}

/// A weighted distance map whose targets each start at their own cost:
/// the distance from a point is the least, over targets, of the weighted
/// cost of a walk to the target plus the target's weight (#198), as
/// [`crate::distance_map_weighted`] is for plain maps. With every weight
/// zero it is [`weighted_distance_map`].
///
/// # Errors
///
/// As [`weighted_distance_map`], and [`MapError::InvalidWeight`] for a
/// weight that is negative or not finite.
pub fn weighted_distance_map_seeded(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    targets: &[(Point2, f64)],
    costs: &[CostRegion],
    spacing: f64,
) -> Result<WeightedMap, MapError> {
    weighted_distance_map_seeded_within(
        region,
        barriers,
        targets,
        costs,
        spacing,
        MAX_WEIGHTED_NODES,
    )
}

/// [`weighted_distance_map_seeded`] with a caller-chosen node budget.
///
/// # Errors
///
/// As [`weighted_distance_map_seeded`], with `budget` for
/// [`MAX_WEIGHTED_NODES`].
pub fn weighted_distance_map_seeded_within(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    seeded: &[(Point2, f64)],
    costs: &[CostRegion],
    spacing: f64,
    budget: usize,
) -> Result<WeightedMap, MapError> {
    validate_region(region, barriers)?;
    if seeded.is_empty() {
        return Err(MapError::NoTargets);
    }
    let targets: Vec<Point2> = seeded.iter().map(|(t, _)| *t).collect();
    let seeds: Vec<f64> = seeded.iter().map(|(_, w)| *w).collect();
    if !targets.iter().all(|t| t.is_finite()) {
        return Err(RouteError::NonFinitePoint.into());
    }
    if let Some(index) = seeds.iter().position(|w| !(w.is_finite() && *w >= 0.0)) {
        return Err(MapError::InvalidWeight { index });
    }
    if !(spacing.is_finite() && spacing > 0.0) {
        return Err(MapError::InvalidSpacing);
    }
    for (index, cost) in costs.iter().enumerate() {
        if !(cost.factor.is_finite() && cost.factor >= 1.0) {
            return Err(MapError::InvalidFactor { index });
        }
    }
    validate_region(
        &costs.iter().map(|c| c.polygon.clone()).collect::<Vec<_>>(),
        &[],
    )?;
    let reach = touch_reach(region);
    let costs = &touch_walls(costs, region, reach)?;
    let polygons: Vec<Polygon> = costs.iter().map(|c| c.polygon.clone()).collect();
    validate_region(&polygons, &[])?;
    for (index, t) in targets.iter().enumerate() {
        if !contains(region, *t)? {
            return Err(MapError::TargetOutside { index });
        }
    }
    let obstacles = obstacle_segments(region, barriers);
    let walls: Vec<(Point2, Point2)> = barriers
        .iter()
        .flat_map(|b| b.windows(2).map(|w| (w[0], w[1])))
        .collect();
    // Vertices first, then the points along cost edges.
    let mut nodes = targets.clone();
    for polygon in region.iter().chain(polygons.iter()) {
        for ring in core::iter::once(&polygon.outer).chain(polygon.holes.iter()) {
            nodes.extend(ring.points.iter().copied());
        }
    }
    for barrier in barriers {
        nodes.extend(barrier.iter().copied());
    }
    dedup_points(&mut nodes);
    // Cost edges are cut at every vertex on them, so no interval's point
    // is a vertex.
    let weights = Weights::new(costs, &nodes, region)?;
    weights.check_crossings(&obstacles, &walls, reach)?;
    let mut kinds = vec![Kind::Vertex; nodes.len()];
    let mut spans: Vec<(Point2, Point2, u32)> = Vec::new();
    for piece in &weights.pieces {
        let (p, q) = if (piece.p.x, piece.p.y) <= (piece.q.x, piece.q.y) {
            (piece.p, piece.q)
        } else {
            (piece.q, piece.p)
        };
        // An edge two polygons share gets its points once.
        if !spans.iter().any(|(u, v, _)| *u == p && *v == q) {
            spans.push((p, q, piece.line));
        }
    }
    for (piece, (p, q, line)) in spans.into_iter().enumerate() {
        // Along a wall (or a hair beyond it) no walk crosses, and none
        // bends except at a vertex: the one free side has one factor, so a
        // chord is shorter. Every other piece has both sides free.
        if weights.beyond_wall(p, q, reach)? {
            continue;
        }
        let count = ((q - p).length() / spacing).ceil().max(1.0) as usize;
        // Even cuts, graded toward both ends: a walk wrapping the corner
        // loses at most the length of the interval it touches there.
        let mut cuts: Vec<f64> = (0..=count).map(|k| k as f64 / count as f64).collect();
        let first = 1.0 / count as f64;
        for level in 1..=GRADING {
            let t = first / f64::from(1u32 << level);
            cuts.push(t);
            cuts.push(1.0 - t);
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
        let count = cuts.len() - 1;
        if nodes.len() + count > budget {
            return Err(RouteError::TooManyVertices {
                supplied: nodes.len() + count,
                budget,
                lower_bound: 0.0,
            }
            .into());
        }
        let at = |k: usize| {
            if k == 0 {
                p
            } else if k == count {
                q
            } else {
                let t = cuts[k];
                Point2::new(p.x + (q.x - p.x) * t, p.y + (q.y - p.y) * t)
            }
        };
        for k in 0..count {
            let (a, b) = (at(k), at(k + 1));
            nodes.push(Point2::new(0.5 * a.x + 0.5 * b.x, 0.5 * a.y + 0.5 * b.y));
            let (left, right) = weights.sides(p, q, a, b)?;
            let along = left.min(right);
            kinds.push(Kind::Interval {
                a,
                b,
                p,
                q,
                line,
                piece: piece as u32,
                left,
                right,
                along,
                // A piece runs from vertex to vertex; its inner cuts are
                // not vertices.
                a_vertex: k == 0,
                b_vertex: k + 1 == count,
            });
        }
    }
    if nodes.len() > budget {
        return Err(RouteError::TooManyVertices {
            supplied: nodes.len(),
            budget,
            lower_bound: 0.0,
        }
        .into());
    }
    let lines: Vec<Vec<u32>> = nodes
        .iter()
        .zip(&kinds)
        .map(|(v, kind)| match kind {
            Kind::Interval { line, .. } => Ok(vec![*line]),
            Kind::Vertex => weights.lines_through(*v),
        })
        .collect::<Result<_, RouteError>>()?;
    let graph = Graph::build(&nodes, region, barriers, &obstacles)?;
    let scale = nodes
        .iter()
        .fold(1.0f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));

    let mut map = WeightedMap {
        region: region.to_vec(),
        walls,
        obstacles,
        weights,
        nodes,
        kinds,
        lines,
        graph,
        upper: Vec::new(),
        next: Vec::new(),
        target: Vec::new(),
        lower: Vec::new(),
        sites: targets.clone(),
        seeds: seeds.clone(),
        costs: costs.to_vec(),
        scale,
    };
    // Every state of every target at its weight; of targets on one point,
    // the lightest (then the first) seeds it.
    let mut sources = Vec::new();
    let mut seed = vec![usize::MAX; map.graph.adjacency.len()];
    for (i, node) in map.nodes.iter().enumerate() {
        let lightest = (0..targets.len())
            .filter(|&t| targets[t] == *node)
            .min_by(|&a, &b| seeds[a].total_cmp(&seeds[b]).then(a.cmp(&b)));
        if let Some(t) = lightest {
            for state in map.graph.states(i) {
                sources.push((state, seeds[t]));
                seed[state] = t;
            }
        }
    }
    map.upper_bounds(&sources, &seed)?;
    map.lower_bounds(&sources)?;
    Ok(map)
}

impl WeightedMap {
    /// Graph nodes: vertices and points along cost edges.
    #[must_use]
    pub fn graph_vertices(&self) -> usize {
        self.nodes.len()
    }

    /// How many targets the map was built from.
    #[must_use]
    pub fn targets(&self) -> usize {
        self.sites.len()
    }

    /// The nearest target from `point` by weighted distance, the distance
    /// bracketed, and a walk there.
    ///
    /// `Ok(Err(StartOutside))` for a point outside the region,
    /// `Ok(Err(DisconnectedComponents))` for one no target can reach.
    ///
    /// # Errors
    ///
    /// [`RouteError`] for a non-finite point or an undecidable predicate.
    pub fn nearest(&self, point: Point2) -> Result<Result<WeightedReach, Unreachable>, RouteError> {
        if !point.is_finite() {
            return Err(RouteError::NonFinitePoint);
        }
        if !contains(&self.region, point)? {
            return Ok(Err(Unreachable::StartOutside));
        }
        let Some((upper, first)) = self.upper_at(point)? else {
            return Ok(Err(Unreachable::DisconnectedComponents));
        };
        let lower = self.lower_at(point, upper)?.min(upper);
        let mut polyline = vec![point];
        let mut state = first;
        loop {
            let at = self.nodes[self.graph.node(state)];
            if polyline.last() != Some(&at) {
                polyline.push(at);
            }
            if self.next[state] == usize::MAX {
                break;
            }
            state = self.next[state];
        }
        let length = polyline.windows(2).map(|w| (w[1] - w[0]).length()).sum();
        Ok(Ok(WeightedReach {
            target: self.target[first],
            cost: LengthInterval { lower, upper },
            route: Route {
                polyline,
                length,
                graph_vertices: self.nodes.len(),
            },
        }))
    }

    /// The region the map covers.
    pub(crate) fn region(&self) -> &[Polygon] {
        &self.region
    }

    /// Region edges and barrier segments.
    pub(crate) fn obstacles(&self) -> &[(Point2, Point2)] {
        &self.obstacles
    }

    /// Barrier segments.
    pub(crate) fn walls(&self) -> &[(Point2, Point2)] {
        &self.walls
    }

    /// The targets given, and the cost each starts at.
    pub(crate) fn seeded(&self) -> impl Iterator<Item = (Point2, f64)> + '_ {
        self.sites.iter().copied().zip(self.seeds.iter().copied())
    }

    /// Whether two maps cover the same free space at the same costs.
    pub(crate) fn same_space(&self, other: &Self) -> bool {
        self.region == other.region && self.walls == other.walls && self.costs == other.costs
    }

    /// The greatest factor meeting the closed triangle: the most a metre
    /// can cost in it.
    pub(crate) fn steepest(&self, t: &[Point2; 3]) -> Result<f64, RouteError> {
        self.weights.steepest(t)
    }

    /// The factor all through the closed triangle when no cost edge meets
    /// it, else 1: what a straight piece inside it costs a metre at least.
    pub(crate) fn inside_factor(&self, t: &[Point2; 3]) -> Result<f64, RouteError> {
        for piece in &self.weights.pieces {
            if meets_triangle(piece.p, piece.q, t)? {
                return Ok(1.0);
            }
        }
        let mut best = 1.0f64;
        for (polygon, &factor) in self.weights.polygons.iter().zip(&self.weights.factors) {
            if factor > best && in_polygon(polygon, t[0])? {
                best = factor;
            }
        }
        Ok(best)
    }

    /// Every node's span -- a vertex twice, or an interval -- with the
    /// least lower bound a walk from a point of it has, when finite.
    pub(crate) fn spans(&self) -> Vec<((Point2, Point2), f64)> {
        (0..self.nodes.len())
            .filter_map(|i| {
                let least = self
                    .graph
                    .states(i)
                    .flat_map(|s| self.lower[s].iter().map(|(_, d)| *d))
                    .fold(f64::INFINITY, f64::min);
                least.is_finite().then(|| (self.span(i), least))
            })
            .collect()
    }

    /// The bracket at a point known to be inside the region; `None` when
    /// no target is reached.
    pub(crate) fn bracket(&self, point: Point2) -> Result<Option<(f64, f64)>, RouteError> {
        let Some((upper, _)) = self.upper_at(point)? else {
            return Ok(None);
        };
        Ok(Some((self.lower_at(point, upper)?.min(upper), upper)))
    }

    /// Upper-bound Dijkstra over the exact graph, each edge costing its
    /// segment's weighted length rounded up.
    fn upper_bounds(&mut self, sources: &[(usize, f64)], seed: &[usize]) -> Result<(), RouteError> {
        let states = self.graph.adjacency.len();
        let mut cost: HashMap<(usize, usize), f64> = HashMap::new();
        let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); states];
        for (u, edges) in self.graph.adjacency.iter().enumerate() {
            let i = self.graph.node(u);
            for &(w, _) in edges {
                let j = self.graph.node(w);
                let key = (i.min(j), i.max(j));
                let c = match cost.get(&key) {
                    Some(c) => *c,
                    None => {
                        let c = self
                            .weights
                            .segment(self.nodes[i], self.nodes[j], self.scale)?
                            .1;
                        cost.insert(key, c);
                        c
                    }
                };
                adjacency[u].push((w, c));
            }
        }
        let (distance, previous) = dijkstra(&adjacency, sources);
        // Follow each state's chain to its end: a heavy target may itself
        // be reached from a lighter one, so being seeded does not end it.
        let mut target = seed.to_vec();
        for (state, t) in target.iter_mut().enumerate() {
            let mut at = state;
            let mut steps = 0;
            while previous[at] != usize::MAX && steps <= seed.len() {
                at = previous[at];
                steps += 1;
            }
            *t = seed[at];
        }
        self.upper = distance;
        self.next = previous;
        self.target = target;
        Ok(())
    }

    /// The states a hop from node `i` toward the segment `b1`-`b2` may
    /// leave in: the cone's sectors at a vertex, every free state of an
    /// interval.
    fn hop_states(&self, i: usize, b1: Point2, b2: Point2) -> Result<Vec<usize>, RouteError> {
        match self.kinds[i] {
            Kind::Vertex => self.graph.cone_states(i, self.nodes[i], b1, b2),
            Kind::Interval { .. } => Ok(self.graph.states(i).collect()),
        }
    }

    /// The segment a node stands for: its interval, or the vertex twice.
    fn span(&self, i: usize) -> (Point2, Point2) {
        match self.kinds[i] {
            Kind::Vertex => (self.nodes[i], self.nodes[i]),
            Kind::Interval { a, b, .. } => (a, b),
        }
    }

    /// The line both spans lie on, if they share a cost edge's line.
    fn shared_line(&self, i: usize, j: usize) -> u32 {
        self.lines[i]
            .iter()
            .find(|l| self.lines[j].contains(l))
            .copied()
            .unwrap_or(NO_LINE)
    }

    /// A lower bound on the cost of a straight piece from a point of span
    /// `(a1, a2)` to a point of span `(b1, b2)`, or `None` when an obstacle
    /// certainly blocks every such piece. Spans on one line bound the
    /// piece by the weighted length of the gap between them, however many
    /// cells it crosses; any other piece stays in one cell and crosses no
    /// cost edge.
    fn hop(
        &self,
        (a1, a2): (Point2, Point2),
        (b1, b2): (Point2, Point2),
        ka: Kind,
        kb: Kind,
    ) -> Result<Option<f64>, RouteError> {
        let (ea, eb) = (vertex_ends(ka), vertex_ends(kb));
        let blocks = |p: Point2, q: Point2| -> Result<bool, RouteError> {
            for (u, skip_u) in [(a1, ea.0), (a2, ea.1)] {
                for (v, skip_v) in [(b1, eb.0), (b2, eb.1)] {
                    if !crosses(u, v, p, q)? && !((skip_u || skip_v) && through_end(u, v, p, q)?) {
                        return Ok(false);
                    }
                }
            }
            Ok(true)
        };
        for &(p, q) in &self.obstacles {
            if blocks(p, q)? {
                return Ok(None);
            }
        }
        if convex_hull(&[a1, a2, b1, b2])?.is_none() {
            // On one line: at least the weighted length of the gap between
            // the spans, which every such piece covers, whatever it runs
            // along.
            let Some((u, v)) = gap((a1, a2), (b1, b2)) else {
                return Ok(Some(0.0));
            };
            return Ok(Some(self.weights.segment(u, v, self.scale)?.0));
        }
        // A piece on an interval's own line never crosses a piece from it;
        // its interpolated ends may lie a rounding off that line, beyond.
        let own = |kind: Kind| match kind {
            Kind::Interval { line, .. } => line,
            Kind::Vertex => NO_LINE,
        };
        let (la, lb) = (own(ka), own(kb));
        for piece in &self.weights.pieces {
            if piece.line != la && piece.line != lb && blocks(piece.p, piece.q)? {
                return Ok(None);
            }
        }
        // A piece from inside an interval leaves it into the side the
        // other span lies on, into the cell there; one from the interval's
        // end, if a vertex, is the vertex's hop.
        let factor = self
            .weights
            .hull_factor(&[a1, a2, b1, b2])?
            .max(leaving(ka, (b1, b2), eb)?)
            .max(leaving(kb, (a1, a2), ea)?);
        let d = span_distance((a1, a2), (b1, b2))?;
        Ok(Some(round_down(factor * d, self.scale)))
    }

    /// Which side of interval node `k`'s line node `n`'s span lies
    /// strictly on: `LEFT`, `RIGHT`, or `NONE` when it touches the line or
    /// `k` is a vertex. An end of `n`'s interval that is a vertex does not
    /// count: a walk through it has its event at the vertex's own node.
    fn side_of(&self, k: usize, n: usize) -> Result<u8, RouteError> {
        let (u, v) = self.span(n);
        let ends = match self.kinds[n] {
            Kind::Interval {
                a_vertex, b_vertex, ..
            } => (a_vertex, b_vertex),
            Kind::Vertex => (false, false),
        };
        self.side_of_span(k, (u, v), ends)
    }

    fn side_of_span(
        &self,
        k: usize,
        (u, v): (Point2, Point2),
        (skip_u, skip_v): (bool, bool),
    ) -> Result<u8, RouteError> {
        let Kind::Interval { p, q, .. } = self.kinds[k] else {
            return Ok(NONE);
        };
        let (su, sv) = (side(p, q, u)?, side(p, q, v)?);
        let su = if su == Sign::Zero && skip_u { sv } else { su };
        let sv = if sv == Sign::Zero && skip_v { su } else { sv };
        Ok(match (su, sv) {
            (Sign::Positive, Sign::Positive) => LEFT,
            (Sign::Negative, Sign::Negative) => RIGHT,
            _ => NONE,
        })
    }

    /// The edge from node `n` (its span `sn`) into node `k`, with the
    /// sides it leaves `n` and arrives at `k` on.
    fn edge(&self, n: usize, k: usize, from: usize, weight: f64) -> Result<Edge, RouteError> {
        let line = self.shared_line(n, k);
        let (leave, arrive) = if line == NO_LINE {
            (self.side_of(n, k)?, self.side_of(k, n)?)
        } else {
            (NONE, NONE)
        };
        let same_piece = match (self.kinds[n], self.kinds[k]) {
            (Kind::Interval { piece: pn, .. }, Kind::Interval { piece: pk, .. }) => pn == pk,
            _ => false,
        };
        Ok(Edge {
            from,
            weight,
            line,
            leave,
            arrive,
            same_piece,
        })
    }

    /// Lower-bound Dijkstra: exact vertex-to-vertex edges costing their
    /// weighted length rounded down, and hops to and between intervals,
    /// over states that remember the first hop's line and side (see
    /// [`allowed`]).
    fn lower_bounds(&mut self, sources: &[(usize, f64)]) -> Result<(), RouteError> {
        let states = self.graph.adjacency.len();
        let mut incoming: Vec<Vec<Edge>> = vec![Vec::new(); states];
        let is_vertex = |i: usize, kinds: &[Kind]| matches!(kinds[i], Kind::Vertex);
        let mut cost: HashMap<(usize, usize), f64> = HashMap::new();
        for (u, edges) in self.graph.adjacency.iter().enumerate() {
            let i = self.graph.node(u);
            if !is_vertex(i, &self.kinds) {
                continue;
            }
            for &(w, _) in edges {
                let j = self.graph.node(w);
                if !is_vertex(j, &self.kinds) {
                    continue;
                }
                let key = (i.min(j), i.max(j));
                let c = match cost.get(&key) {
                    Some(c) => *c,
                    None => {
                        let c = self
                            .weights
                            .segment(self.nodes[i], self.nodes[j], self.scale)?
                            .0;
                        cost.insert(key, c);
                        c
                    }
                };
                // Travelled from `u` into `w`.
                incoming[w].push(self.edge(i, j, u, c)?);
            }
        }
        let n = self.nodes.len();
        for i in 0..n {
            for j in i + 1..n {
                if is_vertex(i, &self.kinds) && is_vertex(j, &self.kinds) {
                    continue;
                }
                let (si, sj) = (self.span(i), self.span(j));
                let Some(c) = self.hop(si, sj, self.kinds[i], self.kinds[j])? else {
                    continue;
                };
                let from = self.hop_states(i, sj.0, sj.1)?;
                let to = self.hop_states(j, si.0, si.1)?;
                for &u in &from {
                    for &w in &to {
                        incoming[w].push(self.edge(i, j, u, c)?);
                        incoming[u].push(self.edge(j, i, w, c)?);
                    }
                }
            }
        }
        let graph = &self.graph;
        let kinds = &self.kinds;
        self.lower = tagged_dijkstra(&incoming, sources, |s| kinds[graph.node(s)]);
        Ok(())
    }

    /// The least upper bound from `point` and the first state on the way.
    fn upper_at(&self, point: Point2) -> Result<Option<(f64, usize)>, RouteError> {
        let mut best: Option<(f64, usize)> = None;
        let offer = |value: f64, state: usize, best: &mut Option<(f64, usize)>| {
            if best.is_none_or(|(b, s)| value < b || (value == b && state < s)) {
                *best = Some((value, state));
            }
        };
        for (i, node) in self.nodes.iter().enumerate() {
            let nearest = self
                .graph
                .states(i)
                .map(|s| self.upper[s])
                .fold(f64::INFINITY, f64::min);
            if nearest.is_infinite() {
                continue;
            }
            if *node == point {
                for s in self.graph.states(i) {
                    offer(self.upper[s], s, &mut best);
                }
                continue;
            }
            // A leg costs at least its length.
            if best.is_some_and(|(b, _)| (*node - point).length() + nearest > b) {
                continue;
            }
            let ok = crate::graph::sides(
                point,
                *node,
                &self.region,
                &self.obstacles,
                &self.nodes,
                &self.graph.stars,
                &self.graph.rayed,
            )?;
            if ok == [false, false] {
                continue;
            }
            let leg = self.weights.segment(point, *node, self.scale)?.1;
            for (k, on) in [Side::Left, Side::Right].into_iter().enumerate() {
                if !ok[k] {
                    continue;
                }
                if let Some(s) = self.graph.arrival(i, point, on)? {
                    if self.upper[s].is_finite() {
                        offer(leg + self.upper[s], s, &mut best);
                    }
                }
            }
        }
        Ok(best)
    }

    /// The least lower bound from `point`: over its first piece to a node
    /// and the node's bound, the piece exact to a vertex, fattened to an
    /// interval. Starts from `ceiling`, a value no less than the answer.
    fn lower_at(&self, point: Point2, ceiling: f64) -> Result<f64, RouteError> {
        let mut best = ceiling;
        let point_lines = self.weights.lines_through(point)?;
        for (i, node) in self.nodes.iter().enumerate() {
            let least = self
                .graph
                .states(i)
                .flat_map(|s| self.lower[s].iter().map(|(_, d)| *d))
                .fold(f64::INFINITY, f64::min);
            let kind = self.kinds[i];
            if least.is_infinite() {
                continue;
            }
            let line = self.lines[i]
                .iter()
                .find(|l| point_lines.contains(l))
                .copied()
                .unwrap_or(NO_LINE);
            // The first hop, as an edge into the node: the states it may
            // go on from.
            let first = Edge {
                from: usize::MAX,
                weight: 0.0,
                line,
                leave: NONE,
                arrive: if line == NO_LINE {
                    self.side_of_span(i, (point, point), (false, false))?
                } else {
                    NONE
                },
                same_piece: false,
            };
            let from = |states: &mut dyn Iterator<Item = usize>| -> f64 {
                states
                    .flat_map(|s| self.lower[s].iter())
                    .filter(|(tag, _)| allowed(kind, *tag, &first))
                    .map(|(_, d)| *d)
                    .fold(f64::INFINITY, f64::min)
            };
            if *node == point {
                // The walk starts here: any way on.
                let any = self
                    .graph
                    .states(i)
                    .flat_map(|s| self.lower[s].iter().map(|(_, d)| *d))
                    .fold(f64::INFINITY, f64::min);
                best = best.min(any);
                continue;
            }
            let (a, b) = self.span(i);
            if span_distance((point, point), (a, b))? + least >= best {
                continue;
            }
            match self.kinds[i] {
                Kind::Vertex => {
                    let ok = crate::graph::sides(
                        point,
                        *node,
                        &self.region,
                        &self.obstacles,
                        &self.nodes,
                        &self.graph.stars,
                        &self.graph.rayed,
                    )?;
                    if ok == [false, false] {
                        continue;
                    }
                    let leg = self.weights.segment(point, *node, self.scale)?.0;
                    for (k, on) in [Side::Left, Side::Right].into_iter().enumerate() {
                        if !ok[k] {
                            continue;
                        }
                        if let Some(s) = self.graph.arrival(i, point, on)? {
                            best = best.min(leg + from(&mut core::iter::once(s)));
                        }
                    }
                }
                Kind::Interval { .. } => {
                    let Some(leg) =
                        self.hop((point, point), (a, b), Kind::Vertex, self.kinds[i])?
                    else {
                        continue;
                    };
                    best = best.min(leg + from(&mut self.graph.states(i)));
                }
            }
        }
        Ok(best)
    }
}

/// Dijkstra with a binary heap from several sources, each at its start
/// distance: distance and the previous state. Ties break on the lowest
/// state.
fn dijkstra(adjacency: &[Vec<(usize, f64)>], sources: &[(usize, f64)]) -> (Vec<f64>, Vec<usize>) {
    let n = adjacency.len();
    let mut distance = vec![f64::INFINITY; n];
    let mut previous = vec![usize::MAX; n];
    let mut heap = BinaryHeap::new();
    for &(s, start) in sources {
        if start < distance[s] {
            distance[s] = start;
            heap.push(Reverse((Ordered(start), s)));
        }
    }
    while let Some(Reverse((Ordered(d), u))) = heap.pop() {
        if d > distance[u] {
            continue;
        }
        for &(w, c) in &adjacency[u] {
            let candidate = d + c;
            if candidate < distance[w] || (candidate == distance[w] && u < previous[w]) {
                if candidate < distance[w] {
                    heap.push(Reverse((Ordered(candidate), w)));
                }
                distance[w] = candidate;
                previous[w] = u;
            }
        }
    }
    (distance, previous)
}

/// No side, or unknown.
const NONE: u8 = 0;
const LEFT: u8 = 1;
const RIGHT: u8 = 2;
/// A first hop along the interval's piece, the next node's first hop
/// leaving on no known side, or on the left or right.
const ALONG: u8 = 3;

/// How a path leaves a state toward the targets: the line of its first hop
/// if that runs along a cost line, and at an interval the side it leaves
/// on (`LEFT`, `RIGHT`), or `ALONG + side` for a first hop along the
/// interval's own piece followed by one leaving on `side`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Tag {
    line: u32,
    side: u8,
}

/// An edge travelled from state `from` into the state it is listed under.
#[derive(Debug, Clone, Copy)]
struct Edge {
    from: usize,
    weight: f64,
    /// The cost line it runs along, or `NO_LINE`.
    line: u32,
    /// The side of the source interval's line it leaves on.
    leave: u8,
    /// The side of the destination interval's line it arrives from.
    arrive: u8,
    /// Along one cost piece, from one of its intervals to another.
    same_piece: bool,
}

/// Whether a path may arrive at a node of `kind` by `edge` and go on as
/// `tag` says. An optimal walk never takes two pieces along one line in a
/// row (they would be one piece), never touches a cost edge and turns
/// back to the side it came from (the chord is shorter), and never runs
/// along a cost piece between two such touches on one side unless along
/// is cheaper than that side (else the chord is no dearer).
fn allowed(kind: Kind, tag: Tag, edge: &Edge) -> bool {
    if edge.line != NO_LINE && tag.line == edge.line {
        return false;
    }
    if let Kind::Interval {
        left, right, along, ..
    } = kind
    {
        if edge.line == NO_LINE && edge.arrive != NONE {
            let factor = if edge.arrive == LEFT { left } else { right };
            if tag.side == edge.arrive {
                return false;
            }
            if tag.side > ALONG && tag.side - ALONG == edge.arrive && along >= factor {
                return false;
            }
        }
    }
    true
}

/// The tag of the source state of `edge`, whose destination went on as
/// `tag`.
fn tag_before(source: Kind, tag: Tag, edge: &Edge) -> Tag {
    let side = match source {
        Kind::Vertex => NONE,
        Kind::Interval { .. } if edge.line != NO_LINE => {
            if edge.same_piece && (tag.side == LEFT || tag.side == RIGHT) {
                ALONG + tag.side
            } else {
                ALONG
            }
        }
        Kind::Interval { .. } => edge.leave,
    };
    Tag {
        line: edge.line,
        side,
    }
}

/// Dijkstra backward from the targets over (state, tag), taking only the
/// transitions [`allowed`] permits. Per state, the settled distances by
/// tag.
fn tagged_dijkstra(
    incoming: &[Vec<Edge>],
    sources: &[(usize, f64)],
    kind: impl Fn(usize) -> Kind,
) -> Vec<Vec<(Tag, f64)>> {
    let mut settled: Vec<Vec<(Tag, f64)>> = vec![Vec::new(); incoming.len()];
    let mut best: HashMap<(usize, Tag), f64> = HashMap::new();
    let mut heap = BinaryHeap::new();
    let start = Tag {
        line: NO_LINE,
        side: NONE,
    };
    for &(s, weight) in sources {
        if best.get(&(s, start)).is_none_or(|b| weight < *b) {
            best.insert((s, start), weight);
            heap.push(Reverse((Ordered(weight), s, start)));
        }
    }
    while let Some(Reverse((Ordered(d), k, tag))) = heap.pop() {
        if settled[k].iter().any(|(t, _)| *t == tag) {
            continue;
        }
        settled[k].push((tag, d));
        let here = kind(k);
        for edge in &incoming[k] {
            if !allowed(here, tag, edge) {
                continue;
            }
            let before = tag_before(kind(edge.from), tag, edge);
            let candidate = d + edge.weight;
            let key = (edge.from, before);
            if best.get(&key).is_none_or(|b| candidate < *b) {
                best.insert(key, candidate);
                heap.push(Reverse((Ordered(candidate), edge.from, before)));
            }
        }
    }
    settled
}

/// A total order on finite and infinite floats for the heaps.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Ordered(f64);

impl Eq for Ordered {}

impl PartialOrd for Ordered {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Ordered {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// Round a nonnegative length down past the error of the few operations
/// that produced it.
fn round_down(value: f64, scale: f64) -> f64 {
    (value * (1.0 - 16.0 * f64::EPSILON) - 16.0 * f64::EPSILON * scale).max(0.0)
}

fn round_up(value: f64, scale: f64) -> f64 {
    value * (1.0 + 16.0 * f64::EPSILON) + 16.0 * f64::EPSILON * scale
}

/// The factor a piece from inside an interval node's span toward the span
/// `other` must pay: that of the side of the interval `other` lies
/// strictly on, or the lesser when it touches the interval's line. 1 for
/// a vertex.
fn leaving(
    kind: Kind,
    (o1, o2): (Point2, Point2),
    (skip1, skip2): (bool, bool),
) -> Result<f64, RouteError> {
    let Kind::Interval {
        p, q, left, right, ..
    } = kind
    else {
        return Ok(1.0);
    };
    let (s1, s2) = (side(p, q, o1)?, side(p, q, o2)?);
    // An end of the other span that is a vertex does not count.
    let s1 = if s1 == Sign::Zero && skip1 { s2 } else { s1 };
    let s2 = if s2 == Sign::Zero && skip2 { s1 } else { s2 };
    Ok(match (s1, s2) {
        (Sign::Positive, Sign::Positive) => left,
        (Sign::Negative, Sign::Negative) => right,
        _ => left.min(right),
    })
}

/// Which ends of a node's span are vertices, whose walks the vertices'
/// own nodes carry.
fn vertex_ends(kind: Kind) -> (bool, bool) {
    match kind {
        Kind::Interval {
            a_vertex, b_vertex, ..
        } => (a_vertex, b_vertex),
        Kind::Vertex => (false, false),
    }
}

/// Whether the segment `u`-`v` passes through an end of the segment
/// `p`-`q`, with neither of its own ends on that segment's line: a touch
/// that still separates the points just beside it, when every other corner
/// segment crosses properly.
fn through_end(u: Point2, v: Point2, p: Point2, q: Point2) -> Result<bool, RouteError> {
    if side(p, q, u)? == Sign::Zero || side(p, q, v)? == Sign::Zero {
        return Ok(false);
    }
    for e in [p, q] {
        if within(u, v, e) && side(u, v, e)? == Sign::Zero {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The gap between two spans on one line: the facing ends, or `None`
/// when the spans overlap.
fn gap((a1, a2): (Point2, Point2), (b1, b2): (Point2, Point2)) -> Option<(Point2, Point2)> {
    let d = if a1 != a2 { a2 - a1 } else { b2 - b1 };
    let key = |p: Point2| (p - a1).dot(d);
    let (alo, ahi) = if key(a1) <= key(a2) {
        (a1, a2)
    } else {
        (a2, a1)
    };
    let (blo, bhi) = if key(b1) <= key(b2) {
        (b1, b2)
    } else {
        (b2, b1)
    };
    if key(ahi) < key(blo) {
        Some((ahi, blo))
    } else if key(bhi) < key(alo) {
        Some((bhi, alo))
    } else {
        None
    }
}

/// The least distance between two closed segments (either may be a
/// point): zero when they meet, decided exactly.
pub(crate) fn span_distance(
    (a1, a2): (Point2, Point2),
    (b1, b2): (Point2, Point2),
) -> Result<f64, RouteError> {
    if segments_meet(a1, a2, b1, b2)? {
        return Ok(0.0);
    }
    Ok(point_segment(a1, b1, b2)
        .min(point_segment(a2, b1, b2))
        .min(point_segment(b1, a1, a2))
        .min(point_segment(b2, a1, a2)))
}

fn point_segment(v: Point2, a: Point2, b: Point2) -> f64 {
    let d = b - a;
    let length2 = d.dot(d);
    let s = if length2 > 0.0 {
        ((v - a).dot(d) / length2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (v - (a + d * s)).length()
}

/// A piece of cost edge: part of one polygon's ring, cut at every vertex
/// lying on it, with the polygon and which side its inside lies on.
#[derive(Debug, Clone, Copy)]
struct Piece {
    p: Point2,
    q: Point2,
    polygon: usize,
    /// The polygon's inside lies left of `p` to `q`.
    inside_left: bool,
    /// Index of the line the piece lies on, shared by collinear pieces.
    line: u32,
}

/// The cost polygons, their edges and the lines those lie on.
#[derive(Debug, Clone)]
struct Weights {
    polygons: Vec<Polygon>,
    factors: Vec<f64>,
    pieces: Vec<Piece>,
    /// Two points on each line cost pieces lie on.
    lines: Vec<(Point2, Point2)>,
    /// Region boundary edges, with whether the region lies left of them:
    /// along one, only that side is free.
    walls: Vec<(Point2, Point2, bool)>,
    greatest: f64,
}

impl Weights {
    fn new(costs: &[CostRegion], cuts: &[Point2], region: &[Polygon]) -> Result<Self, RouteError> {
        let mut walls = Vec::new();
        for polygon in region {
            for (hole, ring) in core::iter::once((false, &polygon.outer))
                .chain(polygon.holes.iter().map(|h| (true, h)))
            {
                let left = region_left(ring, hole)?;
                for (p, q) in ring_edges(ring) {
                    walls.push((p, q, left));
                }
            }
        }
        let polygons: Vec<Polygon> = costs.iter().map(|c| c.polygon.clone()).collect();
        let factors: Vec<f64> = costs.iter().map(|c| c.factor).collect();
        let mut vertices: Vec<Point2> = cuts.to_vec();
        for polygon in &polygons {
            for ring in core::iter::once(&polygon.outer).chain(polygon.holes.iter()) {
                vertices.extend(ring.points.iter().copied());
            }
        }
        dedup_points(&mut vertices);
        let mut pieces = Vec::new();
        let mut lines: Vec<(Point2, Point2)> = Vec::new();
        for (index, polygon) in polygons.iter().enumerate() {
            for (hole, ring) in core::iter::once((false, &polygon.outer))
                .chain(polygon.holes.iter().map(|h| (true, h)))
            {
                let left = region_left(ring, hole)?;
                for (p, q) in ring_edges(ring) {
                    if p == q {
                        continue;
                    }
                    let d = q - p;
                    let mut cuts = vec![p, q];
                    for &v in &vertices {
                        if v != p && v != q && within(p, q, v) && side(p, q, v)? == Sign::Zero {
                            cuts.push(v);
                        }
                    }
                    cuts.sort_by(|u, v| (*u - p).dot(d).total_cmp(&(*v - p).dot(d)));
                    let mut line = NO_LINE;
                    for (k, &(r, s)) in lines.iter().enumerate() {
                        if side(r, s, p)? == Sign::Zero && side(r, s, q)? == Sign::Zero {
                            line = k as u32;
                            break;
                        }
                    }
                    if line == NO_LINE {
                        line = lines.len() as u32;
                        lines.push((p, q));
                    }
                    for pair in cuts.windows(2) {
                        if pair[0] != pair[1] {
                            pieces.push(Piece {
                                p: pair[0],
                                q: pair[1],
                                polygon: index,
                                inside_left: left,
                                line,
                            });
                        }
                    }
                }
            }
        }
        let greatest = factors.iter().copied().fold(1.0, f64::max);
        Ok(Self {
            polygons,
            factors,
            pieces,
            lines,
            walls,
            greatest,
        })
    }

    /// Refuse a cost edge that properly crosses an obstacle or another
    /// cost edge, or runs along or through a barrier. A crossing of a
    /// region edge within `reach` of an end of either edge is a touch
    /// rounding moved (see [`touch_walls`]), and stands.
    fn check_crossings(
        &self,
        obstacles: &[(Point2, Point2)],
        walls: &[(Point2, Point2)],
        reach: f64,
    ) -> Result<(), MapError> {
        for piece in &self.pieces {
            let index = piece.polygon;
            for &(r, s) in obstacles {
                if crosses(piece.p, piece.q, r, s)? {
                    let region_edge = self
                        .walls
                        .iter()
                        .any(|&(p, q, _)| (p, q) == (r, s) || (p, q) == (s, r));
                    let near = point_segment(piece.p, r, s) <= reach
                        || point_segment(piece.q, r, s) <= reach
                        || point_segment(r, piece.p, piece.q) <= reach
                        || point_segment(s, piece.p, piece.q) <= reach;
                    if region_edge && near {
                        continue;
                    }
                    return Err(MapError::CostCrossing { index });
                }
            }
            for &(r, s) in walls {
                let collinear =
                    side(r, s, piece.p)? == Sign::Zero && side(r, s, piece.q)? == Sign::Zero;
                if collinear && segments_meet(piece.p, piece.q, r, s)? {
                    let overlap = [piece.p, piece.q, r, s]
                        .iter()
                        .filter(|v| within(r, s, **v) && within(piece.p, piece.q, **v))
                        .count();
                    if overlap >= 2 {
                        return Err(MapError::CostCrossing { index });
                    }
                }
            }
            for other in &self.pieces {
                if crosses(piece.p, piece.q, other.p, other.q)? {
                    return Err(MapError::CostCrossing { index });
                }
            }
        }
        Ok(())
    }

    /// The factors just left and right of the stretch `a`-`b` of the cost
    /// piece `p`-`q`, which no vertex interrupts. Which edges the stretch
    /// runs along is decided on the piece's own exact ends: `a` and `b`
    /// are interpolated, and on an edge at an angle rounding puts them off
    /// its line.
    fn sides(&self, p: Point2, q: Point2, a: Point2, b: Point2) -> Result<(f64, f64), RouteError> {
        let m = Point2::new(0.5 * a.x + 0.5 * b.x, 0.5 * a.y + 0.5 * b.y);
        let d = q - p;
        let (a, b) = (p, q);
        let (mut left, mut right) = (1.0f64, 1.0f64);
        for (index, polygon) in self.polygons.iter().enumerate() {
            let factor = self.factors[index];
            let mut on = false;
            for piece in self.pieces.iter().filter(|p| p.polygon == index) {
                if side(piece.p, piece.q, a)? == Sign::Zero
                    && side(piece.p, piece.q, b)? == Sign::Zero
                    && within(piece.p, piece.q, m)
                {
                    on = true;
                    let same = (piece.q - piece.p).dot(d) > 0.0;
                    if piece.inside_left == same {
                        left = left.max(factor);
                    } else {
                        right = right.max(factor);
                    }
                }
            }
            if !on && in_polygon(polygon, m)? {
                left = left.max(factor);
                right = right.max(factor);
            }
        }
        Ok((left, right))
    }

    /// Whether the stretch `p`-`q` lies within `reach` of one region edge,
    /// on or beyond it: out of the free space, a cost edge touching the
    /// wall (see [`touch_walls`]) that no walk crosses, only runs along.
    fn beyond_wall(&self, p: Point2, q: Point2, reach: f64) -> Result<bool, RouteError> {
        for &(r, s, left) in &self.walls {
            let free = if left { Sign::Positive } else { Sign::Negative };
            if point_segment(p, r, s) <= reach
                && point_segment(q, r, s) <= reach
                && side(r, s, p)? != free
                && side(r, s, q)? != free
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// The lines of cost pieces that pass through `v`.
    fn lines_through(&self, v: Point2) -> Result<Vec<u32>, RouteError> {
        let mut out = Vec::new();
        for (k, &(p, q)) in self.lines.iter().enumerate() {
            if side(p, q, v)? == Sign::Zero {
                out.push(k as u32);
            }
        }
        Ok(out)
    }

    /// The greatest factor over the polygons holding the whole convex hull
    /// of `corners`, and 1 if none does; 1 for a hull with no interior.
    fn hull_factor(&self, corners: &[Point2]) -> Result<f64, RouteError> {
        let Some(hull) = convex_hull(corners)? else {
            return Ok(1.0);
        };
        let mut best = 1.0f64;
        for (polygon, &factor) in self.polygons.iter().zip(&self.factors) {
            if factor <= best {
                continue;
            }
            if holds(polygon, &hull)? {
                best = factor;
            }
        }
        Ok(best)
    }

    /// The greatest factor over the polygons meeting the closed triangle,
    /// and 1: the most a metre can cost in it.
    fn steepest(&self, t: &[Point2; 3]) -> Result<f64, RouteError> {
        let mut best = 1.0f64;
        for (polygon, &factor) in self.polygons.iter().zip(&self.factors) {
            if factor <= best {
                continue;
            }
            let meets = t.iter().try_fold(false, |m, c| {
                Ok::<_, RouteError>(m || in_polygon(polygon, *c)?)
            })? || core::iter::once(&polygon.outer)
                .chain(polygon.holes.iter())
                .flat_map(ring_edges)
                .try_fold(false, |m, (p, q)| {
                    Ok::<_, RouteError>(m || meets_triangle(p, q, t)?)
                })?;
            if meets {
                best = factor;
            }
        }
        Ok(best)
    }

    /// The weighted length of the segment `a`-`b`, bracketed: the factor
    /// is decided piece by piece between the points where the segment
    /// meets cost edges. Along a cost edge the cheaper side counts. A piece
    /// whose factor rounding could misjudge -- shorter than a few ulps, or
    /// whose middle lies within rounding of a cost edge -- counts 1 below
    /// and the greatest factor above, with two exceptions for a piece
    /// within rounding of a cost edge all along. Above, it counts as the
    /// walk along that edge, which is off it by no more than rounding:
    /// the upper end bounds a walk that close to the segment, not the
    /// segment itself. Below, a piece exactly along a region edge, the
    /// cost edge on or beyond it, counts its free side's factor.
    fn segment(&self, a: Point2, b: Point2, scale: f64) -> Result<(f64, f64), RouteError> {
        let length = (b - a).length();
        if length == 0.0 {
            return Ok((0.0, 0.0));
        }
        if self.polygons.is_empty() {
            return Ok((round_down(length, scale), round_up(length, scale)));
        }
        let d = b - a;
        let param = |v: Point2| (v - a).dot(d) / d.dot(d);
        let mut breaks = vec![0.0, 1.0];
        // Stretches along cost pieces: (from, to, polygon, inside left of a-b).
        let mut along: Vec<(f64, f64, usize, bool)> = Vec::new();
        // Pieces not along `a`-`b`: (p, q, polygon, inside left of p-q).
        let mut near: Vec<(Point2, Point2, usize, bool)> = Vec::new();
        // Stretches along region edges: (from, to, region left of a-b).
        let mut walled: Vec<(f64, f64, bool)> = Vec::new();
        // Region edges not along `a`-`b`: (p, q, region left of p-q).
        let mut near_walls: Vec<(Point2, Point2, bool)> = Vec::new();
        for &(p, q, left) in &self.walls {
            if p != q && side(a, b, p)? == Sign::Zero && side(a, b, q)? == Sign::Zero {
                let (tp, tq) = (param(p), param(q));
                let (lo, hi) = (tp.min(tq).max(0.0), tp.max(tq).min(1.0));
                if hi > lo {
                    walled.push((lo, hi, left == ((q - p).dot(d) > 0.0)));
                    breaks.extend([lo, hi]);
                }
            } else if p != q {
                near_walls.push((p, q, left));
            }
        }
        for piece in &self.pieces {
            let (sp, sq) = (side(a, b, piece.p)?, side(a, b, piece.q)?);
            if sp == Sign::Zero && sq == Sign::Zero {
                let (tp, tq) = (param(piece.p), param(piece.q));
                let (lo, hi) = (tp.min(tq).max(0.0), tp.max(tq).min(1.0));
                if hi > lo {
                    let same = (piece.q - piece.p).dot(d) > 0.0;
                    along.push((lo, hi, piece.polygon, piece.inside_left == same));
                    breaks.extend([lo, hi]);
                }
                continue;
            }
            near.push((piece.p, piece.q, piece.polygon, piece.inside_left));
            let (sa, sb) = (side(piece.p, piece.q, a)?, side(piece.p, piece.q, b)?);
            if sp != sq && sa != sb {
                // The lines cross between the ends of both: at the
                // parameter where `a`-`b` meets the piece's line.
                let e = piece.q - piece.p;
                let denom = d.perp_dot(e);
                if denom != 0.0 {
                    let t = (piece.p - a).perp_dot(e) / denom;
                    if t > 0.0 && t < 1.0 {
                        breaks.push(t);
                    }
                }
            }
        }
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        let (mut lower, mut upper) = (0.0, 0.0);
        let tiny = 64.0 * f64::EPSILON * scale;
        let mut hugged = false;
        for pair in breaks.windows(2) {
            let (t0, t1) = (pair[0], pair[1]);
            let stretch = (t1 - t0) * length;
            if stretch <= 0.0 {
                continue;
            }
            let tm = 0.5 * t0 + 0.5 * t1;
            let at = |t: f64| Point2::new(a.x + d.x * t, a.y + d.y * t);
            let m = at(tm);
            // Along cost pieces: polygon and whether its inside lies left
            // of `a`-`b`.
            let mut on: Vec<(usize, bool)> = along
                .iter()
                .filter(|(lo, hi, ..)| *lo <= tm && tm <= *hi)
                .map(|(_, _, p, l)| (*p, *l))
                .collect();
            // Along region edges: whether the region lies left of `a`-`b`.
            let mut free: Vec<bool> = walled
                .iter()
                .filter(|(lo, hi, _)| *lo <= tm && tm <= *hi)
                .map(|(_, _, l)| *l)
                .collect();
            let doubtful = stretch <= tiny
                || near
                    .iter()
                    .any(|&(p, q, ..)| point_segment(m, p, q) <= tiny);
            let mut resolved = !doubtful;
            if doubtful {
                // A stretch within rounding of a cost piece all along --
                // `a`-`b` joins two points interpolated on a cost edge at
                // an angle -- may be on either side of it. The walk along
                // the piece itself, off the stretch by at most `tiny` at
                // each end, costs the piece's along factor: take that one
                // instead, above. Otherwise the greatest factor.
                let (s0, s1) = (at(t0), at(t1));
                let hugs = |p: Point2, q: Point2| {
                    point_segment(s0, p, q) <= tiny && point_segment(s1, p, q) <= tiny
                };
                let mut all_hug = stretch > tiny;
                // Along a region edge, with every hugged piece exactly on
                // or beyond it -- as a cost edge touching a wall is left
                // (see [`touch_walls`]) -- the free side lies on the
                // piece's inside side, and its factor is known below too.
                let beyond = match (free.contains(&true), free.contains(&false)) {
                    (true, false) => Some(Sign::Negative),
                    (false, true) => Some(Sign::Positive),
                    _ => None,
                };
                let mut known = beyond.is_some();
                for &(p, q, polygon, inside_left) in &near {
                    if point_segment(m, p, q) > tiny {
                        continue;
                    }
                    if !hugs(p, q) {
                        all_hug = false;
                        break;
                    }
                    if let Some(solid) = beyond {
                        for end in [p, q] {
                            let s = side(a, b, end)?;
                            known &= s == Sign::Zero || s == solid;
                        }
                    }
                    on.push((polygon, inside_left == ((q - p).dot(d) > 0.0)));
                }
                if !all_hug {
                    lower += stretch;
                    upper += self.greatest * stretch;
                    continue;
                }
                resolved = known;
                if !resolved {
                    lower += stretch;
                }
                for &(p, q, left) in &near_walls {
                    if hugs(p, q) {
                        free.push(left == ((q - p).dot(d) > 0.0));
                    }
                }
                hugged = true;
            }
            // Extra factor on each side of the segment at this stretch.
            let (mut left, mut right) = (0.0f64, 0.0f64);
            for (index, polygon) in self.polygons.iter().enumerate() {
                let extra = self.factors[index] - 1.0;
                if extra <= left.min(right) {
                    continue;
                }
                let sides: Vec<bool> = on
                    .iter()
                    .filter(|(p, _)| *p == index)
                    .map(|(_, l)| *l)
                    .collect();
                if sides.is_empty() {
                    if in_polygon(polygon, m)? {
                        left = left.max(extra);
                        right = right.max(extra);
                    }
                } else {
                    // On the polygon's edge: its inside on one side (both,
                    // for an edge the polygon has twice).
                    if sides.contains(&true) {
                        left = left.max(extra);
                    }
                    if sides.contains(&false) {
                        right = right.max(extra);
                    }
                }
            }
            // Along a cost edge the cheaper side counts -- of the sides a
            // walk can move to: along a wall, only the region's.
            let (free_left, free_right) = if free.is_empty() {
                (true, true)
            } else {
                (free.contains(&true), free.contains(&false))
            };
            let factor = 1.0
                + if on.is_empty() {
                    left
                } else {
                    match (free_left, free_right) {
                        (true, true) => left.min(right),
                        (true, false) => left,
                        (false, true) => right,
                        (false, false) => left.max(right),
                    }
                };
            if resolved {
                lower += factor * stretch;
            }
            upper += factor * stretch;
        }
        if hugged {
            // The steps between the stretch and the piece it hugs.
            upper += 4.0 * self.greatest * tiny * breaks.len() as f64;
        }
        // Each break's parameter is rounded: a stretch may be off by a few
        // ulps of the length at each end.
        let slack = (self.greatest - 1.0) * breaks.len() as f64 * 8.0 * f64::EPSILON * length;
        Ok((
            round_down(lower - slack, scale),
            round_up(upper + slack, scale),
        ))
    }
}

/// How far off a region edge a cost vertex may lie and still be taken to
/// touch it: 2^-24 of the region's extent, four steps of the grid overlay
/// results were once rounded to, and a few ulps of its coordinates.
fn touch_reach(region: &[Polygon]) -> f64 {
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    let mut scale = 0.0f64;
    for polygon in region {
        for p in core::iter::once(&polygon.outer)
            .chain(polygon.holes.iter())
            .flat_map(|r| r.points.iter())
        {
            lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
            hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            scale = scale.max(p.x.abs()).max(p.y.abs());
        }
    }
    let extent = (hi.x - lo.x).max(hi.y - lo.y).max(0.0);
    extent / f64::from(1u32 << 24) + 64.0 * f64::EPSILON * scale
}

/// The cost regions with every vertex that lies strictly on the free side
/// of a region edge and within `reach` of it moved just across that
/// edge's line: a cost region clipped to the free region touches its walls
/// only up to rounding, and a vertex left a hair inside would open a
/// sliver along the wall a walk could slip through at factor 1. Beyond the
/// wall the region costs nothing, so moving out is taking it as touching.
/// A vertex in a corner crosses each wall it is that near to, moving by at
/// most `reach` and a few ulps for each.
fn touch_walls(
    costs: &[CostRegion],
    region: &[Polygon],
    reach: f64,
) -> Result<Vec<CostRegion>, RouteError> {
    let mut edges = Vec::new();
    for polygon in region {
        for (hole, ring) in
            core::iter::once((false, &polygon.outer)).chain(polygon.holes.iter().map(|h| (true, h)))
        {
            let left = region_left(ring, hole)?;
            for (p, q) in ring_edges(ring) {
                if p != q {
                    edges.push((p, q, left));
                }
            }
        }
    }
    let mut out = costs.to_vec();
    for cost in &mut out {
        for ring in core::iter::once(&mut cost.polygon.outer).chain(cost.polygon.holes.iter_mut()) {
            for v in &mut ring.points {
                for &(p, q, left) in edges.iter().chain(edges.iter()) {
                    let free = if left { Sign::Positive } else { Sign::Negative };
                    if point_segment(*v, p, q) > reach || side(p, q, *v)? != free {
                        continue;
                    }
                    // The foot on the line, then out along the normal until
                    // the exact side says it is no longer on the free side.
                    let d = q - p;
                    let t = (*v - p).dot(d) / d.dot(d);
                    let foot = Point2::new(p.x + d.x * t, p.y + d.y * t);
                    let out = if left {
                        Point2::new(d.y, -d.x)
                    } else {
                        Point2::new(-d.y, d.x)
                    };
                    let unit = Point2::new(out.x / d.length(), out.y / d.length());
                    let mut step = f64::EPSILON * foot.x.abs().max(foot.y.abs()).max(reach);
                    let mut moved = foot;
                    while side(p, q, moved)? == free {
                        moved = Point2::new(foot.x + unit.x * step, foot.y + unit.y * step);
                        step *= 2.0;
                    }
                    *v = moved;
                }
            }
        }
    }
    Ok(out)
}

/// The corners of the convex hull of `points`, counter-clockwise, or
/// `None` when they are collinear. Decided exactly.
fn convex_hull(points: &[Point2]) -> Result<Option<Vec<Point2>>, RouteError> {
    let mut pts: Vec<Point2> = points.to_vec();
    pts.sort_by(|p, q| p.x.total_cmp(&q.x).then(p.y.total_cmp(&q.y)));
    pts.dedup();
    if pts.len() < 3 {
        return Ok(None);
    }
    // Andrew's monotone chain with exact turns.
    let mut hull: Vec<Point2> = Vec::new();
    for pass in 0..2 {
        let start = hull.len();
        let iter: Box<dyn Iterator<Item = &Point2>> = if pass == 0 {
            Box::new(pts.iter())
        } else {
            Box::new(pts.iter().rev())
        };
        for &p in iter {
            while hull.len() >= start + 2
                && side(hull[hull.len() - 2], hull[hull.len() - 1], p)? != Sign::Positive
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
    }
    Ok((hull.len() >= 3).then_some(hull))
}

/// Whether the closed polygon holds the convex polygon `hull`
/// (counter-clockwise, with interior). No edge of the polygon may meet the
/// hull's interior, which then lies wholly inside the polygon or wholly
/// outside, so one point strictly inside the hull decides -- corners on
/// the polygon's boundary do not: a hull can span a notch of a polygon,
/// touching it only along two sides. A hull too thin to hold a point
/// provably inside it is not held.
fn holds(polygon: &Polygon, hull: &[Point2]) -> Result<bool, RouteError> {
    for &c in hull {
        if !in_polygon(polygon, c)? {
            return Ok(false);
        }
    }
    let n = hull.len();
    for ring in core::iter::once(&polygon.outer).chain(polygon.holes.iter()) {
        for (p, q) in ring_edges(ring) {
            if meets_open_convex(p, q, hull, n)? {
                return Ok(false);
            }
        }
    }
    let inv = 1.0 / n as f64;
    let centre = hull.iter().fold(Point2::ZERO, |sum, c| sum + *c * inv);
    for i in 0..n {
        if side(hull[i], hull[(i + 1) % n], centre)? != Sign::Positive {
            return Ok(false);
        }
    }
    let on_boundary = core::iter::once(&polygon.outer)
        .chain(polygon.holes.iter())
        .flat_map(ring_edges)
        .try_fold(false, |on, (p, q)| {
            Ok::<_, RouteError>(on || (within(p, q, centre) && side(p, q, centre)? == Sign::Zero))
        })?;
    Ok(!on_boundary && in_polygon(polygon, centre)?)
}

/// Whether the closed segment `p`-`q` meets the open interior of the
/// counter-clockwise convex polygon: no edge line of the polygon, nor the
/// segment's line, separates them.
fn meets_open_convex(p: Point2, q: Point2, hull: &[Point2], n: usize) -> Result<bool, RouteError> {
    for i in 0..n {
        let (u, v) = (hull[i], hull[(i + 1) % n]);
        if side(u, v, p)? != Sign::Positive && side(u, v, q)? != Sign::Positive {
            return Ok(false);
        }
    }
    if p != q {
        let mut all_left = true;
        let mut all_right = true;
        for &c in hull {
            let s = side(p, q, c)?;
            all_left &= s != Sign::Negative;
            all_right &= s != Sign::Positive;
        }
        if all_left || all_right {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The greatest weighted distance from a subregion to the nearest target
/// over its points in the free space, bracketed to within `tolerance`
/// where the map's own bracket allows: the result is never narrower than
/// the gap between the map's bounds at the farthest point.
///
/// # Errors
///
/// As [`crate::farthest_point`].
pub fn weighted_farthest_point(
    map: &WeightedMap,
    subregion: &Polygon,
    tolerance: f64,
) -> Result<Farthest, FarthestError> {
    weighted_farthest_point_within(map, subregion, tolerance, MAX_CELLS)
}

#[derive(Debug, Clone, Copy)]
struct Cell {
    corners: [Point2; 3],
    root: usize,
    anchor: (Point2, f64),
    upper: f64,
    order: usize,
    depth: u32,
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
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper
            .total_cmp(&other.upper)
            .then(other.order.cmp(&self.order))
    }
}

/// [`weighted_farthest_point`] with a caller-chosen cell budget.
///
/// # Errors
///
/// As [`weighted_farthest_point`].
pub fn weighted_farthest_point_within(
    map: &WeightedMap,
    subregion: &Polygon,
    tolerance: f64,
    max_cells: usize,
) -> Result<Farthest, FarthestError> {
    if !(tolerance.is_finite() && tolerance >= 0.0) {
        return Err(FarthestError::InvalidTolerance);
    }
    validate_region(core::slice::from_ref(subregion), &[])?;
    let roots = free_triangles_in(&map.region, &map.obstacles)?;
    let steep: Vec<f64> = roots
        .iter()
        .map(|t| map.weights.steepest(t))
        .collect::<Result<_, _>>()?;
    let slack = |depth: u32, radius: f64| {
        f64::from(depth + 4) * 4.0 * f64::EPSILON * map.scale + 4.0 * f64::EPSILON * radius
    };
    let mut evaluated: HashMap<(u64, u64), Option<(f64, f64)>> = HashMap::new();
    let mut lower: Option<(f64, Point2)> = None;
    let mut heap = BinaryHeap::new();
    let mut cells = 0usize;
    let mut cell = |corners: [Point2; 3],
                    root: usize,
                    inherited: Option<(Point2, f64)>,
                    depth: u32,
                    order: usize,
                    lower: &mut Option<(f64, Point2)>|
     -> Result<Option<Cell>, FarthestError> {
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
        let k = steep[root];
        let mut best: Option<(f64, (Point2, f64))> =
            inherited.map(|(a, d)| (d + k * radius(a), (a, d)));
        for a in [corners[0], corners[1], corners[2], centroid] {
            if !in_triangle(&roots[root], a)? {
                continue;
            }
            if map.walls.iter().try_fold(false, |m, &(p, q)| {
                Ok::<_, RouteError>(m || (within(p, q, a) && side(p, q, a)? == Sign::Zero))
            })? {
                continue;
            }
            let key = (a.x.to_bits(), a.y.to_bits());
            let bracket = match evaluated.get(&key) {
                Some(b) => *b,
                None => {
                    let b = map.bracket(a)?;
                    evaluated.insert(key, b);
                    b
                }
            };
            let Some((lo, hi)) = bracket else {
                return Err(FarthestError::Unreachable {
                    triangle: roots[root],
                });
            };
            if in_polygon(subregion, a)? && lower.is_none_or(|(l, _)| lo > l) {
                *lower = Some((lo, a));
            }
            let bound = hi + k * radius(a);
            if best.is_none_or(|(b, _)| bound < b) {
                best = Some((bound, (a, hi)));
            }
        }
        Ok(best.map(|(bound, anchor)| Cell {
            corners,
            root,
            anchor,
            upper: bound + k * slack(depth, radius(anchor.0)),
            order,
            depth,
        }))
    };
    for (root, corners) in roots.iter().enumerate() {
        if outside(subregion, corners)? {
            continue;
        }
        cells += 1;
        let Some(c) = cell(*corners, root, None, 0, cells, &mut lower)? else {
            return Err(FarthestError::Triangulation);
        };
        heap.push(c);
    }
    if cells == 0 {
        return Err(FarthestError::Empty);
    }
    let best = |lower: &Option<(f64, Point2)>| lower.map_or(f64::NEG_INFINITY, |(d, _)| d);
    while let Some(top) = heap.peek() {
        if top.upper <= best(&lower) + tolerance || cells >= max_cells {
            break;
        }
        let parent = heap.pop().expect("peeked");
        let [a, b, c] = parent.corners;
        let mid = |p: Point2, q: Point2| Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y);
        let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
        for corners in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
            if outside(subregion, &corners)? {
                continue;
            }
            cells += 1;
            if let Some(child) = cell(
                corners,
                parent.root,
                Some(parent.anchor),
                parent.depth + 1,
                cells,
                &mut lower,
            )? {
                if child.upper > best(&lower) {
                    heap.push(child);
                }
            }
        }
    }
    let Some((lo, witness)) = lower else {
        if heap.is_empty() {
            return Err(FarthestError::Empty);
        }
        return Ok(Farthest {
            distance: LengthInterval {
                lower: 0.0,
                upper: heap.peek().map_or(0.0, |c| c.upper),
            },
            witness: None,
            converged: false,
            cells,
        });
    };
    let hi = heap.peek().map_or(lo, |c| c.upper.max(lo));
    Ok(Farthest {
        distance: LengthInterval {
            lower: lo,
            upper: hi,
        },
        witness: Some(witness),
        converged: hi - lo <= tolerance,
        cells,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_overlay::Ring;

    fn square() -> Weights {
        let p = |x: f64, y: f64| Point2::new(x, y);
        let polygon = Polygon {
            outer: Ring {
                points: vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)],
            },
            holes: Vec::new(),
        };
        Weights::new(&[CostRegion::new(polygon, 2.0)], &[], &[]).unwrap()
    }

    /// A segment a rounding step off a cost edge cannot be told inside
    /// from outside by its middle: the stretch counts 1 below and the
    /// greatest factor above, whichever side it is on.
    #[test]
    fn a_stretch_within_rounding_of_a_cost_edge_is_bracketed_both_ways() {
        let weights = square();
        let tiny = 1e-17;
        // Just inside the square along its lower edge: 1 + 2 + 1 below.
        // Above, the walk along the edge itself, a rounding off the
        // segment, where the cheaper side counts: 3.
        let (lo, hi) = weights
            .segment(Point2::new(-1.0, tiny), Point2::new(2.0, tiny), 2.0)
            .unwrap();
        assert!(lo <= 3.0 && (hi - 3.0).abs() < 1e-12, "[{lo}, {hi}]");
        // Just outside it: 3 at factor 1.
        let (lo, hi) = weights
            .segment(Point2::new(-1.0, -tiny), Point2::new(2.0, -tiny), 2.0)
            .unwrap();
        assert!(lo <= 3.0 && hi >= 3.0, "[{lo}, {hi}]");
        // Inside a notched polygon, the notch's tip a rounding above the
        // middle: doubtful there, and hugging nothing, so the greatest
        // factor above, which is the true one.
        let p = |x: f64, y: f64| Point2::new(x, y);
        let notched = Polygon {
            outer: Ring {
                points: vec![
                    p(0.0, 0.0),
                    p(4.0, 0.0),
                    p(4.0, 2.0),
                    p(2.5, 2.0),
                    p(2.0, 1.0 + f64::EPSILON),
                    p(1.5, 2.0),
                    p(0.0, 2.0),
                ],
            },
            holes: Vec::new(),
        };
        let weights = Weights::new(&[CostRegion::new(notched, 2.0)], &[], &[]).unwrap();
        let (lo, hi) = weights.segment(p(0.5, 1.0), p(3.5, 1.0), 4.0).unwrap();
        assert!(lo <= 6.0 && hi >= 6.0, "[{lo}, {hi}]");
        let weights = square();
        // Crossing the edge within rounding of it, not along it: the
        // greatest factor above.
        let (lo, hi) = weights
            .segment(Point2::new(0.5, -tiny), Point2::new(0.5, tiny), 2.0)
            .unwrap();
        assert!(lo <= 2.0 * tiny && hi >= 2.0 * 2.0 * tiny, "[{lo}, {hi}]");
        // Clear of every edge the stretches are decided.
        let (lo, hi) = weights
            .segment(Point2::new(-1.0, 0.5), Point2::new(2.0, 0.5), 2.0)
            .unwrap();
        assert!(
            (lo - 4.0).abs() < 1e-12 && (hi - 4.0).abs() < 1e-12,
            "[{lo}, {hi}]"
        );
    }

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
        Polygon {
            outer: Ring {
                points: vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)],
            },
            holes: Vec::new(),
        }
    }

    /// The interval of `map` whose span has `end` as one end and runs
    /// toward `toward`.
    fn interval(map: &WeightedMap, end: Point2, toward: Point2) -> usize {
        (0..map.nodes.len())
            .find(|&i| match map.kinds[i] {
                Kind::Interval { a, b, .. } => {
                    (a == end && (b - a).dot(toward - a) > 0.0)
                        || (b == end && (a - b).dot(toward - b) > 0.0)
                }
                Kind::Vertex => false,
            })
            .expect("an interval")
    }

    /// A barrier's foot on the segment from an interval's inner end: the
    /// pieces from just beside that end pass under the foot, so the hop
    /// stands. From an end that is a vertex, the same touch is left to the
    /// vertex's own node, and the hop is blocked when every other piece
    /// crosses the barrier.
    #[test]
    fn a_touch_blocks_a_hop_only_at_a_vertex_end() {
        let room = [rect(0.0, 0.0, 10.0, 4.0)];
        let stair = CostRegion::new(rect(2.0, 0.0, 4.0, 4.0), 2.0);
        // The stair's east edge x = 4 is cut at y = 1 (spacing 1).
        let barrier = vec![vec![p(5.0, 1.0), p(5.0, 3.0)]];
        let map = weighted_distance_map(&room, &barrier, &[p(0.5, 2.0)], &[stair], 1.0).unwrap();
        let above = interval(&map, p(4.0, 1.0), p(4.0, 2.0));
        let q = p(6.0, 1.0);
        let hop = map
            .hop(map.span(above), (q, q), map.kinds[above], Kind::Vertex)
            .unwrap();
        assert!(hop.is_some(), "the inner end's pieces graze the foot");
        // The same geometry with the foot's end a vertex: the stair's
        // corner (4, 4) and a barrier foot on the line to the query.
        // An end that is a vertex: the stair's corner (4, 4), with a
        // barrier's foot on the line from it to the query.
        let walled = vec![vec![p(5.0, 3.0), p(5.0, 1.0)]];
        let stair = CostRegion::new(rect(2.0, 0.0, 4.0, 4.0), 2.0);
        let map = weighted_distance_map(&room, &walled, &[p(0.5, 2.0)], &[stair], 1.0).unwrap();
        let corner = interval(&map, p(4.0, 4.0), p(4.0, 3.0));
        let q = p(6.0, 2.0);
        // Segments from (4, y), y just under 4, to (6, 2) cross x = 5 just
        // under y = 3: through the barrier; from the corner itself, at its
        // end (5, 3).
        let hop = map
            .hop(map.span(corner), (q, q), map.kinds[corner], Kind::Vertex)
            .unwrap();
        assert!(hop.is_none(), "{hop:?}");
    }

    /// A notch poking into a hull through its side: every corner and the
    /// centre of the hull are inside the polygon, but part of the hull is
    /// not, so it is not held.
    #[test]
    fn a_hull_is_held_only_whole() {
        let square = rect(0.0, 0.0, 4.0, 4.0);
        let hull = [p(1.0, 1.0), p(3.0, 1.0), p(3.0, 3.0)];
        assert!(holds(&square, &hull).unwrap());
        // A notch cut in from the west side, its tip at (2.2, 1.6), inside
        // the hull, below its diagonal.
        let notched = Polygon {
            outer: Ring {
                points: vec![
                    p(0.0, 0.0),
                    p(4.0, 0.0),
                    p(4.0, 4.0),
                    p(0.0, 4.0),
                    p(0.0, 1.9),
                    p(2.2, 1.6),
                    p(0.0, 1.5),
                ],
            },
            holes: Vec::new(),
        };
        for c in hull {
            assert!(in_polygon(&notched, c).unwrap());
        }
        assert!(in_polygon(&notched, p(7.0 / 3.0, 5.0 / 3.0)).unwrap());
        assert!(!holds(&notched, &hull).unwrap());
    }

    /// A segment along a blocker's own line, through its end, runs along
    /// it: no touch that separates anything.
    #[test]
    fn a_segment_along_a_blocker_does_not_pass_through_its_end() {
        assert!(!through_end(p(0.0, 0.0), p(4.0, 0.0), p(2.0, 0.0), p(3.0, 0.0)).unwrap());
        assert!(!through_end(p(0.0, 0.0), p(4.0, 0.0), p(2.0, 0.0), p(1.0, 0.0)).unwrap());
        // Across the blocker's line, through its end.
        assert!(through_end(p(0.0, -1.0), p(0.0, 1.0), p(0.0, 0.0), p(1.0, 0.0)).unwrap());
    }
}
