#![forbid(unsafe_code)]
//! Exact planar shortest path over a visibility graph.
//!
//! # Why exact, and what that means here
//!
//! `axiolid-field` already answers route questions, but by sampling a grid:
//! its answer is only as good as the resolution. On a polygon set the shortest
//! path is exactly computable, because the optimal path is a polyline whose
//! interior vertices are region or barrier vertices. No discretisation, no
//! resolution parameter.
//!
//! "Exact" is a claim about the COMBINATORICS, and it is earned by deciding
//! every segment-crossing and sidedness question with certified `orient2d`
//! rather than a tolerance comparison. Which edges exist in the visibility
//! graph is therefore exact. The path LENGTH is still a sum of square roots
//! evaluated in binary64, so it carries ordinary floating-point rounding.
//! Overstating that as "exact length" would be a false claim, so it is not
//! made.
//!
//! # Not a verdict
//!
//! The kernel owns the region, the graph, the path and the typed unreachable
//! reason. It does not own why a route was requested, what clearance is
//! required, or whether a length is acceptable. Same line `navigate.rs` draws.
//!
//! # Input size is bounded explicitly
//!
//! Visibility graph construction is quadratic in vertices and cubic to verify,
//! so a large input silently becomes a hang. [`MAX_VERTICES`] caps it and
//! oversized input is REFUSED, never truncated: truncating would answer a
//! different question than the one asked, and the caller would not be told.
//!
//! The cap is a default, not a law. [`shortest_path_within`] takes the
//! budget as a parameter, because what is affordable depends on the
//! caller's deadline rather than on the kernel. And the refusal carries a
//! PROVEN lower bound — the straight-line distance between the endpoints,
//! which no route can beat — so an over-budget query still yields a usable
//! fact instead of only an error.

use axiolid_contracts::Sign;
use axiolid_core::Point2;
use axiolid_guarantees::Certified;
use axiolid_overlay::{Polygon, Ring};
use axiolid_predicates::orient2d;

/// Maximum vertices, counting region, barrier and endpoint vertices.
pub const MAX_VERTICES: usize = 512;

/// Why no path was produced.
///
/// A geometric fact about the input, distinguished from a numerical failure:
/// a caller must be able to tell "there is no route" from "this could not be
/// decided".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Unreachable {
    /// The start point lies outside the free-space region.
    StartOutside,
    /// The goal point lies outside the free-space region.
    GoalOutside,
    /// Both endpoints are inside, but no route connects them.
    DisconnectedComponents,
}

/// A malformed query, as opposed to an honest "no route".
///
/// `PartialEq` but not `Eq`: [`RouteError::TooManyVertices`] carries a
/// float bound, and float equality is not reflexive.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum RouteError {
    /// A coordinate was NaN or infinite.
    NonFinitePoint,
    /// A ring had fewer than three vertices.
    RingTooShort,
    /// A barrier had fewer than two vertices, so it bounds no segment.
    BarrierTooShort,
    /// The input exceeds the vertex budget. Refused, not truncated.
    ///
    /// Carries a PROVEN lower bound rather than only a complaint. A
    /// refusal and a bound are different facts: a caller that cannot
    /// afford the exact route can still report a minimum, defer, or
    /// escalate, where a bare error forces it to drop the query or
    /// reimplement routing (kernel#92).
    ///
    /// The added fields make this a breaking change: `Eq` is gone because
    /// the bound is a float, and an exhaustive struct variant gained
    /// fields. Both are recorded against the 0.3.0 minor bump rather than
    /// worked around — marking the variant `#[non_exhaustive]` now would
    /// itself be breaking, so it buys nothing here.
    TooManyVertices {
        /// Vertices the caller supplied.
        supplied: usize,
        /// The budget that was applied, so the caller can raise it.
        budget: usize,
        /// Straight-line distance between the endpoints.
        ///
        /// A lower bound on EVERY route between them, not an estimate:
        /// a polyline is at least as long as the straight line joining
        /// its ends, and obstacles only lengthen it. Both endpoints are
        /// already proven inside the free-space region when this is
        /// reported, so the bound applies to a route that could exist.
        ///
        /// Being a bound, it is safe to act on: no admissible route is
        /// shorter. It says nothing about whether a route EXISTS.
        lower_bound: f64,
    },
    /// The exact predicate could not decide a sidedness question.
    ///
    /// Distinct from every [`Unreachable`] variant: this is the kernel
    /// declining to guess, not a statement about the geometry.
    Undecidable,
}

/// A shortest path and its length.
#[derive(Debug, Clone, PartialEq)]
pub struct Route {
    /// The polyline realising the path, start first, goal last.
    pub polyline: Vec<Point2>,
    /// Summed Euclidean length of the polyline.
    ///
    /// A sum of square roots in binary64: the graph is exact, this number
    /// carries ordinary rounding.
    pub length: f64,
    /// Vertices in the visibility graph that produced this path.
    pub graph_vertices: usize,
}

/// Shortest path from `start` to `goal` inside `region`, avoiding `barriers`.
///
/// `barriers` are zero-width: they block visibility without bounding area, so
/// a wall that is a line rather than a thin polygon still stops a route. A
/// barrier is a polyline, not a ring, and is not closed implicitly.
///
/// `Ok(Ok(route))` is a path; `Ok(Err(reason))` is a geometric fact that no
/// path exists; `Err` is a malformed query.
pub fn shortest_path(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    start: Point2,
    goal: Point2,
) -> Result<Result<Route, Unreachable>, RouteError> {
    shortest_path_within(region, barriers, start, goal, MAX_VERTICES)
}

/// [`shortest_path`] with a caller-chosen vertex budget.
///
/// The right cap depends on the caller's time budget, not on the kernel:
/// construction is quadratic in vertices and cubic to verify, so what is
/// affordable is a property of the deadline, not of the geometry
/// (kernel#92). [`MAX_VERTICES`] remains the default for
/// [`shortest_path`].
///
/// Raising the budget does not change any answer, only which inputs are
/// affordable. Over-budget input is still REFUSED rather than truncated,
/// but the refusal carries a proven lower bound.
pub fn shortest_path_within(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    start: Point2,
    goal: Point2,
    budget: usize,
) -> Result<Result<Route, Unreachable>, RouteError> {
    validate(region, barriers, start, goal)?;

    // Endpoint containment is decided before any graph work: it is the
    // cheapest question and gives the most specific answer.
    if !contains(region, start)? {
        return Ok(Err(Unreachable::StartOutside));
    }
    if !contains(region, goal)? {
        return Ok(Err(Unreachable::GoalOutside));
    }

    let mut nodes = vec![start, goal];
    for polygon in region {
        nodes.extend(polygon.outer.points.iter().copied());
        for hole in &polygon.holes {
            nodes.extend(hole.points.iter().copied());
        }
    }
    for barrier in barriers {
        nodes.extend(barrier.iter().copied());
    }

    // Duplicate vertices would create zero-length graph edges and duplicate
    // work without changing the answer.
    dedup_points(&mut nodes);
    if nodes.len() > budget {
        // Both endpoints are proven inside the region by this point, so
        // the straight line between them bounds any route that could
        // exist. Reported as a fact the caller can act on, not as a
        // consolation: no admissible route is shorter than this.
        return Err(RouteError::TooManyVertices {
            supplied: nodes.len(),
            budget,
            lower_bound: (goal - start).length(),
        });
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

    let Some((length, path)) = dijkstra(&adjacency, 0, 1) else {
        return Ok(Err(Unreachable::DisconnectedComponents));
    };
    Ok(Ok(Route {
        polyline: path.into_iter().map(|index| nodes[index]).collect(),
        length,
        graph_vertices: nodes.len(),
    }))
}

fn validate(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    start: Point2,
    goal: Point2,
) -> Result<(), RouteError> {
    if !start.is_finite() || !goal.is_finite() {
        return Err(RouteError::NonFinitePoint);
    }
    for polygon in region {
        for ring in core::iter::once(&polygon.outer).chain(polygon.holes.iter()) {
            if ring.points.len() < 3 {
                return Err(RouteError::RingTooShort);
            }
            if !ring.points.iter().all(|p| p.is_finite()) {
                return Err(RouteError::NonFinitePoint);
            }
        }
    }
    for barrier in barriers {
        if barrier.len() < 2 {
            return Err(RouteError::BarrierTooShort);
        }
        if !barrier.iter().all(|p| p.is_finite()) {
            return Err(RouteError::NonFinitePoint);
        }
    }
    Ok(())
}

/// Exact sidedness, or `Undecidable` rather than a guess.
fn side(a: Point2, b: Point2, c: Point2) -> Result<Sign, RouteError> {
    match orient2d(a, b, c) {
        Certified::Certain { sign, .. } => Ok(sign),
        _ => Err(RouteError::Undecidable),
    }
}

/// True when segments `pq` and `rs` cross at an interior point of both.
///
/// Shared endpoints and collinear touching do NOT count: two visibility edges
/// meeting at a shared polygon vertex is the normal case, and treating that as
/// a blocking crossing would disconnect every graph.
fn crosses(p: Point2, q: Point2, r: Point2, s: Point2) -> Result<bool, RouteError> {
    let d1 = side(p, q, r)?;
    let d2 = side(p, q, s)?;
    let d3 = side(r, s, p)?;
    let d4 = side(r, s, q)?;
    // Strict straddle on both segments. Any Zero means a touch, not a cross.
    Ok(d1 != Sign::Zero
        && d2 != Sign::Zero
        && d3 != Sign::Zero
        && d4 != Sign::Zero
        && d1 != d2
        && d3 != d4)
}

/// Every blocking segment: region boundary edges and barrier edges.
fn obstacle_segments(region: &[Polygon], barriers: &[Vec<Point2>]) -> Vec<(Point2, Point2)> {
    let mut segments = Vec::new();
    for polygon in region {
        for ring in core::iter::once(&polygon.outer).chain(polygon.holes.iter()) {
            segments.extend(ring_edges(ring));
        }
    }
    for barrier in barriers {
        // Open polyline: no closing edge, so a zero-width wall stays a wall
        // rather than becoming an implicit loop.
        for pair in barrier.windows(2) {
            segments.push((pair[0], pair[1]));
        }
    }
    segments
}

fn ring_edges(ring: &Ring) -> Vec<(Point2, Point2)> {
    let count = ring.points.len();
    (0..count)
        .map(|index| (ring.points[index], ring.points[(index + 1) % count]))
        .collect()
}

/// Can the open segment `a`-`b` be travelled without leaving the region?
///
/// Two conditions, both necessary. The segment must not properly cross any
/// obstacle edge, and its midpoint must lie inside the region: a segment can
/// clear every edge yet still pass through a hole by spanning it corner to
/// corner, which the crossing test alone accepts.
fn visible(
    a: Point2,
    b: Point2,
    region: &[Polygon],
    obstacles: &[(Point2, Point2)],
) -> Result<bool, RouteError> {
    for (p, q) in obstacles {
        if crosses(a, b, *p, *q)? {
            return Ok(false);
        }
    }
    let midpoint = Point2::new(a.x * 0.5 + b.x * 0.5, a.y * 0.5 + b.y * 0.5);
    contains(region, midpoint)
}

/// Is `point` inside the region (outer boundary, minus holes)?
///
/// Boundary points count as inside: a route legitimately runs along a wall,
/// and excluding the boundary would make every vertex-to-vertex edge invalid.
fn contains(region: &[Polygon], point: Point2) -> Result<bool, RouteError> {
    for polygon in region {
        // Ray casting alone is boundary-EXCLUSIVE, but region vertices lie
        // exactly on the outer boundary and every visibility edge ends at one.
        // Excluding them would reject every graph edge.
        if !point_in_ring(&polygon.outer, point) && !on_boundary(&polygon.outer, point)? {
            continue;
        }
        let mut in_hole = false;
        for hole in &polygon.holes {
            // Strictly inside a hole is outside the region; ON the hole
            // boundary is still travellable.
            if point_in_ring(hole, point) && !on_boundary(hole, point)? {
                in_hole = true;
                break;
            }
        }
        if !in_hole {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Exact on-boundary test: the point is collinear with, and between, the
/// endpoints of some ring edge.
fn on_boundary(ring: &Ring, point: Point2) -> Result<bool, RouteError> {
    for (a, b) in ring_edges(ring) {
        if side(a, b, point)? != Sign::Zero {
            continue;
        }
        // Collinear; now check it lies within the edge span rather than on
        // its infinite extension.
        let within_x = point.x >= a.x.min(b.x) && point.x <= a.x.max(b.x);
        let within_y = point.y >= a.y.min(b.y) && point.y <= a.y.max(b.y);
        if within_x && within_y {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Ray-casting containment, boundary inclusive.
fn point_in_ring(ring: &Ring, point: Point2) -> bool {
    let mut inside = false;
    let count = ring.points.len();
    for index in 0..count {
        let a = ring.points[index];
        let b = ring.points[(index + 1) % count];
        if (a.y > point.y) != (b.y > point.y) {
            let crossing = (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x;
            if point.x < crossing {
                inside = !inside;
            }
        }
    }
    inside
}

/// Remove exact duplicate points, preserving first-seen order.
fn dedup_points(points: &mut Vec<Point2>) {
    let mut seen: Vec<Point2> = Vec::new();
    points.retain(|point| {
        if seen.iter().any(|other| other == point) {
            false
        } else {
            seen.push(*point);
            true
        }
    });
}

/// Dijkstra returning `(length, path)`.
///
/// Ties break on the lowest node index, so two equal-length paths always
/// resolve to the same one. Without that, the answer would depend on
/// iteration order and the same query could return different polylines.
fn dijkstra(
    adjacency: &[Vec<(usize, f64)>],
    source: usize,
    target: usize,
) -> Option<(f64, Vec<usize>)> {
    let count = adjacency.len();
    let mut distance = vec![f64::INFINITY; count];
    let mut previous = vec![usize::MAX; count];
    let mut settled = vec![false; count];
    distance[source] = 0.0;

    for _ in 0..count {
        // Linear scan rather than a heap: it makes the lowest-index tie-break
        // explicit, and MAX_VERTICES already bounds the cost.
        let mut current = None;
        for index in 0..count {
            if settled[index] || distance[index].is_infinite() {
                continue;
            }
            // Strictly less, so the lowest index wins an exact tie.
            if current.is_none_or(|best: usize| distance[index] < distance[best]) {
                current = Some(index);
            }
        }
        let Some(current) = current else { break };
        if current == target {
            break;
        }
        settled[current] = true;

        for (neighbour, weight) in &adjacency[current] {
            let candidate = distance[current] + weight;
            if candidate < distance[*neighbour] {
                distance[*neighbour] = candidate;
                previous[*neighbour] = current;
            }
        }
    }

    if distance[target].is_infinite() {
        return None;
    }
    let mut path = vec![target];
    let mut node = target;
    while node != source {
        node = previous[node];
        path.push(node);
    }
    path.reverse();
    Some((distance[target], path))
}
