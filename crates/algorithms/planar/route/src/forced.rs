//! The shortest walk forced through a region (#196).
//!
//! A walk from an origin to a target that visits a point `p` is at least
//! `d_origin(p) + d_targets(p)` long, and the best such walk is exactly
//! that long. So the shortest walk that enters a polygon has length
//!
//! ```text
//! W = min over p in the polygon (and the free space) of  d_origin(p) + d_targets(p).
//! ```
//!
//! If `W` exceeds the shortest walk overall, no shortest walk touches the
//! polygon; if `W` equals it, one does.
//!
//! # The bracket
//!
//! The same branch and bound as [`crate::farthest_point`], turned round to
//! find a minimum. Both distances are 1-Lipschitz along segments in the free
//! space, so their sum `f` is 2-Lipschitz: on a triangle of the free space
//! with an anchor `a`,
//!
//! ```text
//! min over T of f  >=  f(a) - 2 (greatest distance from a to a corner of T).
//! ```
//!
//! That bound is weak where `f` is least along a whole stretch of path, as
//! it is on every straight run of the shortest walk through the polygon: it
//! would need cells as small as the tolerance all along the run. So a cell
//! is also bounded through the maps themselves. A distance at `y` is
//! `|y - v| + D(v)` for some graph vertex `v` that `y` sees, so for any
//! vertices `u` of one map and `v` of the other,
//!
//! ```text
//! f(y) >= D(u) + D'(v) + max(|u - v|, dist(T, u) + dist(T, v)),
//! ```
//!
//! and the least of these over the vertices the cell might see bounds `f`
//! on the cell. A vertex is left out only when it certainly sees no point
//! of the cell: one obstacle edge crosses every segment from it to the
//! cell, or the cell lies strictly inside a solid corner at the vertex. On a cell the walk passes straight
//! through, the pair of vertices it runs between gives `|u - v|` exactly.
//!
//! Every value of `f` at a point of the polygon is an upper bound on `W`,
//! and no walk from an origin to a target is shorter than the shortest one,
//! `L`, so `W >= L` bounds every cell from below as well. Cells are split,
//! the one with the least lower bound first, and dropped once their lower
//! bound reaches the best upper one. Rounding widens the result as in
//! [`crate::farthest_point`].

use core::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};

use axiolid_contracts::Sign;
use axiolid_core::Point2;
use axiolid_overlay::Polygon;

use crate::map::{admissible, free_triangles, meets_triangle, outside};
use crate::{
    crosses, side, validate_region, DistanceMap, FarthestError, LengthInterval, RouteError,
    MAX_CELLS,
};

/// The shortest walk from an origin to a target that enters a polygon.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct ForcedWalk {
    /// Contains the length of the shortest walk that enters the polygon.
    /// Both ends are infinite when no walk from an origin to a target
    /// reaches the polygon.
    pub length: LengthInterval,
    /// The shortest walk from any origin to any target, ignoring the
    /// polygon: `length.lower` is never below it, up to rounding. Infinite
    /// when no target is reachable from an origin.
    pub shortest: f64,
    /// A point of the polygon, in the free space, through which a walk of
    /// at most `length.upper` passes; `None` when no walk reaches the
    /// polygon.
    pub witness: Option<Point2>,
    /// Whether the interval is no wider than the tolerance asked for. When
    /// the cell budget runs out first it is still sound, only wider.
    pub converged: bool,
    /// Cells examined.
    pub cells: usize,
}

/// The shortest walk from an origin of `from` to a target of `to` that
/// enters `through`, bracketed to within `tolerance`.
///
/// `from` is a distance map whose targets are the walk's origins, `to` one
/// whose targets are its destinations; both must cover the same region and
/// barriers.
///
/// # Errors
///
/// [`FarthestError::MismatchedMaps`] when the maps cover different free
/// space; otherwise as [`crate::farthest_point`], except that a part of the
/// polygon no walk reaches is not an error: it only contributes nothing.
pub fn forced_walk(
    from: &DistanceMap,
    to: &DistanceMap,
    through: &Polygon,
    tolerance: f64,
) -> Result<ForcedWalk, FarthestError> {
    forced_walk_within(from, to, through, tolerance, MAX_CELLS)
}

/// A cell: a triangle inside the free-space triangle `root`, with the best
/// anchor known for it and the lower bound that gives.
#[derive(Debug, Clone, Copy)]
struct Cell {
    corners: [Point2; 3],
    root: usize,
    anchor: (Point2, f64),
    lower: f64,
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
    /// Least lower bound first (the heap is a max-heap); the earlier cell
    /// on a tie.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .lower
            .total_cmp(&self.lower)
            .then(other.order.cmp(&self.order))
    }
}

/// [`forced_walk`] with a caller-chosen cell budget.
///
/// # Errors
///
/// As [`forced_walk`].
pub fn forced_walk_within(
    from: &DistanceMap,
    to: &DistanceMap,
    through: &Polygon,
    tolerance: f64,
    max_cells: usize,
) -> Result<ForcedWalk, FarthestError> {
    if !(tolerance.is_finite() && tolerance >= 0.0) {
        return Err(FarthestError::InvalidTolerance);
    }
    if !from.same_space(to) {
        return Err(FarthestError::MismatchedMaps);
    }
    validate_region(core::slice::from_ref(through), &[])?;
    let roots = free_triangles(from)?;
    // Each sum carries the roundings of both maps' lengths.
    let relative = (from.hops() as f64 + to.hops() as f64 + 16.0) * 2.0 * f64::EPSILON;
    let scale = from
        .nodes()
        .iter()
        .chain(through.outer.points.iter())
        .fold(0.0f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));
    // As for the farthest point, doubled: two distances move with a point.
    let slack = |depth: u32, radius: f64| {
        2.0 * (f64::from(depth + 4) * 4.0 * f64::EPSILON * scale + 4.0 * f64::EPSILON * radius)
    };
    // The shortest walk overall: no walk through the polygon beats it.
    let mut shortest = f64::INFINITY;
    for (&origin, &weight) in from.sites().iter().zip(from.weights()) {
        if let Some(d) = to.at(origin)? {
            shortest = shortest.min(weight + d);
        }
    }
    let floor = shortest * (1.0 - relative);
    let vertices = |map: &DistanceMap| -> Result<Vec<Vertex>, RouteError> {
        map.nodes()
            .iter()
            .copied()
            .zip(map.vertex_distances())
            .filter(|(_, d)| d.is_finite())
            .map(|(at, distance)| {
                Ok(Vertex {
                    at,
                    distance,
                    solid: solid_corner(map.region(), map.obstacles(), at)?,
                })
            })
            .collect()
    };
    let mut search = Search {
        from,
        to,
        through,
        from_vertices: vertices(from)?,
        to_vertices: vertices(to)?,
        roots: &roots,
        evaluated: HashMap::new(),
        upper: None,
    };
    let mut heap = BinaryHeap::new();
    let mut cells = 0usize;
    let mut met = false;
    for (root, corners) in roots.iter().enumerate() {
        if outside(through, corners)? {
            continue;
        }
        met = true;
        cells += 1;
        match search.cell(*corners, root, None, 0, cells, floor, &slack)? {
            Anchored::Cell(cell) => heap.push(cell),
            // No walk reaches this part of the free space.
            Anchored::Unreached => {}
            Anchored::None => return Err(FarthestError::Triangulation),
        }
    }
    if !met {
        return Err(FarthestError::Empty);
    }
    let best = |search: &Search| search.upper.map_or(f64::INFINITY, |(d, _)| d);
    while let Some(top) = heap.peek() {
        if top.lower >= best(&search) - tolerance || cells >= max_cells {
            break;
        }
        let cell = heap.pop().expect("peeked");
        let [a, b, c] = cell.corners;
        let mid = |p: Point2, q: Point2| Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y);
        let (ab, bc, ca) = (mid(a, b), mid(b, c), mid(c, a));
        for corners in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
            if outside(through, &corners)? {
                continue;
            }
            cells += 1;
            let child = search.cell(
                corners,
                cell.root,
                Some(cell.anchor),
                cell.depth + 1,
                cells,
                floor,
                &slack,
            )?;
            if let Anchored::Cell(child) = child {
                if child.lower < best(&search) {
                    heap.push(child);
                }
            }
        }
    }
    let Some((upper, witness)) = search.upper else {
        if heap.is_empty() {
            // Every part of the polygon in the free space is unreached.
            return Ok(ForcedWalk {
                length: LengthInterval {
                    lower: f64::INFINITY,
                    upper: f64::INFINITY,
                },
                shortest,
                witness: None,
                converged: true,
                cells,
            });
        }
        let lower = heap.peek().map_or(floor, |c| c.lower.max(floor));
        return Ok(ForcedWalk {
            length: LengthInterval {
                lower: lower * (1.0 - relative),
                upper: f64::INFINITY,
            },
            shortest,
            witness: None,
            converged: false,
            cells,
        });
    };
    let lower = heap.peek().map_or(upper, |c| c.lower.min(upper)).max(floor);
    Ok(ForcedWalk {
        length: LengthInterval {
            lower: (lower * (1.0 - relative)).min(upper),
            upper: upper * (1.0 + relative),
        },
        shortest,
        witness: Some(witness),
        converged: upper - lower <= tolerance,
        cells,
    })
}

enum Anchored {
    Cell(Cell),
    /// The cell's root triangle is unreached from an origin or a target.
    Unreached,
    /// No anchor lies in the cell (never for a root, unless the triangle
    /// holds no representable interior point).
    None,
}

struct Search<'a> {
    from: &'a DistanceMap,
    to: &'a DistanceMap,
    through: &'a Polygon,
    /// Graph vertices a target is reached from.
    from_vertices: Vec<Vertex>,
    to_vertices: Vec<Vertex>,
    roots: &'a [[Point2; 3]],
    evaluated: HashMap<(u64, u64), Option<f64>>,
    /// The least walk length at a point of the polygon, and the point.
    upper: Option<(f64, Point2)>,
}

impl Search<'_> {
    /// The least of `D(u) + D'(v) + max(|u - v|, dist(T, u) + dist(T, v))`
    /// over the vertices `u` of `from` and `v` of `to` that the cell might
    /// see; infinite when it sees none of one map's.
    fn pair_bound(&self, t: &[Point2; 3]) -> Result<f64, RouteError> {
        let origins = candidates(&self.from_vertices, t, self.from.obstacles())?;
        let targets = candidates(&self.to_vertices, t, self.to.obstacles())?;
        let Some(&(least, ..)) = targets.first() else {
            return Ok(f64::INFINITY);
        };
        let mut best = f64::INFINITY;
        for &(s, u, du, gu) in &origins {
            if s + least >= best {
                break;
            }
            for &(r, v, dv, gv) in &targets {
                if s + r >= best {
                    break;
                }
                best = best.min(du + dv + (u - v).length().max(gu + gv));
            }
        }
        Ok(best)
    }

    /// `d_origin(p) + d_targets(p)`, or `None` when either is unreachable.
    fn walk(&mut self, p: Point2) -> Result<Option<f64>, FarthestError> {
        let key = (p.x.to_bits(), p.y.to_bits());
        if let Some(d) = self.evaluated.get(&key) {
            return Ok(*d);
        }
        let d = match (self.from.at(p)?, self.to.at(p)?) {
            (Some(a), Some(b)) => Some(a + b),
            _ => None,
        };
        self.evaluated.insert(key, d);
        Ok(d)
    }

    #[allow(clippy::too_many_arguments)]
    fn cell(
        &mut self,
        corners: [Point2; 3],
        root: usize,
        inherited: Option<(Point2, f64)>,
        depth: u32,
        order: usize,
        floor: f64,
        slack: &dyn Fn(u32, f64) -> f64,
    ) -> Result<Anchored, FarthestError> {
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
            inherited.map(|(a, f)| (lipschitz_lower(f, radius(a)), (a, f)));
        for a in [corners[0], corners[1], corners[2], centroid] {
            if !admissible(self.from, &self.roots[root], a)? {
                continue;
            }
            // Each map reaches all of a free-space triangle or none of it.
            let Some(f) = self.walk(a)? else {
                return Ok(Anchored::Unreached);
            };
            if crate::map::in_polygon(self.through, a)? && self.upper.is_none_or(|(u, _)| f < u) {
                self.upper = Some((f, a));
            }
            let bound = lipschitz_lower(f, radius(a));
            if best.is_none_or(|(b, _)| bound > b) {
                best = Some((bound, (a, f)));
            }
        }
        let pairs = self.pair_bound(&corners)?;
        Ok(match best {
            Some((bound, anchor)) => Anchored::Cell(Cell {
                corners,
                root,
                anchor,
                lower: (bound.max(pairs) - slack(depth, radius(anchor.0))).max(floor),
                depth,
                order,
            }),
            None => Anchored::None,
        })
    }
}

/// The vertices a cell might see, as `(dist + D, vertex, D, dist)` in
/// increasing order, where `dist` is the vertex's distance to the cell. A
/// vertex is left out only when one obstacle edge properly crosses the
/// segment from it to every corner of the cell, which then blocks it from
/// every point of the cell.
fn candidates(
    vertices: &[Vertex],
    t: &[Point2; 3],
    obstacles: &[(Point2, Point2)],
) -> Result<Vec<(f64, Point2, f64, f64)>, RouteError> {
    let mut out = Vec::new();
    'vertex: for vertex in vertices {
        let (v, d) = (vertex.at, vertex.distance);
        if let Some((a, b)) = vertex.solid {
            // Strictly inside the convex wedge from `v` between `a` and `b`.
            let inside = |c: Point2| -> Result<bool, RouteError> {
                let (ab, ac) = (side(v, a, b)?, side(v, a, c)?);
                let (ba, bc) = (side(v, b, a)?, side(v, b, c)?);
                Ok(ac != Sign::Zero && ac == ab && bc != Sign::Zero && bc == ba)
            };
            if inside(t[0])? && inside(t[1])? && inside(t[2])? {
                continue 'vertex;
            }
        }
        for &(p, q) in obstacles {
            if crosses(v, t[0], p, q)? && crosses(v, t[1], p, q)? && crosses(v, t[2], p, q)? {
                continue 'vertex;
            }
        }
        let g = triangle_distance(t, v)?;
        out.push((g + d, v, d, g));
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(out)
}

/// The least `f` can be within `radius` of a point where it is `f`: both
/// distances are 1-Lipschitz, so their sum falls at most twice as fast.
fn lipschitz_lower(f: f64, radius: f64) -> f64 {
    f - 2.0 * radius
}

/// A graph vertex a target is reached from.
struct Vertex {
    at: Point2,
    distance: f64,
    /// Neighbours `a`, `b` on its ring such that the open convex wedge from
    /// the vertex between them lies outside the free space: the inside of
    /// a hole at a convex corner, or the outside of the region at a reflex
    /// corner. Nothing strictly inside that wedge is seen from the vertex.
    solid: Option<(Point2, Point2)>,
}

/// The solid convex wedge at `v`, when `v` is a corner of exactly one ring
/// and no barrier, and the wedge less than a half turn lies outside the
/// free space. Decided exactly.
fn solid_corner(
    region: &[Polygon],
    obstacles: &[(Point2, Point2)],
    v: Point2,
) -> Result<Option<(Point2, Point2)>, RouteError> {
    if obstacles.iter().filter(|(p, q)| *p == v || *q == v).count() != 2 {
        return Ok(None);
    }
    for polygon in region {
        for (hole, ring) in
            core::iter::once((false, &polygon.outer)).chain(polygon.holes.iter().map(|h| (true, h)))
        {
            let n = ring.points.len();
            let Some(i) = ring.points.iter().position(|p| *p == v) else {
                continue;
            };
            let (prev, next) = (ring.points[(i + n - 1) % n], ring.points[(i + 1) % n]);
            let turn = side(prev, v, next)?;
            if turn == Sign::Zero {
                return Ok(None);
            }
            // The ring's orientation, from its lowest-leftmost corner,
            // which is convex.
            let k = (0..n)
                .min_by(|&a, &b| {
                    let (p, q) = (ring.points[a], ring.points[b]);
                    p.y.total_cmp(&q.y).then(p.x.total_cmp(&q.x))
                })
                .expect("a ring has corners");
            let orientation = side(
                ring.points[(k + n - 1) % n],
                ring.points[k],
                ring.points[(k + 1) % n],
            )?;
            // The ring's inside is the small wedge at a corner that turns
            // with the ring; solid means inside a hole, outside the outer.
            let small_is_inside = turn == orientation;
            return Ok((small_is_inside == hole).then_some((prev, next)));
        }
    }
    Ok(None)
}

/// The distance from `v` to the closed triangle: zero inside it, decided
/// exactly, else the least distance to an edge.
fn triangle_distance(t: &[Point2; 3], v: Point2) -> Result<f64, RouteError> {
    if meets_triangle(v, v, t)? {
        return Ok(0.0);
    }
    Ok((0..3)
        .map(|i| segment_distance(t[i], t[(i + 1) % 3], v))
        .fold(f64::INFINITY, f64::min))
}

fn segment_distance(a: Point2, b: Point2, v: Point2) -> f64 {
    let d = b - a;
    let length2 = d.dot(d);
    let s = if length2 > 0.0 {
        ((v - a).dot(d) / length2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (v - (a + d * s)).length()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance_map;
    use axiolid_overlay::Ring;

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    fn ring(points: &[(f64, f64)]) -> Ring {
        Ring {
            points: points.iter().map(|(x, y)| p(*x, *y)).collect(),
        }
    }

    fn two_holes() -> Polygon {
        Polygon {
            outer: ring(&[(0.0, 0.0), (12.0, 0.0), (12.0, 8.0), (0.0, 8.0)]),
            holes: vec![
                ring(&[(2.0, 2.0), (5.0, 2.0), (5.0, 6.0), (2.0, 6.0)]),
                ring(&[(7.0, 2.0), (10.0, 2.0), (10.0, 6.0), (7.0, 6.0)]),
            ],
        }
    }

    fn vertices(map: &DistanceMap) -> Vec<Vertex> {
        map.nodes()
            .iter()
            .copied()
            .zip(map.vertex_distances())
            .filter(|(_, d)| d.is_finite())
            .map(|(at, distance)| Vertex {
                at,
                distance,
                solid: solid_corner(map.region(), map.obstacles(), at).unwrap(),
            })
            .collect()
    }

    /// `f(y) = 2 |y|` when origin and target coincide: it falls at twice
    /// the rate of either distance, and the bound must allow for that.
    #[test]
    fn the_lipschitz_bound_allows_both_distances_to_fall_together() {
        let room = [Polygon {
            outer: ring(&[(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0)]),
            holes: Vec::new(),
        }];
        let map = distance_map(&room, &[], &[p(0.0, 0.0)]).unwrap();
        let a = p(5.0, 0.0);
        let f = 2.0 * map.at(a).unwrap().unwrap();
        // The nearest point within radius 1 of `a` to the origin.
        let y = p(4.0, 0.0);
        let least = 2.0 * map.at(y).unwrap().unwrap();
        assert!(lipschitz_lower(f, 1.0) <= least);
        assert_eq!(lipschitz_lower(f, 1.0), least);
    }

    /// The pair search stops early only when no later pair can win: its
    /// result is the brute-force least over every pair of candidates.
    #[test]
    fn the_pair_bound_is_the_least_over_every_pair() {
        let region = [two_holes()];
        let from = distance_map(&region, &[], &[p(1.0, 1.0), p(6.0, 7.5)]).unwrap();
        let to = distance_map(&region, &[], &[p(11.0, 7.0), p(6.0, 0.5)]).unwrap();
        let search = Search {
            from: &from,
            to: &to,
            through: &region[0],
            from_vertices: vertices(&from),
            to_vertices: vertices(&to),
            roots: &[],
            evaluated: HashMap::new(),
            upper: None,
        };
        let mut checked = 0;
        for i in 0..24 {
            for j in 0..16 {
                let (x, y) = (0.25 + 0.5 * f64::from(i), 0.25 + 0.5 * f64::from(j));
                let t = [p(x, y), p(x + 0.4, y), p(x, y + 0.4)];
                let origins = candidates(&search.from_vertices, &t, from.obstacles()).unwrap();
                let targets = candidates(&search.to_vertices, &t, to.obstacles()).unwrap();
                let mut brute = f64::INFINITY;
                for &(_, u, du, gu) in &origins {
                    for &(_, v, dv, gv) in &targets {
                        brute = brute.min(du + dv + (u - v).length().max(gu + gv));
                    }
                }
                assert_eq!(search.pair_bound(&t).unwrap(), brute, "{t:?}");
                checked += 1;
            }
        }
        assert_eq!(checked, 384);
    }

    /// A vertex is hidden only from a cell it sees no point of.
    #[test]
    fn a_vertex_is_hidden_only_from_the_whole_cell() {
        let region = [two_holes()];
        let map = distance_map(&region, &[], &[p(1.0, 1.0)]).unwrap();
        let find = |at: Point2| {
            vertices(&map)
                .into_iter()
                .find(|v| v.at == at)
                .expect("a vertex")
        };
        let hidden = |v: &Vertex, t: [Point2; 3]| {
            !candidates(core::slice::from_ref(v), &t, map.obstacles())
                .unwrap()
                .iter()
                .any(|c| c.1 == v.at)
        };
        // The origin, behind the first hole from a cell right of it.
        let origin = find(p(1.0, 1.0));
        let behind = [p(6.0, 4.5), p(6.5, 4.5), p(6.0, 5.0)];
        assert!(hidden(&origin, behind));
        // One corner of this cell pokes out below the hole's shadow.
        let peeking = [p(6.0, 4.5), p(6.5, 1.5), p(6.0, 5.0)];
        assert!(!hidden(&origin, peeking));
        // The hole's corner (5, 2) sees nothing strictly inside the hole's
        // corner wedge -- the cell above-left of it, inside the hole's
        // interior directions -- but does see a cell with one corner out.
        let corner = find(p(5.0, 2.0));
        assert!(corner.solid.is_some());
        let inside = [p(4.5, 2.5), p(4.8, 2.5), p(4.5, 2.8)];
        assert!(hidden(&corner, inside));
        let straddling = [p(4.5, 2.5), p(5.5, 2.5), p(4.5, 2.8)];
        assert!(!hidden(&corner, straddling));
        // The room's convex corner has no solid wedge.
        assert!(find(p(0.0, 0.0)).solid.is_none());
    }
}
