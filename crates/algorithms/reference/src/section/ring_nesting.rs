//! Nesting of cut rings into regions with holes (ADR 0092, #377).
//!
//! Pure plane geometry on the exact cut points (`PlanePoints`): source
//! vertices in the plane and exact crossings of source edges with it. Every
//! decision is an exact sign on those points; the returned coordinates only
//! label the output. The steps, in order:
//!
//! 1. A crossing sweep over the ring edges, by bounding boxes. A proper
//!    crossing, or a collinear overlap of positive length whose edges run the
//!    same way (the material then lies on the same side), refuses the cut. So
//!    does any contact between two non-adjacent edges of one ring (a repeated
//!    vertex), and a spike, where two adjacent edges fold back on each other.
//!    A contact between two rings is a local test at the contact point `p`:
//!    the directions in which one ring leaves `p` are classified against the
//!    material side of the other (its two edges at a vertex, its half-plane on
//!    an edge). Directions on both sides refuse the cut: the rings cross
//!    through `p`. Rings that only touch, as neighbours sharing an edge with
//!    opposite directions do, pass.
//! 2. A ring is an outer ring or a hole by its turn at its lexicographically
//!    smallest vertex. The sweep has proven every ring simple, so that vertex
//!    is never straight.
//! 3. Containment: ring `r` lies in ring `s` when every vertex of `r` that is
//!    not on `s` is strictly inside it. A mixed answer, or one with no
//!    decisive vertex, refuses the cut.
//! 4. The innermost container of a ring decides its role. An outer ring has
//!    no container or lies in a hole; a hole lies in an outer ring, which owns
//!    it. Anything else, or a tie between containers, refuses the cut.
//!
//! Candidate edge pairs of the sweep and ring pairs of the containment test
//! share one budget, `CutRegionLimits::max_pairs`.

use std::cmp::Ordering;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation, Sign};
use axiolid_core::Point2;
use axiolid_mesh_section_contract::{CutRegion, CutRegionLimits, REFUSED_CROSSING_CUT_RINGS};

use super::plane_points::PlanePoints;
use super::regions::canonical_rotation;
use super::{sequence_order, ScalarSection};

/// A traced ring: its exact points, and the returned coordinates of the same
/// vertices.
pub(super) struct TracedRing {
    pub(super) nodes: Vec<usize>,
    pub(super) points: Vec<Point2>,
}

struct Ring {
    nodes: Vec<usize>,
    points: Vec<Point2>,
    /// Decided after the sweep has proven the ring simple.
    hole: bool,
    /// `[min_u, min_v, max_u, max_v]` of the exact points, with every vertex
    /// box at its outer limit: contains the ring.
    bounds: [f64; 4],
    /// The same with every vertex box at its inner limit: the extremes of the
    /// ring reach at least this far.
    core: [f64; 4],
}

impl Ring {
    fn new(exact: &PlanePoints, traced: TracedRing) -> Self {
        let mut bounds = [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ];
        let mut core = [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ];
        for &node in &traced.nodes {
            let [u_lo, v_lo, u_hi, v_hi] = exact.bounds(node);
            bounds = [
                bounds[0].min(u_lo),
                bounds[1].min(v_lo),
                bounds[2].max(u_hi),
                bounds[3].max(v_hi),
            ];
            core = [
                core[0].min(u_hi),
                core[1].min(v_hi),
                core[2].max(u_lo),
                core[3].max(v_lo),
            ];
        }
        Self {
            nodes: traced.nodes,
            points: traced.points,
            hole: false,
            bounds,
            core,
        }
    }

    /// Whether this ring can lie in `other` by their boxes: its smallest
    /// coordinate cannot be below the box of `other`, nor its largest above.
    fn possibly_within(&self, other: &Self) -> bool {
        self.core[0] >= other.bounds[0]
            && self.core[1] >= other.bounds[1]
            && self.core[2] <= other.bounds[2]
            && self.core[3] <= other.bounds[3]
    }
}

fn refused() -> GeomError {
    GeomError::UnsupportedInput {
        backend: ScalarSection::ID,
        operation: Operation::MeshPlaneSection,
        input: REFUSED_CROSSING_CUT_RINGS,
    }
}

/// Count one candidate pair against the budget.
fn count_pair(pairs: &mut usize, limits: CutRegionLimits) -> GeomResult<()> {
    *pairs += 1;
    if *pairs > limits.max_pairs {
        return Err(GeomError::BudgetExceeded {
            resource: "section cut pairs",
        });
    }
    Ok(())
}

/// Whether the ring is clockwise: its turn at the lexicographically smallest
/// vertex, which is a corner of the convex hull and, in a simple ring, never
/// straight.
fn is_hole(exact: &PlanePoints, nodes: &[usize]) -> GeomResult<bool> {
    let mut smallest = 0;
    for index in 1..nodes.len() {
        if exact.lexicographic(nodes[index], nodes[smallest])? == Ordering::Less {
            smallest = index;
        }
    }
    let count = nodes.len();
    let turn = exact.orient(
        nodes[(smallest + count - 1) % count],
        nodes[smallest],
        nodes[(smallest + 1) % count],
    )?;
    match turn {
        Sign::Positive => Ok(false),
        Sign::Negative => Ok(true),
        // Not reachable after the sweep. Kept as a guard.
        _ => Err(GeomError::Degenerate(
            "cut region ring has no signed area".into(),
        )),
    }
}

/// Group the rings into regions, outer rings counter-clockwise and holes
/// clockwise; also return the number of candidate pairs compared.
pub(super) fn nest_rings(
    exact: &PlanePoints,
    rings: Vec<TracedRing>,
    limits: CutRegionLimits,
    options: &ExecutionOptions,
) -> GeomResult<(Vec<CutRegion>, usize)> {
    let vertices: usize = rings.iter().map(|ring| ring.nodes.len()).sum();
    if vertices > limits.max_region_vertices {
        return Err(GeomError::BudgetExceeded {
            resource: "section cut region vertices",
        });
    }
    let mut rings: Vec<Ring> = rings
        .into_iter()
        .map(|traced| Ring::new(exact, traced))
        .collect();

    let mut pairs = 0usize;
    refuse_crossings(exact, &rings, limits, &mut pairs, options)?;
    for ring in &mut rings {
        ring.hole = is_hole(exact, &ring.nodes)?;
    }
    let containers = containers(exact, &rings, limits, &mut pairs, options)?;

    // The depth of a ring is the number of rings around it.
    let mut owned: Vec<Vec<Vec<Point2>>> = vec![Vec::new(); rings.len()];
    for (index, ring) in rings.iter().enumerate() {
        let depth = |container: usize| containers[container].len();
        let deepest = containers[index].iter().map(|&c| depth(c)).max();
        let innermost = match deepest {
            None => None,
            Some(deepest) => {
                let mut tied = containers[index].iter().filter(|&&c| depth(c) == deepest);
                match (tied.next(), tied.next()) {
                    (Some(&only), None) => Some(only),
                    // Two unrelated containers of equal depth would have to
                    // overlap each other, and the sweep and the containment
                    // test refuse that first. Kept as a guard, not reachable
                    // by a test.
                    _ => return Err(refused()),
                }
            }
        };
        match (ring.hole, innermost) {
            (false, None) => {}
            (false, Some(container)) if rings[container].hole => {}
            (true, Some(container)) if !rings[container].hole => {
                owned[container].push(canonical_rotation(ring.points.clone()));
            }
            _ => return Err(refused()),
        }
    }

    let outer_count = rings.iter().filter(|ring| !ring.hole).count();
    if outer_count > limits.max_regions {
        return Err(GeomError::BudgetExceeded {
            resource: "section cut regions",
        });
    }
    let mut regions: Vec<CutRegion> = rings
        .iter()
        .zip(owned)
        .filter(|(ring, _)| !ring.hole)
        .map(|(ring, mut holes)| {
            holes.sort_by(|left, right| sequence_order(left, right));
            CutRegion::new(canonical_rotation(ring.points.clone()), holes)
        })
        .collect();
    regions.sort_by(|left, right| sequence_order(&left.outer, &right.outer));
    Ok((regions, pairs))
}

struct Edge {
    ring: usize,
    index: usize,
    start: usize,
    end: usize,
    /// `[min_u, min_v, max_u, max_v]`, a box of the exact edge.
    bounds: [f64; 4],
}

/// Refuse a proper crossing or a same-direction collinear overlap of two
/// ring edges.
fn refuse_crossings(
    exact: &PlanePoints,
    rings: &[Ring],
    limits: CutRegionLimits,
    pairs: &mut usize,
    options: &ExecutionOptions,
) -> GeomResult<()> {
    for ring in rings {
        refuse_spikes(exact, &ring.nodes)?;
    }
    let mut edges = Vec::new();
    for (ring_index, ring) in rings.iter().enumerate() {
        for index in 0..ring.nodes.len() {
            let (start, end) = (
                ring.nodes[index],
                ring.nodes[(index + 1) % ring.nodes.len()],
            );
            let (first, second) = (exact.bounds(start), exact.bounds(end));
            edges.push(Edge {
                ring: ring_index,
                index,
                start,
                end,
                bounds: [
                    first[0].min(second[0]),
                    first[1].min(second[1]),
                    first[2].max(second[2]),
                    first[3].max(second[3]),
                ],
            });
        }
    }
    edges.sort_by(|left, right| left.bounds[0].total_cmp(&right.bounds[0]));
    for (first_index, first) in edges.iter().enumerate() {
        options.check_cancelled()?;
        for second in &edges[first_index + 1..] {
            if second.bounds[0] > first.bounds[2] {
                break;
            }
            if second.bounds[1] > first.bounds[3] || first.bounds[1] > second.bounds[3] {
                continue;
            }
            if first.ring == second.ring {
                let length = rings[first.ring].nodes.len();
                let gap = first.index.abs_diff(second.index);
                if gap == 1 || gap == length - 1 {
                    continue;
                }
            }
            count_pair(pairs, limits)?;
            if edges_conflict(exact, first, second)? {
                return Err(refused());
            }
            for contact in contacts(exact, first, second)? {
                if first.ring == second.ring
                    || crosses_through(exact, &rings[first.ring], &rings[second.ring], contact)?
                    || crosses_through(exact, &rings[second.ring], &rings[first.ring], contact)?
                {
                    return Err(refused());
                }
            }
        }
    }
    Ok(())
}

/// Refuse a vertex whose two edges are collinear and fold back: the middle
/// point must lie strictly between the outer two.
fn refuse_spikes(exact: &PlanePoints, ring: &[usize]) -> GeomResult<()> {
    for index in 0..ring.len() {
        let (before, at, after) = (
            ring[(index + ring.len() - 1) % ring.len()],
            ring[index],
            ring[(index + 1) % ring.len()],
        );
        if exact.orient(before, at, after)? == Sign::Zero
            && !strictly_between(exact, before, at, after)?
        {
            return Err(refused());
        }
    }
    Ok(())
}

/// The axis of the projection along which two distinct points of a line
/// differ.
fn dominant(exact: &PlanePoints, a: usize, b: usize) -> GeomResult<usize> {
    Ok(if exact.compare(a, b, 0)? == Ordering::Equal {
        1
    } else {
        0
    })
}

/// `(low, high)` of two points along `axis`.
fn ordered(exact: &PlanePoints, a: usize, b: usize, axis: usize) -> GeomResult<(usize, usize)> {
    Ok(if exact.compare(a, b, axis)? == Ordering::Greater {
        (b, a)
    } else {
        (a, b)
    })
}

/// Whether collinear `at` lies strictly between `before` and `after`.
fn strictly_between(
    exact: &PlanePoints,
    before: usize,
    at: usize,
    after: usize,
) -> GeomResult<bool> {
    let axis = dominant(exact, before, after)?;
    let (low, high) = ordered(exact, before, after, axis)?;
    Ok(exact.compare(low, at, axis)? == Ordering::Less
        && exact.compare(at, high, axis)? == Ordering::Less)
}

/// Whether collinear `point` lies on the closed segment `start`-`end`.
fn within_segment(exact: &PlanePoints, point: usize, start: usize, end: usize) -> GeomResult<bool> {
    let axis = dominant(exact, start, end)?;
    let (low, high) = ordered(exact, start, end, axis)?;
    Ok(exact.compare(low, point, axis)? != Ordering::Greater
        && exact.compare(point, high, axis)? != Ordering::Greater)
}

/// Whether `point` lies on the closed segment of `edge`.
fn on_edge(exact: &PlanePoints, point: usize, edge: &Edge) -> GeomResult<bool> {
    Ok(exact.orient(edge.start, edge.end, point)? == Sign::Zero
        && within_segment(exact, point, edge.start, edge.end)?)
}

/// The points where two closed edges touch: an end of one lying on the other.
fn contacts(exact: &PlanePoints, first: &Edge, second: &Edge) -> GeomResult<Vec<usize>> {
    let mut found = Vec::new();
    for point in [first.start, first.end] {
        if on_edge(exact, point, second)? {
            found.push(point);
        }
    }
    for point in [second.start, second.end] {
        if on_edge(exact, point, first)? {
            found.push(point);
        }
    }
    Ok(found)
}

/// How a ring passes through a point on its boundary.
enum Passage {
    /// At a vertex: the previous and the next vertex.
    Vertex(usize, usize),
    /// Inside an edge: its start and end.
    Edge(usize, usize),
}

fn passages(exact: &PlanePoints, ring: &[usize], point: usize) -> GeomResult<Vec<Passage>> {
    let mut found = Vec::new();
    for index in 0..ring.len() {
        let (start, end) = (ring[index], ring[(index + 1) % ring.len()]);
        if exact.same(start, point)? {
            found.push(Passage::Vertex(
                ring[(index + ring.len() - 1) % ring.len()],
                end,
            ));
        } else if !exact.same(end, point)?
            && exact.orient(start, end, point)? == Sign::Zero
            && strictly_between(exact, start, point, end)?
        {
            found.push(Passage::Edge(start, end));
        }
    }
    Ok(found)
}

/// Whether `towards` lies strictly on the material side of `passage` at
/// `point`. Material is on the left of every ring edge. A direction along an
/// edge of the passage is not strictly inside.
fn in_material(
    exact: &PlanePoints,
    passage: &Passage,
    point: usize,
    towards: usize,
) -> GeomResult<bool> {
    Ok(match *passage {
        Passage::Edge(start, end) => exact.orient(start, end, towards)? == Sign::Positive,
        Passage::Vertex(before, after) => match exact.orient(point, after, before)? {
            // Convex: strictly between the next edge and the reversed previous one.
            Sign::Positive => {
                exact.orient(point, after, towards)? == Sign::Positive
                    && exact.orient(point, towards, before)? == Sign::Positive
            }
            // Reflex: anything but the closed wedge outside the material.
            Sign::Negative => {
                !(exact.orient(point, before, towards)? != Sign::Negative
                    && exact.orient(point, towards, after)? != Sign::Negative)
            }
            _ => exact.orient(point, after, towards)? == Sign::Positive,
        },
    })
}

/// Whether `mover` leaves the contact `point` on both sides of `container`'s
/// material: some direction strictly inside it and some not.
fn crosses_through(
    exact: &PlanePoints,
    mover: &Ring,
    container: &Ring,
    point: usize,
) -> GeomResult<bool> {
    let mut directions = Vec::new();
    for passage in passages(exact, &mover.nodes, point)? {
        match passage {
            Passage::Vertex(before, after) | Passage::Edge(before, after) => {
                directions.extend([before, after]);
            }
        }
    }
    let Some(passage) = passages(exact, &container.nodes, point)?.into_iter().next() else {
        return Ok(false);
    };
    let (mut inside, mut outside) = (false, false);
    for towards in directions {
        if in_material(exact, &passage, point, towards)? {
            inside = true;
        } else {
            outside = true;
        }
    }
    Ok(inside && outside)
}

fn edges_conflict(exact: &PlanePoints, first: &Edge, second: &Edge) -> GeomResult<bool> {
    let (a, b, c, d) = (first.start, first.end, second.start, second.end);
    let around_first = [exact.orient(a, b, c)?, exact.orient(a, b, d)?];
    let around_second = [exact.orient(c, d, a)?, exact.orient(c, d, b)?];
    let signs = [
        around_first[0],
        around_first[1],
        around_second[0],
        around_second[1],
    ];
    if signs.iter().all(|&sign| sign != Sign::Zero) {
        return Ok(around_first[0] != around_first[1] && around_second[0] != around_second[1]);
    }
    if signs.iter().any(|&sign| sign != Sign::Zero) {
        return Ok(false);
    }
    // Collinear: compare along the line's dominant axis.
    let axis = dominant(exact, a, b)?;
    let (low, high) = ordered(exact, a, b, axis)?;
    let (other_low, other_high) = ordered(exact, c, d, axis)?;
    let upper_low = if exact.compare(low, other_low, axis)? == Ordering::Less {
        other_low
    } else {
        low
    };
    let lower_high = if exact.compare(high, other_high, axis)? == Ordering::Greater {
        other_high
    } else {
        high
    };
    let overlap = exact.compare(upper_low, lower_high, axis)? == Ordering::Less;
    let same_way = (exact.compare(b, a, axis)? == Ordering::Greater)
        == (exact.compare(d, c, axis)? == Ordering::Greater);
    Ok(overlap && same_way)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Location {
    Inside,
    Outside,
    Boundary,
}

/// Where `point` lies relative to `ring`: a crossing-number walk with exact
/// orientation tests.
///
/// `orient` is in the frame's sense and the `v` comparisons are in the
/// projection's. When the projection mirrors the frame the walk counts the
/// crossings of the opposite ray, which is the winding number negated; only
/// zero or nonzero is used.
fn locate(exact: &PlanePoints, point: usize, ring: &Ring) -> GeomResult<Location> {
    let mut winding = 0i32;
    for index in 0..ring.nodes.len() {
        let (a, b) = (
            ring.nodes[index],
            ring.nodes[(index + 1) % ring.nodes.len()],
        );
        let side = exact.orient(a, b, point)?;
        if side == Sign::Zero && within_segment(exact, point, a, b)? {
            return Ok(Location::Boundary);
        }
        if exact.compare(a, point, 1)? != Ordering::Greater {
            if exact.compare(b, point, 1)? == Ordering::Greater && side == Sign::Positive {
                winding += 1;
            }
        } else if exact.compare(b, point, 1)? != Ordering::Greater && side == Sign::Negative {
            winding -= 1;
        }
    }
    Ok(if winding == 0 {
        Location::Outside
    } else {
        Location::Inside
    })
}

/// For each ring, the rings that contain it.
fn containers(
    exact: &PlanePoints,
    rings: &[Ring],
    limits: CutRegionLimits,
    pairs: &mut usize,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Vec<usize>>> {
    let mut containers = vec![Vec::new(); rings.len()];
    for (inner_index, inner) in rings.iter().enumerate() {
        for (outer_index, outer) in rings.iter().enumerate() {
            if inner_index == outer_index || !inner.possibly_within(outer) {
                continue;
            }
            options.check_cancelled()?;
            count_pair(pairs, limits)?;
            let (mut inside, mut outside) = (false, false);
            for &node in &inner.nodes {
                match locate(exact, node, outer)? {
                    Location::Inside => inside = true,
                    Location::Outside => outside = true,
                    Location::Boundary => {}
                }
            }
            match (inside, outside) {
                (true, false) => containers[inner_index].push(outer_index),
                (false, true) => {}
                // Mixed: some vertices strictly inside and some strictly
                // outside needs a crossing or a contact, which the sweep
                // refuses first. Kept as a guard.
                _ => return Err(refused()),
            }
        }
    }
    Ok(containers)
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Frame3, Point3, Tolerance, Vec3};

    use super::super::regions::Projection;
    use super::super::{orient_sign, project_local, ExactSectionPlane};
    use super::*;

    fn square(low: f64, high: f64, clockwise: bool) -> Vec<Point2> {
        let mut points = vec![
            Point2::new(low, low),
            Point2::new(high, low),
            Point2::new(high, high),
            Point2::new(low, high),
        ];
        if clockwise {
            points[1..].reverse();
        }
        points
    }

    fn rectangle(x: (f64, f64), y: (f64, f64)) -> Vec<Point2> {
        vec![
            Point2::new(x.0, y.0),
            Point2::new(x.1, y.0),
            Point2::new(x.1, y.1),
            Point2::new(x.0, y.1),
        ]
    }

    /// An exact plane z = 0 with the identity frame.
    fn flat_plane() -> (Frame3, PlanePoints) {
        let frame = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let plane = ExactSectionPlane::new(frame).expect("plane");
        let projection = Projection::new(&plane).expect("projection");
        (frame, PlanePoints::new(&plane, projection).expect("points"))
    }

    /// Each point as an exact source vertex of the plane z = 0, in its own
    /// node, so equal points only compare equal through the exact test.
    fn traced(exact: &mut PlanePoints, rings: Vec<Vec<Point2>>) -> Vec<TracedRing> {
        rings
            .into_iter()
            .map(|points| TracedRing {
                nodes: points
                    .iter()
                    .map(|point| {
                        exact
                            .vertex(Point3::new(point.x, point.y, 0.0))
                            .expect("vertex")
                    })
                    .collect(),
                points,
            })
            .collect()
    }

    fn nest_with(
        rings: Vec<Vec<Point2>>,
        limits: CutRegionLimits,
    ) -> GeomResult<(Vec<CutRegion>, usize)> {
        let (_, mut exact) = flat_plane();
        let rings = traced(&mut exact, rings);
        nest_rings(
            &exact,
            rings,
            limits,
            &ExecutionOptions::new(Tolerance::METRE),
        )
    }

    /// The role of a ring follows its orientation, not the list it is in.
    fn nest(
        outers: Vec<Vec<Point2>>,
        holes: Vec<Vec<Point2>>,
        max_pairs: usize,
    ) -> GeomResult<(Vec<CutRegion>, usize)> {
        nest_with(
            outers.into_iter().chain(holes).collect(),
            CutRegionLimits::new(10, 1_000, max_pairs),
        )
    }

    fn is_refused(result: GeomResult<(Vec<CutRegion>, usize)>) -> bool {
        matches!(
            result,
            Err(GeomError::UnsupportedInput {
                input: REFUSED_CROSSING_CUT_RINGS,
                ..
            })
        )
    }

    #[test]
    fn nested_squares_alternate_outer_and_hole_by_the_innermost_container() {
        // A [0,10] holds hole B [1,9], which holds outer C [2,8] with hole D
        // [3,7]. C also lies in A, but B is its innermost container.
        let (regions, pairs) = nest(
            vec![square(0.0, 10.0, false), square(2.0, 8.0, false)],
            vec![square(1.0, 9.0, true), square(3.0, 7.0, true)],
            1_000,
        )
        .expect("nested");
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].outer, square(0.0, 10.0, false));
        assert_eq!(regions[0].holes, [square(1.0, 9.0, true)]);
        assert_eq!(regions[1].outer, square(2.0, 8.0, false));
        assert_eq!(regions[1].holes, [square(3.0, 7.0, true)]);
        assert!(pairs > 0);
    }

    #[test]
    fn a_proper_crossing_is_refused_by_name() {
        let plus = nest(
            vec![
                rectangle((0.0, 3.0), (1.0, 2.0)),
                rectangle((1.0, 2.0), (0.0, 3.0)),
            ],
            vec![],
            1_000,
        );
        assert!(is_refused(plus));
    }

    #[test]
    fn a_collinear_overlap_is_refused_only_when_the_edges_run_the_same_way() {
        let overlapping = nest(
            vec![
                rectangle((0.0, 2.0), (0.0, 1.0)),
                rectangle((1.0, 3.0), (0.0, 1.0)),
            ],
            vec![],
            1_000,
        );
        assert!(is_refused(overlapping));
        // Neighbours sharing the edge x = 1 run in opposite directions there.
        let (regions, _) = nest(
            vec![
                rectangle((0.0, 1.0), (0.0, 1.0)),
                rectangle((1.0, 2.0), (0.0, 1.0)),
            ],
            vec![],
            1_000,
        )
        .expect("neighbours");
        assert_eq!(regions.len(), 2);
    }

    #[test]
    fn an_outer_ring_in_an_outer_ring_is_refused() {
        let nested = nest(
            vec![square(0.0, 10.0, false), square(2.0, 3.0, false)],
            vec![],
            1_000,
        );
        assert!(is_refused(nested));
    }

    #[test]
    fn a_hole_without_a_container_is_refused() {
        assert!(is_refused(nest(
            vec![],
            vec![square(0.0, 1.0, true)],
            1_000
        )));
        let beside = nest(
            vec![square(0.0, 1.0, false)],
            vec![rectangle((5.0, 6.0), (5.0, 6.0))
                .into_iter()
                .rev()
                .collect()],
            1_000,
        );
        assert!(is_refused(beside));
    }

    #[test]
    fn a_hole_that_pokes_out_of_its_outer_ring_is_refused() {
        // The hole is inside on three vertices and outside on one.
        let poking = nest(
            vec![square(0.0, 4.0, false)],
            vec![vec![
                Point2::new(1.0, 1.0),
                Point2::new(1.0, 3.0),
                Point2::new(3.0, 3.0),
                Point2::new(3.0, 5.0),
            ]],
            1_000,
        );
        assert!(is_refused(poking));
    }

    #[test]
    fn a_ring_that_leaves_its_container_through_vertices_only_is_refused() {
        // The kite A (3,1), V (3,2), B (3,3), C (2,2) lies in the box of the
        // L-shape. V is on the notch edge and C is the reflex corner, so no
        // edge crosses properly, yet A is inside the L and B is outside.
        let ell = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 2.0),
            Point2::new(2.0, 2.0),
            Point2::new(2.0, 4.0),
            Point2::new(0.0, 4.0),
        ];
        let kite = vec![
            Point2::new(2.0, 2.0),
            Point2::new(3.0, 1.0),
            Point2::new(3.0, 2.0),
            Point2::new(3.0, 3.0),
        ];
        assert!(is_refused(nest(vec![ell, kite], vec![], 1_000)));
    }

    #[test]
    fn a_hole_in_a_hole_is_refused() {
        let nested = nest(
            vec![square(0.0, 10.0, false)],
            vec![square(1.0, 9.0, true), square(3.0, 7.0, true)],
            1_000,
        );
        assert!(is_refused(nested));
    }

    #[test]
    fn a_ring_crossing_itself_is_refused() {
        let bow_tie = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(0.0, 2.0),
            Point2::new(4.0, 4.0),
            Point2::new(5.0, 2.0),
        ];
        assert!(is_refused(nest(vec![bow_tie], vec![], 1_000)));
    }

    fn points(list: &[(f64, f64)]) -> Vec<Point2> {
        list.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    #[test]
    fn rings_crossing_through_shared_vertices_are_refused() {
        // R's edge (4,0)-(0,4) runs inside S while its other edges leave S,
        // and no edge crosses properly.
        let first = nest(
            vec![
                points(&[(0.0, 0.0), (4.0, 0.0), (5.0, 5.0), (0.0, 4.0)]),
                points(&[(4.0, 0.0), (0.0, 4.0), (-1.0, -1.0)]),
            ],
            vec![],
            1_000,
        );
        assert!(is_refused(first));
        // A hexagon and a kite share four vertices and overlap in a lens.
        let second = nest(
            vec![
                points(&[
                    (0.0, 0.0),
                    (1.0, 0.0),
                    (2.0, 0.5),
                    (1.0, 1.0),
                    (0.0, 1.0),
                    (-1.0, 0.5),
                ]),
                points(&[
                    (0.0, 0.0),
                    (0.5, -1.0),
                    (1.0, 0.0),
                    (1.0, 1.0),
                    (0.5, 2.0),
                    (0.0, 1.0),
                ]),
            ],
            vec![],
            1_000,
        );
        assert!(is_refused(second));
    }

    #[test]
    fn a_ring_that_repeats_a_vertex_is_refused() {
        let figure = points(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (2.0, 2.0),
            (1.0, 3.0),
            (0.0, 2.0),
            (2.0, 2.0),
        ]);
        assert!(is_refused(nest(vec![figure], vec![], 1_000)));
    }

    #[test]
    fn a_spike_is_refused_and_a_straight_vertex_is_not() {
        let spike = points(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (4.0, 5.0),
            (4.0, 4.0),
            (0.0, 4.0),
        ]);
        assert!(is_refused(nest(vec![spike], vec![], 1_000)));
        let straight = points(&[(0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]);
        assert!(nest(vec![straight], vec![], 1_000).is_ok());
        // Three collinear points fold back with all edges adjacent, so only
        // the spike rule sees them.
        let flat = points(&[(0.0, 0.0), (2.0, 0.0), (1.0, 0.0)]);
        assert!(is_refused(nest(vec![flat], vec![], 1_000)));
    }

    #[test]
    fn a_hole_touching_its_outer_ring_from_inside_is_accepted() {
        // The hole's vertex (2,0) lies on the outer bottom edge; both of its
        // edges lead into the material.
        let touching = points(&[(2.0, 0.0), (1.0, 1.0), (3.0, 1.0)]);
        let (regions, _) = nest(vec![square(0.0, 4.0, false)], vec![touching], 1_000)
            .expect("a hole touching the outer ring");
        assert_eq!(regions.len(), 1);
    }

    #[test]
    fn the_pair_budget_covers_the_sweep_and_the_containment_test() {
        let rings = || (vec![square(0.0, 10.0, false)], vec![square(1.0, 9.0, true)]);
        let (outers, holes) = rings();
        let (_, pairs) = nest(outers, holes, 1_000).expect("nested");
        assert_eq!(pairs, 1, "one containment pair, no candidate edges");
        let (outers, holes) = rings();
        let error = nest(outers, holes, pairs - 1).expect_err("one pair too few");
        assert!(matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "section cut pairs"
            }
        ));
    }

    #[test]
    fn candidate_edge_pairs_of_the_sweep_count_against_the_budget() {
        let neighbours = || {
            vec![
                rectangle((0.0, 1.0), (0.0, 1.0)),
                rectangle((1.0, 2.0), (0.0, 1.0)),
            ]
        };
        let (_, pairs) = nest(neighbours(), vec![], 1_000).expect("neighbours");
        assert!(pairs > 0, "the shared edge is a candidate pair");
        let error = nest(neighbours(), vec![], 0).expect_err("no pair allowed");
        assert!(matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "section cut pairs"
            }
        ));
    }

    #[test]
    fn region_and_vertex_budgets_are_enforced() {
        let two = || vec![square(0.0, 1.0, false), square(5.0, 6.0, false)];
        let regions = nest_with(two(), CutRegionLimits::new(1, 100, 100));
        assert!(matches!(
            regions,
            Err(GeomError::BudgetExceeded {
                resource: "section cut regions"
            })
        ));
        let vertices = nest_with(two(), CutRegionLimits::new(10, 7, 100));
        assert!(matches!(
            vertices,
            Err(GeomError::BudgetExceeded {
                resource: "section cut region vertices"
            })
        ));
    }

    #[test]
    fn a_vertex_one_ulp_outside_a_distant_ring_is_not_inside() {
        // Far from the origin, translated by 2^27, the point (2,1) lies on
        // the edge from (0,0) to (4,2); two ulps lower it is outside.
        let shift = 134_217_728.0;
        let triangle = vec![
            Point2::new(shift, shift),
            Point2::new(shift + 4.0, shift + 2.0),
            Point2::new(shift, shift + 4.0),
        ];
        let on_edge = Point2::new(shift + 2.0, shift + 1.0);
        let past = Point2::new(shift + 2.0, shift + 1.0 - 2.0 * f64::EPSILON * shift);
        let (_, mut exact) = flat_plane();
        let ring = traced(&mut exact, vec![triangle]).remove(0);
        let ring = Ring::new(&exact, ring);
        let probes = traced(&mut exact, vec![vec![on_edge], vec![past]]);
        assert_eq!(
            locate(&exact, probes[0].nodes[0], &ring).expect("locate"),
            Location::Boundary
        );
        assert_eq!(
            locate(&exact, probes[1].nodes[0], &ring).expect("locate"),
            Location::Outside
        );
    }

    #[test]
    fn rings_touching_on_a_line_pass_where_their_rounded_points_cross() {
        // Two source rings meet along y = 0.1 on the plane z = 0. Under a
        // turned frame the returned points are rounded, and the corners of B
        // can fall inside A there; the exact points never do.
        let ring_a = [(0.0, 0.0), (2.0, 0.0), (2.0, 0.1), (0.0, 0.1)];
        let ring_b = [(0.9, 0.1), (1.1, 0.1), (1.1, 1.0), (0.9, 1.0)];
        let mut penetrating = false;
        for k in 1..=200 {
            let (s, c) = (std::f64::consts::PI * f64::from(k) / 201.0).sin_cos();
            let frame = Frame3 {
                origin: Point3::ZERO,
                x: Vec3::new(c, s, 0.0),
                y: Vec3::new(-s, c, 0.0),
                z: Vec3::Z,
            };
            let plane = ExactSectionPlane::new(frame).expect("plane");
            let projection = Projection::new(&plane).expect("projection");
            let mut exact = PlanePoints::new(&plane, projection).expect("points");
            let rings: Vec<TracedRing> = [ring_a.as_slice(), ring_b.as_slice()]
                .into_iter()
                .map(|ring| {
                    let world: Vec<Point3> =
                        ring.iter().map(|&(x, y)| Point3::new(x, y, 0.0)).collect();
                    TracedRing {
                        nodes: world
                            .iter()
                            .map(|&point| exact.vertex(point).expect("vertex"))
                            .collect(),
                        points: world
                            .iter()
                            .map(|&point| project_local(frame, point).expect("local"))
                            .collect(),
                    }
                })
                .collect();
            let (a, b) = (&rings[0].points, &rings[1].points);
            // The top edge of A runs a[2] -> a[3]; material is on its left.
            penetrating |= [b[0], b[1]]
                .into_iter()
                .any(|corner| orient_sign(a[2], a[3], corner) == Ok(Sign::Positive));
            let options = ExecutionOptions::new(Tolerance::METRE);
            let (regions, _) = nest_rings(
                &exact,
                rings,
                CutRegionLimits::new(10, 1_000, 1_000),
                &options,
            )
            .unwrap_or_else(|error| panic!("frame {k} refused: {error:?}"));
            assert_eq!(regions.len(), 2, "frame {k}");
        }
        assert!(penetrating, "the rounded points never overlap");
    }

    #[test]
    fn the_lexicographically_smallest_vertex_decides_outer_or_hole() {
        // The L starts at its reflex corner (2,2), where it turns clockwise
        // although the ring is counter-clockwise.
        let ell = points(&[
            (2.0, 2.0),
            (2.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 2.0),
        ]);
        let (alone, _) = nest(vec![ell.clone()], vec![], 1_000).expect("an outer ring");
        assert_eq!(alone.len(), 1);
        assert!(alone[0].holes.is_empty());
        // Reversed, it is a hole of the big square around it.
        let reversed: Vec<Point2> = std::iter::once(ell[0])
            .chain(ell[1..].iter().rev().copied())
            .collect();
        let big = rectangle((-1.0, 5.0), (-1.0, 5.0));
        let (holed, _) = nest(vec![big], vec![reversed], 1_000).expect("a hole");
        assert_eq!(holed.len(), 1);
        assert_eq!(holed[0].holes.len(), 1);
    }
}
