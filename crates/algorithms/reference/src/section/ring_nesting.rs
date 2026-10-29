//! Nesting of cut rings into regions with holes (ADR 0092).
//!
//! Pure plane geometry on the returned coordinates; every decision is an
//! exact `orient2d`. The steps, in order:
//!
//! 1. A crossing sweep over the ring edges, by bounding boxes. A proper
//!    crossing, or a collinear overlap of positive length whose edges run the
//!    same way (the material then lies on the same side), refuses the cut.
//! 2. Containment: ring `r` lies in ring `s` when every vertex of `r` that is
//!    not on `s` is strictly inside it. A mixed answer, or one with no
//!    decisive vertex, refuses the cut.
//! 3. The innermost container of a ring decides its role. An outer ring has
//!    no container or lies in a hole; a hole lies in an outer ring, which owns
//!    it. Anything else, or a tie between containers, refuses the cut.
//!
//! Candidate edge pairs of the sweep and ring pairs of the containment test
//! share one budget, `CutRegionLimits::max_pairs`.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation, Sign};
use axiolid_core::Point2;
use axiolid_mesh_section_contract::{CutRegion, CutRegionLimits, REFUSED_CROSSING_CUT_RINGS};

use super::{orient_sign, sequence_order, ScalarSection};

struct Ring {
    points: Vec<Point2>,
    hole: bool,
    /// `[min_x, min_y, max_x, max_y]`.
    bounds: [f64; 4],
}

impl Ring {
    fn new(points: Vec<Point2>, hole: bool) -> Self {
        let mut bounds = [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ];
        for point in &points {
            bounds[0] = bounds[0].min(point.x);
            bounds[1] = bounds[1].min(point.y);
            bounds[2] = bounds[2].max(point.x);
            bounds[3] = bounds[3].max(point.y);
        }
        Self {
            points,
            hole,
            bounds,
        }
    }

    fn within(&self, other: &Self) -> bool {
        self.bounds[0] >= other.bounds[0]
            && self.bounds[1] >= other.bounds[1]
            && self.bounds[2] <= other.bounds[2]
            && self.bounds[3] <= other.bounds[3]
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

/// Group `outers` (counter-clockwise) and `holes` (clockwise) into regions;
/// also return the number of candidate pairs compared.
pub(super) fn nest_rings(
    outers: Vec<Vec<Point2>>,
    holes: Vec<Vec<Point2>>,
    limits: CutRegionLimits,
    options: &ExecutionOptions,
) -> GeomResult<(Vec<CutRegion>, usize)> {
    let vertices: usize = outers.iter().chain(&holes).map(Vec::len).sum();
    if vertices > limits.max_region_vertices {
        return Err(GeomError::BudgetExceeded {
            resource: "section cut region vertices",
        });
    }
    let rings: Vec<Ring> = outers
        .into_iter()
        .map(|points| Ring::new(points, false))
        .chain(holes.into_iter().map(|points| Ring::new(points, true)))
        .collect();

    let mut pairs = 0usize;
    refuse_crossings(&rings, limits, &mut pairs, options)?;
    let containers = containers(&rings, limits, &mut pairs, options)?;

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
                    _ => return Err(refused()),
                }
            }
        };
        match (ring.hole, innermost) {
            (false, None) => {}
            (false, Some(container)) if rings[container].hole => {}
            (true, Some(container)) if !rings[container].hole => {
                owned[container].push(ring.points.clone());
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
            CutRegion::new(ring.points.clone(), holes)
        })
        .collect();
    regions.sort_by(|left, right| sequence_order(&left.outer, &right.outer));
    Ok((regions, pairs))
}

struct Edge {
    ring: usize,
    index: usize,
    start: Point2,
    end: Point2,
    /// `[min_x, min_y, max_x, max_y]`.
    bounds: [f64; 4],
}

/// Refuse a proper crossing or a same-direction collinear overlap of two
/// ring edges.
fn refuse_crossings(
    rings: &[Ring],
    limits: CutRegionLimits,
    pairs: &mut usize,
    options: &ExecutionOptions,
) -> GeomResult<()> {
    let mut edges = Vec::new();
    for (ring_index, ring) in rings.iter().enumerate() {
        for index in 0..ring.points.len() {
            let (start, end) = (
                ring.points[index],
                ring.points[(index + 1) % ring.points.len()],
            );
            edges.push(Edge {
                ring: ring_index,
                index,
                start,
                end,
                bounds: [
                    start.x.min(end.x),
                    start.y.min(end.y),
                    start.x.max(end.x),
                    start.y.max(end.y),
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
                let length = rings[first.ring].points.len();
                let gap = first.index.abs_diff(second.index);
                if gap == 1 || gap == length - 1 {
                    continue;
                }
            }
            count_pair(pairs, limits)?;
            if edges_conflict(first, second)? {
                return Err(refused());
            }
        }
    }
    Ok(())
}

fn edges_conflict(first: &Edge, second: &Edge) -> GeomResult<bool> {
    let (a, b, c, d) = (first.start, first.end, second.start, second.end);
    let around_first = [orient_sign(a, b, c)?, orient_sign(a, b, d)?];
    let around_second = [orient_sign(c, d, a)?, orient_sign(c, d, b)?];
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
    let along = |p: Point2| if a.x == b.x { p.y } else { p.x };
    let (low, high) = (along(a).min(along(b)), along(a).max(along(b)));
    let (other_low, other_high) = (along(c).min(along(d)), along(c).max(along(d)));
    let overlap = low.max(other_low) < high.min(other_high);
    let same_way = (along(b) > along(a)) == (along(d) > along(c));
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
fn locate(point: Point2, ring: &[Point2]) -> GeomResult<Location> {
    let mut winding = 0i32;
    for index in 0..ring.len() {
        let (a, b) = (ring[index], ring[(index + 1) % ring.len()]);
        let side = orient_sign(a, b, point)?;
        if side == Sign::Zero
            && a.x.min(b.x) <= point.x
            && point.x <= a.x.max(b.x)
            && a.y.min(b.y) <= point.y
            && point.y <= a.y.max(b.y)
        {
            return Ok(Location::Boundary);
        }
        if a.y <= point.y {
            if b.y > point.y && side == Sign::Positive {
                winding += 1;
            }
        } else if b.y <= point.y && side == Sign::Negative {
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
    rings: &[Ring],
    limits: CutRegionLimits,
    pairs: &mut usize,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Vec<usize>>> {
    let mut containers = vec![Vec::new(); rings.len()];
    for (inner_index, inner) in rings.iter().enumerate() {
        for (outer_index, outer) in rings.iter().enumerate() {
            if inner_index == outer_index || !inner.within(outer) {
                continue;
            }
            options.check_cancelled()?;
            count_pair(pairs, limits)?;
            let (mut inside, mut outside) = (false, false);
            for &point in &inner.points {
                match locate(point, &outer.points)? {
                    Location::Inside => inside = true,
                    Location::Outside => outside = true,
                    Location::Boundary => {}
                }
            }
            match (inside, outside) {
                (true, false) => containers[inner_index].push(outer_index),
                (false, true) => {}
                _ => return Err(refused()),
            }
        }
    }
    Ok(containers)
}

#[cfg(test)]
mod tests {
    use axiolid_core::Tolerance;

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

    fn nest(
        outers: Vec<Vec<Point2>>,
        holes: Vec<Vec<Point2>>,
        max_pairs: usize,
    ) -> GeomResult<(Vec<CutRegion>, usize)> {
        nest_rings(
            outers,
            holes,
            CutRegionLimits::new(10, 1_000, max_pairs),
            &ExecutionOptions::new(Tolerance::METRE),
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
        let options = ExecutionOptions::new(Tolerance::METRE);
        let two = vec![square(0.0, 1.0, false), square(5.0, 6.0, false)];
        let regions = nest_rings(
            two.clone(),
            vec![],
            CutRegionLimits::new(1, 100, 100),
            &options,
        );
        assert!(matches!(
            regions,
            Err(GeomError::BudgetExceeded {
                resource: "section cut regions"
            })
        ));
        let vertices = nest_rings(two, vec![], CutRegionLimits::new(10, 7, 100), &options);
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
        let triangle = [
            Point2::new(shift, shift),
            Point2::new(shift + 4.0, shift + 2.0),
            Point2::new(shift, shift + 4.0),
        ];
        let on_edge = Point2::new(shift + 2.0, shift + 1.0);
        assert_eq!(
            locate(on_edge, &triangle).expect("locate"),
            Location::Boundary
        );
        let past = Point2::new(shift + 2.0, shift + 1.0 - 2.0 * f64::EPSILON * shift);
        assert_eq!(locate(past, &triangle).expect("locate"), Location::Outside);
    }
}
