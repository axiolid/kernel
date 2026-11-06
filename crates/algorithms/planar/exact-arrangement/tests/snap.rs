//! Snap rounding through the public API only (ADR 0100).

use std::collections::BTreeMap;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Point2, Tolerance};
use axiolid_exact::{Arith, Dyadic};
use axiolid_exact_arrangement::{
    arrange, canonical_rotation, check_winding, snap_round, trace, ExactPoints, HalfEdge, Location,
    PairBudget, PointId, Segment, Snapped, TracedRing,
};
use num_bigint::BigInt;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn budget() -> PairBudget {
    PairBudget::new(10_000_000, "test pairs")
}

fn d(value: f64) -> Dyadic {
    Dyadic::try_from_f64(value).expect("finite")
}

/// `multiple * 2^exponent`.
fn scaled(multiple: i64, exponent: i64) -> Dyadic {
    Dyadic::from_parts(BigInt::from(multiple), exponent)
}

fn flat(points: &mut ExactPoints, x: f64, y: f64) -> PointId {
    points.push_flat(Point2::new(x, y)).unwrap()
}

fn snap(points: &ExactPoints, edges: &[Segment]) -> Result<Snapped, GeomError> {
    snap_round(points, edges, &mut budget(), 0, &options())
}

fn ring_of(ids: &[PointId]) -> Vec<Segment> {
    (0..ids.len())
        .map(|i| Segment {
            ends: [ids[i], ids[(i + 1) % ids.len()]],
        })
        .collect()
}

/// The rings of the class that runs along `edges` in order, as positions.
fn traced(snapped: &Snapped, chain: &[(usize, bool)]) -> Vec<(Vec<Point2>, bool)> {
    let boundary = snapped.class_boundary(chain).unwrap();
    let (rings, regions) = trace(
        &snapped.points,
        &snapped.arrangement,
        &boundary.half_edges,
        &|_| false,
        &mut budget(),
        &options(),
    )
    .unwrap();
    check_winding(
        &snapped.points,
        &snapped.arrangement,
        &rings,
        &regions,
        &mut budget(),
        &options(),
    )
    .unwrap();
    rings
        .iter()
        .map(|ring| {
            let at: Vec<Point2> = ring
                .vertices
                .iter()
                .map(|&v| snapped.point(v).unwrap())
                .collect();
            (canonical_rotation(at), ring.hole)
        })
        .collect()
}

fn forward(count: usize) -> Vec<(usize, bool)> {
    (0..count).map(|i| (i, true)).collect()
}

#[test]
fn flat_edges_snap_to_themselves() {
    let mut points = ExactPoints::new();
    let ids: Vec<PointId> = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
        .iter()
        .map(|&(x, y)| flat(&mut points, x, y))
        .collect();
    let edges = ring_of(&ids);
    let snapped = snap(&points, &edges).unwrap();
    assert_eq!(snapped.displacement.get(), 0.0);
    assert_eq!(snapped.rerouted, 0);
    assert!(snapped.routes.iter().all(|route| route.len() == 1));
    let rings = traced(&snapped, &forward(4));
    assert_eq!(rings.len(), 1);
    assert_eq!(rings[0].0.len(), 4);
    assert!(!rings[0].1);
}

#[test]
fn a_thin_triangle_that_rounding_would_flip_collapses() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 1.0, 0.0);
    // The apex lies 2^-1080 above the base's midpoint and rounds onto it.
    let c = points.push([d(1.5), scaled(3, -1080), d(3.0)]).unwrap();
    assert!(points.flat(c).is_none());
    let edges = ring_of(&[a, b, c]);
    let snapped = snap(&points, &edges).unwrap();
    assert_eq!(snapped.rerouted, 1, "the base runs through the apex's cell");
    let boundary = snapped.class_boundary(&forward(3)).unwrap();
    assert!(boundary.half_edges.is_empty());
    assert_eq!(boundary.survives, [false, false, false]);
}

#[test]
fn a_vertex_rounding_onto_a_non_adjacent_edge_becomes_a_shared_vertex() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 4.0, 0.0);
    let c = flat(&mut points, 4.0, 4.0);
    let v = points.push([d(6.0), scaled(3, -1080), d(3.0)]).unwrap();
    let edges = ring_of(&[a, b, c, v]);
    let snapped = snap(&points, &edges).unwrap();
    assert_eq!(snapped.rerouted, 1);
    assert_eq!(
        snapped.routes[0].len(),
        2,
        "the base is split at the vertex"
    );
    let boundary = snapped.class_boundary(&forward(4)).unwrap();
    assert_eq!(boundary.survives, [true, true, true, false]);
    let rings = traced(&snapped, &forward(4));
    assert_eq!(rings.len(), 1);
    assert_eq!(
        rings[0].0,
        [(2.0, 0.0), (4.0, 0.0), (4.0, 4.0)].map(|(x, y)| Point2::new(x, y))
    );
}

#[test]
fn an_island_that_rounds_to_one_point_collapses_inside_a_surviving_ring() {
    let mut points = ExactPoints::new();
    let corners: Vec<PointId> = [(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (0.0, 8.0)]
        .iter()
        .map(|&(x, y)| flat(&mut points, x, y))
        .collect();
    let mut edges = ring_of(&corners);
    // A clockwise hole of three points within 2^-60 of (1, 1).
    let hole: Vec<PointId> = [(0, 0), (0, 1), (1, 0)]
        .iter()
        .map(|&(i, j)| {
            points
                .push([
                    d(3.0).add(&scaled(i, -60)),
                    d(3.0).add(&scaled(j, -60)),
                    d(3.0),
                ])
                .unwrap()
        })
        .collect();
    edges.extend(ring_of(&hole));
    let snapped = snap(&points, &edges).unwrap();
    let chain = forward(7);
    let boundary = snapped.class_boundary(&chain).unwrap();
    assert_eq!(
        boundary.survives,
        [true, true, true, true, false, false, false]
    );
    let rings = traced(&snapped, &chain);
    assert_eq!(rings.len(), 1, "the hole is gone");
    assert!(snapped.displacement.get() > 0.0);
}

#[test]
fn the_displacement_of_a_reroute_is_the_half_width_of_the_cell() {
    let mut points = ExactPoints::new();
    // Every end is a binary64 point, so only the reroute moves anything. The
    // first edge meets (1, 1) at (1, 1 + 2^-53), a tie that the even centre
    // 1.0 owns, where the spacing changes: 2^-53 below 1.0, 2^-52 above.
    let a = flat(&mut points, 0.5, 0.5);
    let b = flat(&mut points, 1.5, 1.5 + 2f64.powi(-52));
    let f = flat(&mut points, 1.0, 1.0);
    let w = flat(&mut points, 1.0, -1.0);
    let edges = [Segment { ends: [a, b] }, Segment { ends: [f, w] }];
    let snapped = snap(&points, &edges).unwrap();
    assert_eq!(snapped.rerouted, 1);
    assert_eq!(snapped.routes[0].len(), 2);
    assert_eq!(snapped.displacement.get(), 2f64.powi(-53));
}

#[test]
fn cells_and_fragments_count_against_the_memory_budget() {
    let mut points = ExactPoints::new();
    let ids: Vec<PointId> = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
        .iter()
        .map(|&(x, y)| flat(&mut points, x, y))
        .collect();
    let edges = ring_of(&ids);
    // The documented charges: 256 bytes per cell, 64 per fragment.
    let (cells, fragments) = (4 * 256, 4 * 64);
    let run = |bytes: usize| {
        snap_round(
            &points,
            &edges,
            &mut budget(),
            0,
            &options().with_memory_budget(bytes),
        )
    };
    run(cells + fragments).unwrap();
    for short in [cells + fragments - 1, cells - 1] {
        assert!(
            matches!(
                run(short),
                Err(GeomError::BudgetExceeded { resource: "memory" })
            ),
            "{short}"
        );
    }
}

#[test]
fn an_edge_visits_its_cells_in_order() {
    let mut points = ExactPoints::new();
    let a = flat(&mut points, 0.0, 0.0);
    let b = flat(&mut points, 3.0, 3.0);
    let mut edges = vec![Segment { ends: [a, b] }];
    // Two vertices just above the diagonal, whose cells it passes through.
    for (x, top) in [(1.0, 3.0), (2.0, 4.0)] {
        let v = points
            .push([d(3.0 * x), d(3.0 * x).add(&scaled(1, -60)), d(3.0)])
            .unwrap();
        let w = flat(&mut points, x, top);
        edges.push(Segment { ends: [v, w] });
    }
    let snapped = snap(&points, &edges).unwrap();
    assert_eq!(snapped.rerouted, 1);
    let mut at = vec![];
    let mut previous = None;
    for &(piece, forward) in &snapped.routes[0] {
        let ends = snapped.arrangement.pieces[piece.0].ends;
        let (from, to) = if forward {
            (ends[0], ends[1])
        } else {
            (ends[1], ends[0])
        };
        assert!(previous.is_none_or(|p| p == from), "the route is connected");
        at.push(snapped.point(from).unwrap());
        previous = Some(to);
    }
    at.push(snapped.point(previous.unwrap()).unwrap());
    assert_eq!(
        at,
        [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0)].map(|(x, y)| Point2::new(x, y))
    );
}

#[test]
fn crossing_edges_and_doubled_chains_are_refused() {
    let mut points = ExactPoints::new();
    let ids: Vec<PointId> = [(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)]
        .iter()
        .map(|&(x, y)| flat(&mut points, x, y))
        .collect();
    let crossing = [
        Segment {
            ends: [ids[0], ids[1]],
        },
        Segment {
            ends: [ids[2], ids[3]],
        },
    ];
    match snap(&points, &crossing) {
        Err(GeomError::Degenerate(m)) => {
            assert_eq!(m, "snap rounding left two fragments crossing");
        }
        other => panic!("{other:?}"),
    }
    let single = [crossing[0]];
    let snapped = snap(&points, &single).unwrap();
    match snapped.class_boundary(&[(0, true), (0, true)]) {
        Err(GeomError::Degenerate(m)) => {
            assert_eq!(m, "snap rounding folded a boundary onto itself");
        }
        other => panic!("{other:?}"),
    }
    let cancelled = snapped.class_boundary(&[(0, true), (0, false)]).unwrap();
    assert!(cancelled.half_edges.is_empty());
    assert_eq!(cancelled.survives, [false, false]);
}

/// Rings given as lists of positions, each traversed as listed.
fn winding_of(rings: &[&[(f64, f64)]]) -> Result<(), GeomError> {
    let mut points = ExactPoints::new();
    let mut segments = Vec::new();
    for ring in rings {
        let ids: Vec<PointId> = ring.iter().map(|&(x, y)| flat(&mut points, x, y)).collect();
        segments.extend(ring_of(&ids));
    }
    let arrangement = arrange(&mut points, &segments, &mut budget(), 0, &options()).unwrap();
    let half_edges: Vec<HalfEdge> = arrangement
        .paths
        .iter()
        .flatten()
        .map(|&(piece, forward)| HalfEdge { piece, forward })
        .collect();
    let (rings, regions) = trace(
        &points,
        &arrangement,
        &half_edges,
        &|_| false,
        &mut budget(),
        &options(),
    )
    .unwrap();
    check_winding(
        &points,
        &arrangement,
        &rings,
        &regions,
        &mut budget(),
        &options(),
    )
}

#[test]
fn a_hole_inside_a_hole_is_refused_and_an_island_is_not() {
    let outer: &[(f64, f64)] = &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
    let clockwise =
        |lo: f64, hi: f64| -> Vec<(f64, f64)> { vec![(lo, lo), (lo, hi), (hi, hi), (hi, lo)] };
    let (hole, inner_hole) = (clockwise(2.0, 8.0), clockwise(4.0, 6.0));
    winding_of(&[outer]).unwrap();
    winding_of(&[outer, &hole]).unwrap();
    match winding_of(&[outer, &hole, &inner_hole]) {
        Err(GeomError::Degenerate(m)) => {
            assert_eq!(m, "a class winds twice after snap rounding");
        }
        other => panic!("{other:?}"),
    }
    let mut island = clockwise(4.0, 6.0);
    island.reverse();
    winding_of(&[outer, &hole, &island]).unwrap();
}

/// A small seeded generator; the tests need spread, not quality.
struct Lcg(u64);

impl Lcg {
    fn below(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) % n
    }
}

/// Whether the sample is farther than `delta` from every exact edge box.
fn clear_of(points: &ExactPoints, edges: &[Segment], at: Point2, delta: f64) -> bool {
    edges.iter().all(|edge| {
        let [a, b] = edge.ends.map(|end| points.bounds(end));
        let (x_lo, y_lo) = (a[0].min(b[0]), a[1].min(b[1]));
        let (x_hi, y_hi) = (a[2].max(b[2]), a[3].max(b[3]));
        // `2 delta` of slack: the sample's own box and the rounding.
        let slack = 2.0 * delta;
        at.x < x_lo - slack || at.x > x_hi + slack || at.y < y_lo - slack || at.y > y_hi + slack
    })
}

fn inside_rings(points: &mut ExactPoints, vertices: &[Vec<PointId>], at: Point2) -> Option<bool> {
    let q = points.push_flat(at).unwrap();
    let mut winding = 0;
    for ring in vertices {
        match points.locate(q, ring).unwrap() {
            Location::Boundary => return None,
            Location::Inside => winding += 1,
            Location::Outside => {}
        }
    }
    // An outer ring and its holes: inside an odd number of rings.
    Some(winding % 2 == 1)
}

#[test]
fn snapping_small_triangles_keeps_the_topology_and_the_winding() {
    let mut rng = Lcg(0xfeed);
    let bases = [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (0.5, 0.5), (-1.0, 2.0)];
    let (mut compared, mut collapsed, mut rerouted) = (0usize, 0usize, 0usize);
    for case in 0..2000 {
        let (bx, by) = bases[rng.below(5) as usize];
        let (bx, by) = (
            if rng.below(2) == 0 && bx == 0.0 {
                -0.0
            } else {
                bx
            },
            by,
        );
        let sign = [1 - 2 * rng.below(2) as i64, 1 - 2 * rng.below(2) as i64];
        let exponent = 40 + rng.below(31) as i64;
        let mut points = ExactPoints::new();
        let mut known: BTreeMap<(i64, i64), PointId> = BTreeMap::new();
        let mut vertex = |points: &mut ExactPoints, g: (i64, i64), rng: &mut Lcg| {
            *known.entry(g).or_insert_with(|| {
                let w = [3i64, 5, 7][rng.below(3) as usize];
                let jitter = |rng: &mut Lcg| rng.below(3) as i64 - 1;
                let coordinate = |base: f64, grid: i64, sign: i64, j: i64| {
                    d(base)
                        .mul(&d(w as f64))
                        .add(&scaled(w * grid * sign + j, -exponent))
                };
                let (jx, jy) = (jitter(rng), jitter(rng));
                points
                    .push([
                        coordinate(bx, g.0, sign[0], jx),
                        coordinate(by, g.1, sign[1], jy),
                        d(w as f64),
                    ])
                    .unwrap()
            })
        };
        // Triangles in a row: apart, or touching at the previous one's
        // second vertex.
        let mut triangles: Vec<[PointId; 3]> = Vec::new();
        let mut origin = (0i64, 0i64);
        for index in 0..2 + rng.below(3) as usize {
            if index > 0 {
                origin = if rng.below(2) == 0 {
                    (origin.0 + 4, origin.1)
                } else {
                    (origin.0 + 10, origin.1)
                };
            }
            let mut corners = [origin, (origin.0 + 4, origin.1), (origin.0, origin.1 + 4)];
            if sign[0] * sign[1] < 0 {
                corners.swap(1, 2);
            }
            let ids = corners.map(|g| vertex(&mut points, g, &mut rng));
            triangles.push(ids);
        }
        let mut edge_of: BTreeMap<(PointId, PointId), usize> = BTreeMap::new();
        let mut edges: Vec<Segment> = Vec::new();
        let mut chains: Vec<Vec<(usize, bool)>> = Vec::new();
        for triangle in &triangles {
            let mut chain = Vec::new();
            for i in 0..3 {
                let (u, v) = (triangle[i], triangle[(i + 1) % 3]);
                let key = (u.min(v), u.max(v));
                let index = *edge_of.entry(key).or_insert_with(|| {
                    edges.push(Segment {
                        ends: [key.0, key.1],
                    });
                    edges.len() - 1
                });
                chain.push((index, u < v));
            }
            chains.push(chain);
        }
        let snapped = snap_round(&points, &edges, &mut budget(), 0, &options())
            .unwrap_or_else(|e| panic!("case {case}: {e}"));
        rerouted += usize::from(snapped.rerouted > 0);
        let delta = snapped.displacement.get();

        for (triangle, chain) in triangles.iter().zip(&chains) {
            let boundary = snapped
                .class_boundary(chain)
                .unwrap_or_else(|e| panic!("case {case}: {e}"));
            collapsed += usize::from(boundary.half_edges.is_empty());
            let (rings, regions) = trace(
                &snapped.points,
                &snapped.arrangement,
                &boundary.half_edges,
                &|_| false,
                &mut budget(),
                &options(),
            )
            .unwrap_or_else(|e| panic!("case {case}: {e}"));
            check_winding(
                &snapped.points,
                &snapped.arrangement,
                &rings,
                &regions,
                &mut budget(),
                &options(),
            )
            .unwrap_or_else(|e| panic!("case {case}: {e}"));

            // Samples at half-grid positions, rounded to binary64.
            let mut after = snapped.points.clone();
            let mut before = points.clone();
            let vertices = |ring: &TracedRing| -> Vec<PointId> {
                ring.vertices
                    .iter()
                    .map(|&v| snapped.arrangement.vertices[v.0].point)
                    .collect()
            };
            let rings_after: Vec<Vec<PointId>> = rings.iter().map(vertices).collect();
            let exact_ring: Vec<PointId> = triangle.to_vec();
            for gx in -2..=44 {
                for gy in -2..=8 {
                    let sample = {
                        let coordinate = |base: f64, grid: i64, sign: i64| {
                            d(base).add(&scaled(2 * grid * sign + sign, -exponent - 1))
                        };
                        let mut store = ExactPoints::new();
                        let id = store
                            .push([
                                coordinate(bx, gx, sign[0]),
                                coordinate(by, gy, sign[1]),
                                d(1.0),
                            ])
                            .unwrap();
                        store.round(id).unwrap()
                    };
                    if !clear_of(&points, &edges, sample, delta) {
                        continue;
                    }
                    let want = inside_rings(&mut before, std::slice::from_ref(&exact_ring), sample);
                    let got = inside_rings(&mut after, &rings_after, sample);
                    if let (Some(want), Some(got)) = (want, got) {
                        compared += 1;
                        assert_eq!(got, want, "case {case}: sample {sample:?}, delta {delta}");
                    }
                }
            }
        }
    }
    assert!(compared > 20_000, "only {compared} samples compared");
    assert!(collapsed > 100, "only {collapsed} collapsed triangles");
    assert!(rerouted > 40, "only {rerouted} rerouted scenes");
}
