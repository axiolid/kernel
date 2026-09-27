//! Region skeletons (#139): corridor ends and junctions, the wall ahead of
//! each end, and certified clearances.

use axiolid_core::Point2;
use axiolid_guarantees::Sign;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{skeleton, NodeKind, Skeleton, SkeletonError, Wall};

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
    }
}

fn polygon(outer: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Vec<Polygon> {
    vec![Polygon {
        outer: ring(outer),
        holes: holes.iter().map(|h| ring(h)).collect(),
    }]
}

const PRUNE: f64 = 1.5;

/// Every node's clearance interval holds its true distance to the walls,
/// decided exactly (squared, in dyadics), and tightly; every edge joins
/// two nodes.
fn check(region: &[Polygon], s: &Skeleton) {
    use axiolid_exact::{Arith, Dyadic};
    let x = Dyadic::from_f64;
    let below = |a: &Dyadic, b: &Dyadic| a.sub(b).sign() == Some(Sign::Negative);
    for node in &s.nodes {
        let p = node.point;
        // Exact squared distance to each wall as num / den, the least kept.
        let mut best: Option<(Dyadic, Dyadic)> = None;
        for poly in region {
            for r in std::iter::once(&poly.outer).chain(&poly.holes) {
                let n = r.points.len();
                for i in 0..n {
                    let (a, b) = (r.points[i], r.points[(i + 1) % n]);
                    let (ex, ey) = (x(b.x).sub(&x(a.x)), x(b.y).sub(&x(a.y)));
                    let (wx, wy) = (x(p.x).sub(&x(a.x)), x(p.y).sub(&x(a.y)));
                    let (vx, vy) = (x(p.x).sub(&x(b.x)), x(p.y).sub(&x(b.y)));
                    let dot = ex.mul(&wx).add(&ey.mul(&wy));
                    let len2 = ex.mul(&ex).add(&ey.mul(&ey));
                    let one = x(1.0);
                    let d = if dot.sign() != Some(Sign::Positive) {
                        (wx.mul(&wx).add(&wy.mul(&wy)), one)
                    } else if !below(&dot, &len2) {
                        (vx.mul(&vx).add(&vy.mul(&vy)), one)
                    } else {
                        let c = ex.mul(&wy).sub(&ey.mul(&wx));
                        (c.mul(&c), len2)
                    };
                    let smaller = best
                        .as_ref()
                        .is_none_or(|(bn, bd)| below(&d.0.mul(bd), &bn.mul(&d.1)));
                    if smaller {
                        best = Some(d);
                    }
                }
            }
        }
        let (num, den) = best.unwrap();
        let (lo, hi) = node.clearance;
        assert!(
            !below(&num, &x(lo).mul(&x(lo)).mul(&den)),
            "lower bound {lo} above the distance at {p:?}"
        );
        assert!(
            !below(&x(hi).mul(&x(hi)).mul(&den), &num),
            "upper bound {hi} below the distance at {p:?}"
        );
        assert!(hi - lo < 1e-12 * (1.0 + hi), "{lo} {hi}");
    }
    for &(a, b) in &s.edges {
        assert!(a < b && b < s.nodes.len());
    }
}

#[test]
fn a_dead_end_corridor_ends_at_its_end_walls() {
    // 10 long, 1 wide. Its axis runs lengthwise; each end sits about half
    // a width from the end wall, which is the wall ahead of it.
    let region = polygon(&[(0.0, 0.0), (10.0, 0.0), (10.0, 1.0), (0.0, 1.0)], &[]);
    let s = skeleton(&region, 0.1, PRUNE).unwrap();
    check(&region, &s);
    let ends = s.ends();
    assert_eq!(ends.len(), 2, "{s:?}");
    assert!(s.junctions().is_empty());
    let mut xs: Vec<(f64, Wall)> = ends
        .iter()
        .map(|&e| (s.nodes[e].point.x, s.nodes[e].ahead.unwrap()))
        .collect();
    xs.sort_by(|a, b| a.0.total_cmp(&b.0));
    // Edge 3 runs (0, 1) -> (0, 0): the left end wall; edge 1 the right.
    assert!(xs[0].0 < 1.0 && xs[1].0 > 9.0, "{xs:?}");
    assert_eq!(
        xs[0].1,
        Wall {
            polygon: 0,
            ring: 0,
            edge: 3
        }
    );
    assert_eq!(
        xs[1].1,
        Wall {
            polygon: 0,
            ring: 0,
            edge: 1
        }
    );
    // Along the middle the clearance is half the width.
    let mid = s
        .nodes
        .iter()
        .filter(|n| (n.point.x - 5.0).abs() < 1.0)
        .map(|n| n.clearance.1)
        .fold(0.0, f64::max);
    assert!((mid - 0.5).abs() < 0.05, "{mid}");
}

#[test]
fn an_l_corridor_has_two_ends() {
    let region = polygon(
        &[
            (0.0, 0.0),
            (8.0, 0.0),
            (8.0, 6.0),
            (7.0, 6.0),
            (7.0, 1.0),
            (0.0, 1.0),
        ],
        &[],
    );
    let s = skeleton(&region, 0.1, PRUNE).unwrap();
    check(&region, &s);
    let ends = s.ends();
    assert_eq!(ends.len(), 2, "{s:?}");
    let walls: Vec<Wall> = ends.iter().map(|&e| s.nodes[e].ahead.unwrap()).collect();
    // The far ends: the left wall (edge 5) and the top wall (edge 2).
    assert!(
        walls.contains(&Wall {
            polygon: 0,
            ring: 0,
            edge: 5
        }),
        "{walls:?}"
    );
    assert!(
        walls.contains(&Wall {
            polygon: 0,
            ring: 0,
            edge: 2
        }),
        "{walls:?}"
    );
}

#[test]
fn a_t_junction_has_three_ends_and_meets_once() {
    let region = polygon(
        &[
            (0.0, 0.0),
            (9.0, 0.0),
            (9.0, 1.0),
            (5.0, 1.0),
            (5.0, 6.0),
            (4.0, 6.0),
            (4.0, 1.0),
            (0.0, 1.0),
        ],
        &[],
    );
    let s = skeleton(&region, 0.1, PRUNE).unwrap();
    check(&region, &s);
    assert_eq!(s.ends().len(), 3, "{s:?}");
    let junctions = s.junctions();
    assert!(!junctions.is_empty());
    for j in junctions {
        let p = s.nodes[j].point;
        assert!((p.x - 4.5).abs() < 0.6 && p.y < 1.5, "{p:?}");
    }
}

#[test]
fn a_loop_round_a_pillar_has_no_ends() {
    // The corners' spurs go; the corridor round the pillar is a cycle.
    let region = polygon(
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        &[&[(4.0, 4.0), (4.0, 6.0), (6.0, 6.0), (6.0, 4.0)]],
    );
    let s = skeleton(&region, 0.2, PRUNE).unwrap();
    check(&region, &s);
    assert!(
        s.ends().is_empty(),
        "{:?}",
        s.ends()
            .iter()
            .map(|&e| s.nodes[e].point)
            .collect::<Vec<_>>()
    );
    assert!(s.nodes.iter().all(|n| n.kind != NodeKind::Isolated));
    // Without pruning, the four corners are ends.
    let raw = skeleton(&region, 0.2, 0.0).unwrap();
    assert!(raw.ends().len() >= 4);
}

#[test]
fn nodes_lie_inside() {
    let region = polygon(
        &[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 1.0),
            (1.0, 1.0),
            (1.0, 5.0),
            (0.0, 5.0),
        ],
        &[],
    );
    let s = skeleton(&region, 0.05, PRUNE).unwrap();
    check(&region, &s);
    for n in &s.nodes {
        let p = n.point;
        let in_l = (0.0..=6.0).contains(&p.x) && (0.0..=1.0).contains(&p.y)
            || (0.0..=1.0).contains(&p.x) && (0.0..=5.0).contains(&p.y);
        assert!(in_l, "{p:?}");
    }
}

#[test]
fn parameters_are_checked() {
    let region = polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], &[]);
    for (spacing, prune) in [(0.0, 1.0), (-1.0, 1.0), (0.1, -1.0), (f64::NAN, 1.0)] {
        assert_eq!(
            skeleton(&region, spacing, prune),
            Err(SkeletonError::InvalidParameter)
        );
    }
}

#[test]
fn the_wall_ahead_is_the_nearest_one() {
    // A corridor ending at x = 5, and a separate room beyond it: the path
    // runs into the corridor's own end wall, not the room's.
    let mut region = polygon(&[(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (0.0, 1.0)], &[]);
    region.extend(polygon(
        &[(7.0, 0.0), (8.0, 0.0), (8.0, 1.0), (7.0, 1.0)],
        &[],
    ));
    let s = skeleton(&region, 0.1, PRUNE).unwrap();
    check(&region, &s);
    let right = s
        .ends()
        .into_iter()
        .filter(|&e| (4.0..5.0).contains(&s.nodes[e].point.x))
        .map(|e| s.nodes[e].ahead.unwrap())
        .collect::<Vec<_>>();
    assert_eq!(
        right,
        vec![Wall {
            polygon: 0,
            ring: 0,
            edge: 1
        }]
    );
}
