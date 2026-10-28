// SPDX-License-Identifier: MPL-2.0

//! `triangulate_domain` keeps exactly the triangles inside the outline and
//! outside the holes (#353).
//!
//! Every accepted case is held to an independent integer oracle
//! (`common::assert_domain`): the triangles' centroids are inside, their
//! area is the signed ring area, and each point is a corner exactly when it
//! belongs to the domain.

mod common;

use axiolid_core::Point2;
use axiolid_triangulate::{
    triangulate_domain, Constraint, Quality, Ring, Triangulation, TriangulationError as E,
};
use common::{assert_constrained_delaunay, assert_domain, triangle_set, Rng};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

/// Append the corners of an axis-aligned rectangle, counter-clockwise, and
/// return their indices.
fn rect(pts: &mut Vec<Point2>, x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<u32> {
    let base = pts.len() as u32;
    pts.extend([p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)]);
    (base..base + 4).collect()
}

fn reversed(v: &[u32]) -> Vec<u32> {
    v.iter().rev().copied().collect()
}

fn domain(
    name: &str,
    points: &[Point2],
    constraints: &[Constraint],
    rings: &[Ring],
) -> Triangulation {
    let tri =
        triangulate_domain(points, constraints, rings).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_domain(name, points, rings, &tri);
    tri
}

#[test]
fn a_concave_outline_keeps_only_its_inside() {
    // An L, with a point in the notch that the hull covers.
    let pts = vec![
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 2.0),
        p(2.0, 2.0),
        p(2.0, 4.0),
        p(0.0, 4.0),
        p(3.0, 3.0),
    ];
    let tri = domain("L", &pts, &[], &[Ring::outer([0, 1, 2, 3, 4, 5])]);
    assert_eq!(tri.triangle_count(), 4);
    assert_eq!(tri.outside_points(), [6]);
    assert_eq!(tri.points().len(), 7);
    assert_constrained_delaunay(&tri);
}

#[test]
fn a_hole_is_left_out() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    let hole = rect(&mut pts, 3.0, 3.0, 7.0, 7.0);
    pts.extend([p(5.0, 5.0), p(1.0, 1.0)]);
    let tri = domain("hole", &pts, &[], &[Ring::outer(outer), Ring::hole(hole)]);
    assert_eq!(tri.outside_points(), [8]);
}

#[test]
fn an_island_in_a_hole_is_kept() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 20.0, 20.0);
    let hole = rect(&mut pts, 4.0, 4.0, 16.0, 16.0);
    let island = rect(&mut pts, 8.0, 8.0, 12.0, 12.0);
    pts.extend([p(2.0, 2.0), p(6.0, 6.0), p(10.0, 10.0)]);
    let tri = domain(
        "island",
        &pts,
        &[],
        &[Ring::outer(outer), Ring::hole(hole), Ring::outer(island)],
    );
    // The moat point (6, 6) is in the hole; the island centre is kept.
    assert_eq!(tri.outside_points(), [13]);
}

#[test]
fn rings_may_come_either_way_round() {
    // A reflex first vertex, an outer ring and a hole, all four orientation
    // combinations: the same triangles each time.
    let pts = vec![
        p(5.0, 2.0),
        p(10.0, 0.0),
        p(10.0, 10.0),
        p(0.0, 10.0),
        p(0.0, 0.0),
        p(4.0, 5.0),
        p(6.0, 5.0),
        p(5.0, 7.0),
    ];
    let outer = vec![0, 1, 2, 3, 4];
    let hole = vec![5, 6, 7];
    let mut sets = Vec::new();
    for flip_outer in [false, true] {
        for flip_hole in [false, true] {
            let o = if flip_outer {
                reversed(&outer)
            } else {
                outer.clone()
            };
            let h = if flip_hole {
                reversed(&hole)
            } else {
                hole.clone()
            };
            let tri = domain("orientation", &pts, &[], &[Ring::outer(o), Ring::hole(h)]);
            sets.push(triangle_set(&tri));
        }
    }
    assert!(sets.windows(2).all(|w| w[0] == w[1]));
}

#[test]
fn a_point_on_a_ring_edge_splits_it_and_is_kept() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 4.0, 4.0);
    pts.push(p(2.0, 0.0));
    let tri = domain("on edge", &pts, &[], &[Ring::outer(outer)]);
    assert!(tri.constraints().contains(&Constraint::new(0, 4)));
    assert!(tri.constraints().contains(&Constraint::new(1, 4)));
    assert!(tri.outside_points().is_empty());
}

#[test]
fn a_hole_touching_the_outline_at_a_vertex() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 8.0, 8.0);
    pts.extend([p(4.0, 1.0), p(1.0, 4.0)]);
    domain(
        "touch",
        &pts,
        &[],
        &[Ring::outer(outer), Ring::hole([0, 4, 5])],
    );
}

#[test]
fn a_hole_vertex_on_the_outline_edge() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 8.0, 8.0);
    pts.extend([p(4.0, 0.0), p(2.0, 3.0), p(6.0, 3.0)]);
    let tri = domain(
        "vertex on edge",
        &pts,
        &[],
        &[Ring::outer(outer), Ring::hole([4, 6, 5])],
    );
    assert!(tri.constraints().contains(&Constraint::new(0, 4)));
    assert!(tri.constraints().contains(&Constraint::new(1, 4)));
}

#[test]
fn a_hole_sharing_an_edge_with_the_outline() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    pts.extend([p(4.0, 0.0), p(4.0, 3.0), p(0.0, 3.0)]);
    domain(
        "shared edge",
        &pts,
        &[],
        &[Ring::outer(outer), Ring::hole([0, 4, 5, 6])],
    );
}

#[test]
fn two_holes_sharing_an_edge_are_one_void() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    pts.extend([
        p(2.0, 2.0),
        p(4.0, 2.0),
        p(4.0, 6.0),
        p(2.0, 6.0),
        p(6.0, 2.0),
        p(6.0, 6.0),
    ]);
    domain(
        "two holes",
        &pts,
        &[],
        &[
            Ring::outer(outer),
            Ring::hole([4, 5, 6, 7]),
            Ring::hole([5, 8, 9, 6]),
        ],
    );
}

#[test]
fn outer_rings_sharing_an_edge_are_one_domain() {
    let mut pts = Vec::new();
    let left = rect(&mut pts, 0.0, 0.0, 4.0, 4.0);
    pts.extend([p(8.0, 0.0), p(8.0, 4.0)]);
    let tri = domain(
        "two outers",
        &pts,
        &[],
        &[Ring::outer(left), Ring::outer([1, 4, 5, 2])],
    );
    assert_eq!(tri.triangle_count(), 4);
}

#[test]
fn disjoint_outer_rings_are_both_kept() {
    let mut pts = Vec::new();
    let a = rect(&mut pts, 0.0, 0.0, 4.0, 4.0);
    let b = rect(&mut pts, 10.0, 0.0, 14.0, 4.0);
    let tri = domain("disjoint", &pts, &[], &[Ring::outer(a), Ring::outer(b)]);
    assert_eq!(tri.triangle_count(), 4);
}

#[test]
fn a_breakline_does_not_change_the_domain() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    let hole = rect(&mut pts, 4.0, 4.0, 6.0, 6.0);
    pts.extend([p(1.0, 5.0), p(9.0, 5.0)]);
    let rings = [Ring::outer(outer), Ring::hole(hole)];
    let plain = domain("plain", &pts, &[], &rings);
    let line = Constraint::new(8, 9);
    let err = triangulate_domain(&pts, &[line], &rings).unwrap_err();
    assert!(matches!(err, E::CrossingConstraints { .. }), "{err:?}");
    // A breakline that stays in the domain only adds an edge.
    let side = Constraint::new(8, 0);
    let with = domain("side", &pts, &[side], &rings);
    assert!(with.constraints().contains(&side));
    assert_eq!(with.triangle_count(), plain.triangle_count());
}

#[test]
fn a_breakline_crossing_a_hole_boundary_is_refused() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    let hole = rect(&mut pts, 4.0, 4.0, 6.0, 6.0);
    pts.extend([p(2.0, 5.0), p(5.0, 5.0)]);
    let err = triangulate_domain(
        &pts,
        &[Constraint::new(8, 9)],
        &[Ring::outer(outer), Ring::hole(hole)],
    )
    .unwrap_err();
    assert!(matches!(err, E::CrossingConstraints { .. }), "{err:?}");
}

#[test]
fn a_breakline_through_a_ring_vertex_is_split_there() {
    let mut pts = Vec::new();
    rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    pts.extend([p(5.0, 0.0), p(5.0, 5.0)]);
    // The breakline (0, 0)-(10, 0) runs along the bottom edge through the
    // ring vertex (5, 0), and another leaves that vertex into the domain.
    let tri = domain(
        "through a vertex",
        &pts,
        &[Constraint::new(4, 5), Constraint::new(0, 1)],
        &[Ring::outer([0, 4, 1, 2, 3])],
    );
    assert!(tri.constraints().contains(&Constraint::new(4, 5)));
    assert!(tri.constraints().contains(&Constraint::new(0, 4)));
}

#[test]
fn a_breakline_outside_the_domain_is_dropped() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 4.0, 4.0);
    pts.extend([p(8.0, 0.0), p(8.0, 4.0), p(12.0, 2.0)]);
    let dropped = Constraint::new(4, 6);
    let tri = domain(
        "dropped",
        &pts,
        &[dropped, Constraint::new(0, 2)],
        &[Ring::outer(outer)],
    );
    assert!(!tri.constraints().contains(&dropped));
    assert!(tri.constraints().contains(&Constraint::new(0, 2)));
    assert_eq!(tri.outside_points(), [4, 5, 6]);
}

#[test]
fn a_self_crossing_ring_is_refused() {
    let pts = vec![p(0.0, 0.0), p(4.0, 4.0), p(4.0, 0.0), p(0.0, 4.0)];
    let err = triangulate_domain(&pts, &[], &[Ring::outer([0, 1, 2, 3])]).unwrap_err();
    assert!(matches!(err, E::CrossingConstraints { .. }), "{err:?}");
}

#[test]
fn crossing_rings_are_refused() {
    let mut pts = Vec::new();
    let a = rect(&mut pts, 0.0, 0.0, 6.0, 6.0);
    let b = rect(&mut pts, 3.0, 3.0, 9.0, 9.0);
    let err = triangulate_domain(&pts, &[], &[Ring::outer(a), Ring::outer(b)]).unwrap_err();
    assert!(matches!(err, E::CrossingConstraints { .. }), "{err:?}");
}

#[test]
fn a_hole_outside_every_outer_ring_is_refused() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 4.0, 4.0);
    let hole = rect(&mut pts, 8.0, 0.0, 12.0, 4.0);
    assert_eq!(
        triangulate_domain(&pts, &[], &[Ring::outer(outer), Ring::hole(hole)]).unwrap_err(),
        E::MisnestedRing {
            ring: 1,
            winding: -1
        }
    );
}

#[test]
fn a_hole_inside_a_hole_is_refused() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 20.0, 20.0);
    let hole = rect(&mut pts, 4.0, 4.0, 16.0, 16.0);
    let inner = rect(&mut pts, 8.0, 8.0, 12.0, 12.0);
    assert_eq!(
        triangulate_domain(
            &pts,
            &[],
            &[Ring::outer(outer), Ring::hole(hole), Ring::hole(inner)]
        )
        .unwrap_err(),
        E::MisnestedRing {
            ring: 2,
            winding: -1
        }
    );
}

#[test]
fn a_hole_listed_twice_is_refused() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    let hole = rect(&mut pts, 3.0, 3.0, 7.0, 7.0);
    assert_eq!(
        triangulate_domain(
            &pts,
            &[],
            &[
                Ring::outer(outer),
                Ring::hole(hole.clone()),
                Ring::hole(hole)
            ]
        )
        .unwrap_err(),
        E::MisnestedRing {
            ring: 1,
            winding: -1
        }
    );
}

#[test]
fn an_outer_ring_inside_an_outer_ring_is_refused() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    let inner = rect(&mut pts, 3.0, 3.0, 7.0, 7.0);
    assert_eq!(
        triangulate_domain(&pts, &[], &[Ring::outer(outer), Ring::outer(inner)]).unwrap_err(),
        E::MisnestedRing {
            ring: 1,
            winding: 2
        }
    );
}

#[test]
fn no_rings_is_an_empty_domain() {
    let mut pts = Vec::new();
    rect(&mut pts, 0.0, 0.0, 4.0, 4.0);
    assert_eq!(
        triangulate_domain(&pts, &[], &[]).unwrap_err(),
        E::EmptyDomain
    );
}

#[test]
fn a_hole_equal_to_its_outline_is_an_empty_domain() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 4.0, 4.0);
    assert_eq!(
        triangulate_domain(&pts, &[], &[Ring::outer(outer.clone()), Ring::hole(outer)])
            .unwrap_err(),
        E::EmptyDomain
    );
}

#[test]
fn coincident_points_in_rings_are_merged() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 6.0, 6.0);
    // A second copy of corner 2, which the ring uses instead.
    pts.push(p(6.0, 6.0));
    let ring = vec![outer[0], outer[1], 4, outer[3]];
    let tri = domain("merged", &pts, &[], &[Ring::outer(ring)]);
    assert_eq!(tri.merged_points().len(), 1);
    assert!(tri.outside_points().is_empty());
}

#[test]
fn the_domain_is_constrained_delaunay_inside() {
    let mut r = Rng(0x5eed_1234_abcd_ef01);
    let mut pts = vec![
        p(0.0, 0.0),
        p(40.0, 0.0),
        p(40.0, 20.0),
        p(20.0, 20.0),
        p(20.0, 40.0),
        p(0.0, 40.0),
    ];
    for _ in 0..60 {
        pts.push(p((r.unit() * 40.0).floor(), (r.unit() * 40.0).floor()));
    }
    let tri = domain("delaunay", &pts, &[], &[Ring::outer([0, 1, 2, 3, 4, 5])]);
    assert_constrained_delaunay(&tri);
}

#[test]
fn the_output_is_deterministic() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 10.0, 10.0);
    let hole = rect(&mut pts, 3.0, 3.0, 7.0, 6.0);
    pts.extend([p(1.0, 8.0), p(9.0, 2.0), p(5.0, 9.0)]);
    let rings = [Ring::outer(outer.clone()), Ring::hole(hole.clone())];
    let a = domain("first", &pts, &[], &rings);
    let b = domain("second", &pts, &[], &rings);
    assert_eq!(a.triangles(), b.triangles());
    let rotate = |v: &[u32], k: usize| [&v[k..], &v[..k]].concat();
    let turned = [Ring::outer(rotate(&outer, 2)), Ring::hole(rotate(&hole, 1))];
    let c = domain("rotated", &pts, &[], &turned);
    assert_eq!(triangle_set(&a), triangle_set(&c));
}

/// A star-shaped polygon about the origin: `count` jittered directions at
/// radii in `[lo, hi]`, rounded to integers and sorted by exact angle.
fn star(r: &mut Rng, count: usize, lo: f64, hi: f64) -> Vec<(i64, i64)> {
    let mut v: Vec<(i64, i64)> = (0..count)
        .map(|k| {
            let angle = (k as f64 + 0.5 * r.unit()) * core::f64::consts::TAU / count as f64;
            let radius = lo + (hi - lo) * r.unit();
            (
                (radius * angle.cos()).round() as i64,
                (radius * angle.sin()).round() as i64,
            )
        })
        .collect();
    let half = |q: &(i64, i64)| i32::from(q.1 < 0 || (q.1 == 0 && q.0 < 0));
    v.sort_by(|a, b| {
        half(a)
            .cmp(&half(b))
            .then_with(|| (a.0 * b.1 - a.1 * b.0).cmp(&0).reverse())
    });
    v
}

/// Whether consecutive vertices turn strictly left about the origin, so
/// the origin is strictly inside the star.
fn is_star(v: &[(i64, i64)]) -> bool {
    (0..v.len()).all(|i| {
        let (a, b) = (v[i], v[(i + 1) % v.len()]);
        a.0 * b.1 - a.1 * b.0 > 0
    })
}

#[test]
fn random_star_domains_match_the_oracle() {
    let mut checked = 0;
    for seed in 1..=60u64 {
        let mut r = Rng(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1);
        let outer = star(&mut r, 11, 30.0, 40.0);
        let hole = star(&mut r, 8, 9.0, 14.0);
        let island = star(&mut r, 5, 2.0, 4.0);
        if ![&outer, &hole, &island].iter().all(|s| is_star(s)) {
            continue;
        }
        let mut pts = Vec::new();
        let mut ring = |poly: &[(i64, i64)]| -> Vec<u32> {
            let base = pts.len() as u32;
            pts.extend(poly.iter().map(|&(x, y)| p(x as f64, y as f64)));
            (base..pts.len() as u32).collect()
        };
        let (o, h, i) = (ring(&outer), ring(&hole), ring(&island));
        let flip = |v: Vec<u32>, on: bool| if on { reversed(&v) } else { v };
        let rings = [
            Ring::outer(flip(o, r.unit() < 0.5)),
            Ring::hole(flip(h, r.unit() < 0.5)),
            Ring::outer(flip(i, r.unit() < 0.5)),
        ];
        for _ in 0..40 {
            pts.push(p(
                (r.unit() * 100.0).floor() - 50.0,
                (r.unit() * 100.0).floor() - 50.0,
            ));
        }
        let tri = domain(&format!("seed {seed}"), &pts, &[], &rings);
        assert_constrained_delaunay(&tri);
        checked += 1;
    }
    assert!(checked >= 40, "only {checked} seeds gave valid stars");
}

#[test]
fn refinement_stays_inside_the_domain() {
    let mut pts = Vec::new();
    let outer = rect(&mut pts, 0.0, 0.0, 12.0, 12.0);
    let hole = rect(&mut pts, 8.0, 8.0, 11.0, 11.0);
    // Two interior points close together make a fixable sliver.
    pts.extend([p(5.8, 6.0), p(6.2, 6.05), p(9.5, 9.5)]);
    let rings = [Ring::outer(outer), Ring::hole(hole)];
    let mut tri = domain("before", &pts, &[], &rings);
    assert_eq!(tri.outside_points(), [10]);
    let original = tri.points().len();
    let outcome = axiolid_triangulate::refine(
        &mut tri,
        Quality {
            min_angle_degrees: 20.0,
            max_steiner_points: 256,
        },
    )
    .unwrap();
    assert!(
        tri.points().len() > original,
        "nothing was inserted: {outcome:?}"
    );
    assert_eq!(tri.outside_points(), [10]);
    // The area is still the domain's: refinement did not re-expand to the
    // hull or fill the hole.
    let area: f64 = tri
        .triangles()
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let (a, b, c) = (
                tri.points()[t[0] as usize],
                tri.points()[t[1] as usize],
                tri.points()[t[2] as usize],
            );
            ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)) / 2.0
        })
        .sum();
    assert!((area - (144.0 - 9.0)).abs() < 1e-9, "area {area}");
    for s in &tri.points()[original..] {
        let in_hole = s.x > 8.0 && s.x < 11.0 && s.y > 8.0 && s.y < 11.0;
        assert!(!in_hole, "Steiner point {s:?} inside the hole");
    }
    assert_constrained_delaunay(&tri);
}
