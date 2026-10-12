// SPDX-License-Identifier: MPL-2.0

//! The unconstrained triangulation covers the exact convex hull (#299).
//!
//! The build used to bracket the input with a finite super triangle and
//! lost the thin hull triangles whose circumcircles reached its corners:
//! 3,977 instead of 3,979 triangles for 2,000 random points, with nothing
//! reported. Every case here is held to the count a triangulation of `n`
//! points with `h` of them on the hull boundary must have, `2n - 2 - h`,
//! with `h` counted by an exact convex hull, and to a boundary that is that
//! hull.

use std::collections::{HashMap, HashSet};

use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{incircle, orient2d};
use axiolid_triangulate::{triangulate, Triangulation};

fn sign(c: Certified) -> Sign {
    match c {
        Certified::Certain { sign, .. } => sign,
        other => panic!("undecided predicate: {other:?}"),
    }
}

/// A seeded xorshift stream of uniform values in `[0, 1)`.
struct Rng(u64);

impl Rng {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn uniform(count: usize, seed: u64) -> Vec<Point2> {
    let mut r = Rng(seed);
    (0..count)
        .map(|_| Point2::new(r.unit(), r.unit()))
        .collect()
}

/// Every input point on the boundary of the exact convex hull, collinear
/// ones included.
fn hull_boundary(points: &[Point2]) -> HashSet<(u64, u64)> {
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    sorted.dedup();
    let mut corners: Vec<Point2> = Vec::new();
    for pass in 0..2 {
        let start = corners.len();
        let iter: Box<dyn Iterator<Item = &Point2>> = if pass == 0 {
            Box::new(sorted.iter())
        } else {
            Box::new(sorted.iter().rev())
        };
        for &p in iter {
            while corners.len() >= start + 2
                && sign(orient2d(
                    corners[corners.len() - 2],
                    corners[corners.len() - 1],
                    p,
                )) != Sign::Positive
            {
                corners.pop();
            }
            corners.push(p);
        }
        corners.pop();
    }
    let within = |a: f64, b: f64, x: f64| a.min(b) <= x && x <= a.max(b);
    let mut on = HashSet::new();
    for k in 0..corners.len() {
        let (a, b) = (corners[k], corners[(k + 1) % corners.len()]);
        for &p in &sorted {
            if sign(orient2d(a, b, p)) == Sign::Zero
                && within(a.x, b.x, p.x)
                && within(a.y, b.y, p.y)
            {
                on.insert((p.x.to_bits(), p.y.to_bits()));
            }
        }
    }
    on
}

/// Fail unless `tri` is a Delaunay triangulation of the distinct points
/// of `points` whose boundary is their exact convex hull.
fn assert_covers_hull(name: &str, points: &[Point2], tri: &Triangulation) {
    let at = |k: u32| points[k as usize];
    let mut edges: HashMap<(u32, u32), u32> = HashMap::new();
    let mut used = vec![false; points.len()];
    for t in tri.triangles().chunks_exact(3) {
        assert_eq!(
            sign(orient2d(at(t[0]), at(t[1]), at(t[2]))),
            Sign::Positive,
            "{name}: triangle {t:?} is not counter-clockwise"
        );
        for i in 0..3 {
            used[t[i] as usize] = true;
            let before = edges.insert((t[i], t[(i + 1) % 3]), t[(i + 2) % 3]);
            assert!(before.is_none(), "{name}: a directed edge occurs twice");
        }
    }
    let distinct: HashSet<(u64, u64)> = points
        .iter()
        .map(|p| (p.x.to_bits(), p.y.to_bits()))
        .collect();
    let corners: HashSet<(u64, u64)> = (0..points.len())
        .filter(|&k| used[k])
        .map(|k| (points[k].x.to_bits(), points[k].y.to_bits()))
        .collect();
    assert_eq!(corners, distinct, "{name}: a point is not a corner");

    let mut boundary = HashSet::new();
    let mut out = HashMap::new();
    for (&(a, b), &apex) in &edges {
        match edges.get(&(b, a)) {
            Some(&other) => assert_ne!(
                sign(incircle(at(a), at(b), at(apex), at(other))),
                Sign::Positive,
                "{name}: edge ({a}, {b}) is not locally Delaunay"
            ),
            None => {
                assert!(out.insert(a, b).is_none(), "{name}: boundary pinches");
                boundary.insert((at(a).x.to_bits(), at(a).y.to_bits()));
                // A hull edge has no point strictly outside it.
                assert!(
                    points
                        .iter()
                        .all(|&p| sign(orient2d(at(a), at(b), p)) != Sign::Negative),
                    "{name}: boundary edge ({a}, {b}) is not on the hull"
                );
            }
        }
    }
    // One closed walk through every boundary vertex.
    let first = *out.keys().next().expect("a boundary");
    let (mut v, mut steps) = (first, 0);
    loop {
        v = out[&v];
        steps += 1;
        if v == first {
            break;
        }
    }
    assert_eq!(steps, out.len(), "{name}: the boundary is not one loop");

    let hull = hull_boundary(points);
    assert_eq!(boundary, hull, "{name}: boundary is not the exact hull");
    let (n, h) = (distinct.len(), hull.len());
    assert_eq!(
        tri.triangle_count(),
        2 * n - 2 - h,
        "{name}: {n} points with {h} on the hull"
    );
}

fn check(name: &str, points: &[Point2]) {
    let tri = triangulate(points, &[]).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_covers_hull(name, points, &tri);
}

#[test]
fn random_sets_from_three_to_sixteen_thousand_points_cover_the_hull() {
    for (k, count) in [3usize, 4, 5, 7, 10, 30, 100, 500, 2000, 4000, 16_000]
        .into_iter()
        .enumerate()
    {
        let points = uniform(count, 0x2990 + k as u64);
        check(&format!("{count} random"), &points);
    }
    // Many small sets: the hull is most of the input.
    for seed in 0..300 {
        let count = 3 + (seed as usize % 12);
        check(&format!("small {seed}"), &uniform(count, 0x29_9000 + seed));
    }
}

#[test]
fn clustered_sets_cover_the_hull() {
    let mut r = Rng(0x2991);
    for round in 0..6 {
        let mut points = Vec::new();
        let clusters = 2 + round;
        for _ in 0..clusters {
            let (cx, cy) = (r.unit() * 100.0, r.unit() * 100.0);
            let spread = 10f64.powi(-round - 1);
            for _ in 0..300 {
                // A sum of uniforms peaks at the centre.
                let dx = (r.unit() + r.unit() + r.unit() - 1.5) * spread;
                let dy = (r.unit() + r.unit() + r.unit() - 1.5) * spread;
                points.push(Point2::new(cx + dx, cy + dy));
            }
        }
        check(&format!("clusters {round}"), &points);
    }
}

#[test]
fn nearly_collinear_hull_runs_cover_the_hull() {
    let mut r = Rng(0x2992);
    // Long runs along a barely curved hull: a flat parabola below, a long
    // turned line above, and random points between. These thin hull
    // triangles are the ones a finite super triangle loses.
    for (round, bend) in [1e-3, 1e-6, 1e-9, 1e-12, 0.0].into_iter().enumerate() {
        let mut points = Vec::new();
        for k in 0..200 {
            let x = k as f64 / 199.0 * 10.0 - 5.0;
            points.push(Point2::new(x, bend * x * x));
        }
        let (c, s) = (0.3f64.cos(), 0.3f64.sin());
        for k in 0..150 {
            let t = k as f64 / 149.0 * 10.0;
            // Not exactly on one line after rounding: a nearly collinear run.
            points.push(Point2::new(-5.0 + c * t, 4.0 + s * t));
        }
        for _ in 0..300 {
            points.push(Point2::new(r.unit() * 8.0 - 4.0, 0.5 + r.unit() * 3.0));
        }
        check(&format!("bend {bend}"), &points);
        // Inserted along the run first, so the build starts on it.
        let mut reversed = points.clone();
        reversed.reverse();
        check(&format!("bend {bend}, reversed {round}"), &reversed);
    }
    // Exactly collinear hull points: a grid, in row order and shuffled.
    for side in [2usize, 3, 10, 40] {
        let mut grid: Vec<Point2> = (0..side * side)
            .map(|k| Point2::new((k % side) as f64, (k / side) as f64))
            .collect();
        check(&format!("grid {side}"), &grid);
        for k in (1..grid.len()).rev() {
            let j = (r.unit() * (k + 1) as f64) as usize;
            grid.swap(k, j);
        }
        check(&format!("shuffled grid {side}"), &grid);
    }
    // The first points collinear, so the first triangle comes late and the
    // collinear ones are inserted on and beyond its hull edges.
    let mut line: Vec<Point2> = (0..50).map(|k| Point2::new(f64::from(k), 0.0)).collect();
    line.push(Point2::new(20.5, 3.0));
    line.extend((0..20).map(|k| Point2::new(-10.0 - f64::from(k), 0.0)));
    line.push(Point2::new(5.0, -1.0));
    check("collinear prefix", &line);
}

#[test]
fn sets_far_from_the_origin_cover_the_hull() {
    for (k, offset) in [1e5, 5e6].into_iter().enumerate() {
        for (j, span) in [1.0, 1e-3, 100.0].into_iter().enumerate() {
            let points: Vec<Point2> = uniform(2000, 0x2993 + 10 * k as u64 + j as u64)
                .into_iter()
                .map(|p| Point2::new(offset + span * p.x, offset - span * p.y))
                .collect();
            check(&format!("offset {offset}, span {span}"), &points);
        }
    }
}

#[test]
fn the_reported_two_thousand_point_case_keeps_its_hull_triangles() {
    // #299: the benchmark's 2,000 points (its splitmix stream and seed)
    // gave 3,977 triangles; the hull holds 19 points, so 3,979 are due.
    let mut state = 0x2680_0000u64 + 2000;
    let mut unit = || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    };
    let points: Vec<Point2> = (0..2000).map(|_| Point2::new(unit(), unit())).collect();
    let tri = triangulate(&points, &[]).expect("triangulates");
    assert_eq!(hull_boundary(&points).len(), 19);
    assert_eq!(tri.triangle_count(), 3979);
    assert_covers_hull("reported case", &points, &tri);
}
