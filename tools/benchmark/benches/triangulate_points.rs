//! Wall-clock cost of the planar Delaunay triangulation against point
//! count (#268).
//!
//! The incremental build used to recompute every twin pointer after each
//! insertion and after every Lawson flip, so n points cost O(n^2) or worse.
//! It now maintains adjacency locally and starts each point-location walk
//! from a grid hint, and the cost from 1k to 16k uniform random points
//! grows close to n log n: n log n is 22x across the range, the build
//! measured 29x (0.61 ms to 17.5 ms), and the old one over 300x (0.21 s
//! to 72 s).
//!
//! The points are seeded and identical on every machine. Each case is
//! validated before it is timed: every point is a corner of some triangle,
//! every triangle is counter-clockwise, no directed edge occurs twice, the
//! boundary is one closed walk through exactly the points on the exact
//! convex hull, the triangle count is `2n - 2 - h` for those `h` points
//! (Euler's formula for a triangulated disc, so no triangle is missing
//! inside or along the hull, or doubled), and every interior edge is
//! locally Delaunay, decided by the certified predicates.
//!
//! Until #299 the boundary could not be held to the hull: the build's
//! finite super triangle lost a few thin hull triangles (2 of 3979 at 2000
//! points). The vertex at infinity that replaced it keeps them.

use std::collections::HashMap;
use std::hint::black_box;

use axiolid_benchmark::Rng;
use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{incircle, orient2d};
use axiolid_triangulate::{triangulate, Triangulation};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

/// The point counts measured.
const COUNTS: [usize; 5] = [1_000, 2_000, 4_000, 8_000, 16_000];

/// `count` uniform random points in the unit square, from a fixed seed.
fn points(count: usize) -> Vec<Point2> {
    let mut rng = Rng::new(0x2680_0000 + count as u64);
    (0..count)
        .map(|_| Point2::new(rng.next_unit(), rng.next_unit()))
        .collect()
}

/// The proven sign, or `Zero` when the predicate did not decide.
fn sign(certified: Certified) -> Sign {
    match certified {
        Certified::Certain { sign, .. } => sign,
        _ => Sign::Zero,
    }
}

/// Indices of the points on the boundary of the exact convex hull,
/// collinear ones included.
fn hull_points(points: &[Point2]) -> Vec<bool> {
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    sorted.dedup();
    let chain = |iter: &mut dyn Iterator<Item = Point2>| {
        let mut hull: Vec<Point2> = Vec::new();
        for p in iter {
            while hull.len() >= 2
                && sign(orient2d(hull[hull.len() - 2], hull[hull.len() - 1], p)) != Sign::Positive
            {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
        hull
    };
    let mut corners = chain(&mut sorted.iter().copied());
    corners.extend(chain(&mut sorted.iter().rev().copied()));
    let within = |a: f64, b: f64, x: f64| a.min(b) <= x && x <= a.max(b);
    points
        .iter()
        .map(|&p| {
            (0..corners.len()).any(|k| {
                let (a, b) = (corners[k], corners[(k + 1) % corners.len()]);
                sign(orient2d(a, b, p)) == Sign::Zero
                    && within(a.x, b.x, p.x)
                    && within(a.y, b.y, p.y)
            })
        })
        .collect()
}

/// Fail unless `tri` is a complete Delaunay triangulation of `points`.
fn validate(points: &[Point2], tri: &Triangulation) {
    let n = points.len();
    let corner = |k: u32| points[k as usize];
    let mut used = vec![false; n];
    let mut edges: HashMap<(u32, u32), u32> = HashMap::new();
    for t in tri.triangles().as_chunks::<3>().0 {
        assert_eq!(
            sign(orient2d(corner(t[0]), corner(t[1]), corner(t[2]))),
            Sign::Positive,
            "{n} points: triangle {t:?} is not counter-clockwise"
        );
        for i in 0..3 {
            used[t[i] as usize] = true;
            let previous = edges.insert((t[i], t[(i + 1) % 3]), t[(i + 2) % 3]);
            assert!(previous.is_none(), "{n} points: edge used twice");
        }
    }
    assert!(used.iter().all(|&u| u), "{n} points: a point was dropped");
    let mut leaving = vec![0u32; n];
    let mut entering = vec![0u32; n];
    let mut next = HashMap::new();
    for (&(a, b), &apex) in &edges {
        match edges.get(&(b, a)) {
            Some(&other) => assert_ne!(
                sign(incircle(corner(a), corner(b), corner(apex), corner(other))),
                Sign::Positive,
                "{n} points: edge ({a}, {b}) is not locally Delaunay"
            ),
            None => {
                next.insert(a, b);
                leaving[a as usize] += 1;
                entering[b as usize] += 1;
            }
        }
    }
    assert!(
        leaving
            .iter()
            .zip(&entering)
            .all(|(l, e)| l == e && *l <= 1),
        "{n} points: the boundary is not one simple closed walk per vertex"
    );
    let first = *next.keys().next().expect("a boundary");
    let (mut v, mut steps) = (first, 1);
    while next[&v] != first {
        v = next[&v];
        steps += 1;
    }
    assert_eq!(
        steps,
        next.len(),
        "{n} points: the boundary is not one loop"
    );
    let hull = hull_points(points);
    assert!(
        (0..n).all(|k| hull[k] == (leaving[k] == 1)),
        "{n} points: the boundary is not the exact convex hull"
    );
    let h = hull.iter().filter(|&&on| on).count();
    assert_eq!(
        tri.triangle_count(),
        2 * n - 2 - h,
        "{n} points: triangle count does not match the hull"
    );
}

fn random_points(c: &mut Criterion) {
    let mut group = c.benchmark_group("triangulate/random_points");
    group.sample_size(10);
    for count in COUNTS {
        let input = points(count);
        // Validated inside the routine, once, so a filtered run does not
        // pay for the sizes it skips.
        let mut validated = false;
        group.bench_with_input(
            BenchmarkId::from_parameter(count),
            &input,
            |bench, input| {
                if !validated {
                    let tri = triangulate(input, &[]).expect("random points triangulate");
                    validate(input, &tri);
                    validated = true;
                }
                bench.iter(|| black_box(triangulate(black_box(input), &[]).expect("triangulates")))
            },
        );
    }
    group.finish();
}

criterion_group!(benches, random_points);
criterion_main!(benches);
