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
//! boundary is closed and has at least the `h` hull vertices, the triangle
//! count is `2n - 2 - b` for the `b` boundary vertices (Euler's formula for
//! a triangulated disc, so no triangle is missing inside or doubled), and
//! every interior edge is locally Delaunay, decided by the certified
//! predicates.
//!
//! The boundary is not required to be the hull: the build's super triangle
//! is finite, and a few thin triangles along the hull can be lost with it
//! (2 of 3979 at 2000 points). That is the build's behaviour before and
//! after #268 alike.

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

/// Number of vertices on the convex hull, strictly convex corners only.
fn hull_vertices(points: &[Point2]) -> usize {
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
        hull.len() - 1
    };
    chain(&mut sorted.iter().copied()) + chain(&mut sorted.iter().rev().copied())
}

/// Fail unless `tri` is a complete Delaunay triangulation of `points`.
fn validate(points: &[Point2], tri: &Triangulation) {
    let n = points.len();
    let corner = |k: u32| points[k as usize];
    let mut used = vec![false; n];
    let mut edges: HashMap<(u32, u32), u32> = HashMap::new();
    for t in tri.triangles().chunks_exact(3) {
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
    let mut boundary = 0;
    for (&(a, b), &apex) in &edges {
        match edges.get(&(b, a)) {
            Some(&other) => assert_ne!(
                sign(incircle(corner(a), corner(b), corner(apex), corner(other))),
                Sign::Positive,
                "{n} points: edge ({a}, {b}) is not locally Delaunay"
            ),
            None => {
                boundary += 1;
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
    assert!(
        boundary >= hull_vertices(points),
        "{n} points: fewer boundary vertices than hull vertices"
    );
    assert_eq!(
        tri.triangle_count(),
        2 * n - 2 - boundary,
        "{n} points: triangle count does not match the boundary"
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
