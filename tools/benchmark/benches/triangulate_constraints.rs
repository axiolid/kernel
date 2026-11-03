//! Wall-clock cost of constraint recovery in the planar triangulation
//! against segment count (#355).
//!
//! Recovery used to scan the whole triangle array for every missing
//! constraint (to find its crossing edges, to look an edge up, and to find
//! the vertices lying on it), so points x segments dominated: a jittered
//! grid of 49,729 points with 12,432 zigzag segments took 96 s in
//! rheonda/cad's measurement and 168 s (best of three) in this benchmark on
//! the same code. It now
//! reaches the edges of a constraint through a vertex-to-triangle index and
//! walks from one end to the other with the exact orientation predicate, so
//! the cost follows the edges a constraint crosses.
//!
//! `cargo bench -p axiolid-benchmark --bench triangulate_constraints`
//!
//! The input is a square grid with every point jittered by up to 0.3 cells
//! from a fixed seed, so it is identical on every machine. Breaklines are
//! zigzags: line `k` joins the grid points `(i, 2k + i mod 2)` for
//! consecutive `i`, which a Delaunay triangulation of a jittered grid
//! mostly lacks. Each case is validated before it is timed: every
//! constraint is an edge of the result, and the result covers the convex
//! hull (`2n - 2 - h` triangles for `h` hull points).
//!
//! Each case is timed `RUNS` times, whole, with `Instant`; the best and the
//! median are printed. The cases are too long on the old code for a
//! sampling harness.

use std::collections::HashSet;
use std::hint::black_box;
use std::time::Instant;

use axiolid_benchmark::Rng;
use axiolid_core::Point2;
use axiolid_triangulate::{triangulate, Constraint};

/// Timed runs per case.
const RUNS: usize = 3;

/// `(grid side, zigzag lines)`: 10,000 points with 4,950 segments, and
/// 49,729 points with 12,432 segments.
const CASES: [(usize, usize); 2] = [(100, 50), (223, 56)];

/// A jittered `side` x `side` grid and the zigzag breaklines over it.
fn input(side: usize, lines: usize) -> (Vec<Point2>, Vec<Constraint>) {
    let mut rng = Rng::new(0x3550_0000 + side as u64);
    let mut points = Vec::with_capacity(side * side);
    for row in 0..side {
        for col in 0..side {
            let x = col as f64 + (rng.next_unit() - 0.5) * 0.6;
            let y = row as f64 + (rng.next_unit() - 0.5) * 0.6;
            points.push(Point2::new(x, y));
        }
    }
    let at = |col: usize, row: usize| (row * side + col) as u32;
    let mut constraints = Vec::new();
    for line in 0..lines {
        let row = |col: usize| 2 * line + col % 2;
        for col in 0..side - 1 {
            constraints.push(Constraint::new(
                at(col, row(col)),
                at(col + 1, row(col + 1)),
            ));
        }
    }
    (points, constraints)
}

/// Fail unless every constraint is an edge of the triangulation and the
/// triangle count matches a triangulated disc over the points.
fn validate(points: &[Point2], constraints: &[Constraint]) {
    let tri = triangulate(points, constraints).expect("the breaklines are recoverable");
    let mut edges = HashSet::new();
    for t in tri.triangles().as_chunks::<3>().0 {
        for i in 0..3 {
            edges.insert((t[i], t[(i + 1) % 3]));
        }
    }
    let hull = edges
        .iter()
        .filter(|&&(a, b)| !edges.contains(&(b, a)))
        .count();
    assert_eq!(
        tri.triangle_count(),
        2 * points.len() - 2 - hull,
        "{} points: triangle count does not match the hull",
        points.len()
    );
    for c in constraints {
        assert!(
            edges.contains(&(c.a, c.b)) || edges.contains(&(c.b, c.a)),
            "{} points: constraint {c:?} is not an edge",
            points.len()
        );
    }
}

fn main() {
    for (side, lines) in CASES {
        let (points, constraints) = input(side, lines);
        validate(&points, &constraints);
        let mut times: Vec<f64> = (0..RUNS)
            .map(|_| {
                let start = Instant::now();
                black_box(triangulate(black_box(&points), black_box(&constraints)))
                    .expect("triangulates");
                start.elapsed().as_secs_f64()
            })
            .collect();
        times.sort_by(f64::total_cmp);
        println!(
            "triangulate/constraints/{side}x{side}: {} points, {} segments: best {:.3} s, median {:.3} s",
            points.len(),
            constraints.len(),
            times[0],
            times[RUNS / 2]
        );
    }
}
