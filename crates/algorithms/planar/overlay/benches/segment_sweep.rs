//! Cost growth of the segment sweep (#146) with segment count.
//!
//! `cargo bench -p axiolid-overlay --bench segment_sweep`
//!
//! Short random segments on an integer grid (so shared ends, overlaps and
//! verticals occur), and random long segments in a unit square (many
//! crossings, general position). Beside each sweep, the time a pairwise
//! scan spends only testing every pair's bounding boxes, which is a lower
//! bound on any `O(n^2)` method.

use std::hint::black_box;
use std::time::Instant;

use axiolid_core::Point2;
use axiolid_overlay::segment_intersections;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn short(n: usize, rng: &mut Rng) -> Vec<[Point2; 2]> {
    let side = (n as f64).sqrt().ceil() * 4.0;
    (0..n)
        .map(|_| {
            let x = (rng.unit() * side).floor();
            let y = (rng.unit() * side).floor();
            let dx = (rng.unit() * 9.0).floor() - 4.0;
            let dy = (rng.unit() * 9.0).floor() - 4.0;
            [Point2::new(x, y), Point2::new(x + dx, y + dy)]
        })
        .collect()
}

fn long(n: usize, rng: &mut Rng) -> Vec<[Point2; 2]> {
    (0..n)
        .map(|_| {
            let a = Point2::new(rng.unit(), rng.unit());
            let length = 4.0 / (n as f64).sqrt();
            let b = Point2::new(
                a.x + (rng.unit() - 0.5) * length,
                a.y + (rng.unit() - 0.5) * length,
            );
            [a, b]
        })
        .collect()
}

fn pairwise_boxes(segments: &[[Point2; 2]]) -> usize {
    let mut hits = 0;
    for (i, [a, b]) in segments.iter().enumerate() {
        for [c, d] in &segments[i + 1..] {
            if a.x.min(b.x) <= c.x.max(d.x)
                && c.x.min(d.x) <= a.x.max(b.x)
                && a.y.min(b.y) <= c.y.max(d.y)
                && c.y.min(d.y) <= a.y.max(b.y)
            {
                hits += 1;
            }
        }
    }
    hits
}

fn main() {
    for (name, make) in [
        (
            "short grid",
            short as fn(usize, &mut Rng) -> Vec<[Point2; 2]>,
        ),
        ("long float", long),
    ] {
        for n in [1_000usize, 4_000, 16_000, 64_000] {
            let segments = make(n, &mut Rng(0x9e37_79b9_7f4a_7c15));
            let start = Instant::now();
            let result = black_box(segment_intersections(black_box(&segments)).unwrap());
            let sweep = start.elapsed().as_secs_f64() * 1e3;
            let incidences: usize = result.points.iter().map(|p| p.incidences.len()).sum();
            let pairwise = if n <= 16_000 {
                let start = Instant::now();
                black_box(pairwise_boxes(black_box(&segments)));
                format!("{:9.1} ms", start.elapsed().as_secs_f64() * 1e3)
            } else {
                "        -".to_owned()
            };
            println!(
                "{name:<10} n={n:>6} points={:>7} incidences={incidences:>7} overlaps={:>5} \
                 sweep {sweep:9.1} ms  pairwise boxes {pairwise}",
                result.points.len(),
                result.overlaps.len(),
            );
        }
    }
}
