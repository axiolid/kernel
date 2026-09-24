//! Filter throughput, escalation rate and exact-tier cost, side by side.
//!
//! Run with: `cargo bench -p axiolid-exact`
//!
//! A throughput number is only interpretable next to the escalation rate
//! that explains it (the same pairing as the predicates bench). Crossing
//! workloads: random coordinates (almost never escalate) and crossings
//! placed exactly on the query line (always escalate: the answer is an
//! exact zero, and no interval can prove zero).

use std::hint::black_box;
use std::time::Instant;

use axiolid_core::Point2;
use axiolid_exact::{
    compare_along, crossing_orientation, crossing_orientation_filter, line_circle_hits, Circle,
    HitCount, Line,
};

const SAMPLES: usize = 100_000;

struct Rng(u64);

impl Rng {
    fn f(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2000.0 - 1000.0
    }

    fn point(&mut self) -> Point2 {
        Point2::new(self.f(), self.f())
    }
}

type Case = (Line, Line, Point2, Point2);

fn random_case(rng: &mut Rng) -> Case {
    loop {
        let first = Line::new(rng.point(), rng.point());
        let second = Line::new(rng.point(), rng.point());
        if let (Ok(first), Ok(second)) = (first, second) {
            return (first, second, rng.point(), rng.point());
        }
    }
}

/// Lines crossing at an integer point, and a query line through it.
fn degenerate_case(rng: &mut Rng) -> Case {
    let (cx, cy) = (rng.f().round(), rng.f().round());
    let at = |dx: f64, dy: f64| Point2::new(cx + dx, cy + dy);
    let first = Line::new(at(-3.0, -1.0), at(3.0, 1.0)).expect("distinct points");
    let second = Line::new(at(-1.0, 2.0), at(1.0, -2.0)).expect("distinct points");
    (first, second, at(-5.0, -5.0), at(7.0, 7.0))
}

fn row(label: &str, nanos: f64, escalated: Option<usize>) {
    let per = nanos / SAMPLES as f64;
    let rate = escalated.map_or_else(
        || "        -".to_owned(),
        |n| format!("{:>8.2}%", n as f64 * 100.0 / SAMPLES as f64),
    );
    println!("  {label:<36} {per:>9.1} ns/call   escalated {rate}");
}

fn bench_crossings(label: &str, make: fn(&mut Rng) -> Case) {
    let mut rng = Rng(0x5EED_1234_ABCD_0001);
    let cases: Vec<Case> = (0..SAMPLES).map(|_| make(&mut rng)).collect();
    let escalated = cases
        .iter()
        .filter(|(f, s, a, b)| !crossing_orientation_filter(*f, *s, *a, *b).is_certain())
        .count();
    for (f, s, a, b) in &cases {
        black_box(crossing_orientation(*f, *s, *a, *b).ok());
    }
    let start = Instant::now();
    for (f, s, a, b) in &cases {
        black_box(crossing_orientation(black_box(*f), *s, *a, *b).ok());
    }
    row(label, start.elapsed().as_nanos() as f64, Some(escalated));
}

fn bench_line_circle() {
    let mut rng = Rng(0x5EED_1234_ABCD_0002);
    let cases: Vec<(Line, Circle)> = (0..SAMPLES)
        .map(|_| loop {
            let line = Line::new(rng.point(), rng.point());
            let circle = Circle::new(rng.point(), rng.f().abs());
            if let (Ok(line), Ok(circle)) = (line, circle) {
                break (line, circle);
            }
        })
        .collect();
    let mut secants = 0;
    let start = Instant::now();
    for (line, circle) in &cases {
        if let Ok(HitCount::Secant(first, second)) = line_circle_hits(*line, *circle) {
            secants += 1;
            black_box(compare_along(first, second).ok());
        }
    }
    row(
        "classify + order both hits",
        start.elapsed().as_nanos() as f64,
        None,
    );
    println!("  ({secants} of {SAMPLES} lines were secants)");
}

fn main() {
    println!("axiolid-exact: filtered exact signs, {SAMPLES} cases per row");
    println!();
    println!("crossing_orientation (two lines crossed, side of a third)");
    bench_crossings("random coordinates", random_case);
    bench_crossings("crossing exactly on the query line", degenerate_case);
    println!();
    println!("line_circle_hits (random line, random circle)");
    bench_line_circle();
}
