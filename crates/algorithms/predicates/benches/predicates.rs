//! Predicate throughput and escalation rate across degeneracy tiers.
//!
//! Run with: `cargo bench -p axiolid-predicates`
//!
//! No criterion dependency: this reports wall-clock throughput and the
//! escalation rate side by side, which is the pairing that makes the numbers
//! interpretable. A throughput drop is only meaningful next to the escalation
//! rate that caused it.

use std::hint::black_box;
use std::time::Instant;

use axiolid_core::Point2;
use axiolid_predicates::scene::{orient2_scene, orient3_scene, DegeneracyRate};
use axiolid_predicates::{
    incircle, incircle_filter, orient2d, orient2d_filter, orient3d, orient3d_filter, StaticFilter,
};

const SAMPLES: usize = 200_000;

fn main() {
    println!("predicate throughput and escalation by degeneracy rate");
    println!("samples per tier: {SAMPLES}");
    println!();
    bench_orient2d();
    println!();
    bench_orient2d_collinear_grid();
    println!();
    bench_orient3d();
    println!();
    bench_incircle();
    println!();
    bench_static_filter();
}

/// Report one row: throughput plus the escalation rate that explains it.
fn row(label: &str, elapsed_ns: f64, escalated: usize) {
    let per_call = elapsed_ns / SAMPLES as f64;
    let millions_per_second = 1_000.0 / per_call;
    println!(
        "  {label:>7}  {per_call:>8.2} ns/call  {millions_per_second:>8.2} M/s  \
         escalated {:>7.4}%",
        escalated as f64 * 100.0 / SAMPLES as f64
    );
}

fn bench_orient2d() {
    println!("orient2d (filtered cascade)");
    for rate in DegeneracyRate::ALL {
        let scene = orient2_scene(SAMPLES, rate, 0x51ED_2701);
        let escalated = scene
            .iter()
            .filter(|c| !orient2d_filter(c[0], c[1], c[2]).is_certain())
            .count();

        // Warm the cache so the first tier is not penalised for cold data.
        for c in &scene {
            black_box(orient2d(c[0], c[1], c[2]));
        }
        let start = Instant::now();
        for c in &scene {
            black_box(orient2d(black_box(c[0]), black_box(c[1]), black_box(c[2])));
        }
        row(rate.label(), start.elapsed().as_nanos() as f64, escalated);
    }
}

/// Exactly collinear triples with full-mantissa coordinates on the line
/// `y = 2x`. Every call escalates past the filter, so this prices the exact
/// stage that the common degenerate case lands in. Two scenes: a millimetre
/// grid in [100, 190] m, whose coordinate differences are all exact, and one
/// in [0, 200) m, which crosses binades so that many differences round.
fn bench_orient2d_collinear_grid() {
    println!("orient2d (exactly collinear, full-mantissa millimetre grid)");
    for (label, offset, count) in [("exact", 100.0, 90_000_u64), ("binades", 0.0, 200_000)] {
        let mut state = 0x51ED_2702_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let scene: Vec<[axiolid_core::Point2; 3]> = (0..SAMPLES)
            .map(|_| {
                let mut vertex = || {
                    let x = offset + (next() % count) as f64 * 0.001;
                    axiolid_core::Point2::new(x, 2.0 * x)
                };
                [vertex(), vertex(), vertex()]
            })
            .collect();
        let escalated = scene
            .iter()
            .filter(|c| !orient2d_filter(c[0], c[1], c[2]).is_certain())
            .count();
        for c in &scene {
            black_box(orient2d(c[0], c[1], c[2]));
        }
        let start = Instant::now();
        for c in &scene {
            black_box(orient2d(black_box(c[0]), black_box(c[1]), black_box(c[2])));
        }
        row(label, start.elapsed().as_nanos() as f64, escalated);
    }
}

fn bench_orient3d() {
    println!("orient3d (filtered cascade)");
    for rate in DegeneracyRate::ALL {
        let scene = orient3_scene(SAMPLES, rate, 0x51ED_2701);
        let escalated = scene
            .iter()
            .filter(|c| !orient3d_filter(c[0], c[1], c[2], c[3]).is_certain())
            .count();

        for c in &scene {
            black_box(orient3d(c[0], c[1], c[2], c[3]));
        }
        let start = Instant::now();
        for c in &scene {
            black_box(orient3d(
                black_box(c[0]),
                black_box(c[1]),
                black_box(c[2]),
                black_box(c[3]),
            ));
        }
        row(rate.label(), start.elapsed().as_nanos() as f64, escalated);
    }
}

/// `incircle` on random points (the filter settles every call) and on exactly
/// cocircular lattice points, where every call escalates past the filter. The
/// lattice points of the circle of radius 5525 are scaled by a power of two,
/// so they stay exactly cocircular: 2^-10 keeps every difference exact and in
/// range (the expansion stage), 2^266 (about 1e80) puts the differences out
/// of the filter's range (the dyadic stage).
fn bench_incircle() {
    println!("incircle");
    let mut state = 0x51ED_2703_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let random: Vec<[Point2; 4]> = (0..SAMPLES)
        .map(|_| {
            let mut vertex = || {
                Point2::new(
                    (next() % 2_000_001) as f64 - 1_000_000.0,
                    (next() % 2_000_001) as f64 - 1_000_000.0,
                )
            };
            [vertex(), vertex(), vertex(), vertex()]
        })
        .collect();
    time_incircle("random", &random);

    const RADIUS: i64 = 5525;
    let lattice: Vec<(f64, f64)> = (-RADIUS..=RADIUS)
        .flat_map(|x| (-RADIUS..=RADIUS).map(move |y| (x, y)))
        .filter(|(x, y)| x * x + y * y == RADIUS * RADIUS)
        .map(|(x, y)| (x as f64, y as f64))
        .collect();
    for (label, scale) in [("grid", 2f64.powi(-10)), ("huge", 2f64.powi(266))] {
        let scene: Vec<[Point2; 4]> = (0..SAMPLES)
            .map(|_| {
                let mut vertex = || {
                    let (x, y) = lattice[(next() % lattice.len() as u64) as usize];
                    Point2::new(x * scale, y * scale)
                };
                [vertex(), vertex(), vertex(), vertex()]
            })
            .collect();
        time_incircle(label, &scene);
    }
}

fn time_incircle(label: &str, scene: &[[Point2; 4]]) {
    let escalated = scene
        .iter()
        .filter(|c| !incircle_filter(c[0], c[1], c[2], c[3]).is_certain())
        .count();
    for c in scene {
        black_box(incircle(c[0], c[1], c[2], c[3]));
    }
    let start = Instant::now();
    for c in scene {
        black_box(incircle(
            black_box(c[0]),
            black_box(c[1]),
            black_box(c[2]),
            black_box(c[3]),
        ));
    }
    row(label, start.elapsed().as_nanos() as f64, escalated);
}

/// The static filter's whole justification is skipping the per-call permanent.
/// If it is not measurably faster on clean data, it should not exist.
fn bench_static_filter() {
    println!("orient2d: dynamic filter vs static filter (0% degenerate)");
    let scene = orient2_scene(SAMPLES, DegeneracyRate::None, 0x51ED_2701);
    let filter = StaticFilter::new(1_000.0).expect("valid range");

    let start = Instant::now();
    for c in &scene {
        black_box(orient2d_filter(
            black_box(c[0]),
            black_box(c[1]),
            black_box(c[2]),
        ));
    }
    let dynamic = start.elapsed().as_nanos() as f64 / SAMPLES as f64;

    let mut declined = 0usize;
    let start = Instant::now();
    for c in &scene {
        if black_box(filter.orient2d(black_box(c[0]), black_box(c[1]), black_box(c[2]))).is_none() {
            declined += 1;
        }
    }
    let stat = start.elapsed().as_nanos() as f64 / SAMPLES as f64;

    println!("  dynamic  {dynamic:>8.2} ns/call");
    println!(
        "   static  {stat:>8.2} ns/call  declined {:>7.4}%  ({:.2}x)",
        declined as f64 * 100.0 / SAMPLES as f64,
        dynamic / stat
    );
}
