//! Per-call cost of the arc boolean on typical BIM-sized scenes.
//!
//! `cargo bench -p axiolid-overlay --bench arc_overlay`

use std::hint::black_box;
use std::time::Instant;

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{arc_overlay, ArcRing, ArcVertex, OverlayOperation};

fn wall_with_round_end() -> ArcRing {
    ArcRing::new(vec![
        ArcVertex::bulged(Point2::new(0.0, 0.0), 0.0),
        ArcVertex::bulged(Point2::new(5.0, 0.0), 1.0),
        ArcVertex::bulged(Point2::new(5.0, 0.3), 0.0),
        ArcVertex::bulged(Point2::new(0.0, 0.3), 0.0),
    ])
}

fn scenes() -> Vec<(&'static str, ArcRing, ArcRing)> {
    vec![
        (
            "disc crossing a rectangle",
            ArcRing::from_points(&[
                Point2::new(0.0, 0.0),
                Point2::new(4.0, 0.0),
                Point2::new(4.0, 3.0),
                Point2::new(0.0, 3.0),
            ]),
            ArcRing::circle(Point2::new(4.0, 1.5), 1.0),
        ),
        (
            "two overlapping discs",
            ArcRing::circle(Point2::new(0.0, 0.0), 1.0),
            ArcRing::circle(Point2::new(1.2, 0.3), 0.8),
        ),
        (
            "round opening in a rounded wall",
            wall_with_round_end(),
            ArcRing::circle(Point2::new(2.5, 0.15), 0.1),
        ),
        (
            "rounded wall end against a disc",
            wall_with_round_end(),
            ArcRing::circle(Point2::new(5.2, 0.15), 0.25),
        ),
    ]
}

/// A wavy outline with `n` arc edges on a circle of radius 2, against a
/// disc crossing it: cost growth with edge count.
fn scaling() {
    for n in [8usize, 32, 128, 512] {
        let step = std::f64::consts::TAU / n as f64;
        let ring = ArcRing::new(
            (0..n)
                .map(|i| {
                    let t = i as f64 * step;
                    let bulge = if i % 2 == 0 { 0.2 } else { -0.2 };
                    ArcVertex::bulged(Point2::new(2.0 * t.cos(), 2.0 * t.sin()), bulge)
                })
                .collect(),
        );
        let disc = ArcRing::circle(Point2::new(1.9, 0.1), 0.7);
        let iterations = (20_000 / n).max(4) as u32;
        let start = Instant::now();
        for _ in 0..iterations {
            black_box(
                arc_overlay(&ring, &disc, OverlayOperation::Union, Tolerance::METRE).unwrap(),
            );
        }
        let per = start.elapsed().as_secs_f64() / f64::from(iterations) * 1e6;
        println!("{n:>4} arc edges vs disc                 Union: {per:8.1} us/call");
    }
}

/// A wavy ring of `n` alternating-bulge arcs around `centre`.
fn wavy(n: usize, centre: Point2, radius: f64, bulge: f64) -> ArcRing {
    let step = std::f64::consts::TAU / n as f64;
    ArcRing::new(
        (0..n)
            .map(|i| {
                let t = i as f64 * step;
                let b = if i % 2 == 0 { bulge } else { -bulge };
                ArcVertex::bulged(
                    Point2::new(centre.x + radius * t.cos(), centre.y + radius * t.sin()),
                    b,
                )
            })
            .collect(),
    )
}

/// Two `n`-edge wavy rings overlapping by half: every edge of one is near
/// many edges of the other, so this is the all-pairs worst case.
fn scaling_both() {
    for n in [16usize, 64, 256, 1024] {
        let a = wavy(n, Point2::new(0.0, 0.0), 2.0, 0.2);
        let b = wavy(n, Point2::new(1.3, 0.4), 2.0, 0.15);
        let iterations = (4_000 / n).max(2) as u32;
        let start = Instant::now();
        for _ in 0..iterations {
            black_box(arc_overlay(&a, &b, OverlayOperation::Union, Tolerance::METRE).unwrap());
        }
        let per = start.elapsed().as_secs_f64() / f64::from(iterations) * 1e6;
        println!("{n:>4} x {n:<4} wavy rings                Union: {per:8.1} us/call");
    }
}

/// A plate with many small round openings punched by one Difference each
/// is the BIM-typical large case; here: a polyline plate edge count grows
/// while the tool stays one small disc far from most edges.
fn scaling_local() {
    for n in [64usize, 256, 1024, 4096] {
        // A long thin comb: n straight edges, only a few near the disc.
        let mut pts = Vec::with_capacity(n);
        let teeth = n / 4;
        for i in 0..teeth {
            let x = i as f64;
            pts.push(Point2::new(x, 0.0));
            pts.push(Point2::new(x + 0.4, 0.0));
            pts.push(Point2::new(x + 0.4, 1.0));
            pts.push(Point2::new(x + 0.6, 1.0));
        }
        pts.push(Point2::new(teeth as f64, 0.0));
        pts.push(Point2::new(teeth as f64, -1.0));
        pts.push(Point2::new(0.0, -1.0));
        let comb = ArcRing::from_points(&pts);
        let disc = ArcRing::circle(Point2::new(2.5, -0.5), 0.3);
        let iterations = (20_000 / n).max(4) as u32;
        let start = Instant::now();
        for _ in 0..iterations {
            black_box(
                arc_overlay(&comb, &disc, OverlayOperation::Difference, Tolerance::METRE).unwrap(),
            );
        }
        let per = start.elapsed().as_secs_f64() / f64::from(iterations) * 1e6;
        println!("{n:>4} comb edges vs a local disc     Diff:  {per:8.1} us/call");
    }
}

fn main() {
    if std::env::var("SCALE").is_ok() {
        scaling();
        scaling_both();
        scaling_local();
        return;
    }
    scaling();
    let iterations = 2_000;
    let only = std::env::var("ONLY").ok();
    for (name, a, b) in scenes() {
        if only.as_deref().is_some_and(|o| !name.contains(o)) {
            continue;
        }
        for op in [OverlayOperation::Union, OverlayOperation::Difference] {
            let start = Instant::now();
            for _ in 0..iterations {
                black_box(arc_overlay(black_box(&a), black_box(&b), op, Tolerance::METRE).unwrap());
            }
            let per = start.elapsed().as_secs_f64() / f64::from(iterations) * 1e6;
            println!("{name:<34} {op:?}: {per:8.1} us/call");
        }
    }
}
