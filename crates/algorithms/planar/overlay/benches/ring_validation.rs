//! Wall-clock time of validating a large ring, and of booleans on regions
//! built from one, against corner count (#341).
//!
//! The ring is a floor of 40 m by 30 m whose two long walls carry teeth
//! 0.5 mm deep, as overlay output along tessellated bodies does: valid,
//! with most edges near collinear neighbours. Cases:
//!
//! - `new`: `Region::new` of the floor, which validates its ring.
//! - `overlay`: `overlay` of the floor's polygon with a 1 m square inside
//!   it (`Difference`), which validates both operands every call.
//! - `region`: `Region::difference` of the floor region and the square's,
//!   five times in a row on each result, as a disc search subtracts one
//!   obstacle after another. Validated regions are not validated again.
//!
//! Before #341 a ring's validation compared every pair of its edges: at
//! 50,000 corners about 4 s per validation, at 100,000 about 16 s.
//!
//! Each figure is the median of `RUNS` timed runs
//! (`RING_VALIDATION_RUNS` overrides the count). Set
//! `RING_VALIDATION_CORNERS` (for example `50000`) to run other sizes.
//!
//! ```bash
//! cargo bench -p axiolid-overlay --bench ring_validation
//! ```

use std::hint::black_box;
use std::time::Instant;

use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::{overlay, FillRule, OverlayInput, OverlayOperation, Polygon, Region, Ring};

/// Corner counts measured.
const CORNERS: [usize; 3] = [10_000, 50_000, 100_000];
/// Timed runs per case; the median is reported.
const RUNS: usize = 3;
/// The floor's size.
const FLOOR: (f64, f64) = (40.0, 30.0);

/// The floor with about `corners` corners: teeth along both long walls.
fn floor(corners: usize) -> Vec<Point2> {
    let teeth = (corners - 4) / 4;
    let (w, h) = FLOOR;
    let step = w / (teeth as f64 + 1.0);
    let depth = 5.0e-4;
    let mut points = vec![Point2::new(0.0, 0.0)];
    for k in 0..teeth {
        let x = step * (k as f64 + 1.0);
        points.push(Point2::new(x - 0.25 * step, 0.0));
        points.push(Point2::new(x, depth));
    }
    points.push(Point2::new(w, 0.0));
    points.push(Point2::new(w, h));
    for k in 0..teeth {
        let x = w - step * (k as f64 + 1.0);
        points.push(Point2::new(x + 0.25 * step, h));
        points.push(Point2::new(x, h - depth));
    }
    points.push(Point2::new(0.0, h));
    points
}

fn polygon(points: Vec<Point2>) -> Polygon {
    Polygon {
        outer: Ring { points },
        holes: Vec::new(),
    }
}

fn square(x: f64, y: f64, side: f64) -> Polygon {
    polygon(vec![
        Point2::new(x, y),
        Point2::new(x + side, y),
        Point2::new(x + side, y + side),
        Point2::new(x, y + side),
    ])
}

fn input(polygons: Vec<Polygon>) -> OverlayInput {
    OverlayInput {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Point2::X,
            y: Point2::Y,
        },
        polygons,
    }
}

fn median(runs: usize, mut run: impl FnMut()) -> f64 {
    let mut times: Vec<f64> = (0..runs)
        .map(|_| {
            let start = Instant::now();
            run();
            start.elapsed().as_secs_f64()
        })
        .collect();
    times.sort_by(f64::total_cmp);
    times[(runs - 1) / 2]
}

fn main() {
    let runs: usize = std::env::var("RING_VALIDATION_RUNS")
        .ok()
        .and_then(|r| r.parse().ok())
        .unwrap_or(RUNS);
    let corners: Vec<usize> = std::env::var("RING_VALIDATION_CORNERS")
        .ok()
        .map(|list| {
            list.split(',')
                .filter(|n| !n.trim().is_empty())
                .map(|n| n.trim().parse().expect("a corner count"))
                .collect()
        })
        .unwrap_or_else(|| CORNERS.to_vec());
    let t = Tolerance::new(1.0e-9, 1.0e-9).expect("a tolerance");
    println!(
        "{:<8} {:>8} {:>10} {:>10} {:>12}",
        "case", "corners", "median_s", "per_call_s", "result"
    );
    for n in corners {
        let points = floor(n);
        let count = points.len();
        let obstacles: Vec<Polygon> = (0..5)
            .map(|k| square(5.0 + 6.0 * k as f64, 10.0, 1.0))
            .collect();

        let new = median(runs, || {
            black_box(Region::new(vec![polygon(points.clone())], t).expect("a valid floor"));
        });
        println!(
            "{:<8} {:>8} {:>10.3} {:>10.3} {:>12}",
            "new", count, new, new, "valid"
        );

        let (subject, clip) = (
            input(vec![polygon(points.clone())]),
            input(vec![obstacles[0].clone()]),
        );
        let mut corners_out = 0;
        let raw = median(runs, || {
            let out = overlay(
                &subject,
                &clip,
                OverlayOperation::Difference,
                FillRule::NonZero,
                t,
            )
            .expect("the difference");
            corners_out = out.polygons.iter().map(|p| p.outer.points.len()).sum();
            black_box(out);
        });
        println!(
            "{:<8} {:>8} {:>10.3} {:>10.3} {:>12}",
            "overlay", count, raw, raw, corners_out
        );

        let room = Region::new(vec![polygon(points.clone())], t).expect("a valid floor");
        let cuts: Vec<Region> = obstacles
            .iter()
            .map(|o| Region::new(vec![o.clone()], t).expect("a valid square"))
            .collect();
        let mut holes = 0;
        let chain = median(runs, || {
            let mut free = room.clone();
            for cut in &cuts {
                free = free.difference(cut, t).expect("the difference");
            }
            holes = free.evidence().holes;
            black_box(free);
        });
        println!(
            "{:<8} {:>8} {:>10.3} {:>10.3} {:>12}",
            "region",
            count,
            chain,
            chain / cuts.len() as f64,
            format!("{holes} holes")
        );
    }
}
