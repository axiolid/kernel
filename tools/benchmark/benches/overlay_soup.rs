//! Wall-clock time of the exact straight-edge boolean on soups whose rings
//! share no edges, against the grid backend it replaced (#218).
//!
//! Since #173, `overlay` and `union_soup` run on the exact subdivision of
//! `axiolid-overlay`. A mesh given as a soup of its triangles shrinks to its
//! outline first (shared edges cancel), so it is close to the grid backend;
//! rings that share no edges get no such help. Three cases:
//!
//! - `random triangles`: `union_soup` of 2,000 triangles with random double
//!   vertices, each within 4 m of a random centre in a 100 m square, so most
//!   overlap several others and no two share an edge.
//! - `two 1000-gons`: `overlay` union of two star-shaped 1000-gons with
//!   random radii, offset so their boundaries cross about a thousand times.
//! - `two meshes`: `union_soup` of two overlapping 60 x 60 triangle meshes
//!   (14,400 triangles), the shared-edge case, which must not regress.
//!
//! The grid column is `i_overlay` (the integer backend `overlay` used before
//! #173) on the same rings and fill rule, as the reference the issue
//! measured against. Each case is validated before its time counts: the
//! exact result's area must match the grid's to `1e-6` relative, and for
//! the meshes the area known from how they were built.
//!
//! Each figure is the median of `RUNS` timed runs (`OVERLAY_SOUP_RUNS`
//! overrides the count), with the fastest and slowest beside it; set
//! `OVERLAY_SOUP_CASE` to a substring of a case name to run only that case.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench overlay_soup
//! ```

use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::{
    overlay, union_soup, FillRule, OverlayInput, OverlayOperation, Polygon, Ring,
};
use i_overlay::core::fill_rule::FillRule as GridFill;
use i_overlay::core::overlay_rule::OverlayRule;
use i_overlay::float::single::SingleFloatOverlay;
use std::hint::black_box;
use std::time::Instant;

/// Timed runs per case; the median is reported.
const RUNS: usize = 9;

/// A deterministic value in `[0, 1)`.
fn hash(i: u64) -> f64 {
    let mut x = i.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03;
    x ^= x >> 31;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 29;
    (x >> 11) as f64 / (1u64 << 53) as f64
}

fn ring(points: Vec<Point2>) -> Ring {
    Ring { points }
}

/// 2,000 triangles with random vertices near random centres, each
/// counter-clockwise and far from degenerate.
fn random_triangles() -> Vec<Ring> {
    let mut out = Vec::new();
    let mut seed = 0u64;
    let mut next = || {
        seed += 1;
        hash(seed)
    };
    while out.len() < 2_000 {
        let (cx, cy) = (100.0 * next(), 100.0 * next());
        let mut p: Vec<Point2> = (0..3)
            .map(|_| Point2::new(cx + 8.0 * next() - 4.0, cy + 8.0 * next() - 4.0))
            .collect();
        let area2 = (p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[1].y - p[0].y) * (p[2].x - p[0].x);
        if area2.abs() < 1.0 {
            continue;
        }
        if area2 < 0.0 {
            p.swap(1, 2);
        }
        out.push(ring(p));
    }
    out
}

/// A star-shaped `n`-gon around `centre`, radii in `[r, r + 1)`.
fn star(n: usize, centre: Point2, r: f64, seed: u64) -> Ring {
    let step = std::f64::consts::TAU / n as f64;
    ring(
        (0..n)
            .map(|i| {
                let t = i as f64 * step;
                let radius = r + hash(seed + i as u64);
                Point2::new(centre.x + radius * t.cos(), centre.y + radius * t.sin())
            })
            .collect(),
    )
}

/// Two 60 x 60 meshes of unit cells, two triangles per cell, the second
/// shifted by `(20.5, 30.25)`: 14,400 triangles, union area
/// `2 * 3600 - 39.5 * 29.75`.
fn two_meshes() -> (Vec<Ring>, f64) {
    let mut out = Vec::new();
    for (dx, dy) in [(0.0, 0.0), (20.5, 30.25)] {
        for i in 0..60 {
            for j in 0..60 {
                let (x, y) = (dx + f64::from(i), dy + f64::from(j));
                let a = Point2::new(x, y);
                let b = Point2::new(x + 1.0, y);
                let c = Point2::new(x + 1.0, y + 1.0);
                let d = Point2::new(x, y + 1.0);
                out.push(ring(vec![a, b, c]));
                out.push(ring(vec![a, c, d]));
            }
        }
    }
    (out, 2.0 * 3600.0 - 39.5 * 29.75)
}

fn ring_area(points: &[Point2]) -> f64 {
    let n = points.len();
    (0..n)
        .map(|i| {
            let (p, q) = (points[i], points[(i + 1) % n]);
            p.x * q.y - q.x * p.y
        })
        .sum::<f64>()
        / 2.0
}

fn exact_area(polygons: &[Polygon]) -> f64 {
    polygons
        .iter()
        .map(|p| {
            ring_area(&p.outer.points).abs()
                - p.holes
                    .iter()
                    .map(|h| ring_area(&h.points).abs())
                    .sum::<f64>()
        })
        .sum()
}

type Contours = Vec<Vec<[f64; 2]>>;

fn contours(rings: &[Ring]) -> Contours {
    rings
        .iter()
        .map(|r| r.points.iter().map(|p| [p.x, p.y]).collect())
        .collect()
}

/// The grid backend's union, NonZero, and its area.
fn grid_union(subject: &Contours, clip: &Contours) -> f64 {
    let shapes = subject.overlay(clip, OverlayRule::Union, GridFill::NonZero);
    shapes
        .iter()
        .flatten()
        .map(|contour| {
            let n = contour.len();
            (0..n)
                .map(|i| {
                    let (p, q) = (contour[i], contour[(i + 1) % n]);
                    p[0] * q[1] - q[0] * p[1]
                })
                .sum::<f64>()
                / 2.0
        })
        .sum::<f64>()
        .abs()
}

/// Minimum, median and maximum of the timed runs, in milliseconds.
fn median_ms(mut run: impl FnMut()) -> (f64, f64, f64) {
    let runs = std::env::var("OVERLAY_SOUP_RUNS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(RUNS)
        .max(1);
    let mut times: Vec<f64> = (0..runs)
        .map(|_| {
            let start = Instant::now();
            run();
            start.elapsed().as_secs_f64() * 1e3
        })
        .collect();
    times.sort_by(f64::total_cmp);
    (times[0], times[runs / 2], times[runs - 1])
}

fn check(name: &str, exact: f64, reference: f64, what: &str) {
    let error = (exact - reference).abs() / reference.abs().max(1.0);
    assert!(
        error < 1e-6,
        "{name}: exact area {exact} differs from {what} {reference}"
    );
}

fn report(name: &str, exact: (f64, f64, f64), grid: (f64, f64, f64)) {
    println!(
        "{name:<18} exact {:9.2} ms (min {:8.2}, max {:8.2})   grid {:7.2} ms   ratio {:5.1}x",
        exact.1,
        exact.0,
        exact.2,
        grid.1,
        exact.1 / grid.1
    );
}

fn main() {
    let only = std::env::var("OVERLAY_SOUP_CASE").ok();
    let wanted = |name: &str| only.as_deref().is_none_or(|o| name.contains(o));
    let tolerance = Tolerance::METRE;
    let empty: Contours = Vec::new();

    let name = "random triangles";
    if wanted(name) {
        let rings = random_triangles();
        let grid_rings = contours(&rings);
        let reference = grid_union(&grid_rings, &empty);
        check(
            name,
            exact_area(&union_soup(&rings, tolerance).unwrap()),
            reference,
            "the grid's",
        );
        let exact = median_ms(|| {
            black_box(union_soup(black_box(&rings), tolerance).unwrap());
        });
        let grid = median_ms(|| {
            black_box(grid_union(black_box(&grid_rings), &empty));
        });
        report(name, exact, grid);
    }

    let name = "two 1000-gons";
    if wanted(name) {
        let a = star(1000, Point2::new(0.0, 0.0), 10.0, 0);
        let b = star(1000, Point2::new(3.0, 1.7), 10.0, 5_000);
        let frame = Frame2 {
            origin: Point2::ZERO,
            x: Point2::X,
            y: Point2::Y,
        };
        let input = |r: &Ring| OverlayInput {
            frame,
            polygons: vec![Polygon {
                outer: r.clone(),
                holes: Vec::new(),
            }],
        };
        let (sa, sb) = (input(&a), input(&b));
        let run = || {
            overlay(
                &sa,
                &sb,
                OverlayOperation::Union,
                FillRule::NonZero,
                tolerance,
            )
            .unwrap()
        };
        let (ga, gb) = (contours(&[a.clone()]), contours(&[b.clone()]));
        check(
            name,
            exact_area(&run().polygons),
            grid_union(&ga, &gb),
            "the grid's",
        );
        let exact = median_ms(|| {
            black_box(run());
        });
        let grid = median_ms(|| {
            black_box(grid_union(black_box(&ga), &gb));
        });
        report(name, exact, grid);
    }

    let name = "two meshes";
    if wanted(name) {
        let (rings, area) = two_meshes();
        check(
            name,
            exact_area(&union_soup(&rings, tolerance).unwrap()),
            area,
            "the built area",
        );
        let grid_rings = contours(&rings);
        let exact = median_ms(|| {
            black_box(union_soup(black_box(&rings), tolerance).unwrap());
        });
        let grid = median_ms(|| {
            black_box(grid_union(black_box(&grid_rings), &empty));
        });
        report(name, exact, grid);
    }
}
