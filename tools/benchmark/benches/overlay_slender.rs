//! Wall-clock time of the exact straight-edge boolean on the plan footprint
//! of an assembly of overlapping slender parts, against part count (#315).
//!
//! `k` closed boxes, 1 m by 4 cm by 2 cm, at seeded turns and offsets over
//! a square metre, overlap one another in plan as the parts of an assembly
//! do (stair treads and stringers, frames, slats). Each box's twelve
//! triangles are projected to xy, those without area left out and the rest
//! wound counter-clockwise: 4 shadows per box. Two operations, as
//! axioval/engine measures a footprint:
//!
//! - `union`: `overlay` of the soup with itself, `Union`, `NonZero`.
//! - `square`: `overlay` of a 5 cm square at the origin's corner with the
//!   soup, `Intersection`, `NonZero`.
//!
//! The boxes cross about `k^2` times, so the arrangement grows that fast
//! while every union comes out as one or a few polygons. Before #315 the
//! time grew about `k^2.65`.
//!
//! Each case is validated before its time counts: its area must match the
//! grid backend's (`i_overlay`, the integer backend `overlay` used before
//! #173) to `1e-6` relative.
//!
//! Each figure is the median of `RUNS` timed runs (`OVERLAY_SLENDER_RUNS`
//! overrides the count), with the fastest and slowest beside it. Set
//! `OVERLAY_SLENDER_K` (for example `100,200`) to run fewer sizes.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench overlay_slender
//! ```

use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::{overlay, FillRule, OverlayInput, OverlayOperation, Polygon, Ring};
use i_overlay::core::fill_rule::FillRule as GridFill;
use i_overlay::core::overlay_rule::OverlayRule;
use i_overlay::float::single::SingleFloatOverlay;
use std::hint::black_box;
use std::time::Instant;

/// Part counts measured.
const KS: [usize; 3] = [100, 200, 400];
/// Timed runs per case; the median is reported.
const RUNS: usize = 3;

type P3 = [f64; 3];

/// `k` closed slender boxes as triangles, from a fixed xorshift seed.
fn assembly(k: usize) -> Vec<[P3; 3]> {
    let mut out = Vec::new();
    let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..k {
        let (cx, cy, cz, turn) = (next(), next(), next(), next() * std::f64::consts::PI);
        let (c, s) = (turn.cos(), turn.sin());
        let corner = |i: usize| {
            let x = if i & 1 == 0 { -0.5 } else { 0.5 };
            let y = if i & 2 == 0 { -0.02 } else { 0.02 };
            let z = if i & 4 == 0 { 0.0 } else { 0.02 };
            [cx + c * x - s * y, cy + s * x + c * y, cz + z]
        };
        let v: Vec<P3> = (0..8).map(corner).collect();
        for [a, b, c2, d] in [
            [0, 1, 3, 2],
            [4, 6, 7, 5],
            [0, 4, 5, 1],
            [2, 3, 7, 6],
            [0, 2, 6, 4],
            [1, 5, 7, 3],
        ] {
            out.push([v[a], v[b], v[c2]]);
            out.push([v[a], v[c2], v[d]]);
        }
    }
    out
}

/// Each triangle projected to xy and wound counter-clockwise, those
/// without area left out.
fn shadows(triangles: &[[P3; 3]]) -> Vec<Polygon> {
    triangles
        .iter()
        .filter_map(|t| {
            let mut points: Vec<Point2> = t.iter().map(|p| Point2::new(p[0], p[1])).collect();
            let (a, b, c) = (points[0], points[1], points[2]);
            let area = 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
            if area.abs() <= 1e-12 {
                return None;
            }
            if area < 0.0 {
                points.reverse();
            }
            Some(Polygon {
                outer: Ring { points },
                holes: Vec::new(),
            })
        })
        .collect()
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

fn square(x: f64, y: f64, side: f64) -> Polygon {
    Polygon {
        outer: Ring {
            points: vec![
                Point2::new(x, y),
                Point2::new(x + side, y),
                Point2::new(x + side, y + side),
                Point2::new(x, y + side),
            ],
        },
        holes: Vec::new(),
    }
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

fn contours(polygons: &[Polygon]) -> Contours {
    polygons
        .iter()
        .map(|p| p.outer.points.iter().map(|q| [q.x, q.y]).collect())
        .collect()
}

/// The grid backend's area of an operation, `NonZero`.
fn grid_area(subject: &Contours, clip: &Contours, rule: OverlayRule) -> f64 {
    subject
        .overlay(clip, rule, GridFill::NonZero)
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

/// Minimum, median and maximum of the timed runs, in seconds.
fn median_s(mut run: impl FnMut()) -> (f64, f64, f64) {
    let runs = std::env::var("OVERLAY_SLENDER_RUNS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(RUNS)
        .max(1);
    let mut times: Vec<f64> = (0..runs)
        .map(|_| {
            let start = Instant::now();
            run();
            start.elapsed().as_secs_f64()
        })
        .collect();
    times.sort_by(f64::total_cmp);
    (times[0], times[runs / 2], times[runs - 1])
}

fn check(name: &str, exact: f64, reference: f64) {
    let error = (exact - reference).abs() / reference.abs().max(1e-12);
    assert!(
        error < 1e-6,
        "{name}: exact area {exact} differs from the grid's {reference}"
    );
}

fn main() {
    let ks: Vec<usize> = std::env::var("OVERLAY_SLENDER_K")
        .ok()
        .map(|v| v.split(',').filter_map(|k| k.trim().parse().ok()).collect())
        .unwrap_or_else(|| KS.to_vec());
    let tolerance = Tolerance::new(1.0e-9, 1.0e-9).expect("a valid tolerance");
    for k in ks {
        let kept = shadows(&assembly(k));
        let count = kept.len();
        let grid_body = contours(&kept);
        let body = input(kept);
        let small = input(vec![square(0.0, 0.0, 0.05)]);
        let grid_small = contours(&small.polygons);

        let union = || {
            overlay(
                &body,
                &body,
                OverlayOperation::Union,
                FillRule::NonZero,
                tolerance,
            )
            .expect("the union of valid shadows")
        };
        let meet = || {
            overlay(
                &small,
                &body,
                OverlayOperation::Intersection,
                FillRule::NonZero,
                tolerance,
            )
            .expect("the intersection of valid shadows")
        };
        let united = union();
        check(
            &format!("union k = {k}"),
            exact_area(&united.polygons),
            grid_area(&grid_body, &grid_body, OverlayRule::Union),
        );
        let met = meet();
        check(
            &format!("square k = {k}"),
            exact_area(&met.polygons),
            grid_area(&grid_small, &grid_body, OverlayRule::Intersect),
        );
        let u = median_s(|| {
            black_box(union());
        });
        let m = median_s(|| {
            black_box(meet());
        });
        println!(
            "k = {k:>3} ({count:>4} shadows)  union {:8.3} s (min {:7.3}, max {:7.3}, {} polygons)   \
             square {:8.3} s (min {:7.3}, max {:7.3}, {} polygons)",
            u.1,
            u.0,
            u.2,
            united.polygons.len(),
            m.1,
            m.0,
            m.2,
            met.polygons.len(),
        );
    }
}
