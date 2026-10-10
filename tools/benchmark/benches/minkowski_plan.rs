//! Wall-clock time and peak memory of the exact Minkowski sum of a floor's
//! walls with a 0.9 m square, against corner count (#292).
//!
//! The walls of a floor, merged into one region with a hole per room, are
//! summed with the square a person turns in; the free part of the floor
//! (the rooms) is eroded by it. A floor of 58 rooms, about 3,400 corners,
//! needed more than 12 GB before #292, because every piece of the
//! subdivision kept a flag for every ring. Rooms are rectangles or L-shapes
//! on a grid, walls 0.1 to 0.3 m thick, the floor about 130 m by 70 m,
//! turned 2.3 degrees and placed near (640 000, 5 650 000). Two wall
//! finishes are measured: `pilasters`, a 0.4 m pilaster every 3 m along
//! each room wall, and `ribs`, the same corners clustered as fine ribs on
//! one wall per room, where many edge hulls overlap (a 124-corner piece
//! costs about what the consumer's 114-corner piece did).
//!
//! Each case is validated before its time counts: every vertex of a sum
//! lies at Chebyshev distance 0.45 m from the walls (the square centred
//! there touches them without overlapping), and every vertex of an
//! erosion at 0.45 m from outside the rooms.
//!
//! Peak memory is the child process's resident high-water mark (`VmHWM`),
//! so each case runs in a fresh process: this binary re-runs itself with
//! the case in `MINKOWSKI_PLAN_CASE`. Set `MINKOWSKI_PLAN_ROOMS` (for
//! example `2,8`) to run fewer sizes.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench minkowski_plan
//! ```

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{Polygon, Region, Ring};
use std::hint::black_box;
use std::process::Command;
use std::time::Instant;

/// Rooms per measured size: about 114, 500, 1,000, 2,000 and 3,400 corners.
const ROOMS: [usize; 5] = [2, 8, 16, 33, 58];
/// Timed runs per case; the median is reported.
const RUNS: usize = 3;
/// Half the side of the square.
const HALF: f64 = 0.45;
/// Slack of the distance check: above the rounding of the vertex sums at
/// these coordinates and the tolerance by which settling may move a vertex.
const SLACK: f64 = 1e-5;

#[derive(Clone, Copy, PartialEq)]
enum Finish {
    Pilasters,
    Ribs,
}

/// A deterministic value in `[0, 1)`.
fn hash(i: u64) -> f64 {
    let mut x = i.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03;
    x ^= x >> 31;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 29;
    (x >> 11) as f64 / (1u64 << 53) as f64
}

/// The floor's outline, then one ring per room, counter-clockwise, in
/// local metres: four rows of 8.6 m by 17.3 m cells.
fn floor(rooms: usize, finish: Finish) -> Vec<Vec<(f64, f64)>> {
    let rows = rooms.min(4);
    let cols = rooms.div_ceil(rows);
    let (w, h) = (8.6, 17.3);
    let thickness = |i: usize, last: usize, seed: u64| {
        if i == 0 || i == last {
            0.3
        } else {
            0.1 + 0.2 * hash(seed + i as u64)
        }
    };
    let tv: Vec<f64> = (0..=cols).map(|i| thickness(i, cols, 0)).collect();
    let th: Vec<f64> = (0..=rows).map(|i| thickness(i, rows, 1000)).collect();
    let (x0, x1) = (-tv[0] / 2.0, cols as f64 * w + tv[cols] / 2.0);
    let (y0, y1) = (-th[0] / 2.0, rows as f64 * h + th[rows] / 2.0);
    let mut out = vec![vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]];
    for n in 1..=rooms {
        let (c, r) = ((n - 1) / rows, (n - 1) % rows);
        let ax = c as f64 * w + tv[c] / 2.0;
        let bx = (c + 1) as f64 * w - tv[c + 1] / 2.0;
        let ay = r as f64 * h + th[r] / 2.0;
        let by = (r + 1) as f64 * h - th[r + 1] / 2.0;
        let corners = if n % 3 == 0 {
            // An L: a shaft taken from the room's far corner.
            let (sx, sy) = (bx - 2.5, by - 3.1);
            vec![(ax, ay), (bx, ay), (bx, sy), (sx, sy), (sx, by), (ax, by)]
        } else {
            vec![(ax, ay), (bx, ay), (bx, by), (ax, by)]
        };
        out.push(match finish {
            Finish::Pilasters => pilasters(&corners),
            Finish::Ribs => ribs(&corners),
        });
    }
    out
}

/// A point `along` the edge `a -> b` and `inward` to its left.
fn at(a: (f64, f64), b: (f64, f64), along: f64, inward: f64) -> (f64, f64) {
    let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
    let d = ((b.0 - a.0) / len, (b.1 - a.1) / len);
    (
        a.0 + d.0 * along - d.1 * inward,
        a.1 + d.1 * along + d.0 * inward,
    )
}

/// A 0.4 m by 0.25 m pilaster into the room about every 3 m of wall.
fn pilasters(corners: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for i in 0..corners.len() {
        let (a, b) = (corners[i], corners[(i + 1) % corners.len()]);
        out.push(a);
        let len = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
        let count = ((len - 1.0) / 3.0).floor().max(0.0) as usize;
        let step = len / (count as f64 + 1.0);
        for j in 1..=count {
            let s = step * j as f64 - 0.2;
            out.extend([
                at(a, b, s, 0.0),
                at(a, b, s, 0.25),
                at(a, b, s + 0.4, 0.25),
                at(a, b, s + 0.4, 0.0),
            ]);
        }
    }
    out
}

/// Fourteen ribs 0.04 m wide and 0.06 m deep, 0.08 m apart, on the
/// room's first wall only.
fn ribs(corners: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let (a, b) = (corners[0], corners[1]);
    let mut out = vec![a];
    for j in 0..14 {
        let s = 1.0 + 0.08 * j as f64;
        out.extend([
            at(a, b, s, 0.0),
            at(a, b, s, 0.06),
            at(a, b, s + 0.04, 0.06),
            at(a, b, s + 0.04, 0.0),
        ]);
    }
    out.extend_from_slice(&corners[1..]);
    out
}

/// The floor turned 2.3 degrees and placed at survey coordinates.
fn placed(rooms: usize, finish: Finish) -> Vec<Vec<Point2>> {
    let (s, c) = 2.3f64.to_radians().sin_cos();
    floor(rooms, finish)
        .into_iter()
        .map(|ring| {
            ring.into_iter()
                .map(|(x, y)| Point2::new(640_000.0 + c * x - s * y, 5_650_000.0 + s * x + c * y))
                .collect()
        })
        .collect()
}

/// Whether the segment `a -> b` meets the closed box `lo..hi`
/// (Liang-Barsky clipping).
fn segment_meets_box(a: Point2, b: Point2, lo: Point2, hi: Point2) -> bool {
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for (p, q) in [
        (a.x - b.x, a.x - lo.x),
        (b.x - a.x, hi.x - a.x),
        (a.y - b.y, a.y - lo.y),
        (b.y - a.y, hi.y - a.y),
    ] {
        if p == 0.0 {
            if q < 0.0 {
                return false;
            }
        } else if p < 0.0 {
            t0 = t0.max(q / p);
        } else {
            t1 = t1.min(q / p);
        }
    }
    t0 <= t1
}

/// Even-odd membership over all rings.
fn inside(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut odd = false;
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if (a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x {
                odd = !odd;
            }
        }
    }
    odd
}

/// Whether the square of half side `r` about `c` meets the boundary.
fn box_meets_boundary(rings: &[Vec<Point2>], c: Point2, r: f64) -> bool {
    let (lo, hi) = (Point2::new(c.x - r, c.y - r), Point2::new(c.x + r, c.y + r));
    rings.iter().any(|ring| {
        (0..ring.len()).any(|i| segment_meets_box(ring[i], ring[(i + 1) % ring.len()], lo, hi))
    })
}

/// Every vertex of `out` lies at Chebyshev distance `HALF` from the set
/// `rings` bound (a sum), or from its complement (an erosion).
fn validate(rings: &[Vec<Point2>], out: &Region, erosion: bool) {
    assert!(!out.is_empty(), "an empty result");
    for polygon in out.polygons() {
        for ring in std::iter::once(&polygon.outer).chain(&polygon.holes) {
            for &v in &ring.points {
                let (near, far) = (HALF - SLACK, HALF + SLACK);
                let clear = !box_meets_boundary(rings, v, near);
                let touches = box_meets_boundary(rings, v, far);
                let ok = if erosion {
                    // The square fits inside, and a slightly larger one
                    // reaches the boundary.
                    clear && touches && inside(rings, v)
                } else {
                    // The square does not reach into the set, and a
                    // slightly larger one meets it.
                    clear && touches && !inside(rings, v)
                };
                assert!(ok, "vertex {v:?} is not at {HALF} m");
            }
        }
    }
}

/// One case in this process: prints the median time and peak memory.
fn case(spec: &str) {
    let parts: Vec<&str> = spec.split(':').collect();
    let finish = if parts[0] == "ribs" {
        Finish::Ribs
    } else {
        Finish::Pilasters
    };
    let erosion = parts[1] == "erosion";
    let rooms: usize = parts[2].parse().expect("a room count");
    let tolerance = Tolerance::METRE;
    let rings = placed(rooms, finish);
    let ring = |points: &Vec<Point2>| Ring {
        points: points.clone(),
    };
    let (set, polygons): (Vec<Vec<Point2>>, Vec<Polygon>) = if erosion {
        // The free part: the rooms.
        (
            rings[1..].to_vec(),
            rings[1..]
                .iter()
                .map(|r| Polygon {
                    outer: ring(r),
                    holes: Vec::new(),
                })
                .collect(),
        )
    } else {
        (
            rings.clone(),
            vec![Polygon {
                outer: ring(&rings[0]),
                holes: rings[1..].iter().map(ring).collect(),
            }],
        )
    };
    let corners: usize = set.iter().map(Vec::len).sum();
    let region = Region::new(polygons, tolerance).expect("a valid floor");
    let square = Ring {
        points: [(-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF)]
            .iter()
            .map(|&(x, y)| Point2::new(x, y))
            .collect(),
    };
    let run = || {
        if erosion {
            region.minkowski_erosion(&square, tolerance)
        } else {
            region.minkowski_sum(&square, tolerance)
        }
        .expect("the morphology")
    };
    let out = run();
    validate(&set, &out, erosion);
    let out_corners: usize = out
        .polygons()
        .iter()
        .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
        .map(|r| r.points.len())
        .sum();
    let mut times: Vec<f64> = (0..RUNS)
        .map(|_| {
            let start = Instant::now();
            black_box(run());
            start.elapsed().as_secs_f64()
        })
        .collect();
    times.sort_by(f64::total_cmp);
    let peak = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status
                .lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|kb| kb.parse::<f64>().ok())
        })
        .map_or_else(|| "n/a".to_string(), |kb| format!("{:.0}", kb / 1024.0));
    println!(
        "{:<10} {:<8} {:>5} {:>7} {:>7} {:>9.3} {:>8}",
        parts[0],
        parts[1],
        rooms,
        corners,
        out_corners,
        times[RUNS / 2],
        peak
    );
}

fn main() {
    if let Ok(spec) = std::env::var("MINKOWSKI_PLAN_CASE") {
        case(&spec);
        return;
    }
    let rooms: Vec<usize> = std::env::var("MINKOWSKI_PLAN_ROOMS")
        .ok()
        .map(|list| {
            list.split(',')
                .map(|n| n.trim().parse().expect("a room count"))
                .collect()
        })
        .unwrap_or_else(|| ROOMS.to_vec());
    println!(
        "{:<10} {:<8} {:>5} {:>7} {:>7} {:>9} {:>8}",
        "finish", "op", "rooms", "corners", "result", "median_s", "peak_mb"
    );
    let exe = std::env::current_exe().expect("this benchmark");
    for finish in ["pilasters", "ribs"] {
        for op in ["sum", "erosion"] {
            for &n in &rooms {
                let status = Command::new(&exe)
                    .env("MINKOWSKI_PLAN_CASE", format!("{finish}:{op}:{n}"))
                    .status()
                    .expect("a child run");
                assert!(status.success(), "{finish} {op} {n} rooms failed");
            }
        }
    }
}
