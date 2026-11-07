//! Wall-clock time and peak memory of the exact Minkowski sum and erosion
//! of a cluster of furniture outlines with a 1.7 m by 0.9 m rectangle,
//! against corner count (#344).
//!
//! A group is a C-shaped sofa (an annular sector of radii 1.0 and 0.6 over
//! 240 degrees, opening towards +x), a round table in its notch and an
//! L-shaped desk beside it. With 48 points per arc a group has 152
//! corners. Groups sit on a grid 4.0 m by 2.4 m, five to a row, so the
//! rectangle is larger than the gaps between outlines and every edge hull
//! overlaps most others: the case a rule engine meets when it searches a
//! furnished room for a place a rectangle fits. 1, 2, 7 and 20 groups
//! give 152, 304, 1,064 and 3,040 corners. The scene is turned 2.3
//! degrees and placed near (640 000, 5 650 000), as in `minkowski_plan`.
//!
//! The sum takes the outlines as one polygon each, with the rectangle; the
//! erosion takes the free floor (the scene's box grown by 1.2 m, the
//! outlines as holes) and the rectangle. The rectangle is measured
//! axis-aligned (`axis`) and turned 30 degrees (`turned`).
//!
//! Each case is validated before its time counts: a copy of the rectangle
//! centred on a result vertex and shrunk about it by 1e-5 m meets no
//! outline and lies outside the furniture (inside the free floor for an
//! erosion), and a copy grown by 1e-5 m meets an outline.
//!
//! The last column is an FNV-1a hash of every vertex of the result, with
//! the same digest as `tests/minkowski.rs`: equal hashes before and after
//! a change show the result unchanged bit for bit. Peak memory is the
//! child process's resident high-water mark (`VmHWM`), so each case runs in
//! a fresh process: this binary re-runs itself with the case in
//! `MINKOWSKI_FURNITURE_CASE`. Set `MINKOWSKI_FURNITURE_CORNERS` (for
//! example `152,304`) to run fewer sizes and `MINKOWSKI_FURNITURE_RUNS` for
//! fewer timed runs.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench minkowski_furniture
//! ```

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{Polygon, Region, Ring};
use std::hint::black_box;
use std::process::Command;
use std::time::Instant;

/// Groups per measured size: 152, 304, 1,064 and 3,040 corners.
const GROUPS: [usize; 4] = [1, 2, 7, 20];
/// Timed runs per case; the median is reported.
const RUNS: usize = 3;
/// Points per arc, and sides of the table.
const ARC: usize = 48;
/// Groups per row of the grid.
const PER_ROW: usize = 5;
/// Slack of the distance check, as in `minkowski_plan`.
const SLACK: f64 = 1e-5;
/// The rectangle's sides.
const WIDTH: f64 = 1.7;
const HEIGHT: f64 = 0.9;

/// One group's three outlines about the sofa's centre `(cx, cy)`,
/// counter-clockwise, in local metres.
fn group(cx: f64, cy: f64) -> [Vec<(f64, f64)>; 3] {
    let arc = |radius: f64, from: f64, to: f64| -> Vec<(f64, f64)> {
        (0..=ARC)
            .map(|i| {
                let a = (from + (to - from) * i as f64 / ARC as f64).to_radians();
                (cx + radius * a.cos(), cy + radius * a.sin())
            })
            .collect()
    };
    // The gap faces +x: the sector spans 60..300 degrees.
    let mut sofa = arc(1.0, 60.0, 300.0);
    sofa.extend(arc(0.6, 300.0, 60.0));
    let table = (0..ARC)
        .map(|i| {
            let a = (360.0 * i as f64 / ARC as f64).to_radians();
            (cx + 0.35 * a.cos(), cy + 0.35 * a.sin())
        })
        .collect();
    let desk = [
        (1.2, -0.7),
        (2.8, -0.7),
        (2.8, 0.0),
        (1.9, 0.0),
        (1.9, 0.7),
        (1.2, 0.7),
    ]
    .iter()
    .map(|&(x, y)| (cx + x, cy + y))
    .collect();
    [sofa, table, desk]
}

/// The scene's outlines turned 2.3 degrees and placed at survey
/// coordinates, and the corners of its box grown by 1.2 m.
fn placed(groups: usize) -> (Vec<Vec<Point2>>, Vec<Point2>) {
    let (s, c) = 2.3f64.to_radians().sin_cos();
    let at =
        |(x, y): (f64, f64)| Point2::new(640_000.0 + c * x - s * y, 5_650_000.0 + s * x + c * y);
    let mut local = Vec::new();
    for g in 0..groups {
        local.extend(group(
            4.0 * (g % PER_ROW) as f64,
            2.4 * (g / PER_ROW) as f64,
        ));
    }
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(x, y) in local.iter().flatten() {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    let m = 1.2;
    let floor = [
        (lo.0 - m, lo.1 - m),
        (hi.0 + m, lo.1 - m),
        (hi.0 + m, hi.1 + m),
        (lo.0 - m, hi.1 + m),
    ]
    .map(at)
    .to_vec();
    let outlines = local
        .into_iter()
        .map(|ring| ring.into_iter().map(at).collect())
        .collect();
    (outlines, floor)
}

/// The rectangle about the origin, turned by `angle` radians.
fn rectangle(angle: f64) -> Vec<Point2> {
    let (s, c) = angle.sin_cos();
    let (hw, hh) = (0.5 * WIDTH, 0.5 * HEIGHT);
    [(-hw, -hh), (hw, -hh), (hw, hh), (-hw, hh)]
        .iter()
        .map(|&(x, y)| Point2::new(c * x - s * y, s * x + c * y))
        .collect()
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

/// Whether the segment `a -> b` meets the convex counter-clockwise polygon
/// (Cyrus-Beck clipping against its sides).
fn segment_meets_convex(a: Point2, b: Point2, poly: &[Point2]) -> bool {
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for i in 0..poly.len() {
        let (u, v) = (poly[i], poly[(i + 1) % poly.len()]);
        // Inside is the left of u -> v.
        let side = |p: Point2| (v.x - u.x) * (p.y - u.y) - (v.y - u.y) * (p.x - u.x);
        let (da, db) = (side(a), side(b));
        if da < 0.0 && db < 0.0 {
            return false;
        }
        if da < 0.0 {
            t0 = t0.max(da / (da - db));
        } else if db < 0.0 {
            t1 = t1.min(da / (da - db));
        }
    }
    t0 <= t1
}

/// The shape `k` centred on `v` and scaled about `v` by `scale`.
fn around(k: &[Point2], v: Point2, scale: f64) -> Vec<Point2> {
    k.iter()
        .map(|p| Point2::new(v.x + scale * p.x, v.y + scale * p.y))
        .collect()
}

fn meets_boundary(rings: &[Vec<Point2>], poly: &[Point2]) -> bool {
    rings.iter().any(|ring| {
        (0..ring.len()).any(|i| segment_meets_convex(ring[i], ring[(i + 1) % ring.len()], poly))
    })
}

/// Every vertex of `out` is a position where the (centrally symmetric)
/// rectangle `k` touches the boundary of the set `rings` without
/// overlapping it (a sum), or fits in the free floor `rings` bound (an
/// erosion).
fn validate(rings: &[Vec<Point2>], k: &[Point2], out: &Region, erosion: bool) {
    assert!(!out.is_empty(), "an empty result");
    for polygon in out.polygons() {
        for ring in std::iter::once(&polygon.outer).chain(&polygon.holes) {
            for &v in &ring.points {
                let shrunk = around(k, v, 1.0 - SLACK / (0.5 * HEIGHT));
                let grown = around(k, v, 1.0 + SLACK / (0.5 * HEIGHT));
                let clear = !meets_boundary(rings, &shrunk);
                let touches = meets_boundary(rings, &grown);
                // `rings` bound the furniture in a sum and the free floor
                // in an erosion.
                let side = inside(rings, v) == erosion;
                assert!(clear && touches && side, "vertex {v:?} is not at the shape");
            }
        }
    }
}

/// An FNV-1a digest of every vertex of the result, as `tests/minkowski.rs`.
fn digest(out: &Region) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for p in out.boundary_rings().iter().flat_map(|r| r.points.clone()) {
        for b in
            p.x.to_bits()
                .to_le_bytes()
                .into_iter()
                .chain(p.y.to_bits().to_le_bytes())
        {
            h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
        }
    }
    (h ^ 0xff).wrapping_mul(0x0100_0000_01b3)
}

/// One case in this process, given as `op:shape:groups`: prints the median
/// time and peak memory.
fn case(spec: &str) {
    let parts: Vec<&str> = spec.split(':').collect();
    let erosion = parts[0] == "erosion";
    let turned = parts[1] == "turned";
    let groups: usize = parts[2].parse().expect("a group count");
    let runs: usize = std::env::var("MINKOWSKI_FURNITURE_RUNS")
        .ok()
        .map_or(RUNS, |n| n.parse().expect("a run count"));
    let tolerance = Tolerance::METRE;
    let (outlines, floor) = placed(groups);
    let ring = |points: &Vec<Point2>| Ring {
        points: points.clone(),
    };
    let (set, polygons): (Vec<Vec<Point2>>, Vec<Polygon>) = if erosion {
        let mut set = vec![floor.clone()];
        set.extend(outlines.iter().cloned());
        let polygon = Polygon {
            outer: ring(&floor),
            holes: outlines.iter().map(ring).collect(),
        };
        (set, vec![polygon])
    } else {
        let polygons = outlines
            .iter()
            .map(|r| Polygon {
                outer: ring(r),
                holes: Vec::new(),
            })
            .collect();
        (outlines.clone(), polygons)
    };
    let corners: usize = outlines.iter().map(Vec::len).sum();
    let region = Region::new(polygons, tolerance).expect("a valid scene");
    let points = rectangle(if turned { 30f64.to_radians() } else { 0.0 });
    let shape = Ring {
        points: points.clone(),
    };
    let run = || {
        if erosion {
            region.minkowski_erosion(&shape, tolerance)
        } else {
            region.minkowski_sum(&shape, tolerance)
        }
        .expect("the morphology")
    };
    // The first run also feeds the validation; later runs are timed alone.
    let start = Instant::now();
    let out = run();
    let first = start.elapsed().as_secs_f64();
    validate(&set, &points, &out, erosion);
    let out_corners: usize = out
        .polygons()
        .iter()
        .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
        .map(|r| r.points.len())
        .sum();
    let mut times: Vec<f64> = (1..runs)
        .map(|_| {
            let start = Instant::now();
            black_box(run());
            start.elapsed().as_secs_f64()
        })
        .collect();
    times.push(first);
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
        "{:<8} {:<7} {:>7} {:>7} {:>9.3} {:>8}  {:016x}",
        parts[0],
        parts[1],
        corners,
        out_corners,
        times[times.len() / 2],
        peak,
        digest(&out)
    );
}

fn main() {
    if let Ok(spec) = std::env::var("MINKOWSKI_FURNITURE_CASE") {
        case(&spec);
        return;
    }
    let groups: Vec<usize> = std::env::var("MINKOWSKI_FURNITURE_CORNERS")
        .ok()
        .map(|list| {
            list.split(',')
                .map(|n| n.trim().parse::<usize>().expect("a corner count") / (3 * ARC + 8))
                .collect()
        })
        .unwrap_or_else(|| GROUPS.to_vec());
    println!(
        "{:<8} {:<7} {:>7} {:>7} {:>9} {:>8}  {:<16}",
        "op", "shape", "corners", "result", "median_s", "peak_mb", "hash"
    );
    let exe = std::env::current_exe().expect("this benchmark");
    for op in ["sum", "erosion"] {
        for shape in ["axis", "turned"] {
            for &n in &groups {
                let status = Command::new(&exe)
                    .env("MINKOWSKI_FURNITURE_CASE", format!("{op}:{shape}:{n}"))
                    .status()
                    .expect("a child run");
                assert!(status.success(), "{op} {shape} {n} groups failed");
            }
        }
    }
}
