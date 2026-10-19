//! Wall-clock time and peak memory of the disc morphology
//! `Region::erode_inner`, `erode_outer`, `dilate_inner` and `dilate_outer`
//! against corner count (#323, #330).
//!
//! The room is #323's reproducer: 4.5 m by 5 m, its two long walls
//! carrying teeth 0.5 mm deep. Overlay output along tessellated or flush
//! bodies looks like this. Before #323 its erosion by a disc of 0.75 m
//! grew with the square of the corners, 1,202 corners taking 5 s and
//! 1.6 GB. Measured at 102, 602 and 1,202 corners, with teeth (`teeth`)
//! and with the same corners exactly collinear (`flat`, teeth of depth
//! 0), beside the plain rectangle (`plain`, 4 corners).
//!
//! The cases of #330 take a disc of 0.4 m and 2,402, 5,002 and 10,002
//! corners: the toothed room eroded and dilated, and the same toothed
//! ring as a hole (`hole`) in a plain room 1.5 m larger all round, eroded
//! and dilated. A hole's boundary and a dilation sum on the other side of
//! the teeth than a room's erosion does.
//!
//! Each case is validated before its time counts: the result is not
//! empty, every vertex lies inside the region for an erosion and outside
//! it for a dilation, an inner erosion's vertices and an outer dilation's
//! at least the radius from the walls, an outer erosion's and an inner
//! dilation's within the radius plus the reported deviation (all within a
//! small slack for the margin and the rounding). The hash shows the result
//! unchanged bit for bit across a change.
//!
//! Peak memory is the child process's resident high-water mark (`VmHWM`),
//! so each case runs in a fresh process: this binary re-runs itself with
//! the case in `DISC_EROSION_CASE`. Set `DISC_EROSION_CORNERS` (for example
//! `102,602`) or `DISC_EROSION_LARGE` (for example `2402`) to run other
//! sizes of #323's and #330's cases, an empty list to run none, and
//! `DISC_EROSION_RUNS` for fewer timed runs.
//!
//! ```bash
//! cargo bench -p axiolid-benchmark --bench disc_erosion
//! ```

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{Polygon, Region, Ring};
use std::hint::black_box;
use std::process::Command;
use std::time::Instant;

/// Corners per measured size: 25, 150 and 300 teeth per wall.
const CORNERS: [usize; 3] = [102, 602, 1202];
/// Corners of #330's cases: 600, 1,250 and 2,500 teeth per wall.
const LARGE: [usize; 3] = [2402, 5002, 10002];
/// Timed runs per case; the median is reported.
const RUNS: usize = 3;
/// The disc's radius in #323's cases.
const RADIUS: f64 = 0.75;
/// The disc's radius in #330's cases.
const LARGE_RADIUS: f64 = 0.4;
/// How far the room around a hole reaches past it.
const AROUND: f64 = 1.5;
/// The room's width (along the toothed walls) and depth.
const ROOM: (f64, f64) = (4.5, 5.0);
/// Slack of the distance check: the margin kept for rounding and
/// presentation, and the tolerance by which settling may move a vertex.
const SLACK: f64 = 1e-6;

/// The room's corners, counter-clockwise: the issue's reproducer with
/// `teeth` teeth of `depth` on each long wall.
fn room(teeth: usize, depth: f64) -> Vec<Point2> {
    let (w, h) = ROOM;
    if teeth == 0 {
        return vec![
            Point2::new(0.0, 0.0),
            Point2::new(w, 0.0),
            Point2::new(w, h),
            Point2::new(0.0, h),
        ];
    }
    let n = 2 * teeth as u32;
    let at = |k: u32| w * f64::from(k) / f64::from(n);
    let tooth = |k: u32| if k % 2 == 1 { depth } else { 0.0 };
    let mut points: Vec<Point2> = (0..=n).map(|k| Point2::new(at(k), tooth(k))).collect();
    points.extend((0..=n).map(|k| Point2::new(w - at(k), h - tooth(k))));
    points
}

/// The distance from `p` to the segment `a -> b`.
fn to_segment(p: Point2, a: Point2, b: Point2) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
    ((p.x - a.x - t * dx).powi(2) + (p.y - a.y - t * dy).powi(2)).sqrt()
}

/// Even-odd membership in the rings.
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

/// Every vertex of `out` lies inside the region for an erosion, outside it
/// for a dilation, and on the operation's side of the radius.
fn validate(rings: &[Vec<Point2>], out: &Region, erode: bool, outer: bool, radius: f64) {
    assert!(!out.is_empty(), "an empty result");
    let deviation = out.bound().expect("a disc morphology's bound").deviation;
    // An outer erosion and an inner dilation reach no farther than the
    // radius plus the deviation; the other two no nearer than the radius.
    let within = outer == erode;
    for polygon in out.polygons() {
        for r in std::iter::once(&polygon.outer).chain(&polygon.holes) {
            for &v in &r.points {
                let d = rings
                    .iter()
                    .flat_map(|ring| {
                        (0..ring.len()).map(|i| to_segment(v, ring[i], ring[(i + 1) % ring.len()]))
                    })
                    .fold(f64::INFINITY, f64::min);
                assert_eq!(inside(rings, v), erode, "vertex {v:?} on the wrong side");
                if within {
                    assert!(d <= radius + deviation + SLACK, "vertex {v:?} at {d}");
                } else {
                    assert!(d >= radius - SLACK, "vertex {v:?} at {d}");
                }
            }
        }
    }
}

/// An FNV-1a hash of the result's vertex bits, in order: equal hashes
/// before and after a change show the result unchanged bit for bit.
fn fingerprint(out: &Region) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |bits: u64| {
        for byte in bits.to_le_bytes() {
            h = (h ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
        }
    };
    for polygon in out.polygons() {
        for r in std::iter::once(&polygon.outer).chain(&polygon.holes) {
            eat(r.points.len() as u64);
            for p in &r.points {
                eat(p.x.to_bits());
                eat(p.y.to_bits());
            }
        }
    }
    h
}

/// One case in this process, given as `walls:op:side:corners:radius`:
/// prints the median time and peak memory.
fn case(spec: &str) {
    let parts: Vec<&str> = spec.split(':').collect();
    let corners: usize = parts[3].parse().expect("a corner count");
    let radius: f64 = parts[4].parse().expect("a radius");
    let (teeth, depth) = match parts[0] {
        "plain" => (0, 0.0),
        "flat" => ((corners - 2) / 4, 0.0),
        _ => ((corners - 2) / 4, 5.0e-4),
    };
    let erode = parts[1] == "erode";
    let outer = parts[2] == "outer";
    let runs: usize = std::env::var("DISC_EROSION_RUNS")
        .ok()
        .and_then(|r| r.parse().ok())
        .unwrap_or(RUNS);
    let tolerance = Tolerance::new(1.0e-9, 1.0e-9).expect("a tolerance");
    let ring = room(teeth, depth);
    let rings: Vec<Vec<Point2>> = if parts[0] == "hole" {
        let (w, h) = ROOM;
        let around = vec![
            Point2::new(-AROUND, -AROUND),
            Point2::new(w + AROUND, -AROUND),
            Point2::new(w + AROUND, h + AROUND),
            Point2::new(-AROUND, h + AROUND),
        ];
        vec![around, ring.iter().rev().copied().collect()]
    } else {
        vec![ring.clone()]
    };
    let region = Region::new(
        vec![Polygon {
            outer: Ring {
                points: rings[0].clone(),
            },
            holes: rings[1..]
                .iter()
                .map(|r| Ring { points: r.clone() })
                .collect(),
        }],
        tolerance,
    )
    .expect("a valid room");
    let run = || {
        match (erode, outer) {
            (true, true) => region.erode_outer(radius, tolerance),
            (true, false) => region.erode_inner(radius, tolerance),
            (false, true) => region.dilate_outer(radius, tolerance),
            (false, false) => region.dilate_inner(radius, tolerance),
        }
        .expect("the disc morphology")
    };
    let out = run();
    validate(&rings, &out, erode, outer, radius);
    let out_corners: usize = out
        .polygons()
        .iter()
        .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
        .map(|r| r.points.len())
        .sum();
    let mut times: Vec<f64> = (0..runs)
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
        "{:<6} {:<6} {:<6} {:>6} {:>7} {:>7} {:>9.3} {:>8} {:016x}",
        parts[0],
        parts[1],
        parts[2],
        radius,
        ring.len(),
        out_corners,
        times[runs / 2],
        peak,
        fingerprint(&out)
    );
}

fn main() {
    if let Ok(spec) = std::env::var("DISC_EROSION_CASE") {
        case(&spec);
        return;
    }
    let sizes = |name: &str, default: &[usize]| -> Vec<usize> {
        std::env::var(name)
            .ok()
            .map(|list| {
                list.split(',')
                    .filter(|n| !n.trim().is_empty())
                    .map(|n| n.trim().parse().expect("a corner count"))
                    .collect()
            })
            .unwrap_or_else(|| default.to_vec())
    };
    let corners = sizes("DISC_EROSION_CORNERS", &CORNERS);
    let large = sizes("DISC_EROSION_LARGE", &LARGE);
    println!(
        "{:<6} {:<6} {:<6} {:>6} {:>7} {:>7} {:>9} {:>8} {:>16}",
        "walls", "op", "side", "radius", "corners", "result", "median_s", "peak_mb", "hash"
    );
    let exe = std::env::current_exe().expect("this benchmark");
    let mut cases = Vec::new();
    if !corners.is_empty() {
        for side in ["inner", "outer"] {
            cases.push(format!("plain:erode:{side}:4:{RADIUS}"));
            for walls in ["teeth", "flat"] {
                for &n in &corners {
                    cases.push(format!("{walls}:erode:{side}:{n}:{RADIUS}"));
                }
            }
        }
    }
    // #330: the teeth seen from the other side, by a hole or a dilation.
    for walls in ["teeth", "hole"] {
        for op in ["erode", "dilate"] {
            for side in ["inner", "outer"] {
                for &n in &large {
                    cases.push(format!("{walls}:{op}:{side}:{n}:{LARGE_RADIUS}"));
                }
            }
        }
    }
    for spec in cases {
        let status = Command::new(&exe)
            .env("DISC_EROSION_CASE", &spec)
            .status()
            .expect("a child run");
        assert!(status.success(), "{spec} failed");
    }
}
