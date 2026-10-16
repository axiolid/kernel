//! Wall-clock time and peak memory of the disc erosions
//! `Region::erode_inner` and `Region::erode_outer` against corner count
//! (#323).
//!
//! The room is the issue's reproducer: 4.5 m by 5 m, its two long walls
//! carrying teeth 0.5 mm deep, eroded by a disc of 0.75 m. Overlay output
//! along tessellated or flush bodies looks like this. Before #323 the cost
//! grew with the square of the corners, 1,202 corners taking 5 s and
//! 1.6 GB. Measured at 102, 602 and 1,202 corners, with teeth (`teeth`)
//! and with the same corners exactly collinear (`flat`, teeth of depth
//! 0), beside the plain rectangle (`plain`, 4 corners).
//!
//! Each case is validated before its time counts: the result is not
//! empty, every vertex lies inside the room, an inner erosion's vertices
//! at least the radius from the walls and an outer erosion's within the
//! radius plus the reported deviation (both within a small slack for the
//! margin and the rounding).
//!
//! Peak memory is the child process's resident high-water mark (`VmHWM`),
//! so each case runs in a fresh process: this binary re-runs itself with
//! the case in `DISC_EROSION_CASE`. Set `DISC_EROSION_CORNERS` (for example
//! `102,602`) to run fewer sizes.
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
/// Timed runs per case; the median is reported.
const RUNS: usize = 3;
/// The disc's radius.
const RADIUS: f64 = 0.75;
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

/// Even-odd membership in the ring.
fn inside(ring: &[Point2], p: Point2) -> bool {
    let mut odd = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if (a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x {
            odd = !odd;
        }
    }
    odd
}

/// Every vertex of `out` lies inside the room and on the erosion's side.
fn validate(ring: &[Point2], out: &Region, outer: bool) {
    assert!(!out.is_empty(), "an empty erosion");
    let deviation = out.bound().expect("a disc erosion's bound").deviation;
    for polygon in out.polygons() {
        for r in std::iter::once(&polygon.outer).chain(&polygon.holes) {
            for &v in &r.points {
                let d = (0..ring.len())
                    .map(|i| to_segment(v, ring[i], ring[(i + 1) % ring.len()]))
                    .fold(f64::INFINITY, f64::min);
                assert!(inside(ring, v), "vertex {v:?} outside the room");
                if outer {
                    assert!(d <= RADIUS + deviation + SLACK, "vertex {v:?} at {d}");
                } else {
                    assert!(d >= RADIUS - SLACK, "vertex {v:?} at {d}");
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

/// One case in this process: prints the median time and peak memory.
fn case(spec: &str) {
    let parts: Vec<&str> = spec.split(':').collect();
    let corners: usize = parts[2].parse().expect("a corner count");
    let (teeth, depth) = match parts[0] {
        "plain" => (0, 0.0),
        "flat" => ((corners - 2) / 4, 0.0),
        _ => ((corners - 2) / 4, 5.0e-4),
    };
    let outer = parts[1] == "outer";
    let tolerance = Tolerance::new(1.0e-9, 1.0e-9).expect("a tolerance");
    let ring = room(teeth, depth);
    let region = Region::new(
        vec![Polygon {
            outer: Ring {
                points: ring.clone(),
            },
            holes: Vec::new(),
        }],
        tolerance,
    )
    .expect("a valid room");
    let run = || {
        if outer {
            region.erode_outer(RADIUS, tolerance)
        } else {
            region.erode_inner(RADIUS, tolerance)
        }
        .expect("the erosion")
    };
    let out = run();
    validate(&ring, &out, outer);
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
        "{:<6} {:<6} {:>7} {:>7} {:>9.3} {:>8} {:016x}",
        parts[0],
        parts[1],
        ring.len(),
        out_corners,
        times[RUNS / 2],
        peak,
        fingerprint(&out)
    );
}

fn main() {
    if let Ok(spec) = std::env::var("DISC_EROSION_CASE") {
        case(&spec);
        return;
    }
    let corners: Vec<usize> = std::env::var("DISC_EROSION_CORNERS")
        .ok()
        .map(|list| {
            list.split(',')
                .map(|n| n.trim().parse().expect("a corner count"))
                .collect()
        })
        .unwrap_or_else(|| CORNERS.to_vec());
    println!(
        "{:<6} {:<6} {:>7} {:>7} {:>9} {:>8} {:>16}",
        "walls", "side", "corners", "result", "median_s", "peak_mb", "hash"
    );
    let exe = std::env::current_exe().expect("this benchmark");
    let mut cases = Vec::new();
    for side in ["inner", "outer"] {
        cases.push(format!("plain:{side}:4"));
        for walls in ["teeth", "flat"] {
            for &n in &corners {
                cases.push(format!("{walls}:{side}:{n}"));
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
