//! Wall-clock time and peak memory of `union_soup` and `overlay` on a soup
//! of grown triangles whose union is a single polygon (#338).
//!
//! The scene is axioval/engine's, built from coordinates only:
//!
//! 1. **Strips.** For `i in 0..n`, `t_i = 2 pi i / n + 0.1` and
//!    `y_i = y0 + R cos(t_i)`. Each strip gives two counter-clockwise
//!    triangles, `[(0, y_i), (0, y_i+1), (6, y_i+1)]` and
//!    `[(0, y_i), (6, y_i+1), (6, y_i)]`. Clip: `y0 = 0.01, R = 0.12`.
//!    Subject: `y0 = 0, R = 0.1`.
//! 2. **Grow** each clip triangle's corners by a 16-gon: `inscribed`
//!    takes the corners at angles `k pi / 8` and radius 0.049,
//!    `circumscribed` at `pi / 16 + k pi / 8` and radius
//!    `0.051 / cos(pi / 16)`.
//! 3. **Hull** the grown corners (monotone chain, collinear points
//!    dropped) and fan it from its first vertex `h0`.
//!
//! Operations: `soup` is `union_soup` of the grown clip triangles, `inter`
//! `overlay` of the subject strips with them (`Intersection`, `NonZero`),
//! `self` `overlay` of the grown triangles with themselves (`Union`,
//! `NonZero`). Tolerance `1e-9` throughout.
//!
//! With circumscribed growth the hulls' flat sides line up across strips,
//! and before #338 the arrangement grew far beyond the output: 64 facets
//! took 7.2 s and 0.98 GB in `union_soup`, and 128 facets failed to
//! allocate 2.8 GB, against 0.39 s and 53 MB inscribed.
//!
//! Peak memory is the child process's resident high-water mark (`VmHWM`),
//! so each case runs in a fresh process: this binary re-runs itself with
//! the case in `GROWN_SOUP_CASE`. The hash is over every output
//! coordinate's bits and shows a result unchanged across a change. Set
//! `GROWN_SOUP_FACETS` (for example `64,128`), `GROWN_SOUP_OPS` (for
//! example `soup,inter`), `GROWN_SOUP_GROWTH` (`inscribed`,
//! `circumscribed` or both) and `GROWN_SOUP_RUNS` to run fewer cases.
//! Memory-heavy sizes are best run under `ulimit -v`.
//!
//! `GROWN_SOUP_BUDGET` (bytes) runs every case through `union_soup_within`
//! and `overlay_within` under that memory budget; a refused case prints
//! `refused` with the bytes charged when the budget was passed, its time
//! and its peak, which stays within the budget.
//!
//! ```bash
//! cargo bench -p axiolid-overlay --bench grown_soup
//! ```

use std::hint::black_box;
use std::process::Command;
use std::time::Instant;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::{
    overlay_within, union_soup_within, FillRule, OverlayInput, OverlayOperation, OverlayOptions,
    OverlayRefusal, Polygon, Ring,
};

/// Facet counts measured.
const FACETS: [usize; 4] = [32, 64, 128, 256];
/// Timed runs per case; the median is reported.
const RUNS: usize = 3;

fn ccw(mut p: Vec<(f64, f64)>) -> Polygon {
    let a = (p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[1].1 - p[0].1) * (p[2].0 - p[0].0);
    if a < 0.0 {
        p.reverse();
    }
    Polygon {
        outer: Ring {
            points: p.into_iter().map(|(x, y)| Point2::new(x, y)).collect(),
        },
        holes: vec![],
    }
}

/// The plan shadows of a bar of `n` facets along x (0..6 m), section
/// radius `radius` about `y = y0`: each facet's two triangles.
fn strips(n: usize, y0: f64, radius: f64) -> Vec<[(f64, f64); 3]> {
    let y = |i: usize| {
        let t = std::f64::consts::TAU * (i % n) as f64 / n as f64 + 0.1;
        y0 + radius * t.cos()
    };
    let mut out = Vec::new();
    for i in 0..n {
        let (a, b) = (y(i), y(i + 1));
        if a == b {
            continue;
        }
        out.push([(0.0, a), (0.0, b), (6.0, b)]);
        out.push([(0.0, a), (6.0, b), (6.0, a)]);
    }
    out
}

/// The convex hull, counter-clockwise, collinear points dropped.
fn hull(mut pts: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    pts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    pts.dedup();
    let mut h: Vec<(f64, f64)> = Vec::new();
    for pass in 0..2 {
        let floor = h.len() + 1;
        let it: Vec<_> = if pass == 1 {
            pts.iter().rev().skip(1).copied().collect()
        } else {
            pts.clone()
        };
        for p in it {
            while h.len() >= floor.max(2) && cross(h[h.len() - 2], h[h.len() - 1], p) <= 0.0 {
                h.pop();
            }
            h.push(p);
        }
    }
    h.pop();
    h
}

/// The clip strips of `n` facets, each corner grown by the 16-gon, hulled
/// and fanned from the hull's first vertex.
fn grown(n: usize, circumscribed: bool) -> Vec<Polygon> {
    let step = std::f64::consts::PI / 8.0;
    let (offset, r) = if circumscribed {
        (step / 2.0, 0.051 / (step / 2.0).cos())
    } else {
        (0.0, 0.049)
    };
    let disc: Vec<(f64, f64)> = (0..16)
        .map(|k| {
            let a = offset + step * k as f64;
            (r * a.cos(), r * a.sin())
        })
        .collect();
    let mut out = Vec::new();
    for t in strips(n, 0.01, 0.12) {
        let h = hull(
            t.iter()
                .flat_map(|p| disc.iter().map(move |d| (p.0 + d.0, p.1 + d.1)))
                .collect(),
        );
        for w in h[1..].windows(2) {
            out.push(ccw(vec![h[0], w[0], w[1]]));
        }
    }
    out
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

fn fingerprint(polygons: &[Polygon]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    let mut eat = |v: u64| {
        h ^= v;
        h = h.wrapping_mul(0x0100_0000_01b3);
    };
    for p in polygons {
        for r in std::iter::once(&p.outer).chain(&p.holes) {
            eat(r.points.len() as u64);
            for q in &r.points {
                eat(q.x.to_bits());
                eat(q.y.to_bits());
            }
        }
    }
    h
}

/// One case in this process, given as `op:growth:facets`: prints the
/// median time, peak memory and hash.
fn case(spec: &str) {
    let parts: Vec<&str> = spec.split(':').collect();
    let circumscribed = parts[1] == "circumscribed";
    let n: usize = parts[2].parse().expect("a facet count");
    let runs: usize = std::env::var("GROWN_SOUP_RUNS")
        .ok()
        .and_then(|r| r.parse().ok())
        .unwrap_or(RUNS);
    let tolerance = Tolerance::new(1.0e-9, 1.0e-9).expect("a tolerance");
    let clip = grown(n, circumscribed);
    let subject: Vec<Polygon> = strips(n, 0.0, 0.1)
        .into_iter()
        .map(|t| ccw(t.to_vec()))
        .collect();
    let rings: Vec<Ring> = clip.iter().map(|p| p.outer.clone()).collect();
    let (clip, subject) = (input(clip), input(subject));
    let mut execution = ExecutionOptions::new(tolerance);
    if let Some(bytes) = std::env::var("GROWN_SOUP_BUDGET")
        .ok()
        .and_then(|b| b.parse().ok())
    {
        execution = execution.with_memory_budget(bytes);
    }
    let options = OverlayOptions::new();
    let attempt = || match parts[0] {
        "soup" => union_soup_within(&rings, options, &execution).map(|r| r.polygons),
        "inter" => overlay_within(
            &subject,
            &clip,
            OverlayOperation::Intersection,
            FillRule::NonZero,
            options,
            &execution,
        )
        .map(|r| r.result.polygons),
        _ => overlay_within(
            &clip,
            &clip,
            OverlayOperation::Union,
            FillRule::NonZero,
            options,
            &execution,
        )
        .map(|r| r.result.polygons),
    };
    let start = Instant::now();
    let out = match attempt() {
        Ok(out) => out,
        Err(OverlayRefusal::BudgetExceeded { needed, .. }) => {
            println!(
                "{:<6} {:<14} {:>6} {:>9} refused at {} bytes charged, {:.3} s, {} MB",
                parts[0],
                parts[1],
                n,
                rings.len(),
                needed,
                start.elapsed().as_secs_f64(),
                peak_mb()
            );
            return;
        }
        Err(other) => panic!("{spec}: {other:?}"),
    };
    let run = || attempt().expect("admitted once, admitted again");
    let corners: usize = out
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
    let peak = peak_mb();
    println!(
        "{:<6} {:<14} {:>6} {:>9} {:>8} {:>8} {:>9.3} {:>8} {:016x}",
        parts[0],
        parts[1],
        n,
        rings.len(),
        out.len(),
        corners,
        times[runs.saturating_sub(1) / 2],
        peak,
        fingerprint(&out)
    );
}

/// This process's resident high-water mark in MB.
fn peak_mb() -> String {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status
                .lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|kb| kb.parse::<f64>().ok())
        })
        .map_or_else(|| "n/a".to_string(), |kb| format!("{:.0}", kb / 1024.0))
}

fn main() {
    if let Ok(spec) = std::env::var("GROWN_SOUP_CASE") {
        case(&spec);
        return;
    }
    let list = |name: &str, default: &str| -> Vec<String> {
        std::env::var(name)
            .unwrap_or_else(|_| default.to_string())
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    };
    let facets = list(
        "GROWN_SOUP_FACETS",
        &FACETS.map(|n| n.to_string()).join(","),
    );
    let ops = list("GROWN_SOUP_OPS", "soup,inter,self");
    let growths = list("GROWN_SOUP_GROWTH", "inscribed,circumscribed");
    println!(
        "{:<6} {:<14} {:>6} {:>9} {:>8} {:>8} {:>9} {:>8} {:>16}",
        "op", "growth", "facets", "triangles", "polygons", "corners", "median_s", "peak_mb", "hash"
    );
    let exe = std::env::current_exe().expect("this benchmark");
    for op in &ops {
        for growth in &growths {
            for n in &facets {
                let spec = format!("{op}:{growth}:{n}");
                let status = Command::new(&exe)
                    .env("GROWN_SOUP_CASE", &spec)
                    .status()
                    .expect("a child run");
                if !status.success() {
                    println!("{op:<6} {growth:<14} {n:>6} failed: {status}");
                }
            }
        }
    }
}
