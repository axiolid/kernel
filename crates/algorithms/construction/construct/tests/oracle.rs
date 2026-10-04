//! Differential gate: the adopted triangulator is audited, not trusted.
//!
//! ADR 0012 says the scalar reference exists to validate other implementations.
//! `axiolid_reference::triangulate_simple` is certified (its orientation decisions go
//! through exact predicates) but has no hole support. That is exactly enough to
//! audit the profile triangulator on the hole-free case, where a silent regression
//! would otherwise go unnoticed.

use axiolid_construct::profile::{triangulate, Rings};
use axiolid_core::Point2;
use axiolid_reference::{signed_area2, triangulate_simple};

fn area_of(points: &[Point2], tris: &[[u32; 3]]) -> f64 {
    tris.iter()
        .map(|t| {
            let (a, b, c) = (
                points[t[0] as usize],
                points[t[1] as usize],
                points[t[2] as usize],
            );
            ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) / 2.0
        })
        .sum::<f64>()
        .abs()
}

/// Deterministic pseudo-random convex-ish rings, so the comparison covers more
/// than the shapes I happened to think of.
fn wobbly_ring(seed: u64, vertices: usize) -> Vec<Point2> {
    let mut state = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) as f64) / (u32::MAX as f64)
    };
    (0..vertices)
        .map(|i| {
            let t = (i as f64) * core::f64::consts::TAU / (vertices as f64);
            let r = 1.0 + next() * 0.6;
            Point2::new(r * t.cos(), r * t.sin())
        })
        .collect()
}

/// Both triangulators must cover the same area on hole-free polygons.
///
/// This is the check that would catch a triangulator regression, a wrong
/// hole-index convention, or an orientation flip.
#[test]
fn the_triangulator_and_the_certified_oracle_agree_on_covered_area() {
    let mut compared = 0;
    for seed in 0..200u64 {
        for vertices in [3usize, 5, 8, 17] {
            let ring = wobbly_ring(seed, vertices);
            if signed_area2(&ring) <= 0.0 {
                continue;
            }
            let Ok(oracle) = triangulate_simple(&ring) else {
                continue;
            };
            let rings = Rings {
                outer: ring.clone(),
                holes: Vec::new(),
            };
            let (pts, adopted) = triangulate(&rings).expect("triangulates");

            let a = area_of(&ring, &oracle);
            let b = area_of(&pts, &adopted);
            assert!(
                (a - b).abs() < 1e-9 * a.max(1.0),
                "seed {seed}/{vertices}: oracle {a} vs triangulator {b}"
            );
            assert_eq!(oracle.len(), vertices - 2, "a simple polygon has n-2 ears");
            assert_eq!(adopted.len(), vertices - 2);
            compared += 1;
        }
    }
    assert!(compared > 100, "only {compared} polygons compared");
}
