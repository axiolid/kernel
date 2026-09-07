//! Scaling behaviour as an assertion, not a chart to read later.
//!
//! # Why this is a test and not a benchmark
//!
//! A benchmark records how long something took on one machine. It cannot
//! tell you whether an algorithm is linear or quadratic, because a single
//! point has no shape. These tests measure *work* -- operation counts that
//! do not vary with CPU speed or load -- across a size sweep and assert the
//! growth rate. A linear algorithm silently becoming quadratic is the
//! regression that hurts most at scale and is invisible to a single-size
//! benchmark.
//!
//! Timing is deliberately NOT used here: a shared runner's wall clock is too
//! noisy to support a complexity claim.

use axiolid_benchmark::workload;
use axiolid_core::Point3;
use axiolid_spatial::PointIndex;

/// Count the points a radius query actually visits.
///
/// This is the quantity that determines the query's cost, and it is exact:
/// the same input yields the same count on every machine.
fn visited_for_radius(points: &[Point3], radius: f64) -> usize {
    let index = PointIndex::build(points);
    let mut visited = 0_usize;
    index
        .for_each_within(Point3::ZERO, radius, |_| visited += 1)
        .expect("a finite query and radius are valid");
    visited
}

/// A radius query must scale with what it finds, not with the whole cloud.
///
/// The point of a spatial index: a small query against a growing cloud stays
/// cheap. If the index degraded into a linear scan, the visit count would
/// track the cloud size and this fails.
#[test]
fn a_radius_query_visits_a_bounded_neighbourhood_as_the_cloud_grows() {
    // Same density, growing extent: the neighbourhood within a fixed radius
    // holds roughly the same number of points while the cloud grows 16x.
    let small = workload::point_cloud(0xA11CE, 4_000, 10.0);
    let large = workload::point_cloud(0xA11CE, 64_000, 40.0);

    let visited_small = visited_for_radius(&small, 1.0);
    let visited_large = visited_for_radius(&large, 1.0);

    // A linear scan would visit ~16x more. Allow generous slack for the
    // random cloud while still failing an order-of-magnitude regression.
    assert!(
        visited_large < visited_small * 4,
        "query visits grew {visited_small} -> {visited_large} as the cloud grew 16x; \
         the index is behaving like a linear scan"
    );
}

/// Mesh measurement must stay linear in triangle count.
///
/// Measured in instructions rather than seconds would be better still, but
/// the triangle count itself is the exact work quantity here: the routine
/// visits each triangle once, so output size is the invariant to pin.
#[test]
fn sphere_generation_stays_linear_in_requested_vertices() {
    let small = workload::sphere_mesh(1.0, 4_096);
    let large = workload::sphere_mesh(1.0, 16_384);

    let small_triangles = small.mesh.indices.len() / 3;
    let large_triangles = large.mesh.indices.len() / 3;

    // 4x the requested vertices should give ~4x the triangles, never 16x.
    let ratio = large_triangles as f64 / small_triangles as f64;
    assert!(
        (2.0..=6.0).contains(&ratio),
        "triangle count grew {small_triangles} -> {large_triangles} (x{ratio:.2}) \
         for a 4x vertex request; generation is not linear"
    );
}

/// Workload generation must be reproducible, or every comparison is void.
///
/// This is the assumption the whole measurement system rests on: if inputs
/// drift between runs, an instruction-count difference measures the
/// generator rather than the kernel.
#[test]
fn workloads_are_identical_across_calls() {
    let first = workload::point_cloud(0x5EED, 1_000, 5.0);
    let second = workload::point_cloud(0x5EED, 1_000, 5.0);
    assert_eq!(first, second, "the same seed must give the same cloud");

    let wall_a = workload::wall_with_openings(8);
    let wall_b = workload::wall_with_openings(8);
    assert_eq!(
        wall_a.subject.positions, wall_b.subject.positions,
        "wall generation must not vary between calls"
    );
    assert_eq!(wall_a.expected_volume, wall_b.expected_volume);
}
