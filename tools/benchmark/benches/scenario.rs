//! End-to-end scenarios through the public facade.
//!
//! # Why these differ from `micro.rs`
//!
//! A microbenchmark calls one function. These drive `axiolid::application`,
//! the same entry point an integrator uses, so the number includes dispatch,
//! provider selection, validation, and the operation itself. That is the
//! number a user actually experiences, and it is where a regression in the
//! seams -- not the algorithms -- would show up.
//!
//! Every scenario validates its result against a derived expectation before
//! it is timed. A subtraction that silently drops its cutters would be fast
//! and wrong; the harness refuses to report that as a win.

use axiolid::application::Application;
use axiolid::contracts::ExecutionOptions;
use axiolid_benchmark::workload;
use axiolid_core::{Frame3, Point3, Tolerance, Vec3};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::hint::black_box;

/// Tolerance used across every scenario, so runs stay comparable.
fn tolerance() -> Tolerance {
    Tolerance::new(1.0e-9, 1.0e-9).expect("a positive tolerance is valid")
}

/// Subtracting many openings from a wall, through the facade.
///
/// Scaling is the point: the row count sweeps opening count so the shape of
/// the curve is visible, not just one point on it.
fn wall_subtraction(criterion: &mut Criterion) {
    let application = Application::portable().expect("the portable application builds");
    let options = ExecutionOptions::new(tolerance());
    let mut group = criterion.benchmark_group("scenario/wall-subtraction");

    for openings in [1_usize, 4, 16, 64] {
        let wall = workload::wall_with_openings(openings);

        // Validate BEFORE timing: a run that drops cutters must not be
        // reported as a fast run.
        let outcome = application
            .subtract_many(&wall.subject, &wall.tools, &options)
            .expect("the portable provider handles disjoint cutters");
        let measured = application
            .measure_mesh(&outcome.mesh, tolerance())
            .expect("the result is a closed solid");
        let volume = measured.volume.signed_volume.abs();
        assert!(
            (volume - wall.expected_volume).abs() <= 1.0e-6 * wall.expected_volume,
            "openings={openings}: volume {volume} != expected {}",
            wall.expected_volume
        );

        group.throughput(Throughput::Elements(openings as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(openings),
            &wall,
            |bencher, wall| {
                bencher.iter(|| {
                    black_box(
                        application
                            .subtract_many(
                                black_box(&wall.subject),
                                black_box(&wall.tools),
                                &options,
                            )
                            .expect("subtraction succeeds"),
                    )
                });
            },
        );
    }
    group.finish();
}

/// Sectioning a sphere with a plane, through the facade.
///
/// Sweeps mesh density so the cost curve against triangle count is visible.
fn mesh_section(criterion: &mut Criterion) {
    let application = Application::portable().expect("the portable application builds");
    let options = ExecutionOptions::new(tolerance());
    // Limits are explicit by contract: generous enough that the benchmark
    // measures sectioning rather than a refusal.
    let limits = axiolid::mesh_section::SectionLimits::new(1_000_000, 1_000_000, 1_000_000, 10_000);
    // The section plane is the frame's xy-plane, so z is its normal.
    let frame = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let mut group = criterion.benchmark_group("scenario/mesh-section");

    // `sphere_mesh` takes a target vertex count, so these span roughly
    // 500 to 30k triangles rather than a handful.
    for vertices in [256_usize, 4_096, 16_384] {
        let sphere = workload::sphere_mesh(1.0, vertices);
        let triangles = sphere.mesh.indices.len() / 3;

        // A plane through the centre of a unit sphere must produce a
        // non-empty section. An empty result would time beautifully.
        let outcome = application
            .section_mesh(&sphere.mesh, frame, limits, &options)
            .expect("a central plane sections the sphere");
        assert!(
            !outcome.contours.is_empty(),
            "vertices={vertices}: a central section must not be empty"
        );

        group.throughput(Throughput::Elements(triangles as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(triangles),
            &sphere.mesh,
            |bencher, mesh| {
                bencher.iter(|| {
                    black_box(
                        application
                            .section_mesh(black_box(mesh), frame, limits, &options)
                            .expect("sectioning succeeds"),
                    )
                });
            },
        );
    }
    group.finish();
}

criterion_group!(scenarios, wall_subtraction, mesh_section);
criterion_main!(scenarios);
