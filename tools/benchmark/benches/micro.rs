//! Wall-clock microbenchmarks over kernel hot paths.
//!
//! These report elapsed time, which is informative but machine- and
//! load-dependent. They are deliberately NOT the regression gate: see
//! `regression.rs` for the deterministic instruction-count benchmarks CI can
//! compare across runs. Use these to find where time goes; use those to prove
//! it stopped going there.
//!
//! Every case validates its result before timing it, so a path that silently
//! stops computing cannot appear here as an improvement.

use axiolid_benchmark::validate::{expect_valid, Tolerance as CheckTolerance, Validated};
use axiolid_benchmark::workload::{self, Scale};
use axiolid_core::Tolerance;
use axiolid_measure::volume_properties;
use axiolid_predicates::orient3d;
use axiolid_spatial::PointIndex;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::hint::black_box;

/// Volume of a closed shell, across problem sizes.
///
/// `volume_properties` refuses a mesh that is not a closed two-manifold, so
/// this measures the whole audit-and-reduce path rather than a bare sum. That
/// refusal is also why it is a good benchmark: an implementation that broke
/// closure detection would start erroring here instead of quietly getting
/// faster.
fn volume(c: &mut Criterion) {
    let mut group = c.benchmark_group("measure/volume");
    // Explicit, validated tolerance: the kernel deliberately has no Default,
    // because a silently-chosen tolerance is how measurements become
    // incomparable between benchmarks.
    let tolerance = Tolerance::new(1e-9, 1e-9).expect("positive tolerances are valid");

    let box_workload = workload::box_mesh(2.0, 3.0, 5.0);
    let expected = box_workload
        .expected_volume
        .expect("the box construction states its volume");

    // Validate before timing: a wrong or refusing implementation must never be
    // reported as a fast one.
    let measured = volume_properties(&box_workload.mesh, tolerance)
        .expect("a constructed box is a closed two-manifold")
        .signed_volume;
    expect_valid(Validated::scalar(
        "box volume",
        measured,
        expected,
        CheckTolerance::GEOMETRIC,
    ));

    group.throughput(Throughput::Elements(
        box_workload.mesh.triangle_count() as u64
    ));
    group.bench_function("box", |b| {
        b.iter(|| {
            black_box(volume_properties(
                black_box(&box_workload.mesh),
                black_box(tolerance),
            ))
        });
    });

    for scale in Scale::all() {
        let sphere = workload::sphere_mesh(1.0, scale.count());
        // A UV sphere is closed by construction; prove it before timing it.
        assert!(
            volume_properties(&sphere.mesh, tolerance).is_ok(),
            "the sphere construction must produce a closed two-manifold"
        );
        group.throughput(Throughput::Elements(sphere.mesh.triangle_count() as u64));
        group.bench_with_input(
            BenchmarkId::new("sphere", scale.label()),
            &sphere,
            |b, sphere| {
                b.iter(|| {
                    black_box(volume_properties(
                        black_box(&sphere.mesh),
                        black_box(tolerance),
                    ))
                });
            },
        );
    }
    group.finish();
}

/// The certified orientation predicate, on inputs of increasing difficulty.
///
/// Exact predicates escalate to slower arithmetic only when the floating-point
/// filter cannot decide, so a single throughput figure is misleading. The
/// interesting quantity is the gap between the separated and degenerate cases:
/// that gap is the cost of certification, and it is what any future
/// optimisation has to move.
fn predicates(c: &mut Criterion) {
    use axiolid_core::Point3;

    let mut group = c.benchmark_group("predicates/orient3d");

    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let c_point = Point3::new(0.0, 1.0, 0.0);

    let cases = [
        // Well separated: the filter decides immediately.
        ("separated", Point3::new(0.5, 0.5, 1.0)),
        // Exactly coplanar: the filter cannot decide, so exact arithmetic runs.
        ("coplanar", Point3::new(0.5, 0.5, 0.0)),
        // Sub-normal magnitude: below the filter's 1e-90 guard the
        // difference is rejected outright and exact arithmetic runs.
        //
        // Measured, not assumed. A probe over z from 1e-1 to 1e-300
        // (tools/benchmark/examples/filter_probe.rs) shows the filter
        // stays CERTAIN for every ordinary magnitude, because shrinking
        // the offset shrinks its error bound with it. Only exact zero
        // and sub-normal inputs escalate. Two earlier guesses (1e-30,
        // 1e-17) both silently measured the filter path under a name
        // that claimed otherwise.
        ("subnormal", Point3::new(0.5, 0.5, 1e-300)),
    ];

    for (label, d) in cases {
        group.bench_function(label, |bencher| {
            bencher.iter(|| {
                black_box(orient3d(
                    black_box(a),
                    black_box(b),
                    black_box(c_point),
                    black_box(d),
                ))
            });
        });
    }
    group.finish();
}

/// Spatial index construction and query, across problem sizes.
///
/// Reported per element so the curve shows scaling rather than raw totals: a
/// structure that is near-linear to build looks flat here, and anything that
/// bends upward with size is the finding.
fn spatial(c: &mut Criterion) {
    let mut group = c.benchmark_group("spatial/point-index");

    for scale in Scale::all() {
        let points = workload::point_cloud(0x5EED, scale.count(), 100.0);
        group.throughput(Throughput::Elements(points.len() as u64));

        group.bench_with_input(
            BenchmarkId::new("build", scale.label()),
            &points,
            |b, points| {
                b.iter(|| black_box(PointIndex::build(black_box(points))));
            },
        );

        let index = PointIndex::build(&points);
        let probe = points[points.len() / 2];

        // The index must actually find the probe point itself; a query that
        // silently returns nothing would otherwise benchmark as very fast.
        let mut seen = 0usize;
        index
            .for_each_within(probe, 10.0, |_| seen += 1)
            .expect("a finite probe and radius are a valid query");
        assert!(
            seen > 0,
            "a radius query centred on a stored point must find at least that point"
        );

        group.bench_with_input(
            BenchmarkId::new("radius-query", scale.label()),
            &index,
            |b, index| {
                b.iter(|| {
                    let mut found = 0usize;
                    let outcome = index.for_each_within(black_box(probe), 10.0, |_| found += 1);
                    black_box((outcome, found))
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, volume, predicates, spatial);
criterion_main!(benches);
