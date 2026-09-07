//! Deterministic regression benchmarks, measured in instructions rather than seconds.
//!
//! # Why this file exists alongside `micro.rs`
//!
//! Wall-clock timing cannot gate CI. A shared runner's frequency scaling, noisy
//! neighbours, and thermal state move elapsed time by more than most real
//! regressions, so a threshold loose enough to avoid false alarms is too loose
//! to catch anything. Teams that gate on wall-clock end up ignoring the gate.
//!
//! Callgrind counts executed instructions under emulation. The number does not
//! depend on machine load, and for a deterministic input it is reproducible
//! run-to-run, so a change in it is a real change in work done. That is what
//! makes an automated comparison trustworthy.
//!
//! The trade is that instruction count is not time: it ignores cache behaviour,
//! branch prediction, and memory latency, so an optimisation that improves
//! locality without removing instructions will not show up here. The two files
//! answer different questions and neither replaces the other.
//!
//! # Requires valgrind
//!
//! `cargo bench --bench regression` needs `valgrind` on PATH. It is therefore
//! not part of the default gate; see `scripts/bench-regression.sh`, which skips
//! with a clear message rather than failing when valgrind is absent.

use axiolid_benchmark::workload::{self, Scale};
use axiolid_core::{Point3, Tolerance};
use axiolid_measure::volume_properties;
use axiolid_mesh::try_audit_mesh;
use axiolid_predicates::orient3d;
use axiolid_spatial::PointIndex;
use iai_callgrind::{
    library_benchmark, library_benchmark_group, main, Callgrind, EventKind, LibraryBenchmarkConfig,
};
use std::hint::black_box;

/// Tolerance used by every measured case.
///
/// Fixed rather than defaulted: the kernel has no `Default` for `Tolerance`
/// precisely so a measurement cannot silently depend on an unstated value.
fn tolerance() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).expect("positive tolerances are valid")
}

// Volume of a box: the smallest end-to-end closed-shell measurement.
#[library_benchmark]
fn volume_box() -> f64 {
    let workload = workload::box_mesh(2.0, 3.0, 5.0);
    let properties = volume_properties(black_box(&workload.mesh), tolerance())
        .expect("a constructed box is a closed two-manifold");
    black_box(properties.signed_volume)
}

// Volume of a mid-size sphere: the audit-and-reduce path over many triangles.
#[library_benchmark]
fn volume_sphere() -> f64 {
    let workload = workload::sphere_mesh(1.0, Scale::Medium.count());
    let properties = volume_properties(black_box(&workload.mesh), tolerance())
        .expect("the sphere construction is a closed two-manifold");
    black_box(properties.signed_volume)
}

// The predicate's fast path: the floating-point filter certifies immediately.
#[library_benchmark]
fn orient3d_filtered() -> i32 {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let c = Point3::new(0.0, 1.0, 0.0);
    let d = Point3::new(0.5, 0.5, 1.0);
    black_box(format!("{:?}", orient3d(a, b, c, d)).len() as i32)
}

// The predicate's escalation path: exactly coplanar, so exact arithmetic runs.
//
// Kept as a separate case because the gap between this and the filtered case
// is the cost of certification — the quantity any future optimisation of the
// predicates has to move.
#[library_benchmark]
fn orient3d_exact() -> i32 {
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let c = Point3::new(0.0, 1.0, 0.0);
    let d = Point3::new(0.5, 0.5, 0.0);
    black_box(format!("{:?}", orient3d(a, b, c, d)).len() as i32)
}

// Spatial index construction over a deterministic point set.
#[library_benchmark]
fn point_index_build() -> usize {
    let points = workload::point_cloud(0x5EED, Scale::Medium.count(), 100.0);
    let index = PointIndex::build(black_box(&points));
    black_box(index.len())
}

// The mesh audit's edge pass. Profiling volume_sphere showed 68.8% of its
// instructions inside one sort_unstable_by_key over per-triangle edge
// records, so the audit -- not the arithmetic -- is the real cost centre.
// Benchmarked directly to make that hot path visible on its own.
#[library_benchmark]
fn mesh_audit() -> usize {
    let sphere = workload::sphere_mesh(1.0, Scale::Medium.count());
    let health = try_audit_mesh(black_box(&sphere.mesh), tolerance())
        .expect("the medium sphere fits the audit scratch budget");
    black_box(health.boundary_edges)
}

library_benchmark_group!(
    name = regression;
    // Fail the run when a benchmark executes measurably more
    // instructions than the stored baseline. Without a threshold the
    // harness only PRINTS a delta, which nobody reads in CI output.
    //
    // 5% is chosen to sit far above the measured noise floor -- repeated
    // runs are bit-identical, so any movement at all is real work -- and
    // below a change worth investigating. Tighten once a baseline is
    // committed and the true variance across CI runners is known.
    config = LibraryBenchmarkConfig::default()
        .tool(Callgrind::default().soft_limits([(EventKind::Ir, 5.0)]));
    benchmarks =
        volume_box,
        volume_sphere,
        orient3d_filtered,
        orient3d_exact,
        point_index_build,
        mesh_audit
);

main!(library_benchmark_groups = regression);
