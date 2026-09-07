//! How the mesh audit scales with triangle count.
//!
//! The audit's edge sort dominates its cost, so this sweep is the honest way
//! to see whether an algorithmic change helps at the sizes that matter
//! rather than only at one benchmark size.

use axiolid_benchmark::workload;
use axiolid_core::Tolerance;
use axiolid_mesh::try_audit_mesh;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::hint::black_box;

fn tolerance() -> Tolerance {
    Tolerance::new(1.0e-9, 1.0e-9).expect("finite positive tolerance")
}

fn mesh_audit(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("audit/mesh");
    for &vertices in &[1_000usize, 10_000, 50_000, 200_000] {
        let sphere = workload::sphere_mesh(1.0, vertices);
        let triangles = sphere.mesh.triangle_count();
        group.throughput(Throughput::Elements(triangles as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(triangles),
            &sphere,
            |bencher, sphere| {
                bencher.iter(|| {
                    let health = try_audit_mesh(black_box(&sphere.mesh), tolerance())
                        .expect("sphere fits the audit budget");
                    black_box(health.boundary_edges)
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, mesh_audit);
criterion_main!(benches);
