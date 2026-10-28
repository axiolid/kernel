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
use axiolid_hidden_lines::{ScalarHiddenLines, ScanHiddenLines};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineLimits, HiddenLineMesh, HiddenLineRequest, MeshHiddenLines, MeshId,
    MeshRole, OrthographicView, Visibility,
};
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

/// A closed, outward-wound box.
fn hidden_line_box(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let (x0, y0, z0, x1, y1, z1) = (min[0], min[1], min[2], max[0], max[1], max[2]);
    TriMesh::new(
        vec![
            Point3::new(x0, y0, z0),
            Point3::new(x1, y0, z0),
            Point3::new(x1, y1, z0),
            Point3::new(x0, y1, z0),
            Point3::new(x0, y0, z1),
            Point3::new(x1, y0, z1),
            Point3::new(x1, y1, z1),
            Point3::new(x0, y1, z1),
        ],
        vec![
            0, 3, 2, 0, 2, 1, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ],
    )
}

/// Hidden-line drawing of a 20 x 20 and a 50 x 50 grid of boxes (#138).
///
/// A building-like layout: every cell holds a column of height 1 to 3 with a
/// gap to its neighbours. In the isometric view every box overlaps its
/// neighbours in the image; the cut plan sets the near plane at mid-height
/// (2.0), so the taller columns are capped and the lower ones drawn whole.
/// Each case runs with the BVH broad phase and with the linear scan
/// (`ScanHiddenLines`, identical output), so the broad phase's gain is read
/// off the same group; the provider CHANGELOG records a run.
fn hidden_lines(criterion: &mut Criterion) {
    let bvh = ScalarHiddenLines::new();
    let providers: [(&str, &dyn MeshHiddenLines); 2] = [("bvh", &bvh), ("scan", &ScanHiddenLines)];
    let options = ExecutionOptions::new(tolerance());
    let (s2, s6, s3) = (
        std::f64::consts::FRAC_1_SQRT_2,
        1.0 / 6.0_f64.sqrt(),
        1.0 / 3.0_f64.sqrt(),
    );
    let isometric = OrthographicView::new(Frame3 {
        origin: Point3::ZERO,
        x: Vec3::new(s2, -s2, 0.0),
        y: Vec3::new(s6, s6, -2.0 * s6),
        z: Vec3::new(s3, s3, s3),
    });
    let plan = OrthographicView::new(Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    })
    .with_near(2.0);
    let mut group = criterion.benchmark_group("scenario/hidden-lines");
    group.sample_size(10);

    for n in [20_usize, 50] {
        let meshes: Vec<TriMesh> = (0..n * n)
            .map(|k| {
                let (i, j) = ((k % n) as f64, (k / n) as f64);
                let height = 1.0 + ((k * 7 + (k / n) * 3) % 5) as f64 * 0.5;
                hidden_line_box([i, j, 0.0], [i + 0.75, j + 0.75, height])
            })
            .collect();
        let entries: Vec<HiddenLineMesh<'_>> = meshes
            .iter()
            .enumerate()
            .map(|(k, mesh)| HiddenLineMesh::new(MeshId(k as u64 + 1), MeshRole::Solid, mesh))
            .collect();
        let limits =
            HiddenLineLimits::new(n * n, 8 * n * n, 12 * n * n, 4_000_000_000, 100_000_000);
        for (name, view) in [("isometric", isometric), ("cut-plan", plan)] {
            let request =
                HiddenLineRequest::new(&entries, view, EdgeSelection::new(Some(0.5)), limits);
            // A drawing with nothing visible, or nothing hidden, would time
            // beautifully; require both before timing.
            let outcome = bvh
                .hidden_lines(&request, &options)
                .expect("the grid draws");
            for visibility in [Visibility::Visible, Visibility::Hidden] {
                assert!(
                    outcome.segments.iter().any(|s| s.visibility == visibility),
                    "{name} {n}x{n}: no {visibility:?} segment"
                );
            }
            // The scan is the baseline only if it draws the same thing.
            let scanned = ScanHiddenLines
                .hidden_lines(&request, &options)
                .expect("the grid draws");
            assert_eq!(
                format!("{:?}", outcome.segments),
                format!("{:?}", scanned.segments),
                "{name} {n}x{n}: broad phases disagree"
            );
            group.throughput(Throughput::Elements(outcome.evidence.selected_edges as u64));
            for (phase, provider) in providers {
                group.bench_function(format!("{name}/{n}x{n}/{phase}"), |bencher| {
                    bencher.iter(|| {
                        black_box(
                            provider
                                .hidden_lines(black_box(&request), &options)
                                .expect("the grid draws"),
                        )
                    });
                });
            }
        }
    }
    group.finish();
}

criterion_group!(scenarios, wall_subtraction, mesh_section, hidden_lines);
criterion_main!(scenarios);
