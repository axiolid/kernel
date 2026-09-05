//! The reference provider against the contract's conformance suite and
//! against real point data.
//!
//! Conformance proves the provider is well-behaved. It cannot prove the
//! surface is faithful, so the tests here additionally check reconstructed
//! geometry against shapes whose answer is known in closed form.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::audit_mesh;
use axiolid_pointcloud::PointCloud;
use axiolid_pointcloud_reconstruction_contract::{
    conformance, PointcloudReconstruction, Reconstruction, ReconstructionRefusal,
    ReconstructionRequest,
};
use axiolid_pointcloud_reconstruction_sdf::SdfReconstruction;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// Evenly distributed samples on a sphere, with outward normals.
fn sphere(count: usize, radius: Scalar) -> PointCloud {
    let mut points = Vec::with_capacity(count);
    let mut normals = Vec::with_capacity(count);
    let golden = core::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    for i in 0..count {
        let y = 1.0 - (i as Scalar / (count.max(2) - 1) as Scalar) * 2.0;
        let ring = (1.0 - y * y).max(0.0).sqrt();
        let theta = golden * i as Scalar;
        let direction = Vec3::new(theta.cos() * ring, y, theta.sin() * ring);
        points.push(Point3::ZERO + direction * radius);
        normals.push(direction);
    }
    PointCloud::new(points)
        .expect("finite")
        .with_normals(normals)
        .expect("normals match")
}

/// The suite the contract exports must pass, and must actually exercise the
/// provider rather than reaching conformance by refusing everything.
#[test]
fn the_provider_is_conformant() {
    let report = conformance::run(&SdfReconstruction::new());
    assert!(report.is_conformant(), "{report}");
    assert!(
        report.exercised() >= 4,
        "conformance must actually exercise the provider, not skip its way to a pass:\n{report}"
    );
}

/// The central geometric claim: samples on a sphere reconstruct a surface
/// whose vertices lie on that sphere. Checked against the radius directly,
/// not against a volume or bbox proxy that a wrong-but-plausible surface
/// could also satisfy.
#[test]
fn a_sampled_sphere_reconstructs_onto_the_sphere() {
    let radius = 1.0;
    let cloud = sphere(2000, radius);
    let outcome = SdfReconstruction::new()
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("reconstruction runs");

    let outcome = match outcome {
        Reconstruction::Surface(outcome) => outcome,
        Reconstruction::Refused(reason) => panic!("dense sphere refused: {reason}"),
    };

    assert!(
        outcome.evidence.output_triangles > 100,
        "expected a real surface, got {} triangles",
        outcome.evidence.output_triangles
    );

    let mut worst: Scalar = 0.0;
    for vertex in &outcome.mesh.positions {
        worst = worst.max((vertex.length() - radius).abs());
    }
    // One sample spacing of slack: the extraction grid cannot resolve finer
    // than the samples themselves.
    let allowed = outcome.evidence.sample_spacing * 1.5;
    assert!(
        worst <= allowed,
        "worst radial error {worst} exceeds {allowed} (spacing {})",
        outcome.evidence.sample_spacing
    );
}

/// A closed capture must reconstruct a closed solid, and say so.
#[test]
fn a_fully_sampled_object_reconstructs_closed() {
    let cloud = sphere(2000, 1.0);
    let outcome = SdfReconstruction::new()
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("runs");
    let outcome = outcome.outcome().expect("surface");

    let health = audit_mesh(&outcome.mesh, Tolerance::METRE);
    assert!(
        health.is_closed_two_manifold(),
        "expected a closed solid: boundary={} non_manifold={}",
        health.boundary_edges,
        health.non_manifold_edges
    );
    assert!(
        outcome.evidence.closed,
        "evidence must report the closure it achieved"
    );
}

/// Without normals there is no way to know which side is inside, so the
/// result is a shrink-wrap rather than a fit. It must still be a valid
/// solid, and the evidence must say normals were not used.
#[test]
fn a_positions_only_reconstruction_is_reported_as_such() {
    let cloud = sphere(2000, 1.0);
    let outcome = SdfReconstruction::new()
        .reconstruct(
            &cloud,
            &ReconstructionRequest::default().ignoring_normals(),
            &options(),
        )
        .expect("runs");
    let outcome = outcome.outcome().expect("surface");

    assert!(
        !outcome.evidence.used_normals,
        "a positions-only result must not claim it used normals"
    );
    assert!(
        outcome.evidence.output_triangles > 100,
        "expected a surface without normals too"
    );
}

/// Using normals must actually change the answer. Without this, the
/// `use_normals` flag could be ignored entirely and every other test would
/// still pass.
#[test]
fn normals_change_the_reconstruction() {
    let cloud = sphere(1500, 1.0);
    let provider = SdfReconstruction::new();

    let with = provider
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("runs");
    let without = provider
        .reconstruct(
            &cloud,
            &ReconstructionRequest::default().ignoring_normals(),
            &options(),
        )
        .expect("runs");

    let with = with.outcome().expect("surface");
    let without = without.outcome().expect("surface");

    assert!(
        with.mesh.positions != without.mesh.positions,
        "the normals flag must actually affect the surface"
    );

    // The normal-driven surface passes through the samples; the unsigned
    // one sits offset from them, so it must be measurably larger.
    let radius_of = |mesh: &axiolid_mesh::TriMesh| -> Scalar {
        mesh.positions
            .iter()
            .map(|p| p.length())
            .fold(0.0_f64, |a, b| a.max(b))
    };
    assert!(
        radius_of(&without.mesh) > radius_of(&with.mesh),
        "the positions-only shrink-wrap should sit outside the fitted surface"
    );
}

/// Too few points is a named refusal, never an empty mesh a caller could
/// read as "there is no object here".
#[test]
fn too_few_points_is_refused_by_name() {
    let cloud = PointCloud::new(vec![Point3::ZERO, Point3::new(1.0, 0.0, 0.0)]).expect("finite");
    let result = SdfReconstruction::new()
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("runs");
    assert!(matches!(
        result.refusal(),
        Some(ReconstructionRefusal::TooFewPoints { supplied: 2, .. })
    ));
}

/// Coplanar points bound no volume. Fitting a zero-thickness sheet and
/// calling it a solid would misrepresent the capture.
#[test]
fn a_flat_capture_is_refused_rather_than_given_zero_thickness() {
    let points: Vec<Point3> = (0..400)
        .map(|i| {
            let (x, y) = (i % 20, i / 20);
            Point3::new(x as Scalar * 0.1, y as Scalar * 0.1, 0.0)
        })
        .collect();
    let cloud = PointCloud::new(points).expect("finite");
    let result = SdfReconstruction::new()
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("runs");
    assert!(
        matches!(
            result.refusal(),
            Some(ReconstructionRefusal::DegenerateExtent { .. })
        ),
        "a planar capture must be refused, got {result:?}"
    );
}

/// Asking for detail finer than the samples resolve is refused rather than
/// silently clamped: the caller learns the data's limit instead of
/// receiving interpolation presented as measurement.
#[test]
fn resolution_finer_than_the_data_is_refused() {
    let cloud = sphere(500, 1.0);
    let result = SdfReconstruction::new()
        .reconstruct(
            &cloud,
            &ReconstructionRequest::at_edge_length(1e-4),
            &options(),
        )
        .expect("runs");
    assert!(
        matches!(
            result.refusal(),
            Some(ReconstructionRefusal::ResolutionExceedsData { .. })
        ),
        "expected a resolution refusal, got {result:?}"
    );
}

/// A provider declaring bitwise determinism must actually repeat itself.
#[test]
fn reconstruction_is_reproducible() {
    let cloud = sphere(800, 1.0);
    let provider = SdfReconstruction::new();
    let first = provider
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("runs");
    let second = provider
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("runs");
    assert_eq!(
        first.outcome().expect("surface").mesh.positions,
        second.outcome().expect("surface").mesh.positions
    );
}

/// The evidence must describe the mesh actually returned, so a caller can
/// trust the counters without re-deriving them.
#[test]
fn evidence_describes_the_mesh_returned() {
    let cloud = sphere(1200, 1.0);
    let outcome = SdfReconstruction::new()
        .reconstruct(&cloud, &ReconstructionRequest::default(), &options())
        .expect("runs");
    let outcome = outcome.outcome().expect("surface");

    assert_eq!(
        outcome.evidence.output_triangles,
        outcome.mesh.indices.len() / 3
    );
    assert_eq!(outcome.evidence.input_points, cloud.len());
    assert_eq!(
        outcome.evidence.output_components, 1,
        "a sphere reconstructs as one shell"
    );
    assert!(outcome.evidence.sample_spacing > 0.0);
    assert!(outcome.evidence.achieved_edge_length > 0.0);
}
