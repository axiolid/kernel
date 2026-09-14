//! Refinement: exact on planes, convergent on analytic surfaces.

use axiolid_core::{Frame3, Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_refine::{refine, RefineError, RefineTarget};
use axiolid_surface::{Cylinder, Surface};

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

/// Unit cube, outward wound.
fn cube() -> TriMesh {
    let p = |x: Scalar, y: Scalar, z: Scalar| Point3::new(x, y, z);
    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 0.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.0, 0.0, 1.0),
        p(1.0, 0.0, 1.0),
        p(1.0, 1.0, 1.0),
        p(0.0, 1.0, 1.0),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // bottom
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // front
        1, 2, 6, 1, 6, 5, // right
        2, 3, 7, 2, 7, 6, // back
        3, 0, 4, 3, 4, 7, // left
    ];
    TriMesh::new(positions, indices)
}

fn volume(mesh: &TriMesh) -> Scalar {
    axiolid_measure::volume_properties(mesh, tol())
        .expect("closed two-manifold")
        .signed_volume
}

/// Uniform refinement quadruples triangles and must not move the surface.
///
/// A midpoint of a flat triangle lies in that triangle's plane, so a
/// planar refinement that changes volume has moved a vertex it had no
/// business moving. Volume is the sharpest available witness.
#[test]
fn uniform_refinement_quadruples_triangles_and_preserves_volume_exactly() {
    let cube = cube();
    let before = volume(&cube);

    let (refined, report) =
        refine(&cube, RefineTarget::Uniform { levels: 1 }, None, tol()).expect("refines");

    assert_eq!(report.output_triangles, report.input_triangles * 4);
    assert_eq!(refined.triangle_count(), 48, "12 triangles become 48");
    // Not bit-equality: refining changes the number of terms in the
    // divergence sum, so the last ulp moves even though no vertex did.
    // The bound is tight enough that an actually-moved vertex fails.
    assert!(
        (volume(&refined) - before).abs() < 1e-12,
        "planar refinement must not move the surface: {before} -> {}",
        volume(&refined)
    );
    assert_eq!(
        report.max_deviation, 0.0,
        "a linear midpoint deviates from itself by zero"
    );
}

/// Two passes quadruple twice, and shared edges are not duplicated.
///
/// A cracked mesh would still have the right triangle count, so vertex
/// count is what actually proves the midpoints were shared.
#[test]
fn refinement_shares_midpoints_between_adjacent_triangles() {
    let (refined, report) =
        refine(&cube(), RefineTarget::Uniform { levels: 2 }, None, tol()).expect("refines");

    assert_eq!(refined.triangle_count(), 12 * 16);
    // Euler: a closed surface refined this way has V - E + F = 2. Rather
    // than recompute topology, assert the mesh stays closed, which a
    // cracked mesh would fail.
    assert!(
        axiolid_measure::volume_properties(&refined, tol()).is_ok(),
        "shared midpoints keep the mesh closed and two-manifold"
    );
    assert!(report.vertices_added > 0);
}

/// Surface-aware refinement moves new vertices ONTO the analytic cylinder.
///
/// This is the capability a mesh-only kernel cannot match, so it is
/// asserted directly: every introduced vertex must sit at the cylinder's
/// radius, not on the chord between two existing vertices.
#[test]
fn surface_aware_refinement_places_vertices_on_the_real_cylinder() {
    let radius = 1.0;
    let cylinder = Surface::Cylinder(Cylinder {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    });

    // A coarse 8-sided prism approximating the cylinder.
    let sides = 8;
    let mut positions = Vec::new();
    for i in 0..sides {
        let a = std::f64::consts::TAU * (i as Scalar) / (sides as Scalar);
        positions.push(Point3::new(radius * a.cos(), radius * a.sin(), 0.0));
        positions.push(Point3::new(radius * a.cos(), radius * a.sin(), 1.0));
    }
    let mut indices = Vec::new();
    for i in 0..sides {
        let n = (i + 1) % sides;
        let (a, b, c, d) = (
            (2 * i) as u32,
            (2 * i + 1) as u32,
            (2 * n) as u32,
            (2 * n + 1) as u32,
        );
        indices.extend_from_slice(&[a, c, b, b, c, d]);
    }
    let prism = TriMesh::new(positions, indices);

    let worst_before = max_surface_error(&prism, radius);

    let (refined, report) = refine(
        &prism,
        RefineTarget::Uniform { levels: 1 },
        Some(&cylinder),
        tol(),
    )
    .expect("refines");

    assert!(report.surface_aware, "the report must admit what it did");
    assert!(
        report.max_deviation > 0.0,
        "snapping to a curved surface must move vertices off the chord"
    );

    let worst_after = max_surface_error(&refined, radius);
    assert!(
        worst_after < worst_before,
        "surface-aware refinement must get CLOSER to the analytic surface: \
         {worst_before} -> {worst_after}"
    );
}

/// Largest deviation of the faceted SURFACE from the cylinder.
///
/// Sampled at triangle centroids, not at vertices. A prism's vertices sit
/// exactly on the cylinder by construction, so a vertex metric reports
/// zero error for a visibly faceted mesh and can never improve. The
/// approximation error lives in the middle of each facet, which is
/// precisely what refinement removes.
fn max_surface_error(mesh: &TriMesh, radius: Scalar) -> Scalar {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let c = (mesh.positions[t[0] as usize]
                + mesh.positions[t[1] as usize]
                + mesh.positions[t[2] as usize])
                / 3.0;
            ((c.x * c.x + c.y * c.y).sqrt() - radius).abs()
        })
        .fold(0.0, Scalar::max)
}

/// Without a surface, refinement subdivides but does not improve.
///
/// The counterpart to the test above: it proves the improvement comes
/// from the analytic surface and not from subdivision itself.
#[test]
fn linear_refinement_does_not_improve_the_approximation() {
    let radius = 1.0;
    let sides = 8;
    let mut positions = Vec::new();
    for i in 0..sides {
        let a = std::f64::consts::TAU * (i as Scalar) / (sides as Scalar);
        positions.push(Point3::new(radius * a.cos(), radius * a.sin(), 0.0));
        positions.push(Point3::new(radius * a.cos(), radius * a.sin(), 1.0));
    }
    let mut indices = Vec::new();
    for i in 0..sides {
        let n = (i + 1) % sides;
        indices.extend_from_slice(&[
            (2 * i) as u32,
            (2 * n) as u32,
            (2 * i + 1) as u32,
            (2 * i + 1) as u32,
            (2 * n) as u32,
            (2 * n + 1) as u32,
        ]);
    }
    let prism = TriMesh::new(positions, indices);
    let before = max_surface_error(&prism, radius);

    let (refined, report) =
        refine(&prism, RefineTarget::Uniform { levels: 1 }, None, tol()).expect("refines");

    assert!(!report.surface_aware);
    assert!(
        max_surface_error(&refined, radius) >= before,
        "linear subdivision cannot improve on the chord it subdivides"
    );
}

/// An impossible request is refused with its numbers, not truncated.
#[test]
fn an_over_budget_request_is_refused_by_name() {
    let error = refine(&cube(), RefineTarget::Uniform { levels: 20 }, None, tol())
        .expect_err("20 levels of a cube is over budget");
    assert!(matches!(error, RefineError::BudgetExceeded { .. }));
}

/// A non-positive edge target is refused rather than looping forever.
#[test]
fn a_non_positive_edge_target_is_refused() {
    assert_eq!(
        refine(
            &cube(),
            RefineTarget::EdgeLength { max_edge: 0.0 },
            None,
            tol()
        )
        .unwrap_err(),
        RefineError::InvalidTarget(0.0)
    );
}

/// Edge-length refinement reaches the requested edge length.
#[test]
fn edge_length_refinement_reaches_its_target() {
    let target = 0.3;
    let (refined, _) = refine(
        &cube(),
        RefineTarget::EdgeLength { max_edge: target },
        None,
        tol(),
    )
    .expect("refines");

    let longest = refined
        .indices
        .chunks_exact(3)
        .flat_map(|t| {
            [(0, 1), (1, 2), (2, 0)].map(|(a, b)| {
                (refined.positions[t[b] as usize] - refined.positions[t[a] as usize]).length()
            })
        })
        .fold(0.0, Scalar::max);

    assert!(
        longest <= target + 1e-9,
        "longest edge {longest} must not exceed the {target} target"
    );
}

/// Weld caches moved to a hash map, whose iteration order is unstable.
/// That is safe only because they are never iterated: output ordering
/// comes from the triangle walk. A hash map seeded per-process would
/// break the documented cross-process guarantee, so assert on exact
/// vertex positions and indices, not just counts.
#[test]
fn refine_is_deterministic_across_runs() {
    let mesh = cube();
    let tolerance = tol();
    let target = RefineTarget::Uniform { levels: 2 };
    let (first, _) = refine(&mesh, target, None, tolerance).expect("refine");
    for _ in 0..8 {
        let (again, _) = refine(&mesh, target, None, tolerance).expect("refine");
        assert_eq!(first.positions, again.positions, "vertex order drifted");
        assert_eq!(first.indices, again.indices, "index order drifted");
    }
}
