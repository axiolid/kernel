//! Decomposition: convex parts, honest about whether they reproduce the input.

use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_decompose::{convex_decompose, DecomposeError, Fidelity, Strategy};
use axiolid_measure::volume_properties;
use axiolid_mesh::{audit_mesh, TriMesh};

/// An axis-aligned box as a closed two-manifold solid.
fn boxx(min: Point3, max: Point3) -> TriMesh {
    let p = vec![
        Point3::new(min.x, min.y, min.z),
        Point3::new(max.x, min.y, min.z),
        Point3::new(max.x, max.y, min.z),
        Point3::new(min.x, max.y, min.z),
        Point3::new(min.x, min.y, max.z),
        Point3::new(max.x, min.y, max.z),
        Point3::new(max.x, max.y, max.z),
        Point3::new(min.x, max.y, max.z),
    ];
    let i = vec![
        0, 2, 1, 0, 3, 2, // bottom
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // front
        1, 2, 6, 1, 6, 5, // right
        2, 3, 7, 2, 7, 6, // back
        3, 0, 4, 3, 4, 7, // left
    ];
    TriMesh::new(p, i)
}

/// An L-shaped solid: the canonical non-convex test case.
///
/// Built from two boxes sharing a face, written out as one closed shell so
/// the input is a genuine solid rather than two overlapping ones.
fn l_shape() -> TriMesh {
    // Footprint of the L, extruded in z.
    let footprint = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (1.0, 1.0),
        (1.0, 2.0),
        (0.0, 2.0),
    ];
    let height = 1.0;

    let mut positions = Vec::new();
    for &(x, y) in &footprint {
        positions.push(Point3::new(x, y, 0.0));
    }
    for &(x, y) in &footprint {
        positions.push(Point3::new(x, y, height));
    }

    let n = footprint.len() as u32;
    let mut indices = Vec::new();

    // The footprint is convex-fan-able from vertex 0 except at the reflex
    // corner, so triangulate explicitly rather than fanning.
    let caps = [(0u32, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5)];
    for &(a, b, c) in &caps {
        indices.extend_from_slice(&[a, c, b]); // bottom, wound downward
        indices.extend_from_slice(&[a + n, b + n, c + n]); // top
    }
    for i in 0..n {
        let j = (i + 1) % n;
        indices.extend_from_slice(&[i, j, j + n]);
        indices.extend_from_slice(&[i, j + n, i + n]);
    }
    TriMesh::new(positions, indices)
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::METRE)
        .expect("solid measures")
        .signed_volume
        .abs()
}

/// A convex solid is already its own decomposition.
#[test]
fn a_convex_solid_yields_one_part() {
    let cube = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 1.0));
    let result = convex_decompose(&cube, Strategy::Exact, Tolerance::METRE).expect("decomposes");

    assert!(
        result.is_single_part(),
        "a convex solid must not be split, got {} parts",
        result.parts.len()
    );
    assert_eq!(result.splits, 0);
    assert_eq!(result.fidelity, Fidelity::Exact);
}

/// The non-convex case must actually split.
#[test]
fn a_reflex_solid_is_split_into_convex_parts() {
    let l = l_shape();
    let result = convex_decompose(&l, Strategy::Exact, Tolerance::METRE).expect("decomposes");

    assert!(
        result.parts.len() > 1,
        "an L-shape is not convex and must be split"
    );

    // Each part must itself be convex: decomposing a part again is a no-op.
    for (index, part) in result.parts.iter().enumerate() {
        let again =
            convex_decompose(part, Strategy::Exact, Tolerance::METRE).expect("a part decomposes");
        assert!(
            again.is_single_part(),
            "part {index} is not convex: it split into {} further parts",
            again.parts.len()
        );
    }
}

/// The union of the parts must reproduce the input's volume.
///
/// Checked by measurement rather than by index comparison, because the
/// parts are hulls with their own triangulation: only the enclosed volume
/// is meaningfully comparable. Parts of a convex decomposition meet on
/// shared faces without overlapping, so their volumes sum to the input's.
#[test]
fn the_parts_reproduce_the_input_volume() {
    let l = l_shape();
    let expected = volume(&l);
    let result = convex_decompose(&l, Strategy::Exact, Tolerance::METRE).expect("decomposes");

    let total: Scalar = result.parts.iter().map(volume).sum();
    assert!(
        (total - expected).abs() < 1e-9,
        "parts must sum to the input volume: {total} vs {expected}"
    );
}

/// Every part must be a solid in its own right.
#[test]
fn every_part_is_a_closed_two_manifold_solid() {
    let l = l_shape();
    let result = convex_decompose(&l, Strategy::Exact, Tolerance::METRE).expect("decomposes");

    for (index, part) in result.parts.iter().enumerate() {
        let health = audit_mesh(part, Tolerance::METRE);
        assert!(
            health.is_closed_two_manifold(),
            "part {index} is not a solid: boundary={} non_manifold={}",
            health.boundary_edges,
            health.non_manifold_edges
        );
    }
}

/// The approximate strategy trades parts for fidelity, and says so.
#[test]
fn the_approximate_strategy_reports_what_it_achieved() {
    let l = l_shape();

    let exact = convex_decompose(&l, Strategy::Exact, Tolerance::METRE).expect("decomposes");
    let loose = convex_decompose(
        &l,
        Strategy::Approximate {
            max_concavity: 10.0,
        },
        Tolerance::METRE,
    )
    .expect("decomposes");

    // A bound larger than the solid tolerates the whole thing as one part.
    assert!(
        loose.parts.len() <= exact.parts.len(),
        "a looser bound must not need more parts: {} vs {}",
        loose.parts.len(),
        exact.parts.len()
    );

    match loose.fidelity {
        Fidelity::Approximate {
            requested,
            achieved,
        } => {
            assert_eq!(requested, 10.0);
            assert!(
                achieved <= requested,
                "achieved concavity {achieved} must respect the {requested} bound"
            );
        }
        other => panic!("an approximate request must not report {other:?}"),
    }
}

/// Fidelity is never implied: an approximate result says it is approximate.
#[test]
fn an_approximate_result_never_claims_to_be_exact() {
    let l = l_shape();
    let result = convex_decompose(
        &l,
        Strategy::Approximate { max_concavity: 0.5 },
        Tolerance::METRE,
    )
    .expect("decomposes");

    assert!(
        !matches!(result.fidelity, Fidelity::Exact),
        "an approximate strategy must never report Fidelity::Exact"
    );
}

/// An open surface has no meaningful decomposition.
#[test]
fn an_open_surface_is_refused_by_name() {
    // A single triangle: a surface, not a solid.
    let sheet = TriMesh::new(
        vec![
            Point3::ZERO,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2],
    );

    match convex_decompose(&sheet, Strategy::Exact, Tolerance::METRE) {
        Err(DecomposeError::NotASolid { boundary, .. }) => {
            assert!(boundary > 0, "an open sheet has boundary edges");
        }
        other => panic!("an open surface must be refused by name, got {other:?}"),
    }
}

/// Malformed input is refused before any geometry is attempted.
#[test]
fn malformed_input_is_refused() {
    let cube = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 1.0));

    let mut ragged = cube.clone();
    ragged.indices.push(0);
    assert!(matches!(
        convex_decompose(&ragged, Strategy::Exact, Tolerance::METRE),
        Err(DecomposeError::RaggedIndices(_))
    ));

    let mut out_of_range = cube.clone();
    out_of_range.indices[0] = 99;
    assert!(matches!(
        convex_decompose(&out_of_range, Strategy::Exact, Tolerance::METRE),
        Err(DecomposeError::IndexOutOfRange(..))
    ));

    assert!(matches!(
        convex_decompose(
            &cube,
            Strategy::Approximate {
                max_concavity: -1.0
            },
            Tolerance::METRE
        ),
        Err(DecomposeError::InvalidBound(_))
    ));
}

/// Part ordering must not depend on traversal order.
#[test]
fn decomposition_is_deterministic() {
    let l = l_shape();
    let first = convex_decompose(&l, Strategy::Exact, Tolerance::METRE).expect("decomposes");
    let second = convex_decompose(&l, Strategy::Exact, Tolerance::METRE).expect("decomposes");

    assert_eq!(first.parts.len(), second.parts.len());
    for (a, b) in first.parts.iter().zip(second.parts.iter()) {
        assert_eq!(a.positions, b.positions, "part positions must be stable");
        assert_eq!(a.indices, b.indices, "part indices must be stable");
    }
}
