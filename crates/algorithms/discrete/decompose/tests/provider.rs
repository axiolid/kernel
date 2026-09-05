//! The same decomposition, driven by a boolean provider.
//!
//! These prove the switch is real: identical input, the same
//! `Decomposition` type, and guarantees judged by the same measurements.
//! Where the two paths differ, the difference is measured, not assumed.

use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_decompose::split::Splitter;
use axiolid_decompose::{convex_decompose_with, Strategy};
use axiolid_measure::volume_properties;
use axiolid_mesh::{audit_mesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;

/// An L-shaped solid: the canonical non-convex case.
fn l_shape() -> TriMesh {
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
    for &(a, b, c) in &[(0u32, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5)] {
        indices.extend_from_slice(&[a, c, b]);
        indices.extend_from_slice(&[a + n, b + n, c + n]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        indices.extend_from_slice(&[i, j, j + n]);
        indices.extend_from_slice(&[i, j + n, i + n]);
    }
    TriMesh::new(positions, indices)
}

fn cube() -> TriMesh {
    let p = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
    ];
    let i = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(p, i)
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::METRE)
        .expect("solid measures")
        .signed_volume
        .abs()
}

/// A convex solid needs no splitting, whichever splitter is used.
#[test]
fn a_convex_solid_is_one_part_under_the_provider_too() {
    let provider = BoolmeshBoolean::new();
    let result = convex_decompose_with(
        &cube(),
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::Provider(&provider),
    )
    .expect("decomposes");

    assert!(result.is_single_part(), "a cube must not be split");
    assert_eq!(result.splits, 0);
}

/// The provider path produces parts that are genuine solids.
///
/// This is the guarantee the hand-rolled clipper does not yet meet, so it
/// is the reason the switch exists.
#[test]
fn provider_parts_are_closed_two_manifold_solids() {
    let provider = BoolmeshBoolean::new();
    let result = convex_decompose_with(
        &l_shape(),
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::Provider(&provider),
    )
    .expect("decomposes");

    assert!(
        result.parts.len() > 1,
        "an L-shape must be split, got {} part(s)",
        result.parts.len()
    );
    for (index, part) in result.parts.iter().enumerate() {
        let health = audit_mesh(part, Tolerance::METRE);
        assert!(
            health.is_closed_two_manifold(),
            "part {index} is not a solid: boundary={} non_manifold={} winding={}",
            health.boundary_edges,
            health.non_manifold_edges,
            health.inconsistent_winding_edges
        );
    }
}

/// The parts must add back up to the solid they came from.
#[test]
fn provider_parts_reproduce_the_input_volume() {
    let provider = BoolmeshBoolean::new();
    let input = l_shape();
    let expected = volume(&input);

    let result = convex_decompose_with(
        &input,
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::Provider(&provider),
    )
    .expect("decomposes");

    let total: Scalar = result.parts.iter().map(volume).sum();
    assert!(
        (total - expected).abs() < 1e-6,
        "parts must sum to the input volume: {total} vs {expected}"
    );
}

/// Every part must itself be convex, or the decomposition is not one.
#[test]
fn provider_parts_are_convex() {
    let provider = BoolmeshBoolean::new();
    let result = convex_decompose_with(
        &l_shape(),
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::Provider(&provider),
    )
    .expect("decomposes");

    for (index, part) in result.parts.iter().enumerate() {
        let again = convex_decompose_with(
            part,
            Strategy::Exact,
            Tolerance::METRE,
            &Splitter::Provider(&provider),
        )
        .expect("a part decomposes");
        assert!(
            again.is_single_part(),
            "part {index} is not convex: it split into {} further parts",
            again.parts.len()
        );
    }
}

/// Switching splitter must not change what the result CLAIMS about itself.
#[test]
fn both_splitters_report_the_same_fidelity() {
    let provider = BoolmeshBoolean::new();
    let input = l_shape();

    let hand = convex_decompose_with(
        &input,
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::HandRolled,
    )
    .expect("decomposes");
    let via_provider = convex_decompose_with(
        &input,
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::Provider(&provider),
    )
    .expect("decomposes");

    assert_eq!(
        hand.fidelity, via_provider.fidelity,
        "fidelity is a property of the request, not of the splitter"
    );
}

/// Both splitters must agree on the SOLID, not merely on their own claims.
///
/// This is the point of offering a choice: the two paths are independent
/// implementations of the same contract, so each is evidence about the
/// other. Volume is the comparable quantity -- the triangulations differ
/// by construction, so comparing indices would test nothing.
#[test]
fn both_splitters_produce_the_same_solid() {
    let provider = BoolmeshBoolean::new();
    let input = l_shape();
    let expected = volume(&input);

    let hand = convex_decompose_with(
        &input,
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::HandRolled,
    )
    .expect("decomposes");
    let via_provider = convex_decompose_with(
        &input,
        Strategy::Exact,
        Tolerance::METRE,
        &Splitter::Provider(&provider),
    )
    .expect("decomposes");

    let hand_total: Scalar = hand.parts.iter().map(volume).sum();
    let provider_total: Scalar = via_provider.parts.iter().map(volume).sum();

    assert!(
        (hand_total - expected).abs() < 1e-6,
        "hand-rolled parts must sum to the input: {hand_total} vs {expected}"
    );
    assert!(
        (provider_total - expected).abs() < 1e-6,
        "provider parts must sum to the input: {provider_total} vs {expected}"
    );
    assert!(
        (hand_total - provider_total).abs() < 1e-6,
        "the two splitters must agree: {hand_total} vs {provider_total}"
    );

    // Both must produce solids. A path that returned open shells could
    // still pass a volume test by accident, so this is asserted separately.
    for (label, result) in [("hand-rolled", &hand), ("provider", &via_provider)] {
        for (index, part) in result.parts.iter().enumerate() {
            assert!(
                audit_mesh(part, Tolerance::METRE).is_closed_two_manifold(),
                "{label} part {index} is not a closed two-manifold solid"
            );
        }
    }
}

/// A solid with two separate notches, to test more than one split.
///
/// The L-shape needs a single cut. A shape needing several exercises the
/// recursion, where a part produced by one split is itself split again --
/// the case where a capper that is only accidentally correct falls over.
#[test]
fn a_solid_with_several_notches_decomposes_under_both_splitters() {
    // A plus/cross footprint: four reflex corners, extruded.
    let footprint = [
        (1.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (3.0, 1.0),
        (3.0, 2.0),
        (2.0, 2.0),
        (2.0, 3.0),
        (1.0, 3.0),
        (1.0, 2.0),
        (0.0, 2.0),
        (0.0, 1.0),
        (1.0, 1.0),
    ];
    let mut positions = Vec::new();
    for &(x, y) in &footprint {
        positions.push(Point3::new(x, y, 0.0));
    }
    for &(x, y) in &footprint {
        positions.push(Point3::new(x, y, 1.0));
    }
    let n = footprint.len() as u32;
    let mut indices = Vec::new();
    // A cross is NOT star-shaped about any of its corners, so a fan from
    // vertex 0 produces triangles outside the footprint. Ear clipping is
    // the correct triangulation and is already in-tree.
    let ring: Vec<axiolid_core::Point2> = footprint
        .iter()
        .map(|&(x, y)| axiolid_core::Point2::new(x, y))
        .collect();
    let fan = axiolid_reference::polygon::triangulate_simple(&ring).expect("cross triangulates");
    for triple in &fan {
        indices.extend_from_slice(&[triple[0], triple[2], triple[1]]);
        indices.extend_from_slice(&[n + triple[0], n + triple[1], n + triple[2]]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        indices.extend_from_slice(&[i, j, j + n]);
        indices.extend_from_slice(&[i, j + n, i + n]);
    }
    let cross = TriMesh::new(positions, indices);

    let expected = volume(&cross);
    let provider = BoolmeshBoolean::new();

    for (label, splitter) in [
        ("hand-rolled", Splitter::HandRolled),
        ("provider", Splitter::Provider(&provider)),
    ] {
        let result = convex_decompose_with(&cross, Strategy::Exact, Tolerance::METRE, &splitter)
            .unwrap_or_else(|error| panic!("{label} decomposes: {error}"));

        assert!(
            result.parts.len() > 1,
            "{label}: a cross is not convex and must be split"
        );
        let total: Scalar = result.parts.iter().map(volume).sum();
        assert!(
            (total - expected).abs() < 1e-6,
            "{label}: parts must sum to the input volume: {total} vs {expected}"
        );
        for (index, part) in result.parts.iter().enumerate() {
            assert!(
                audit_mesh(part, Tolerance::METRE).is_closed_two_manifold(),
                "{label} part {index} is not a closed two-manifold solid"
            );
        }
    }
}
