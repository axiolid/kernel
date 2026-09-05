//! Minkowski sum and difference, checked against closed forms.

use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_measure::volume_properties;
use axiolid_mesh::{audit_mesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_minkowski::{
    minkowski_difference_with, minkowski_sum, minkowski_sum_with, MinkowskiError,
};

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
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(p, i)
}

/// An L-shaped solid: the canonical non-convex operand.
fn l_shape() -> TriMesh {
    let footprint = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (1.0, 1.0),
        (1.0, 2.0),
        (0.0, 2.0),
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

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::METRE)
        .expect("solid measures")
        .signed_volume
        .abs()
}

/// The sum of two boxes is a box with summed extents.
///
/// Checked against the closed form rather than against another
/// implementation, so the test is evidence about the operation itself.
#[test]
fn the_sum_of_two_boxes_is_a_box_with_summed_extents() {
    let a = boxx(Point3::ZERO, Point3::new(2.0, 3.0, 4.0));
    let b = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 2.0));

    let sum = minkowski_sum(&a, &b, Tolerance::METRE).expect("sums");

    // Extents add: 2+1, 3+1, 4+2.
    let expected = 3.0 * 4.0 * 6.0;
    let measured = volume(&sum);
    assert!(
        (measured - expected).abs() < 1e-9,
        "summed extents give {expected}, got {measured}"
    );

    let mut min = Point3::new(Scalar::INFINITY, Scalar::INFINITY, Scalar::INFINITY);
    let mut max = Point3::new(
        Scalar::NEG_INFINITY,
        Scalar::NEG_INFINITY,
        Scalar::NEG_INFINITY,
    );
    for p in &sum.positions {
        min = Point3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
        max = Point3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
    }
    let span = max - min;
    assert!(
        (span.x - 3.0).abs() < 1e-9 && (span.y - 4.0).abs() < 1e-9 && (span.z - 6.0).abs() < 1e-9,
        "extents must add, got {span:?}"
    );
}

/// A convex sum agrees with the convex hull of the pairwise vertex sums.
///
/// This is the definition, so it is asserted directly rather than trusted.
#[test]
fn a_convex_sum_is_the_hull_of_pairwise_vertex_sums() {
    let a = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 1.0));
    let b = boxx(Point3::ZERO, Point3::new(0.5, 2.0, 1.0));

    let sum = minkowski_sum(&a, &b, Tolerance::METRE).expect("sums");

    let mut pairwise = Vec::new();
    for p in &a.positions {
        for q in &b.positions {
            pairwise.push(*p + Vec3::new(q.x, q.y, q.z));
        }
    }
    let reference = axiolid_construct::hull::convex_hull(&pairwise).expect("hull");

    assert!(
        (volume(&sum) - volume(&reference)).abs() < 1e-9,
        "a convex sum must equal the hull of pairwise sums"
    );
}

/// Every result must be a solid in its own right.
#[test]
fn results_pass_diagnosis() {
    let a = boxx(Point3::ZERO, Point3::new(2.0, 2.0, 2.0));
    let b = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 1.0));

    let sum = minkowski_sum(&a, &b, Tolerance::METRE).expect("sums");
    let health = audit_mesh(&sum, Tolerance::METRE);
    assert!(
        health.is_closed_two_manifold(),
        "a sum must be closed and manifold: boundary={} non_manifold={}",
        health.boundary_edges,
        health.non_manifold_edges
    );
}

/// A non-convex sum agrees with the union of its per-part sums.
///
/// The whole point of decomposition: growing an L by a box is not the same
/// as growing its convex hull, and the difference is measurable.
#[test]
fn a_non_convex_sum_differs_from_treating_the_operand_as_convex() {
    let provider = BoolmeshBoolean::new();
    let l = l_shape();
    let tool = boxx(Point3::ZERO, Point3::new(0.5, 0.5, 0.5));

    let outcome = minkowski_sum_with(&l, &tool, Tolerance::METRE, &provider).expect("sums");

    assert!(
        outcome.evidence.subject_parts > 1,
        "an L-shape must be decomposed, got {} part(s)",
        outcome.evidence.subject_parts
    );
    assert!(!outcome.evidence.was_convex());

    // Treating the L as convex means summing its hull instead, which grows
    // the filled notch as well and is therefore strictly larger.
    let hull = axiolid_construct::hull::convex_hull(&l.positions).expect("hull");
    let convex_sum = minkowski_sum(&hull, &tool, Tolerance::METRE).expect("sums");

    let decomposed_volume = volume(&outcome.mesh);
    let convex_volume = volume(&convex_sum);
    assert!(
        decomposed_volume < convex_volume - 1e-6,
        "the decomposed sum must be smaller than the convex-hull sum: \
         {decomposed_volume} vs {convex_volume}"
    );

    // And it must still be a solid.
    let health = audit_mesh(&outcome.mesh, Tolerance::METRE);
    assert!(
        health.is_closed_two_manifold(),
        "a decomposed sum must be closed and manifold: boundary={} non_manifold={}",
        health.boundary_edges,
        health.non_manifold_edges
    );
}

/// The sum grows the solid; it never shrinks it.
#[test]
fn a_sum_contains_the_original() {
    let provider = BoolmeshBoolean::new();
    let l = l_shape();
    let tool = boxx(Point3::ZERO, Point3::new(0.25, 0.25, 0.25));

    let outcome = minkowski_sum_with(&l, &tool, Tolerance::METRE, &provider).expect("sums");
    assert!(
        volume(&outcome.mesh) > volume(&l),
        "growing a solid must not shrink it"
    );
}

/// Difference is an erosion, not a hull of pairwise differences.
///
/// Eroding a 4-cube by a unit cube leaves a 3-cube: the translations that
/// keep the tool inside span the box's extent minus the tool's.
#[test]
fn the_difference_of_two_boxes_erodes_by_the_tool_extent() {
    let provider = BoolmeshBoolean::new();
    let subject = boxx(Point3::ZERO, Point3::new(4.0, 4.0, 4.0));
    let tool = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 1.0));

    let outcome =
        minkowski_difference_with(&subject, &tool, Tolerance::METRE, &provider).expect("erodes");

    // 4-1 in every axis.
    let expected = 3.0 * 3.0 * 3.0;
    let measured = volume(&outcome.mesh);
    assert!(
        (measured - expected).abs() < 1e-6,
        "erosion must leave {expected}, got {measured}"
    );
}

/// Erosion must respect WHERE the tool sits, not just its size.
///
/// A tool whose vertices are offset from the origin erodes asymmetrically:
/// the valid translations shift by the tool's own position. A test using a
/// tool anchored at the origin cannot see a sign error, because negating a
/// symmetric offset set leaves it unchanged -- that mutant survived until
/// this test existed.
#[test]
fn erosion_places_the_result_where_the_tool_actually_fits() {
    let provider = BoolmeshBoolean::new();
    let subject = boxx(Point3::ZERO, Point3::new(4.0, 4.0, 4.0));
    // Deliberately away from the origin and not centred.
    let tool = boxx(Point3::new(1.0, 2.0, 0.5), Point3::new(2.0, 3.0, 1.5));

    let outcome =
        minkowski_difference_with(&subject, &tool, Tolerance::METRE, &provider).expect("erodes");

    // x + tool must lie inside the subject, so x ranges over
    // [0 - 1, 4 - 2] x [0 - 2, 4 - 3] x [0 - 0.5, 4 - 1.5].
    let mut min = Point3::new(Scalar::INFINITY, Scalar::INFINITY, Scalar::INFINITY);
    let mut max = Point3::new(
        Scalar::NEG_INFINITY,
        Scalar::NEG_INFINITY,
        Scalar::NEG_INFINITY,
    );
    for p in &outcome.mesh.positions {
        min = Point3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
        max = Point3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
    }

    for (axis, got, want) in [
        ("x.min", min.x, -1.0),
        ("x.max", max.x, 2.0),
        ("y.min", min.y, -2.0),
        ("y.max", max.y, 1.0),
        ("z.min", min.z, -0.5),
        ("z.max", max.z, 2.5),
    ] {
        assert!(
            (got - want).abs() < 1e-6,
            "erosion {axis} must be {want}, got {got}"
        );
    }
}

/// A tool larger than the subject fits nowhere.
#[test]
fn eroding_by_an_oversized_tool_leaves_nothing() {
    let provider = BoolmeshBoolean::new();
    let subject = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 1.0));
    let tool = boxx(Point3::ZERO, Point3::new(2.0, 2.0, 2.0));

    let outcome =
        minkowski_difference_with(&subject, &tool, Tolerance::METRE, &provider).expect("erodes");
    assert!(
        outcome.mesh.indices.is_empty() || volume(&outcome.mesh) < 1e-9,
        "a tool larger than the subject must leave nothing"
    );
}

/// Erosion refuses a non-convex subject rather than overstating the result.
#[test]
fn erosion_refuses_a_non_convex_subject() {
    let provider = BoolmeshBoolean::new();
    let tool = boxx(Point3::ZERO, Point3::new(0.25, 0.25, 0.25));

    let error = minkowski_difference_with(&l_shape(), &tool, Tolerance::METRE, &provider)
        .expect_err("must refuse");
    assert!(
        matches!(error, MinkowskiError::ErosionSubjectNotConvex),
        "a non-convex subject must be refused by name, got {error:?}"
    );
}

/// The convex-only entry point refuses a non-convex operand by name.
#[test]
fn the_convex_entry_point_refuses_a_non_convex_operand() {
    let tool = boxx(Point3::ZERO, Point3::new(0.5, 0.5, 0.5));

    let error = minkowski_sum(&l_shape(), &tool, Tolerance::METRE).expect_err("must refuse");
    assert!(
        matches!(error, MinkowskiError::NotConvex { operand: "subject" }),
        "a non-convex subject must be named, got {error:?}"
    );
}

/// Malformed operands are refused before any geometry is attempted.
#[test]
fn malformed_operands_are_refused() {
    let cube = boxx(Point3::ZERO, Point3::new(1.0, 1.0, 1.0));

    // A single triangle is a surface, not a solid.
    let sheet = TriMesh::new(
        vec![
            Point3::ZERO,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2],
    );
    assert!(matches!(
        minkowski_sum(&sheet, &cube, Tolerance::METRE),
        Err(MinkowskiError::NotASolid {
            operand: "subject",
            ..
        })
    ));
    assert!(matches!(
        minkowski_sum(&cube, &sheet, Tolerance::METRE),
        Err(MinkowskiError::NotASolid {
            operand: "tool",
            ..
        })
    ));

    let empty = TriMesh::new(Vec::new(), Vec::new());
    assert!(matches!(
        minkowski_sum(&empty, &cube, Tolerance::METRE),
        Err(MinkowskiError::EmptyOperand("subject"))
    ));
}

/// The sum is commutative, as the definition requires.
#[test]
fn the_sum_is_commutative() {
    let a = boxx(Point3::ZERO, Point3::new(2.0, 1.0, 1.0));
    let b = boxx(Point3::ZERO, Point3::new(1.0, 3.0, 1.0));

    let forward = minkowski_sum(&a, &b, Tolerance::METRE).expect("sums");
    let backward = minkowski_sum(&b, &a, Tolerance::METRE).expect("sums");

    assert!(
        (volume(&forward) - volume(&backward)).abs() < 1e-9,
        "A + B and B + A must agree"
    );
}
