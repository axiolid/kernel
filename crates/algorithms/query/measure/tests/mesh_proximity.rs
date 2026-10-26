//! Mesh-level proximity: distance, witnesses, components.
//!
//! Distances are checked against values derived from the fixture geometry and
//! against an independent brute-force scan over sampled surface points, so a
//! systematic error in the composition cannot pass unnoticed.

use axiolid_core::Point3;
use axiolid_measure::{mesh_distance, proximity_components, MeshProximityError};
use axiolid_mesh::TriMesh;

/// A unit square in the z = `z` plane, offset by (dx, dy).
fn square(dx: f64, dy: f64, z: f64) -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(dx, dy, z),
            Point3::new(dx + 1.0, dy, z),
            Point3::new(dx + 1.0, dy + 1.0, z),
            Point3::new(dx, dy + 1.0, z),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

/// Brute-force minimum distance over densely sampled surface points.
///
/// Deliberately independent of the implementation: it samples barycentric
/// points on every triangle and scans all pairs. It over-estimates (samples
/// need not include the true witness), so it is an UPPER bound on the true
/// distance, and the exact answer must never exceed it.
fn sampled_minimum(first: &TriMesh, second: &TriMesh, steps: usize) -> f64 {
    let samples = |mesh: &TriMesh| {
        let mut points = Vec::new();
        for corner in mesh.indices.as_chunks::<3>().0 {
            let a = mesh.positions[corner[0] as usize];
            let b = mesh.positions[corner[1] as usize];
            let c = mesh.positions[corner[2] as usize];
            for i in 0..=steps {
                for j in 0..=(steps - i) {
                    let u = i as f64 / steps as f64;
                    let v = j as f64 / steps as f64;
                    let w = 1.0 - u - v;
                    points.push(a * w + b * u + c * v);
                }
            }
        }
        points
    };

    let (points_a, points_b) = (samples(first), samples(second));
    let mut best = f64::INFINITY;
    for a in &points_a {
        for b in &points_b {
            best = best.min((*a - *b).length());
        }
    }
    best
}

#[test]
fn parallel_squares_are_separated_by_their_gap() {
    // Two coincident unit squares 3 apart in z: the distance is exactly 3 and
    // the witnesses are vertically aligned.
    let result =
        mesh_distance(&square(0.0, 0.0, 0.0), &square(0.0, 0.0, 3.0)).expect("two valid meshes");

    let distance = result.distance_squared.sqrt();
    assert!((distance - 3.0).abs() < 1e-12, "expected 3, got {distance}");
    assert!(!result.surfaces_cross);
    // The witness must actually realise the reported distance.
    let witness_gap = (result.point_a - result.point_b).length();
    assert!((witness_gap - distance).abs() < 1e-12);
}

#[test]
fn the_exact_distance_never_exceeds_a_sampled_scan() {
    // Offset diagonally so the closest approach is a corner-to-corner pair,
    // which a naive centroid or vertex-only method gets wrong.
    let first = square(0.0, 0.0, 0.0);
    let second = square(2.5, 1.5, 4.0);

    let exact = mesh_distance(&first, &second)
        .expect("two valid meshes")
        .distance_squared
        .sqrt();
    let sampled = sampled_minimum(&first, &second, 12);

    // The scan only visits sampled points, so it can never beat the exact
    // minimum. If it does, the exact path missed a closer pair.
    assert!(
        exact <= sampled + 1e-12,
        "exact {exact} exceeded sampled upper bound {sampled}"
    );
    // And the sampling is dense enough that the two agree closely, so the
    // bound above is a real check rather than a vacuous one.
    assert!(
        (sampled - exact) < 0.2,
        "sampled {sampled} and exact {exact} disagree too much to be meaningful"
    );
}

#[test]
fn touching_and_crossing_both_report_zero_surface_separation() {
    // Edge-to-edge contact: distance zero, surfaces meet.
    let touching =
        mesh_distance(&square(0.0, 0.0, 0.0), &square(1.0, 0.0, 0.0)).expect("two valid meshes");
    assert!(touching.distance_squared.abs() < 1e-24);
    assert!(touching.surfaces_cross);

    // Heavily overlapping: also zero, because SURFACE separation is zero.
    // The two cases are indistinguishable by distance alone, which is exactly
    // why penetration depth is not inferred from this number.
    let crossing =
        mesh_distance(&square(0.0, 0.0, 0.0), &square(0.25, 0.25, 0.0)).expect("two valid meshes");
    assert!(crossing.surfaces_cross);
}

/// Two coplanar squares far apart in x, as one mesh.
fn two_squares(gap: f64) -> TriMesh {
    let mut mesh = square(0.0, 0.0, 0.0);
    let far = square(1.0 + gap, 0.0, 0.0);
    let base = mesh.positions.len() as u32;
    mesh.positions.extend(far.positions);
    mesh.indices.extend(far.indices.iter().map(|i| i + base));
    mesh
}

#[test]
fn separated_approaches_are_reported_as_distinct_components() {
    // A long bar above two squares that are far apart: the bar approaches
    // each square separately, so there are two proximity regions, not one.
    let pair = two_squares(6.0);
    let bar = TriMesh::new(
        vec![
            Point3::new(0.0, 0.4, 0.5),
            Point3::new(8.0, 0.4, 0.5),
            Point3::new(8.0, 0.6, 0.5),
            Point3::new(0.0, 0.6, 0.5),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );

    // The bar sits 0.5 above both squares; a threshold of 0.6 catches both
    // approaches but nothing in the 6-unit gap between the squares.
    let components = proximity_components(&pair, &bar, 0.6).expect("valid query");

    assert_eq!(components.len(), 2, "one approach per square");
    // Nearest first, and each component must be within the threshold.
    for component in &components {
        assert!(component.witness.distance_squared <= 0.6 * 0.6 + 1e-12);
        assert!(!component.triangles_a.is_empty());
        assert!(!component.triangles_b.is_empty());
    }
    assert!(
        components[0].witness.distance_squared <= components[1].witness.distance_squared,
        "components must be ordered nearest first"
    );
}

#[test]
fn a_threshold_below_the_separation_finds_nothing() {
    let components = proximity_components(&square(0.0, 0.0, 0.0), &square(0.0, 0.0, 3.0), 1.0)
        .expect("valid query");
    assert!(components.is_empty(), "3 apart is outside a 1.0 threshold");
}

#[test]
fn an_empty_mesh_is_refused_rather_than_reported_as_infinitely_far() {
    let empty = TriMesh::new(Vec::new(), Vec::new());
    let error = mesh_distance(&square(0.0, 0.0, 0.0), &empty)
        .expect_err("distance to nothing is undefined");
    assert_eq!(error, MeshProximityError::EmptyMesh);
}

#[test]
fn a_negative_threshold_is_refused() {
    let error = proximity_components(&square(0.0, 0.0, 0.0), &square(0.0, 0.0, 1.0), -1.0)
        .expect_err("a negative threshold is a caller error");
    assert_eq!(error, MeshProximityError::InvalidThreshold);
}

/// Two closed boxes where one passes entirely through the other.
///
/// Reduced from a real beam/duct fixture. The beam spans x 0..8000, y -500..500,
/// z 0..1000. The duct spans x 3900..4100, y -1000..1000, z 400..600 — it is
/// narrower in x and z but WIDER in y, so it enters one side of the beam and
/// exits the other. Their surfaces genuinely cross.
///
/// Every duct vertex lies outside the beam and every beam vertex lies outside
/// the duct, and the crossing is edge-through-face, so a scan over
/// vertex/triangle and edge/edge pairs alone never samples the intersection.
#[test]
fn a_box_pierced_by_a_wider_box_reports_zero_separation() {
    let beam = TriMesh::new(
        vec![
            Point3::new(0.0, 500.0, 0.0),
            Point3::new(0.0, -500.0, 0.0),
            Point3::new(0.0, -500.0, 1000.0),
            Point3::new(0.0, 500.0, 1000.0),
            Point3::new(8000.0, -500.0, 0.0),
            Point3::new(8000.0, 500.0, 0.0),
            Point3::new(8000.0, 500.0, 1000.0),
            Point3::new(8000.0, -500.0, 1000.0),
        ],
        vec![
            0, 1, 2, 2, 3, 0, 4, 5, 6, 6, 7, 4, 5, 0, 3, 3, 6, 5, 6, 3, 2, 2, 7, 6, 7, 2, 1, 1, 4,
            7, 4, 1, 0, 0, 5, 4,
        ],
    );
    let duct = TriMesh::new(
        vec![
            Point3::new(3900.0, -1000.0, 600.0),
            Point3::new(3900.0, -1000.0, 400.0),
            Point3::new(4100.0, -1000.0, 400.0),
            Point3::new(4100.0, -1000.0, 600.0),
            Point3::new(3900.0, 1000.0, 400.0),
            Point3::new(3900.0, 1000.0, 600.0),
            Point3::new(4100.0, 1000.0, 600.0),
            Point3::new(4100.0, 1000.0, 400.0),
        ],
        vec![
            0, 1, 2, 2, 3, 0, 4, 5, 6, 6, 7, 4, 5, 0, 3, 3, 6, 5, 6, 3, 2, 2, 7, 6, 7, 2, 1, 1, 4,
            7, 4, 1, 0, 0, 5, 4,
        ],
    );

    let crossing = mesh_distance(&beam, &duct).expect("two valid meshes");
    assert!(
        crossing.surfaces_cross,
        "a pierced box must report zero separation, got distance_squared {}",
        crossing.distance_squared
    );
}
