//! Certified mesh Hausdorff distance (#148): closed forms, asymmetric cases,
//! and a brute-force dense-sampling cross-check.
//!
//! Every closed form here is derived from the fixture geometry, not from the
//! implementation, and every interval must contain it while being at most the
//! requested accuracy wide.

use axiolid_core::Point3;
use axiolid_measure::{
    closest_point_on_triangle, hausdorff_distance, one_sided_hausdorff, HausdorffBounds,
    HausdorffError,
};
use axiolid_mesh::TriMesh;

const ACCURACY: f64 = 1e-9;

/// An `n` by `n` grid over `[x0, x0 + size] x [y0, y0 + size]` at height `z`.
fn grid(x0: f64, y0: f64, size: f64, z: f64, n: u32) -> TriMesh {
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let (u, v) = (f64::from(i) / f64::from(n), f64::from(j) / f64::from(n));
            positions.push(Point3::new(x0 + size * u, y0 + size * v, z));
        }
    }
    let mut indices = Vec::new();
    let at = |i: u32, j: u32| j * (n + 1) + i;
    for j in 0..n {
        for i in 0..n {
            indices.extend([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            indices.extend([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    TriMesh::new(positions, indices)
}

/// The regular tetrahedron with corners at alternate cube corners, scaled.
///
/// Circumradius `sqrt(3) * scale`, inradius `scale / sqrt(3)`.
fn tetrahedron(scale: f64) -> TriMesh {
    let corners = [
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(1.0, -1.0, -1.0),
        Point3::new(-1.0, 1.0, -1.0),
        Point3::new(-1.0, -1.0, 1.0),
    ];
    TriMesh::new(
        corners.iter().map(|c| *c * scale).collect(),
        vec![0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2],
    )
}

/// An icosphere with every vertex on the unit sphere.
fn icosphere(levels: u32) -> TriMesh {
    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut positions: Vec<Point3> = [
        (-1.0, t, 0.0),
        (1.0, t, 0.0),
        (-1.0, -t, 0.0),
        (1.0, -t, 0.0),
        (0.0, -1.0, t),
        (0.0, 1.0, t),
        (0.0, -1.0, -t),
        (0.0, 1.0, -t),
        (t, 0.0, -1.0),
        (t, 0.0, 1.0),
        (-t, 0.0, -1.0),
        (-t, 0.0, 1.0),
    ]
    .iter()
    .map(|&(x, y, z)| Point3::new(x, y, z).normalize())
    .collect();
    let mut faces: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..levels {
        let mut midpoints = std::collections::HashMap::new();
        let mut midpoint = |a: u32, b: u32, positions: &mut Vec<Point3>| {
            *midpoints.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let m = (positions[a as usize] + positions[b as usize]).normalize();
                positions.push(m);
                (positions.len() - 1) as u32
            })
        };
        let mut next = Vec::new();
        for [a, b, c] in faces {
            let ab = midpoint(a, b, &mut positions);
            let bc = midpoint(b, c, &mut positions);
            let ca = midpoint(c, a, &mut positions);
            next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        faces = next;
    }
    TriMesh::new(positions, faces.into_iter().flatten().collect())
}

fn triangles(mesh: &TriMesh) -> Vec<[Point3; 3]> {
    mesh.indices
        .chunks_exact(3)
        .map(|c| {
            [
                mesh.positions[c[0] as usize],
                mesh.positions[c[1] as usize],
                mesh.positions[c[2] as usize],
            ]
        })
        .collect()
}

/// Smallest distance from the origin to a face plane: the inradius of a
/// convex mesh around the origin.
fn inradius(mesh: &TriMesh) -> f64 {
    triangles(mesh)
        .iter()
        .map(|[a, b, c]| (*b - *a).cross(*c - *a).normalize().dot(*a).abs())
        .fold(f64::INFINITY, f64::min)
}

fn check(bounds: &HausdorffBounds, expected: f64, accuracy: f64) {
    assert!(
        bounds.contains(expected),
        "[{}, {}] misses {expected}",
        bounds.lower,
        bounds.upper
    );
    assert!(
        bounds.width() <= accuracy,
        "[{}, {}] wider than {accuracy}",
        bounds.lower,
        bounds.upper
    );
    witness(bounds);
}

/// The witness pair is a point and a point of the other mesh at least
/// `lower` from it.
fn witness(bounds: &HausdorffBounds) {
    let gap = (bounds.point_from - bounds.point_to).length();
    assert!(
        gap >= bounds.lower - 1e-12,
        "witnesses {gap} apart, below the lower bound {}",
        bounds.lower
    );
}

#[test]
fn a_plane_patch_is_its_offset_away_from_a_parallel_copy() {
    // Differently triangulated copies, so vertices do not line up.
    let a = grid(0.0, 0.0, 2.0, 0.0, 3);
    let b = grid(0.0, 0.0, 2.0, 0.75, 5);
    let result = hausdorff_distance(&a, &b, ACCURACY).expect("valid meshes");
    check(&result.forward, 0.75, ACCURACY);
    check(&result.backward, 0.75, ACCURACY);
    check(&result.distance, 0.75, ACCURACY);
}

#[test]
fn a_finely_triangulated_flat_face_closes_to_the_accuracy() {
    // Every piece of one grid sits over seams and vertices of the other;
    // only the flat patches let those pieces close without splitting each
    // one down to the square root of the accuracy.
    let a = grid(0.0, 0.0, 10.0, 0.0, 40);
    let b = grid(0.0, 0.0, 10.0, 0.5, 30);
    let result = hausdorff_distance(&a, &b, ACCURACY).expect("valid meshes");
    check(&result.forward, 0.5, ACCURACY);
    check(&result.backward, 0.5, ACCURACY);
}

#[test]
fn a_patch_inside_a_larger_coplanar_one_is_asymmetric() {
    // The unit square lies in the larger one: h(small, large) = 0, while the
    // larger one's corner is sqrt(2) from the unit square's nearest corner.
    let small = grid(0.0, 0.0, 1.0, 0.0, 2);
    let large = grid(-1.0, -1.0, 3.0, 0.0, 3);
    let result = hausdorff_distance(&small, &large, ACCURACY).expect("valid meshes");
    check(&result.forward, 0.0, ACCURACY);
    check(&result.backward, 2.0_f64.sqrt(), ACCURACY);
    check(&result.distance, 2.0_f64.sqrt(), ACCURACY);
    // The witness of the larger side sits at a corner of the large square.
    let corner = result.backward.point_from;
    assert!((corner.x.abs() - 1.0).abs() < 1e-6 || (corner.x - 2.0).abs() < 1e-6);
}

#[test]
fn a_tetrahedron_and_its_scaled_copy_differ_by_inradius_and_circumradius() {
    // Scaled by s about the centroid: the outer corners are (s - 1) R from
    // the inner ones, the inner faces (s - 1) r from the outer faces.
    let s = 2.5;
    let small = tetrahedron(1.0);
    let large = tetrahedron(s);
    let (big_r, small_r) = (3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt());
    let result = hausdorff_distance(&small, &large, ACCURACY).expect("valid meshes");
    check(&result.forward, (s - 1.0) * small_r, ACCURACY);
    check(&result.backward, (s - 1.0) * big_r, ACCURACY);
    check(&result.distance, (s - 1.0) * big_r, ACCURACY);
}

#[test]
fn a_sphere_mesh_and_a_finer_one_meet_the_analytic_bounds() {
    let coarse = icosphere(1);
    let fine = icosphere(3);
    let (rc, rf) = (inradius(&coarse), inradius(&fine));
    let result = hausdorff_distance(&coarse, &fine, 1e-8).expect("valid meshes");

    // Coarse points lie between radius rc and 1, the fine surface encloses
    // the ball of radius rf > rc: h(coarse, fine) is in [rf - rc, 1 - rc].
    let forward = result.forward;
    assert!(forward.width() <= 1e-8);
    assert!(forward.upper >= rf - rc && forward.lower <= 1.0 - rc);
    // A fine vertex x on the unit sphere is at least 1 - max x.v from the
    // coarse polytope (its support function), and every fine point is
    // within 1 - min(rc, rf) of the coarse surface along its ray.
    let backward = result.backward;
    let support = fine
        .positions
        .iter()
        .map(|x| {
            1.0 - coarse
                .positions
                .iter()
                .map(|v| x.dot(*v))
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .fold(0.0, f64::max);
    assert!(backward.width() <= 1e-8);
    assert!(backward.upper >= support && backward.lower <= 1.0 - rc.min(rf));
    witness(&forward);
    witness(&backward);
    // The two-sided value is the larger side and meets both sides' bounds.
    let distance = result.distance;
    assert!(distance.upper >= (rf - rc).max(support) && distance.lower <= 1.0 - rc.min(rf));
}

/// The largest sampled distance from `from` to `to`, by exhaustive scan,
/// and the sampling radius: every point of `from` is that close to a sample.
fn sampled(from: &TriMesh, to: &TriMesh, steps: u32) -> (f64, f64) {
    let targets = triangles(to);
    let mut estimate = 0.0_f64;
    let mut radius = 0.0_f64;
    for [a, b, c] in triangles(from) {
        radius = radius.max(((b - a).length() + (c - a).length()) / f64::from(steps));
        for i in 0..=steps {
            for j in 0..=(steps - i) {
                let (u, v) = (
                    f64::from(i) / f64::from(steps),
                    f64::from(j) / f64::from(steps),
                );
                let point = a + (b - a) * u + (c - a) * v;
                let nearest = targets
                    .iter()
                    .map(|t| (point - closest_point_on_triangle(point, *t).unwrap()).length())
                    .fold(f64::INFINITY, f64::min);
                estimate = estimate.max(nearest);
            }
        }
    }
    (estimate, radius)
}

/// A deterministic bumpy height field.
fn terrain(seed: u64, n: u32, offset: Point3) -> TriMesh {
    let mut mesh = grid(0.0, 0.0, 1.0, 0.0, n);
    let mut state = seed;
    for p in &mut mesh.positions {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let noise = ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5;
        *p += offset + Point3::new(0.0, 0.0, 0.3 * noise + 0.2 * (4.0 * p.x).sin());
    }
    mesh
}

#[test]
fn the_interval_holds_every_densely_sampled_estimate() {
    for seed in 1..=6 {
        let a = terrain(seed, 4, Point3::ZERO);
        let b = terrain(seed * 7919, 5, Point3::new(0.2, -0.1, 0.1 * seed as f64));
        for (from, to) in [(&a, &b), (&b, &a)] {
            let bounds = one_sided_hausdorff(from, to, 1e-7).expect("valid meshes");
            let (estimate, radius) = sampled(from, to, 24);
            // Samples are points of `from`: none may exceed the upper bound.
            assert!(
                estimate <= bounds.upper + 1e-12,
                "seed {seed}: sampled {estimate} above upper {}",
                bounds.upper
            );
            // Every point is within `radius` of a sample and the distance is
            // 1-Lipschitz, so the true value is at most `estimate + radius`.
            assert!(
                bounds.lower <= estimate + radius,
                "seed {seed}: lower {} above sampled {estimate} + {radius}",
                bounds.lower
            );
            assert!(bounds.width() <= 1e-7);
            witness(&bounds);
        }
    }
}

#[test]
fn the_two_sided_distance_is_the_larger_side() {
    let a = terrain(3, 4, Point3::ZERO);
    let b = terrain(11, 3, Point3::new(0.5, 0.0, 0.4));
    let result = hausdorff_distance(&a, &b, 1e-8).expect("valid meshes");
    let larger = if result.forward.lower >= result.backward.lower {
        result.forward
    } else {
        result.backward
    };
    assert_eq!(result.distance.lower, larger.lower);
    assert_eq!(
        result.distance.upper,
        result.forward.upper.max(result.backward.upper)
    );
    assert!(result.forward.lower != result.backward.lower);
}

#[test]
fn identical_meshes_are_zero_apart() {
    let a = icosphere(1);
    let result = hausdorff_distance(&a, &a, ACCURACY).expect("valid meshes");
    check(&result.distance, 0.0, ACCURACY);
}

#[test]
fn a_coarse_request_stops_early_and_stays_sound() {
    let coarse = icosphere(0);
    let fine = icosphere(3);
    let loose = one_sided_hausdorff(&coarse, &fine, 0.05).expect("valid meshes");
    let tight = one_sided_hausdorff(&coarse, &fine, 1e-9).expect("valid meshes");
    assert!(loose.width() <= 0.05);
    assert!(loose.lower <= tight.upper && tight.lower <= loose.upper);
}

#[test]
fn unusable_input_is_refused() {
    let a = grid(0.0, 0.0, 1.0, 0.0, 1);
    let empty = TriMesh::new(vec![Point3::ZERO], vec![]);
    assert_eq!(
        one_sided_hausdorff(&a, &empty, 1e-6),
        Err(HausdorffError::EmptyMesh)
    );
    assert_eq!(
        one_sided_hausdorff(&empty, &a, 1e-6),
        Err(HausdorffError::EmptyMesh)
    );
    let broken = TriMesh::new(vec![Point3::ZERO; 3], vec![0, 1, 7]);
    assert_eq!(
        hausdorff_distance(&a, &broken, 1e-6),
        Err(HausdorffError::IndexOutOfRange)
    );
    let nan = TriMesh::new(
        vec![Point3::ZERO, Point3::X, Point3::new(f64::NAN, 0.0, 0.0)],
        vec![0, 1, 2],
    );
    assert_eq!(
        hausdorff_distance(&a, &nan, 1e-6),
        Err(HausdorffError::NonFiniteInput)
    );
    for accuracy in [-1.0, f64::NAN] {
        assert_eq!(
            hausdorff_distance(&a, &a, accuracy),
            Err(HausdorffError::InvalidAccuracy)
        );
    }
}

#[test]
fn degenerate_triangles_are_measured_as_their_edges() {
    // A zero-area sliver standing in for the segment from (0,0,0) to (2,0,0).
    let sliver = TriMesh::new(
        vec![
            Point3::ZERO,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        vec![0, 1, 2],
    );
    // The square [0,2] x [1,3] runs 1 beside the segment along its near
    // edge, and its far edge is 3 from it.
    let square = grid(0.0, 1.0, 2.0, 0.0, 2);
    let result = hausdorff_distance(&sliver, &square, ACCURACY).expect("valid");
    check(&result.forward, 1.0, ACCURACY);
    check(&result.backward, 3.0, ACCURACY);
}
