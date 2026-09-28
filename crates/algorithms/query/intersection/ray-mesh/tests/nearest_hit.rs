//! Narrow-phase ray/mesh behaviour, including the oracle fixtures #41 requires.
//!
//! Fixtures are chosen so every expected answer is checkable by hand: unit
//! triangles at known offsets, axis-aligned rays, and exact edge/vertex hits.
//! Keep covering edge-on, vertex-on, parallel-in-plane, back-face,
//! behind-origin and BVH composition: a suite that never produces an edge-on
//! hit proves nothing about tie-breaking.

use std::ops::ControlFlow;

use axiolid_core::{Aabb, Point3, Ray3, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_ray_mesh::{
    intersect_triangle, nearest_hit, nearest_hit_among, FaceSide, RayMeshError,
};
use axiolid_spatial::{Bvh, SpatialIndex, SpatialItem};

fn ray(origin: [f64; 3], direction: [f64; 3]) -> Ray3 {
    Ray3 {
        origin: Point3::from_array(origin),
        direction: Vec3::from_array(direction),
    }
}

/// Two stacked z-planes: triangle 0 at z = 1, triangle 1 at z = 2.
fn stacked_planes() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(4.0, 0.0, 1.0),
            Point3::new(0.0, 4.0, 1.0),
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(4.0, 0.0, 2.0),
            Point3::new(0.0, 4.0, 2.0),
        ],
        vec![0, 1, 2, 3, 4, 5],
    )
}

/// Two coplanar triangles sharing the diagonal edge from (1,0) to (0,1).
fn shared_edge_pair() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
        ],
        vec![0, 1, 2, 3, 2, 1],
    )
}

#[test]
fn nearest_hit_reports_the_closest_plane_not_the_first() {
    let mesh = stacked_planes();
    let hit = nearest_hit(
        &mesh,
        &ray([0.5, 0.5, -1.0], [0.0, 0.0, 1.0]),
        Tolerance::METRE,
    )
    .expect("valid query")
    .expect("the ray crosses both planes");

    assert_eq!(hit.triangle, 0);
    assert!((hit.t - 2.0).abs() < 1e-12, "t was {}", hit.t);
    assert!((hit.point.z - 1.0).abs() < 1e-12);
    assert_eq!(hit.side, FaceSide::Back);
}

#[test]
fn candidate_order_does_not_change_the_answer() {
    let mesh = stacked_planes();
    let query = ray([0.5, 0.5, -1.0], [0.0, 0.0, 1.0]);

    let forward = nearest_hit_among(&mesh, &query, Tolerance::METRE, [0, 1])
        .unwrap()
        .unwrap();
    let reversed = nearest_hit_among(&mesh, &query, Tolerance::METRE, [1, 0])
        .unwrap()
        .unwrap();
    let repeated = nearest_hit_among(&mesh, &query, Tolerance::METRE, [1, 0, 1, 0])
        .unwrap()
        .unwrap();

    assert_eq!(forward, reversed);
    assert_eq!(forward, repeated);
}

#[test]
fn unnormalised_direction_reports_t_in_direction_units() {
    let mesh = stacked_planes();
    let hit = nearest_hit(
        &mesh,
        &ray([0.5, 0.5, 0.0], [0.0, 0.0, 2.0]),
        Tolerance::METRE,
    )
    .unwrap()
    .unwrap();

    // The plane is 1.0 away, and the direction is twice unit length.
    assert!((hit.t - 0.5).abs() < 1e-12, "t was {}", hit.t);
    assert!((hit.point.z - 1.0).abs() < 1e-12);
}

#[test]
fn front_and_back_sides_are_certified_by_orientation() {
    let mesh = stacked_planes();

    let from_below = nearest_hit(
        &mesh,
        &ray([0.5, 0.5, -1.0], [0.0, 0.0, 1.0]),
        Tolerance::METRE,
    )
    .unwrap()
    .unwrap();
    let from_above = nearest_hit(
        &mesh,
        &ray([0.5, 0.5, 3.0], [0.0, 0.0, -1.0]),
        Tolerance::METRE,
    )
    .unwrap()
    .unwrap();

    assert_eq!(from_below.side, FaceSide::Back);
    assert_eq!(from_above.side, FaceSide::Front);
    assert_eq!(from_above.triangle, 1);
}

#[test]
fn a_ray_starting_in_the_triangle_plane_is_reported_coplanar() {
    let corners = [
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(4.0, 0.0, 1.0),
        Point3::new(0.0, 4.0, 1.0),
    ];
    let hit = intersect_triangle(
        &ray([0.5, 0.5, 1.0], [0.0, 0.0, 1.0]),
        corners,
        Tolerance::METRE,
        7,
    )
    .unwrap()
    .expect("a zero-distance hit is still a hit");

    assert_eq!(hit.side, FaceSide::Coplanar);
    assert_eq!(hit.triangle, 7);
    assert!(hit.t.abs() < 1e-12);
}

#[test]
fn edge_on_hits_do_not_fall_through_the_surface() {
    let mesh = shared_edge_pair();
    // Exactly on the shared diagonal edge.
    let hit = nearest_hit(
        &mesh,
        &ray([0.5, 0.5, 0.0], [0.0, 0.0, 1.0]),
        Tolerance::METRE,
    )
    .unwrap()
    .expect("an edge-on ray must hit the surface");

    // Deterministic tie-break: equal t resolves to the lower triangle index.
    assert_eq!(hit.triangle, 0);
    assert!((hit.t - 1.0).abs() < 1e-12);
}

#[test]
fn vertex_on_hits_are_retained() {
    let mesh = shared_edge_pair();
    let hit = nearest_hit(
        &mesh,
        &ray([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        Tolerance::METRE,
    )
    .unwrap()
    .expect("a vertex-on ray must hit the surface");

    assert_eq!(hit.triangle, 0);
    let [w, u, v] = hit.barycentric;
    assert!((w - 1.0).abs() < 1e-12, "w was {w}");
    assert!(u.abs() < 1e-12 && v.abs() < 1e-12);
}

#[test]
fn a_parallel_in_plane_ray_misses_rather_than_dividing_by_zero() {
    let corners = [
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(4.0, 0.0, 1.0),
        Point3::new(0.0, 4.0, 1.0),
    ];
    let miss = intersect_triangle(
        &ray([-1.0, 0.5, 1.0], [1.0, 0.0, 0.0]),
        corners,
        Tolerance::METRE,
        0,
    )
    .unwrap();

    assert!(miss.is_none());
}

#[test]
fn geometry_behind_the_origin_is_not_a_hit() {
    let mesh = stacked_planes();
    let behind = nearest_hit(
        &mesh,
        &ray([0.5, 0.5, 3.0], [0.0, 0.0, 1.0]),
        Tolerance::METRE,
    )
    .unwrap();

    assert!(behind.is_none());
}

#[test]
fn a_degenerate_triangle_refuses_instead_of_silently_missing() {
    let mesh = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(2.0, 0.0, 1.0),
        ],
        vec![0, 1, 2],
    );

    assert_eq!(
        nearest_hit(
            &mesh,
            &ray([0.5, 0.0, 0.0], [0.0, 0.0, 1.0]),
            Tolerance::METRE
        ),
        Err(RayMeshError::DegenerateTriangle { triangle: 0 })
    );
}

#[test]
fn invalid_rays_and_policies_are_typed_refusals() {
    let mesh = stacked_planes();

    assert_eq!(
        nearest_hit(
            &mesh,
            &ray([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            Tolerance::METRE
        ),
        Err(RayMeshError::ZeroDirection)
    );
    assert_eq!(
        nearest_hit(
            &mesh,
            &ray([f64::NAN, 0.0, 0.0], [0.0, 0.0, 1.0]),
            Tolerance::METRE
        ),
        Err(RayMeshError::NonFiniteInput)
    );
    assert_eq!(
        nearest_hit(
            &mesh,
            &ray([0.0, 0.0, 0.0], [0.0, 0.0, f64::INFINITY]),
            Tolerance::METRE
        ),
        Err(RayMeshError::NonFiniteInput)
    );
}

#[test]
fn an_out_of_range_triangle_index_is_reported_not_ignored() {
    struct BadMesh;
    impl axiolid_mesh::TriangleMeshView for BadMesh {
        fn position_count(&self) -> usize {
            1
        }
        fn position(&self, _index: usize) -> Point3 {
            Point3::ZERO
        }
        fn triangle_count(&self) -> usize {
            1
        }
        fn triangle(&self, _index: usize) -> [u64; 3] {
            [0, 1, 2]
        }
    }

    assert_eq!(
        nearest_hit(
            &BadMesh,
            &ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0]),
            Tolerance::METRE
        ),
        Err(RayMeshError::PositionIndexOutOfRange { triangle: 0 })
    );
}

#[test]
fn broad_phase_candidates_compose_with_the_narrow_phase() {
    let mesh = stacked_planes();
    let index = Bvh::build((0..mesh.triangle_count()).map(|triangle| {
        let mut bounds = Aabb::empty();
        for corner in axiolid_mesh::TriangleMeshView::triangle(&mesh, triangle) {
            bounds.extend(axiolid_mesh::TriangleMeshView::position(
                &mesh,
                corner as usize,
            ));
        }
        SpatialItem::new(triangle, bounds)
    }));

    let query = ray([0.5, 0.5, -1.0], [0.0, 0.0, 1.0]);
    let mut candidates = Vec::new();
    index.visit_ray(&query, &mut |hit| {
        candidates.push(*hit.key);
        ControlFlow::Continue(())
    });
    assert!(
        !candidates.is_empty(),
        "the broad phase must yield candidates"
    );

    let through_index = nearest_hit_among(&mesh, &query, Tolerance::METRE, candidates)
        .unwrap()
        .unwrap();
    let exhaustive = nearest_hit(&mesh, &query, Tolerance::METRE)
        .unwrap()
        .unwrap();

    assert_eq!(through_index, exhaustive);
}
