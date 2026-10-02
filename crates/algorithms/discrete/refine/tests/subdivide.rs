//! Subdivision surfaces: topology of each level, known limit positions,
//! invariance of the limit masks across levels, volume convergence,
//! boundary rules, affine invariance and refusals.

use axiolid_core::{Point3, Vec3};
use axiolid_mesh::halfedge::{HalfedgeMesh, VertexId};
use axiolid_refine::{
    limit_positions, subdivide, SubdivisionError, SubdivisionOptions, SubdivisionScheme,
};

const LOOP: SubdivisionScheme = SubdivisionScheme::Loop;
const CC: SubdivisionScheme = SubdivisionScheme::CatmullClark;

fn levels(scheme: SubdivisionScheme, n: u32) -> SubdivisionOptions {
    SubdivisionOptions::new(scheme, n)
}

/// Regular tetrahedron centred at the origin, outward wound.
fn tetrahedron() -> HalfedgeMesh {
    let positions = vec![
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(1.0, -1.0, -1.0),
        Point3::new(-1.0, 1.0, -1.0),
        Point3::new(-1.0, -1.0, 1.0),
    ];
    HalfedgeMesh::from_faces(positions, &[[0u32, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]])
        .expect("closed")
}

fn icosahedron() -> HalfedgeMesh {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let positions = [
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
    let faces: [[u32; 3]; 20] = [
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
    HalfedgeMesh::from_faces(positions, &faces).expect("closed")
}

/// Cube with corners `(+-1, +-1, +-1)`, outward wound quads.
fn cube() -> HalfedgeMesh {
    let positions = (0..8)
        .map(|k| {
            let s = |bit: u32| if k & bit == 0 { -1.0 } else { 1.0 };
            Point3::new(s(1), s(2), s(4))
        })
        .collect();
    let faces: [[u32; 4]; 6] = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    HalfedgeMesh::from_faces(positions, &faces).expect("closed")
}

/// Open triangle grid over `[0, n]^2` with height `z(x, y)`.
fn grid(n: u32, z: impl Fn(f64, f64) -> f64) -> HalfedgeMesh {
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            let (x, y) = (f64::from(i), f64::from(j));
            positions.push(Point3::new(x, y, z(x, y)));
        }
    }
    let id = |i: u32, j: u32| j * (n + 1) + i;
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            faces.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
            faces.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
        }
    }
    HalfedgeMesh::from_faces(positions, &faces).expect("open grid")
}

fn volume(mesh: &HalfedgeMesh) -> f64 {
    mesh.faces()
        .map(|f| {
            let c: Vec<Point3> = mesh.face_vertices(f).map(|v| mesh.position(v)).collect();
            (1..c.len() - 1)
                .map(|k| c[0].dot(c[k].cross(c[k + 1])) / 6.0)
                .sum::<f64>()
        })
        .sum()
}

fn v(i: u32) -> VertexId {
    VertexId::new(i)
}

#[test]
fn loop_preserves_topology_and_counts() {
    let mut mesh = icosahedron();
    for level in 1..=4 {
        let next = subdivide(&mesh, &levels(LOOP, 1)).expect("triangles");
        next.validate().expect("valid");
        assert_eq!(next.vertex_count(), mesh.vertex_count() + mesh.edge_count());
        assert_eq!(next.face_count(), 4 * mesh.face_count());
        assert_eq!(
            next.edge_count(),
            2 * mesh.edge_count() + 3 * mesh.face_count()
        );
        assert_eq!(next.euler_characteristic(), 2, "level {level}");
        assert!(next.boundary_loops().is_empty());
        mesh = next;
    }
    assert_eq!(mesh.face_count(), 20 * 256);
    assert!(volume(&mesh) > 0.0, "winding kept outward");
}

#[test]
fn loop_converges_to_the_known_limit_of_a_tetrahedron_corner() {
    // Valence 3: beta = 3/16 and the limit mask is 2/5 v + 1/5 sum(nbrs).
    // The other three corners sum to -v, so the limit is v / 5.
    let mesh = tetrahedron();
    let limit = limit_positions(&mesh, LOOP).expect("triangles");
    for k in 0..4 {
        let expected = mesh.position(v(k)) / 5.0;
        assert!((limit[k as usize] - expected).length() <= 1e-15);
    }
    let mut errors = Vec::new();
    let mut current = mesh.clone();
    for _ in 0..8 {
        current = subdivide(&current, &levels(LOOP, 1)).expect("triangles");
        errors.push((current.position(v(0)) - mesh.position(v(0)) / 5.0).length());
    }
    // Geometric convergence: the subdominant eigenvalue at valence 3 is
    // 5/8 - beta = 7/16... well under 1/2 per level.
    for pair in errors.windows(2) {
        assert!(pair[1] <= 0.5 * pair[0], "{errors:?}");
    }
    assert!(errors[7] <= 1e-3, "{errors:?}");
}

#[test]
fn loop_limit_positions_are_invariant_across_levels() {
    let mesh = icosahedron();
    let limit0 = limit_positions(&mesh, LOOP).expect("triangles");
    let mut current = mesh.clone();
    for _ in 0..3 {
        current = subdivide(&current, &levels(LOOP, 1)).expect("triangles");
        let limit = limit_positions(&current, LOOP).expect("triangles");
        for k in 0..12 {
            assert!(
                (limit[k] - limit0[k]).length() <= 1e-14,
                "vertex {k}: {:?} vs {:?}",
                limit[k],
                limit0[k]
            );
        }
    }
    // Valence 5 on the unit icosahedron: the limit is a fixed shrink.
    let ratio = limit0[0].length();
    for p in &limit0 {
        assert!((p.length() - ratio).abs() <= 1e-14);
    }
    // Neighbours of a unit icosahedron vertex v have v . u = 1 / sqrt 5, so
    // the limit is v (1 - 5 gamma + 5 gamma / sqrt 5) with Loop's gamma.
    let c = 3.0 / 8.0 + 0.25 * (2.0 * core::f64::consts::PI / 5.0).cos();
    let beta = (5.0 / 8.0 - c * c) / 5.0;
    let gamma = 1.0 / (3.0 / (8.0 * beta) + 5.0);
    let expected = 1.0 - 5.0 * gamma + 5.0 * gamma / 5f64.sqrt();
    assert!((ratio - expected).abs() <= 1e-14, "{ratio} vs {expected}");
}

#[test]
fn loop_volume_converges_geometrically() {
    let mut mesh = icosahedron();
    let mut volumes = vec![volume(&mesh)];
    for _ in 0..5 {
        mesh = subdivide(&mesh, &levels(LOOP, 1)).expect("triangles");
        volumes.push(volume(&mesh));
    }
    let steps: Vec<f64> = volumes.windows(2).map(|w| w[1] - w[0]).collect();
    // Loop shrinks a convex mesh; each change is about a quarter of the one
    // before (the 1/4 eigenvalue of the curvature terms).
    for pair in steps.windows(2) {
        assert!(pair[0] < 0.0 && pair[1] < 0.0, "{volumes:?}");
        assert!(pair[1].abs() <= 0.3 * pair[0].abs(), "{steps:?}");
    }
    let last = steps[4] / steps[3];
    assert!((last - 0.25).abs() <= 0.01, "{steps:?}");
}

#[test]
fn loop_follows_the_cubic_b_spline_on_a_boundary() {
    // A wavy open grid: boundary vertices must depend only on boundary
    // neighbours, and boundary edge points are midpoints.
    let bumpy = grid(4, |x, y| (x * 1.3).sin() * (y * 0.7).cos());
    let mut flat_inside = bumpy.clone();
    for k in flat_inside.vertices().collect::<Vec<_>>() {
        if !flat_inside.is_boundary_vertex(k) {
            let p = flat_inside.position(k);
            flat_inside.set_position(k, p + Vec3::new(0.3, -0.2, 5.0));
        }
    }
    let a = subdivide(&bumpy, &levels(LOOP, 2)).expect("triangles");
    let b = subdivide(&flat_inside, &levels(LOOP, 2)).expect("triangles");
    let mut boundary = 0;
    for k in a.vertices() {
        if a.is_boundary_vertex(k) {
            assert_eq!(a.position(k), b.position(k), "boundary vertex {k}");
            boundary += 1;
        }
    }
    assert_eq!(boundary, 4 * 4 * 4, "boundary vertices after two levels");
    // One level: old boundary vertex gets 3/4 v + 1/8 (b1 + b2).
    let one = subdivide(&bumpy, &levels(LOOP, 1)).expect("triangles");
    let p = |k: u32| bumpy.position(v(k));
    let expected = 0.75 * p(2) + 0.125 * (p(1) + p(3));
    assert!((one.position(v(2)) - expected).length() <= 1e-15);
    // The boundary converges to the B-spline limit (4 v + b1 + b2) / 6.
    let limit = limit_positions(&bumpy, LOOP).expect("triangles");
    assert!((limit[2] - (4.0 * p(2) + p(1) + p(3)) / 6.0).length() <= 1e-15);
    let deep = subdivide(&bumpy, &levels(LOOP, 7)).expect("triangles");
    assert!((deep.position(v(2)) - limit[2]).length() <= 1e-4);
}

#[test]
fn loop_keeps_a_plane_planar_and_commutes_with_affine_maps() {
    let flat = grid(3, |_, _| 0.0);
    let sub = subdivide(&flat, &levels(LOOP, 3)).expect("triangles");
    assert!(sub.vertices().all(|k| sub.position(k).z == 0.0));

    let mesh = icosahedron();
    let map = |p: Point3| {
        axiolid_core::Mat3::from_cols(
            Vec3::new(2.0, 0.1, 0.0),
            Vec3::new(0.3, 1.0, -0.2),
            Vec3::new(0.0, 0.5, 0.7),
        ) * p
            + Vec3::new(10.0, -3.0, 4.0)
    };
    let mut moved = mesh.clone();
    for k in moved.vertices().collect::<Vec<_>>() {
        moved.set_position(k, map(moved.position(k)));
    }
    let a = subdivide(&mesh, &levels(LOOP, 2)).expect("triangles");
    let b = subdivide(&moved, &levels(LOOP, 2)).expect("triangles");
    for k in a.vertices() {
        assert!((map(a.position(k)) - b.position(k)).length() <= 1e-12);
    }
}

#[test]
fn loop_is_deterministic_and_zero_levels_copy() {
    let mesh = icosahedron();
    assert_eq!(subdivide(&mesh, &levels(LOOP, 0)).expect("copy"), mesh);
    assert_eq!(
        subdivide(&mesh, &levels(LOOP, 2)).expect("ok"),
        subdivide(&mesh, &levels(LOOP, 2)).expect("ok")
    );
}

#[test]
fn loop_refuses_polygons_and_budget_overruns() {
    assert_eq!(
        subdivide(&cube(), &levels(LOOP, 1)),
        Err(SubdivisionError::NotATriangle {
            face: axiolid_mesh::halfedge::FaceId::new(0),
            degree: 4
        })
    );
    assert!(matches!(
        limit_positions(&cube(), LOOP),
        Err(SubdivisionError::NotATriangle { degree: 4, .. })
    ));
    let options = SubdivisionOptions {
        max_faces: 20 * 16 - 1,
        ..levels(LOOP, 2)
    };
    assert_eq!(
        subdivide(&icosahedron(), &options),
        Err(SubdivisionError::BudgetExceeded {
            faces: 320,
            limit: 319
        })
    );
    let at_limit = SubdivisionOptions {
        max_faces: 320,
        ..levels(LOOP, 2)
    };
    assert_eq!(
        subdivide(&icosahedron(), &at_limit)
            .expect("at the limit")
            .face_count(),
        320
    );
}

#[test]
fn catmull_clark_turns_any_mesh_into_quads() {
    let tet = subdivide(&tetrahedron(), &levels(CC, 1)).expect("polygons");
    tet.validate().expect("valid");
    assert_eq!(tet.face_count(), 12);
    assert_eq!(tet.vertex_count(), 4 + 6 + 4);
    assert!(tet.faces().all(|f| tet.face_degree(f) == 4));
    assert_eq!(tet.euler_characteristic(), 2);

    let mut mesh = cube();
    for _ in 0..3 {
        let next = subdivide(&mesh, &levels(CC, 1)).expect("quads");
        assert_eq!(
            next.vertex_count(),
            mesh.vertex_count() + mesh.edge_count() + mesh.face_count()
        );
        assert_eq!(next.face_count(), 4 * mesh.face_count());
        assert_eq!(next.euler_characteristic(), 2);
        assert!(volume(&next) > 0.0);
        mesh = next;
    }
}

#[test]
fn catmull_clark_converges_to_the_known_limit_of_a_cube_corner() {
    // Valence 3: limit = (9 v + 4 sum(edge nbrs) + sum(diagonals)) / 24.
    // For corner (1,1,1) that is (9 + 4 - 1) / 24 = 1/2 per coordinate.
    let mesh = cube();
    let limit = limit_positions(&mesh, CC).expect("quads");
    for k in 0..8u32 {
        let expected = mesh.position(v(k)) / 2.0;
        assert!((limit[k as usize] - expected).length() <= 1e-15, "{k}");
    }
    let mut current = mesh.clone();
    let mut errors = Vec::new();
    for _ in 0..7 {
        current = subdivide(&current, &levels(CC, 1)).expect("quads");
        let corner = current.position(v(7));
        errors.push((corner - Point3::splat(0.5)).length());
        let l = limit_positions(&current, CC).expect("quads");
        assert!(
            (l[7] - Point3::splat(0.5)).length() <= 1e-14,
            "limit invariant"
        );
    }
    for pair in errors.windows(2) {
        assert!(pair[1] <= 0.6 * pair[0], "{errors:?}");
    }
    assert!(errors[6] <= 1e-2, "{errors:?}");
}

#[test]
fn catmull_clark_volume_converges() {
    let mut mesh = cube();
    let mut volumes = vec![volume(&mesh)];
    for _ in 0..5 {
        mesh = subdivide(&mesh, &levels(CC, 1)).expect("quads");
        volumes.push(volume(&mesh));
    }
    let steps: Vec<f64> = volumes.windows(2).map(|w| w[1] - w[0]).collect();
    for pair in steps.windows(2) {
        assert!(pair[1].abs() <= 0.35 * pair[0].abs(), "{steps:?}");
    }
}

#[test]
fn catmull_clark_rules_on_one_level() {
    let mesh = cube();
    let one = subdivide(&mesh, &levels(CC, 1)).expect("quads");
    // Face points come last, in face order: face 0 is z = -1.
    let face0 = one.position(v(8 + 12));
    assert_eq!(face0, Point3::new(0.0, 0.0, -1.0));
    // Corner (1,1,1): Q = mean of face points (1,0,0), (0,1,0), (0,0,1)
    // = 1/3, R = mean of edge midpoints = 2/3: (Q + 2R + 0 v) / 3 = 5/9.
    let corner = one.position(v(7));
    assert!(
        (corner - Point3::splat(5.0 / 9.0)).length() <= 1e-15,
        "{corner:?}"
    );
}

#[test]
fn catmull_clark_boundary_follows_the_b_spline() {
    let quads = |z: f64| {
        let positions: Vec<Point3> = (0..16)
            .map(|k| {
                let (i, j) = (f64::from(k % 4), f64::from(k / 4));
                let inner = (1..3).contains(&(k % 4)) && (1..3).contains(&(k / 4));
                Point3::new(i, j, if inner { z } else { (i * 0.9).sin() })
            })
            .collect();
        let mut faces = Vec::new();
        for j in 0..3u32 {
            for i in 0..3u32 {
                let a = j * 4 + i;
                faces.push([a, a + 1, a + 5, a + 4]);
            }
        }
        HalfedgeMesh::from_faces(positions, &faces).expect("open quads")
    };
    let a = subdivide(&quads(0.0), &levels(CC, 2)).expect("quads");
    let b = subdivide(&quads(3.0), &levels(CC, 2)).expect("quads");
    for k in a.vertices() {
        if a.is_boundary_vertex(k) {
            assert_eq!(a.position(k), b.position(k));
        }
    }
    let flat = subdivide(&quads(0.0), &levels(CC, 1)).expect("quads");
    let base = quads(0.0);
    let p = |k: u32| base.position(v(k));
    assert!((flat.position(v(1)) - (0.75 * p(1) + 0.125 * (p(0) + p(2)))).length() <= 1e-15);
    assert!(matches!(
        limit_positions(&tetrahedron(), CC),
        Err(SubdivisionError::NotAQuad { degree: 3, .. })
    ));
}
