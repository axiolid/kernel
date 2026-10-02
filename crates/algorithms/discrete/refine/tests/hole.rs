//! Hole filling and fairing: planar holes fill planar, a hole cut from a
//! sphere is refilled close to the sphere, the patch closes the mesh with
//! consistent winding, and every refusal leaves the mesh unchanged.

use axiolid_core::{Point3, Vec3};
use axiolid_mesh::halfedge::{HalfedgeId, HalfedgeMesh, VertexId};
use axiolid_refine::{
    fair, fill_hole, FairError, FairOptions, FairingOrder, FairingWeights, HoleFillError,
    HoleFillOptions,
};

/// `n x n` squares of the unit grid in the plane `z = 0`, each split into
/// two counter-clockwise triangles, minus the faces whose centroid lies
/// within `radius` of the centre.
fn holed_grid(n: u32, radius: f64) -> HalfedgeMesh {
    let mut positions = Vec::new();
    for j in 0..=n {
        for i in 0..=n {
            positions.push(Point3::new(f64::from(i), f64::from(j), 0.0));
        }
    }
    let id = |i: u32, j: u32| j * (n + 1) + i;
    let centre = Point3::new(f64::from(n) / 2.0, f64::from(n) / 2.0, 0.0);
    let mut faces = Vec::new();
    for j in 0..n {
        for i in 0..n {
            for t in [
                [id(i, j), id(i + 1, j), id(i + 1, j + 1)],
                [id(i, j), id(i + 1, j + 1), id(i, j + 1)],
            ] {
                let c = t.iter().fold(Vec3::ZERO, |s, &v| s + positions[v as usize]) / 3.0;
                if (c - centre).length() > radius {
                    faces.push(t);
                }
            }
        }
    }
    compact(positions, &faces)
}

/// Drop unused positions and build.
fn compact(positions: Vec<Point3>, faces: &[[u32; 3]]) -> HalfedgeMesh {
    let mut map = vec![u32::MAX; positions.len()];
    let mut kept = Vec::new();
    let faces: Vec<[u32; 3]> = faces
        .iter()
        .map(|t| {
            t.map(|v| {
                if map[v as usize] == u32::MAX {
                    map[v as usize] = kept.len() as u32;
                    kept.push(positions[v as usize]);
                }
                map[v as usize]
            })
        })
        .collect();
    HalfedgeMesh::from_faces(kept, &faces).expect("manifold")
}

/// Icosphere of radius 1: an icosahedron subdivided `levels` times with
/// every new vertex pushed onto the sphere.
fn icosphere(levels: u32) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
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
        let mut midpoint = std::collections::BTreeMap::new();
        let mut next = Vec::new();
        for f in &faces {
            let mut mid = [0u32; 3];
            for k in 0..3 {
                let (a, b) = (f[k], f[(k + 1) % 3]);
                let key = (a.min(b), a.max(b));
                mid[k] = *midpoint.entry(key).or_insert_with(|| {
                    positions
                        .push(((positions[a as usize] + positions[b as usize]) / 2.0).normalize());
                    positions.len() as u32 - 1
                });
            }
            next.push([f[0], mid[0], mid[2]]);
            next.push([f[1], mid[1], mid[0]]);
            next.push([f[2], mid[2], mid[1]]);
            next.push(mid);
        }
        faces = next;
    }
    (positions, faces)
}

/// An icosphere with every face touching a vertex within `angle` radians
/// of the north pole removed.
fn holed_sphere(levels: u32, angle: f64) -> HalfedgeMesh {
    let (positions, faces) = icosphere(levels);
    let pole = Vec3::new(0.3, 0.2, 1.0).normalize();
    let kept: Vec<[u32; 3]> = faces
        .into_iter()
        .filter(|t| {
            t.iter()
                .all(|&v| positions[v as usize].dot(pole) < angle.cos())
        })
        .collect();
    compact(positions, &kept)
}

fn only_hole(mesh: &HalfedgeMesh) -> HalfedgeId {
    let loops = mesh.boundary_loops();
    assert_eq!(loops.len(), 1, "fixture has one hole");
    loops[0]
}

/// The hole of a holed grid: the boundary loop that is not the outer one.
fn inner_hole(mesh: &HalfedgeMesh) -> HalfedgeId {
    let loops = mesh.boundary_loops();
    assert_eq!(loops.len(), 2, "outer boundary and one hole");
    *loops
        .iter()
        .min_by_key(|&&h| mesh.loop_halfedges(h).count())
        .expect("two loops")
}

fn signed_volume(mesh: &HalfedgeMesh) -> f64 {
    mesh.faces()
        .map(|f| {
            let c: Vec<Point3> = mesh.face_vertices(f).map(|v| mesh.position(v)).collect();
            c[0].dot(c[1].cross(c[2])) / 6.0
        })
        .sum()
}

fn face_normal(mesh: &HalfedgeMesh, f: axiolid_mesh::halfedge::FaceId) -> Vec3 {
    let c: Vec<Point3> = mesh.face_vertices(f).map(|v| mesh.position(v)).collect();
    (c[1] - c[0]).cross(c[2] - c[0])
}

#[test]
fn a_planar_hole_fills_planar_exactly_and_consistently_wound() {
    let mut mesh = holed_grid(16, 4.0);
    let hole = inner_hole(&mesh);
    let before = mesh.clone();
    let report = fill_hole(&mut mesh, hole, &HoleFillOptions::default()).expect("fillable");
    mesh.validate().expect("valid halfedge mesh");
    assert_eq!(
        mesh.boundary_loops().len(),
        1,
        "only the outer boundary is left"
    );
    assert_eq!(
        mesh.euler_characteristic(),
        before.euler_characteristic() + 1
    );
    assert!(!report.vertices.is_empty(), "the hole was refined");
    assert!(report.fairing.is_some());
    for &v in &report.vertices {
        assert_eq!(mesh.position(v).z, 0.0, "vertex {v} left the plane");
    }
    for &f in &report.faces {
        let n = face_normal(&mesh, f);
        assert!(n.z > 0.0, "patch face {f} is wound against the grid");
    }
    // The patch covers exactly the removed area.
    let area: f64 = report
        .faces
        .iter()
        .map(|&f| face_normal(&mesh, f).length() / 2.0)
        .sum();
    let removed = 2.0 * 16.0 * 16.0 / 2.0 - before.face_count() as f64 / 2.0;
    assert!((area - removed).abs() <= 1e-9, "{area} vs {removed}");
    // Fixed vertices are untouched.
    for v in before.vertices() {
        assert_eq!(mesh.position(v), before.position(v));
    }
    // Refinement reached the grid's density: no patch edge much longer
    // than the grid's diagonal.
    for &f in &report.faces {
        let c: Vec<Point3> = mesh.face_vertices(f).map(|v| mesh.position(v)).collect();
        for k in 0..3 {
            assert!((c[k] - c[(k + 1) % 3]).length() <= 2.5, "long patch edge");
        }
    }
    assert!(report.max_dihedral <= 1e-6, "planar triangulation");
}

#[test]
fn a_tilted_planar_hole_stays_in_its_plane() {
    let mut mesh = holed_grid(16, 4.0);
    let rotation = axiolid_core::Mat3::from_axis_angle(Vec3::new(1.0, 2.0, 0.5).normalize(), 0.7);
    let offset = Vec3::new(1000.0, -2000.0, 500.0);
    for v in mesh.vertices().collect::<Vec<VertexId>>() {
        let p = rotation * mesh.position(v) + offset;
        mesh.set_position(v, p);
    }
    let normal = rotation * Vec3::Z;
    let hole = inner_hole(&mesh);
    let report = fill_hole(&mut mesh, hole, &HoleFillOptions::default()).expect("fillable");
    let worst = report
        .vertices
        .iter()
        .map(|&v| (mesh.position(v) - offset).dot(normal).abs())
        .fold(0.0, f64::max);
    // Hole diameter 8; coordinates around 2000.
    assert!(worst <= 1e-9, "out-of-plane deviation {worst}");
}

#[test]
fn a_hole_cut_from_a_sphere_is_refilled_close_to_the_sphere() {
    // Level-4 icosphere (edge ~0.08). A cap of angular radius `t` removed
    // from the unit sphere would sag by up to 1 - cos t if filled flat. The
    // biharmonic continuation that matches the sphere's height and slope
    // on the rim is the paraboloid z = a - r^2 / (2 cos t), whose apex
    // overshoots the sphere by a - 1 = cos t + sin^2 t / (2 cos t) - 1
    // (about t^4 / 8). The removed faces reach up to one edge (0.08) past
    // the cap, so the cotangent fill is held to the continuum overshoot of
    // a cap 0.08 wider. Uniform weights see connectivity only and distort
    // the irregular patch tangentially; they get half as much again.
    // Measured: cotangent 0.0177 (cap 0.6) and 0.0017 (cap 0.3), uniform
    // 0.0396 (cap 0.6).
    let continuum = |t: f64| t.cos() + t.sin().powi(2) / (2.0 * t.cos()) - 1.0;
    for (weights, cap, slack) in [
        (FairingWeights::Cotangent, 0.6, 1.0),
        (FairingWeights::Uniform, 0.6, 1.5),
        (FairingWeights::Cotangent, 0.3, 1.0),
    ] {
        let bound = slack * continuum(cap + 0.08);
        let mut mesh = holed_sphere(4, cap);
        let hole = only_hole(&mesh);
        let options = HoleFillOptions {
            fairing: Some(FairOptions {
                weights,
                ..FairOptions::default()
            }),
            ..HoleFillOptions::default()
        };
        let report = fill_hole(&mut mesh, hole, &options).expect("fillable");
        mesh.validate().expect("valid halfedge mesh");
        assert!(mesh.boundary_loops().is_empty(), "closed");
        assert_eq!(mesh.euler_characteristic(), 2);
        assert!(
            report.vertices.len() > 50,
            "{} new vertices",
            report.vertices.len()
        );
        let worst = report
            .vertices
            .iter()
            .map(|&v| (mesh.position(v).length() - 1.0).abs())
            .fold(0.0, f64::max);
        eprintln!(
            "{weights:?} cap {cap}: {} new vertices, max radial error {worst}, bound {bound}",
            report.vertices.len()
        );
        assert!(
            worst <= bound,
            "{weights:?} cap {cap}: max radial error {worst}"
        );
        // Outward winding everywhere: the volume is the sphere's.
        let volume = signed_volume(&mesh);
        let sphere = 4.0 / 3.0 * core::f64::consts::PI;
        assert!((volume - sphere).abs() <= 0.02 * sphere, "volume {volume}");
        for &f in &report.faces {
            let c: Vec<Point3> = mesh.face_vertices(f).map(|v| mesh.position(v)).collect();
            let centroid = (c[0] + c[1] + c[2]) / 3.0;
            assert!(
                face_normal(&mesh, f).dot(centroid) > 0.0,
                "face {f} points inward"
            );
        }
    }
}

#[test]
fn fairing_beats_the_flat_patch_and_harmonic_fairing_by_a_wide_margin() {
    let radial = |options: &HoleFillOptions| {
        let mut mesh = holed_sphere(4, 0.6);
        let hole = only_hole(&mesh);
        let report = fill_hole(&mut mesh, hole, options).expect("fillable");
        report
            .vertices
            .iter()
            .map(|&v| (mesh.position(v).length() - 1.0).abs())
            .fold(0.0, f64::max)
    };
    let flat = radial(&HoleFillOptions {
        fairing: None,
        ..HoleFillOptions::default()
    });
    let harmonic = radial(&HoleFillOptions {
        fairing: Some(FairOptions {
            order: FairingOrder::Harmonic,
            ..FairOptions::default()
        }),
        ..HoleFillOptions::default()
    });
    let biharmonic = radial(&HoleFillOptions::default());
    eprintln!("flat {flat}, harmonic {harmonic}, biharmonic {biharmonic}");
    assert!(flat > 0.1, "flat patch sags {flat}");
    // A membrane is pulled flat in the middle; a thin plate continues the
    // curvature.
    assert!(harmonic > 0.05, "harmonic {harmonic}");
    assert!(
        biharmonic < harmonic / 10.0,
        "biharmonic {biharmonic} vs harmonic {harmonic}"
    );
}

#[test]
fn triangulation_only_adds_no_vertex() {
    let mut mesh = holed_sphere(3, 0.5);
    let hole = only_hole(&mesh);
    let boundary = mesh.loop_halfedges(hole).count();
    let vertices = mesh.vertex_count();
    let report = fill_hole(
        &mut mesh,
        hole,
        &HoleFillOptions {
            refine: false,
            ..HoleFillOptions::default()
        },
    )
    .expect("fillable");
    assert_eq!(mesh.vertex_count(), vertices);
    assert_eq!(report.faces.len(), boundary - 2);
    assert!(report.vertices.is_empty() && report.fairing.is_none());
    assert!(mesh.boundary_loops().is_empty());
    mesh.validate().expect("valid");
    // The minimum-dihedral triangulation of a spherical cap bends gently.
    assert!(
        report.max_dihedral < 0.6,
        "max dihedral {}",
        report.max_dihedral
    );
}

#[test]
fn filling_is_deterministic() {
    let base = holed_sphere(3, 0.7);
    let hole = only_hole(&base);
    let mut first = base.clone();
    let mut second = base.clone();
    let a = fill_hole(&mut first, hole, &HoleFillOptions::default()).expect("fillable");
    let b = fill_hole(&mut second, hole, &HoleFillOptions::default()).expect("fillable");
    assert_eq!(a, b);
    assert_eq!(first, second);
}

#[test]
fn a_triangular_hole_becomes_one_triangle() {
    // Tetrahedron without one face.
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ];
    let mut mesh =
        HalfedgeMesh::from_faces(positions, &[[0u32, 1, 3], [1, 2, 3], [2, 0, 3]]).expect("ok");
    let hole = only_hole(&mesh);
    let report = fill_hole(&mut mesh, hole, &HoleFillOptions::default()).expect("fillable");
    assert_eq!(report.faces.len(), 1);
    assert!(report.vertices.is_empty());
    assert!(mesh.boundary_loops().is_empty());
    assert!(signed_volume(&mesh) > 0.0);
}

fn assert_refused_unchanged(
    mesh: &HalfedgeMesh,
    hole: HalfedgeId,
    options: &HoleFillOptions,
) -> HoleFillError {
    let mut work = mesh.clone();
    let error = fill_hole(&mut work, hole, options).expect_err("refused");
    assert_eq!(&work, mesh, "a refusal must leave the mesh unchanged");
    error
}

#[test]
fn refuses_halfedges_that_bound_no_hole() {
    let mesh = holed_grid(8, 2.0);
    let interior = mesh
        .halfedges()
        .find(|&h| !mesh.is_boundary_halfedge(h))
        .expect("a face halfedge");
    assert_eq!(
        assert_refused_unchanged(&mesh, interior, &HoleFillOptions::default()),
        HoleFillError::NotBoundary { halfedge: interior }
    );
    let missing = HalfedgeId::new(1_000_000);
    assert_eq!(
        assert_refused_unchanged(&mesh, missing, &HoleFillOptions::default()),
        HoleFillError::RemovedHalfedge { halfedge: missing }
    );
}

#[test]
fn refuses_bad_options_and_oversized_holes() {
    let mesh = holed_grid(8, 2.0);
    let hole = inner_hole(&mesh);
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            assert_refused_unchanged(
                &mesh,
                hole,
                &HoleFillOptions {
                    density,
                    ..HoleFillOptions::default()
                }
            ),
            HoleFillError::InvalidOption {
                name: "density",
                ..
            }
        ));
    }
    let n = mesh.loop_halfedges(hole).count();
    assert_eq!(
        assert_refused_unchanged(
            &mesh,
            hole,
            &HoleFillOptions {
                max_boundary_vertices: n - 1,
                ..HoleFillOptions::default()
            }
        ),
        HoleFillError::HoleTooLarge {
            boundary_vertices: n,
            limit: n - 1
        }
    );
    // At the limit the hole is accepted.
    let mut work = mesh.clone();
    fill_hole(
        &mut work,
        hole,
        &HoleFillOptions {
            max_boundary_vertices: n,
            ..HoleFillOptions::default()
        },
    )
    .expect("at the limit");
    assert_eq!(
        assert_refused_unchanged(
            &mesh,
            hole,
            &HoleFillOptions {
                max_new_vertices: 2,
                ..HoleFillOptions::default()
            }
        ),
        HoleFillError::BudgetExceeded { limit: 2 }
    );
}

#[test]
fn refuses_a_boundary_that_overlaps_itself_in_projection() {
    // A helical ramp of 1.5 turns: its boundary winds round the axis
    // one and a half times, so its projection covers itself.
    let steps = 36u32;
    let mut positions = Vec::new();
    for k in 0..=steps {
        let angle = f64::from(k) / f64::from(steps) * 3.0 * core::f64::consts::PI;
        let z = f64::from(k) * 0.05;
        positions.push(Point3::new(angle.cos(), angle.sin(), z));
        positions.push(Point3::new(2.0 * angle.cos(), 2.0 * angle.sin(), z));
    }
    let mut faces = Vec::new();
    for k in 0..steps {
        let (i0, o0, i1, o1) = (2 * k, 2 * k + 1, 2 * k + 2, 2 * k + 3);
        faces.push([i0, o0, o1]);
        faces.push([i0, o1, i1]);
    }
    let mesh = HalfedgeMesh::from_faces(positions, &faces).expect("manifold strip");
    let hole = only_hole(&mesh);
    assert!(matches!(
        assert_refused_unchanged(&mesh, hole, &HoleFillOptions::default()),
        HoleFillError::SelfIntersectingBoundary { .. }
    ));
}

#[test]
fn refuses_a_boundary_with_no_area() {
    // Two triangles folded flat onto each other: the loop encloses no area.
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
    ];
    let mesh = HalfedgeMesh::from_faces(positions, &[[0u32, 1, 2], [1, 0, 3]]).expect("ok");
    let hole = only_hole(&mesh);
    assert_eq!(
        assert_refused_unchanged(&mesh, hole, &HoleFillOptions::default()),
        HoleFillError::DegenerateBoundary
    );
}

#[test]
fn refuses_a_hole_every_triangulation_of_which_breaks_the_mesh() {
    // Quad hole a b c d with a, b, c collinear: diagonal a-c gives a
    // degenerate triangle, and diagonal b-d already exists.
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];
    let mesh = HalfedgeMesh::from_faces(positions, &[[2u32, 1, 3], [1, 0, 3]]).expect("ok");
    let hole = only_hole(&mesh);
    assert_eq!(
        assert_refused_unchanged(&mesh, hole, &HoleFillOptions::default()),
        HoleFillError::NoValidTriangulation
    );
}

#[test]
fn a_fairing_that_does_not_converge_is_refused() {
    let mesh = holed_sphere(4, 0.6);
    let hole = only_hole(&mesh);
    let options = HoleFillOptions {
        fairing: Some(FairOptions {
            max_iterations: Some(2),
            ..FairOptions::default()
        }),
        ..HoleFillOptions::default()
    };
    assert!(matches!(
        assert_refused_unchanged(&mesh, hole, &options),
        HoleFillError::Fair(FairError::NotConverged { coordinate: 0, .. })
    ));
}

// ---- fair on its own ------------------------------------------------------

/// Full `n x n` grid (no hole) and its interior vertices.
fn grid(n: u32) -> (HalfedgeMesh, Vec<VertexId>) {
    let mesh = holed_grid(n, -1.0);
    let interior = mesh
        .vertices()
        .filter(|&v| !mesh.is_boundary_vertex(v))
        .collect();
    (mesh, interior)
}

#[test]
fn fairing_flattens_a_bump_in_a_plane_exactly() {
    for order in [FairingOrder::Harmonic, FairingOrder::Biharmonic] {
        for weights in [FairingWeights::Cotangent, FairingWeights::Uniform] {
            let (mut mesh, interior) = grid(8);
            let before = mesh.clone();
            for &v in &interior {
                let p = mesh.position(v);
                mesh.set_position(v, p + Vec3::new(0.0, 0.0, (p.x * p.y).sin() * 0.3));
            }
            let options = FairOptions {
                order,
                weights,
                ..FairOptions::default()
            };
            let report = fair(&mut mesh, &interior, &options).expect("fairable");
            assert_eq!(report.vertices, interior.len());
            assert!(report.relative_residual <= 1e-12);
            for &v in &interior {
                assert_eq!(mesh.position(v).z, 0.0, "{order:?} {weights:?}");
            }
            // In-plane placement depends on the weights: cotangent weights
            // are taken from the bumpy surface, and the biharmonic system
            // sees the one-sided Laplacian of the border. The uniform
            // harmonic one is exact on the regular grid and puts every
            // vertex back.
            if (order, weights) != (FairingOrder::Harmonic, FairingWeights::Uniform) {
                continue;
            }
            for &v in &interior {
                let d = (mesh.position(v) - before.position(v)).length();
                assert!(d <= 1e-9, "{order:?} {weights:?}: moved {d}");
            }
        }
    }
}

#[test]
fn fairing_reproduces_a_quadratic_height_field_biharmonically() {
    // z = x^2 - y^2 is harmonic and so biharmonic too; with uniform weights
    // the discrete Laplacian of the regular grid is exact on it at every
    // interior vertex. Fair only vertices two rings in, so that every row
    // of the system is such a vertex.
    let (mut mesh, _) = grid(10);
    let interior: Vec<VertexId> = mesh
        .vertices()
        .filter(|&v| {
            let p = mesh.position(v);
            (2.0..=8.0).contains(&p.x) && (2.0..=8.0).contains(&p.y)
        })
        .collect();
    for v in mesh.vertices().collect::<Vec<_>>() {
        let p = mesh.position(v);
        mesh.set_position(
            v,
            Point3::new(p.x, p.y, (p.x - 5.0).powi(2) - (p.y - 5.0).powi(2)),
        );
    }
    let exact = mesh.clone();
    for &v in &interior {
        let p = mesh.position(v);
        mesh.set_position(v, Point3::new(p.x, p.y, 0.0));
    }
    let options = FairOptions {
        weights: FairingWeights::Uniform,
        ..FairOptions::default()
    };
    fair(&mut mesh, &interior, &options).expect("fairable");
    for &v in &interior {
        let d = (mesh.position(v) - exact.position(v)).length();
        assert!(d <= 1e-9, "vertex {v} off by {d}");
    }
}

#[test]
fn fair_refuses_by_name_and_leaves_the_mesh_unchanged() {
    let (mesh, interior) = grid(4);
    let check = |vertices: &[VertexId], options: &FairOptions| {
        let mut work = mesh.clone();
        let error = fair(&mut work, vertices, options).expect_err("refused");
        assert_eq!(work, mesh);
        error
    };
    let everything: Vec<VertexId> = mesh.vertices().collect();
    assert!(matches!(
        check(&everything, &FairOptions::default()),
        FairError::Unanchored { .. }
    ));
    let missing = VertexId::new(9999);
    assert_eq!(
        check(&[missing], &FairOptions::default()),
        FairError::RemovedVertex { vertex: missing }
    );
    for tolerance in [-1.0, f64::NAN] {
        assert!(matches!(
            check(
                &interior,
                &FairOptions {
                    relative_tolerance: tolerance,
                    ..FairOptions::default()
                }
            ),
            FairError::InvalidTolerance(_)
        ));
    }
    // Empty selection is a no-op, not an error.
    let mut work = mesh.clone();
    let report = fair(&mut work, &[], &FairOptions::default()).expect("nothing to do");
    assert_eq!(report.vertices, 0);
    assert_eq!(work, mesh);
}

#[test]
fn cotangent_fairing_refuses_degenerate_and_non_triangular_faces() {
    // A fan around vertex 0 in which triangle 0-5-1 is degenerate: vertex 5
    // lies on the segment from 0 to 1.
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(0.5, 0.0, 0.0),
    ];
    let degenerate = HalfedgeMesh::from_faces(
        positions,
        &[[0u32, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1]],
    )
    .expect("ok");
    let mut work = degenerate.clone();
    assert!(matches!(
        fair(&mut work, &[VertexId::new(0)], &FairOptions::default()),
        Err(FairError::DegenerateTriangle { .. })
    ));
    assert_eq!(work, degenerate);
    // Uniform weights do not look at the shape.
    let mut work = degenerate.clone();
    fair(
        &mut work,
        &[VertexId::new(0)],
        &FairOptions {
            weights: FairingWeights::Uniform,
            ..FairOptions::default()
        },
    )
    .expect("uniform weights ignore shape");

    let quads = HalfedgeMesh::from_faces(
        (0..9)
            .map(|k| Point3::new(f64::from(k % 3), f64::from(k / 3), 0.0))
            .collect(),
        &[[0u32, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]],
    )
    .expect("ok");
    let mut work = quads.clone();
    assert!(matches!(
        fair(&mut work, &[VertexId::new(4)], &FairOptions::default()),
        Err(FairError::NotATriangle { degree: 4, .. })
    ));
    let mut work = quads.clone();
    fair(
        &mut work,
        &[VertexId::new(4)],
        &FairOptions {
            weights: FairingWeights::Uniform,
            ..FairOptions::default()
        },
    )
    .expect("uniform weights work on quads");
    assert_eq!(work.position(VertexId::new(4)), Point3::new(1.0, 1.0, 0.0));
}

#[test]
fn a_refined_planar_patch_is_delaunay() {
    // In the plane the flip criterion is the Delaunay one: after filling,
    // every interior patch edge has opposite angles summing to at most pi.
    let mut mesh = holed_grid(16, 5.0);
    let hole = inner_hole(&mesh);
    let report = fill_hole(&mut mesh, hole, &HoleFillOptions::default()).expect("fillable");
    let in_patch = |f: Option<axiolid_mesh::halfedge::FaceId>| {
        f.is_some_and(|f| report.faces.binary_search(&f).is_ok())
    };
    let mut checked = 0;
    for e in mesh.edges() {
        let h = mesh.edge_halfedge(e, 0);
        let g = mesh.opposite(h);
        if !(in_patch(mesh.face(h)) && in_patch(mesh.face(g))) {
            continue;
        }
        let a = mesh.position(mesh.source(h));
        let b = mesh.position(mesh.target(h));
        let c = mesh.position(mesh.target(mesh.next(h)));
        let d = mesh.position(mesh.target(mesh.next(g)));
        let angle = |apex: Point3| (a - apex).angle_between(b - apex);
        assert!(
            angle(c) + angle(d) <= core::f64::consts::PI + 1e-9,
            "edge {e} is not locally Delaunay"
        );
        checked += 1;
    }
    assert!(checked > 50, "{checked} interior patch edges");
}

#[test]
fn the_vertex_budget_is_exact() {
    let mesh = holed_sphere(3, 0.7);
    let hole = only_hole(&mesh);
    let mut work = mesh.clone();
    let needed = fill_hole(&mut work, hole, &HoleFillOptions::default())
        .expect("fillable")
        .vertices
        .len();
    let mut work = mesh.clone();
    fill_hole(
        &mut work,
        hole,
        &HoleFillOptions {
            max_new_vertices: needed,
            ..HoleFillOptions::default()
        },
    )
    .expect("exactly at the budget");
    assert_eq!(
        assert_refused_unchanged(
            &mesh,
            hole,
            &HoleFillOptions {
                max_new_vertices: needed - 1,
                ..HoleFillOptions::default()
            }
        ),
        HoleFillError::BudgetExceeded { limit: needed - 1 }
    );
}

/// Every triangulation of the polygon `0..n` (indices into its loop).
fn all_triangulations(i: usize, k: usize) -> Vec<Vec<[usize; 3]>> {
    if k - i < 2 {
        return vec![Vec::new()];
    }
    let mut out = Vec::new();
    for m in (i + 1)..k {
        for left in all_triangulations(i, m) {
            for right in all_triangulations(m, k) {
                let mut t = left.clone();
                t.extend(right.iter().copied());
                t.push([i, m, k]);
                out.push(t);
            }
        }
    }
    out
}

#[test]
fn the_triangulation_minimises_the_largest_dihedral_angle() {
    // Brute force over every triangulation of small holes in a bumpy
    // surface: the reported largest dihedral angle is the true minimum, the
    // patch realises it, and on some holes it is strictly smaller than that
    // of the least-area triangulation (so the angle, not the area, decided).
    let mut beat_least_area = false;
    let height = |x: f64, y: f64| 0.6 * (1.1 * x).sin() * (0.7 * y).cos() + 0.08 * x * x;
    for centre in [(3.0, 3.0), (2.5, 3.5), (3.5, 2.0)] {
        let mut mesh = holed_grid(7, -1.0);
        for v in mesh.vertices().collect::<Vec<_>>() {
            let p = mesh.position(v);
            mesh.set_position(v, Point3::new(p.x, p.y, height(p.x, p.y)));
        }
        // Rebuild with the faces near `centre` removed.
        let tri = mesh.to_tri_mesh().expect("triangles");
        let faces: Vec<[u32; 3]> = tri
            .triangles()
            .filter(|t| {
                let c = t
                    .iter()
                    .fold(Vec3::ZERO, |s, &v| s + tri.positions[v as usize])
                    / 3.0;
                (c.x - centre.0).hypot(c.y - centre.1) > 1.0
            })
            .collect();
        let mesh = compact(tri.positions.clone(), &faces);
        let hole = inner_hole(&mesh);
        let ring: Vec<HalfedgeId> = mesh.loop_halfedges(hole).collect();
        let n = ring.len();
        assert!((5..=12).contains(&n), "hole of {n} vertices");
        let corner: Vec<VertexId> = ring.iter().map(|&h| mesh.source(h)).collect();
        let p: Vec<Point3> = corner.iter().map(|&v| mesh.position(v)).collect();
        let outside: Vec<Vec3> = ring
            .iter()
            .map(|&h| {
                let f = mesh.face(mesh.opposite(h)).expect("a face across");
                face_normal(&mesh, f).normalize()
            })
            .collect();
        let normal = |t: [usize; 3]| (p[t[1]] - p[t[0]]).cross(p[t[2]] - p[t[0]]).normalize();
        let mut best = f64::INFINITY;
        let mut least_area = (f64::INFINITY, 0.0);
        'triangulations: for t in all_triangulations(0, n - 1) {
            let mut sides: std::collections::BTreeMap<(usize, usize), Vec3> =
                std::collections::BTreeMap::new();
            let mut worst: f64 = 0.0;
            let mut area = 0.0;
            for &tri in &t {
                area += (p[tri[1]] - p[tri[0]])
                    .cross(p[tri[2]] - p[tri[0]])
                    .length()
                    / 2.0;
                for k in 0..3 {
                    let (a, b) = (tri[k], tri[(k + 1) % 3]);
                    let chord = (a.min(b), a.max(b));
                    let boundary = chord.1 - chord.0 == 1 || chord == (0, n - 1);
                    if !boundary && mesh.find_halfedge(corner[a], corner[b]).is_some() {
                        continue 'triangulations;
                    }
                    if boundary {
                        let side = if chord == (0, n - 1) { n - 1 } else { chord.0 };
                        worst = worst.max(normal(tri).angle_between(outside[side]));
                    } else if let Some(other) = sides.insert(chord, normal(tri)) {
                        worst = worst.max(normal(tri).angle_between(other));
                    }
                }
            }
            best = best.min(worst);
            if area < least_area.0 {
                least_area = (area, worst);
            }
        }
        // The optimum does not depend on where the loop starts, so neither
        // may the answer: start from every boundary halfedge in turn, which
        // moves the closing edge of the recurrence round the hole.
        for &start in &ring {
            let mut work = mesh.clone();
            let report = fill_hole(
                &mut work,
                start,
                &HoleFillOptions {
                    refine: false,
                    ..HoleFillOptions::default()
                },
            )
            .expect("fillable");
            let reported = report.max_dihedral;
            assert!(
                (reported - best).abs() <= 1e-9,
                "centre {centre:?}, start {start}: reported {reported} vs brute force {best}"
            );
            beat_least_area |= reported < least_area.1 - 1e-3;
            // And the patch in the mesh realises it.
            let mut realised: f64 = 0.0;
            for &f in &report.faces {
                for h in work.face_halfedges(f).collect::<Vec<_>>() {
                    let g = work.face(work.opposite(h)).expect("closed hole");
                    realised = realised.max(
                        face_normal(&work, f)
                            .normalize()
                            .angle_between(face_normal(&work, g).normalize()),
                    );
                }
            }
            assert!(
                (realised - reported).abs() <= 1e-9,
                "realised {realised} vs {reported}"
            );
        }
    }
    assert!(
        beat_least_area,
        "no fixture separates the angle from the area"
    );
}
