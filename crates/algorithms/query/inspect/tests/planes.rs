//! Planar regions of meshes (#131): every region's certified deviation
//! bounds its corners' actual distances, and stays within the tolerance.

use axiolid_core::{Point3, Vec3};
use axiolid_inspect::{detect_planes, DetectedPlane, PlaneError, PlaneTolerance};
use axiolid_mesh::TriMesh;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn cuboid(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions = vec![
        p(x0, y0, z0),
        p(x1, y0, z0),
        p(x1, y1, z0),
        p(x0, y1, z0),
        p(x0, y0, z1),
        p(x1, y0, z1),
        p(x1, y1, z1),
        p(x0, y1, z1),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Quads as pairs of triangles over a shared vertex list.
fn quads(positions: Vec<Point3>, quads: &[[u32; 4]]) -> TriMesh {
    let indices = quads
        .iter()
        .flat_map(|q| [q[0], q[1], q[2], q[0], q[2], q[3]])
        .collect();
    TriMesh::new(positions, indices)
}

/// A flight of `n` steps, `going` deep and `rise` high, `width` wide:
/// treads and risers as one strip.
fn stair(n: usize, going: f64, rise: f64, width: f64) -> TriMesh {
    let mut positions = Vec::new();
    let mut faces = Vec::new();
    // Profile corners along x-z: up a riser, along a tread, ...
    let mut profile = vec![(0.0, 0.0)];
    for i in 0..n {
        let (x, z) = (i as f64 * going, (i + 1) as f64 * rise);
        profile.push((x, z));
        profile.push((x + going, z));
    }
    for &(x, z) in &profile {
        positions.push(p(x, 0.0, z));
        positions.push(p(x, width, z));
    }
    for k in 0..profile.len() as u32 - 1 {
        let (a, b) = (2 * k, 2 * k + 2);
        faces.push([a, b, b + 1, a + 1]);
    }
    quads(positions, &faces)
}

/// Every member corner lies within the region's certified deviation of
/// its plane, the deviation within the tolerance, and each triangle of
/// area in exactly one region.
fn check(mesh: &TriMesh, planes: &[DetectedPlane], distance: f64) {
    let mut seen = vec![0usize; mesh.indices.len() / 3];
    for plane in planes {
        let n = plane.normal / plane.normal.length();
        let reach = mesh
            .positions
            .iter()
            .map(|q| q.abs().max_element())
            .fold(0.0, f64::max);
        assert!(
            plane.deviation <= distance + 64.0 * f64::EPSILON * reach,
            "{plane:?}"
        );
        for &t in &plane.triangles {
            seen[t as usize] += 1;
            for k in 0..3 {
                let v = mesh.positions[mesh.indices[3 * t as usize + k] as usize];
                let d = n.dot(v - plane.point).abs();
                assert!(
                    d <= plane.deviation * (1.0 + 1e-12),
                    "{d} > {}",
                    plane.deviation
                );
            }
        }
    }
    assert!(seen.iter().all(|&c| c == 1), "{seen:?}");
}

#[test]
fn a_box_has_six_flat_faces() {
    let mesh = cuboid([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 0.0,
            angle: 0.01,
        },
    )
    .unwrap();
    assert_eq!(planes.len(), 6);
    check(&mesh, &planes, 0.0);
    for plane in &planes {
        assert_eq!(plane.triangles.len(), 2);
        assert!(plane.coplanar);
        let axis = [Vec3::X, Vec3::Y, Vec3::Z]
            .iter()
            .any(|a| (plane.normal.dot(*a).abs() - 1.0).abs() < 1e-15);
        assert!(axis, "{plane:?}");
    }
    // Largest first: the 3 x 4 faces.
    assert!((planes[0].area - 12.0).abs() < 1e-12);
}

#[test]
fn treads_and_risers_are_told_apart_by_angle() {
    // A distance loose enough to take a riser into a tread's plane: only
    // the angle keeps them apart.
    let mesh = stair(4, 0.28, 0.17, 1.2);
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 0.3,
            angle: 0.1,
        },
    )
    .unwrap();
    check(&mesh, &planes, 0.3);
    assert_eq!(planes.len(), 8, "{planes:#?}");
    let treads = planes.iter().filter(|q| q.normal.z.abs() > 0.99).count();
    let risers = planes.iter().filter(|q| q.normal.x.abs() > 0.99).count();
    assert_eq!((treads, risers), (4, 4));
    // Tight: the same.
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 0.0,
            angle: 0.01,
        },
    )
    .unwrap();
    assert_eq!(planes.len(), 8);
}

/// A 20 x 20 grid on `z = 0.1 x - 0.05 y + 3`, each vertex lifted by a
/// deterministic bump of at most `noise`.
fn noisy(noise: f64) -> TriMesh {
    let mut positions = Vec::new();
    let mut s = 0x2545_f491_4f6c_dd1du64;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    for i in 0..=20 {
        for j in 0..=20 {
            let (x, y) = (i as f64 * 0.5, j as f64 * 0.5);
            let bump = (2.0 * next() - 1.0) * noise;
            positions.push(p(x, y, 0.1 * x - 0.05 * y + 3.0 + bump));
        }
    }
    let mut faces = Vec::new();
    for i in 0..20u32 {
        for j in 0..20u32 {
            let a = i * 21 + j;
            faces.push([a, a + 21, a + 22, a + 1]);
        }
    }
    quads(positions, &faces)
}

#[test]
fn a_noisy_plane_is_one_region_within_its_bound() {
    let mesh = noisy(0.001);
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 0.005,
            angle: 0.2,
        },
    )
    .unwrap();
    check(&mesh, &planes, 0.005);
    assert_eq!(planes.len(), 1, "{}", planes.len());
    assert!(!planes[0].coplanar);
    let n = planes[0].normal;
    let expected = Vec3::new(-0.1, 0.05, 1.0).normalize();
    assert!(n.dot(expected) > 0.9999, "{n:?}");
}

#[test]
fn a_tolerance_tighter_than_the_noise_splits_but_never_exceeds_it() {
    let mesh = noisy(0.001);
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 0.0004,
            angle: 0.2,
        },
    )
    .unwrap();
    check(&mesh, &planes, 0.0004);
    assert!(planes.len() > 1);
}

#[test]
fn a_faceted_cylinder_keeps_its_facets_at_a_small_angle() {
    // 24 facets, 15 degrees apart.
    let n = 24u32;
    let mut positions = Vec::new();
    for k in 0..n {
        let a = std::f64::consts::TAU * k as f64 / n as f64;
        positions.push(p(a.cos(), a.sin(), 0.0));
        positions.push(p(a.cos(), a.sin(), 2.0));
    }
    let faces: Vec<[u32; 4]> = (0..n)
        .map(|k| {
            let j = (k + 1) % n;
            [2 * k, 2 * j, 2 * j + 1, 2 * k + 1]
        })
        .collect();
    let mesh = quads(positions, &faces);
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 0.01,
            angle: 0.1,
        },
    )
    .unwrap();
    check(&mesh, &planes, 0.01);
    assert_eq!(planes.len(), 24);
    // A loose angle: facets merge only as far as the distance allows
    // (a neighbour's far edge is cos 7.5 - cos 22.5 = 0.068 off a facet's
    // plane).
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 0.08,
            angle: 1.0,
        },
    )
    .unwrap();
    check(&mesh, &planes, 0.08);
    assert!(planes.len() < 24 && planes.len() > 4, "{}", planes.len());
}

#[test]
fn far_from_the_origin() {
    let mesh = cuboid([1.0e6, -2.0e6, 300.0], [1.0e6 + 2.0, -2.0e6 + 3.0, 304.0]);
    let planes = detect_planes(
        &mesh,
        PlaneTolerance {
            distance: 1e-6,
            angle: 0.01,
        },
    )
    .unwrap();
    assert_eq!(planes.len(), 6);
    check(&mesh, &planes, 1e-6);
    assert!(planes.iter().all(|q| q.coplanar && q.deviation < 1e-8));
}

#[test]
fn refusals() {
    let mesh = cuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    for (distance, angle) in [(-1.0, 0.1), (0.1, f64::NAN), (f64::INFINITY, 0.1)] {
        assert_eq!(
            detect_planes(&mesh, PlaneTolerance { distance, angle }),
            Err(PlaneError::InvalidTolerance)
        );
    }
    let mut bad = mesh.clone();
    bad.positions[0] = p(f64::NAN, 0.0, 0.0);
    assert_eq!(
        detect_planes(
            &bad,
            PlaneTolerance {
                distance: 0.1,
                angle: 0.1
            }
        ),
        Err(PlaneError::NonFinite)
    );
}
