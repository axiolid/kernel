//! Authored polygon faces: n-gons and faces with holes (#160).
//!
//! IFC4 `IfcPolygonalFaceSet` and `IfcIndexedPolygonalFaceWithVoids` reach
//! the compiler as `PolygonMesh` faces, deliberately not pre-triangulated.
//! Every oracle here is computed from the authored corners by hand -- a
//! shoelace area, a cross-product normal, a closed-form volume -- never read
//! back from the compiler.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_mesh::{PolygonFace, PolygonMesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode};

fn compile(mesh: PolygonMesh) -> Result<TriMesh, GeomError> {
    compile_at(mesh, Tolerance::METRE)
}

/// Compile at an explicit tolerance. Real exports go through the IFC path
/// at `Tolerance::MILLIMETRE`, so their regression tests use that.
fn compile_at(mesh: PolygonMesh, tolerance: Tolerance) -> Result<TriMesh, GeomError> {
    let mut b = GeometryGraphBuilder::new();
    let node = b.push(GeometryNode::PolygonMesh(mesh)).unwrap();
    let graph = b.finish(vec![node]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(
        &graph,
        node,
        &ExecutionOptions::new(tolerance),
    )
}

fn face(outer: &[u32]) -> PolygonFace {
    PolygonFace {
        outer: outer.to_vec(),
        holes: Vec::new(),
    }
}

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// Per-triangle unnormalised normals, three indices at a time.
fn normals(mesh: &TriMesh) -> Vec<Vec3> {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            (b - a).cross(c - a)
        })
        .collect()
}

fn area(mesh: &TriMesh) -> f64 {
    normals(mesh).iter().map(|n| n.length() / 2.0).sum()
}

fn close(got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-12 * want.abs().max(1.0),
        "got {got}, want {want}"
    );
}

/// Every triangle's normal points along `direction`, none degenerate.
fn all_face(mesh: &TriMesh, direction: Vec3) {
    for n in normals(mesh) {
        assert!(n.length() > 1e-12, "degenerate triangle");
        assert!(
            n.normalize().dot(direction.normalize()) > 1.0 - 1e-12,
            "triangle normal {n:?} is not along {direction:?}"
        );
    }
}

fn is_invalid_naming(error: &GeomError, needle: &str) -> bool {
    matches!(error, GeomError::InvalidInput(text) if text.contains(needle))
}

#[test]
fn a_quad_becomes_two_triangles_with_its_area_and_winding() {
    let mesh = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(3.0, 0.0, 0.0),
            p(3.0, 2.0, 0.0),
            p(0.0, 2.0, 0.0),
        ],
        faces: vec![face(&[0, 1, 2, 3])],
    })
    .expect("a planar quad triangulates");
    assert_eq!(mesh.indices.len(), 6);
    close(area(&mesh), 6.0);
    // Counter-clockwise seen from +z: every triangle faces +z.
    all_face(&mesh, Vec3::Z);
    // Positions are the authored ones, untouched and not duplicated.
    assert_eq!(mesh.positions.len(), 4);
}

#[test]
fn a_clockwise_quad_keeps_facing_the_other_way() {
    let mesh = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(0.0, 2.0, 0.0),
            p(3.0, 2.0, 0.0),
            p(3.0, 0.0, 0.0),
        ],
        faces: vec![face(&[0, 1, 2, 3])],
    })
    .unwrap();
    close(area(&mesh), 6.0);
    all_face(&mesh, -Vec3::Z);
}

#[test]
fn a_concave_pentagon_is_not_fanned_across_its_notch() {
    // An arrow: the notch corner (2, 1) is reflex. A fan from corner 0
    // would cover the notch; the right area proves it is not covered.
    //   shoelace over (0,0) (4,0) (4,2) (2,1) (0,2) = 8 - 2 = 6
    let mesh = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(4.0, 0.0, 0.0),
            p(4.0, 2.0, 0.0),
            p(2.0, 1.0, 0.0),
            p(0.0, 2.0, 0.0),
        ],
        faces: vec![face(&[0, 1, 2, 3, 4])],
    })
    .unwrap();
    assert_eq!(mesh.indices.len(), 9, "a pentagon is three triangles");
    close(area(&mesh), 6.0);
    all_face(&mesh, Vec3::Z);
}

#[test]
fn a_tilted_hexagon_keeps_its_area_and_plane() {
    // A regular hexagon of circumradius 1 in the plane x + y + z = 3,
    // wound counter-clockwise about (1, 1, 1). Area = 3 sqrt(3) / 2.
    let n = Vec3::new(1.0, 1.0, 1.0).normalize();
    let u = Vec3::new(1.0, -1.0, 0.0).normalize();
    let v = n.cross(u);
    let centre = p(1.0, 1.0, 1.0);
    let positions: Vec<Point3> = (0..6)
        .map(|k| {
            let a = std::f64::consts::FRAC_PI_3 * f64::from(k);
            centre + u * a.cos() + v * a.sin()
        })
        .collect();
    let mesh = compile(PolygonMesh {
        positions,
        faces: vec![face(&[0, 1, 2, 3, 4, 5])],
    })
    .unwrap();
    assert_eq!(mesh.indices.len(), 12);
    assert!((area(&mesh) - 1.5 * 3f64.sqrt()).abs() < 1e-12);
    all_face(&mesh, n);
}

#[test]
fn a_face_with_a_hole_covers_only_the_ring() {
    // 4 x 4 square minus a 2 x 2 window: area 12. A polygon with k corners
    // and h holes triangulates into k + 2h - 2 = 8 + 2 - 2 = 8 triangles.
    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(4.0, 0.0, 0.0),
        p(4.0, 4.0, 0.0),
        p(0.0, 4.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(3.0, 1.0, 0.0),
        p(3.0, 3.0, 0.0),
        p(1.0, 3.0, 0.0),
    ];
    // IFC does not fix a hole's winding, so both must give the same face.
    for hole in [vec![4, 7, 6, 5], vec![4, 5, 6, 7]] {
        let mesh = compile(PolygonMesh {
            positions: positions.clone(),
            faces: vec![PolygonFace {
                outer: vec![0, 1, 2, 3],
                holes: vec![hole.clone()],
            }],
        })
        .unwrap();
        assert_eq!(mesh.indices.len(), 24, "hole {hole:?}");
        close(area(&mesh), 12.0);
        all_face(&mesh, Vec3::Z);
        // No triangle covers the window: its centroid is never inside it.
        for t in mesh.indices.chunks_exact(3) {
            let c = t
                .iter()
                .fold(Vec3::ZERO, |s, &i| s + mesh.positions[i as usize])
                / 3.0;
            assert!(
                !(c.x > 1.0 && c.x < 3.0 && c.y > 1.0 && c.y < 3.0),
                "a triangle lies in the window: centroid {c:?}"
            );
        }
    }
}

#[test]
fn a_cube_of_quads_is_a_closed_solid_with_unit_volume() {
    // Six outward quads sharing eight corners. The divergence-theorem
    // volume is only 1 if every face kept its authored winding and the
    // corners stayed shared.
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
    let quads = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ];
    let mesh = compile(PolygonMesh {
        positions,
        faces: quads.iter().map(|q| face(q)).collect(),
    })
    .unwrap();
    assert_eq!(mesh.indices.len(), 36);
    let volume: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0;
    close(volume, 1.0);
    close(area(&mesh), 6.0);
}

/// Is every edge used by exactly two triangles, in opposite directions?
fn is_closed_two_manifold(mesh: &TriMesh) -> bool {
    use std::collections::HashMap;
    let mut uses: HashMap<(u32, u32), i32> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *uses.entry((a.min(b), a.max(b))).or_default() += if a < b { 1 } else { -1 };
        }
    }
    let mut count: HashMap<(u32, u32), usize> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *count.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    count.values().all(|&n| n == 2) && uses.values().all(|&d| d == 0)
}

#[test]
fn a_corner_on_a_straight_run_stays_shared_with_its_neighbour() {
    // A 2 x 1 x 1 box whose front bottom edge carries a midpoint (8): the
    // front face is a pentagon with 8 on its straight bottom run, and the
    // bottom face is a pentagon that uses 8 too. Dropping 8 from either
    // face while the other keeps it leaves a T-junction crack, so the
    // mesh is no longer closed although the authored faces are.
    let positions = vec![
        p(0.0, 0.0, 0.0),
        p(2.0, 0.0, 0.0),
        p(2.0, 1.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.0, 0.0, 1.0),
        p(2.0, 0.0, 1.0),
        p(2.0, 1.0, 1.0),
        p(0.0, 1.0, 1.0),
        p(1.0, 0.0, 0.0),
    ];
    let faces = vec![
        face(&[0, 3, 2, 1, 8]),
        face(&[4, 5, 6, 7]),
        face(&[0, 8, 1, 5, 4]),
        face(&[1, 2, 6, 5]),
        face(&[2, 3, 7, 6]),
        face(&[3, 0, 4, 7]),
    ];
    let mesh = compile(PolygonMesh { positions, faces }).unwrap();
    assert!(
        mesh.indices.contains(&8),
        "the straight-run corner is used by a triangle"
    );
    assert!(is_closed_two_manifold(&mesh), "no T-junction crack");
    let volume: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0;
    close(volume, 2.0);
    close(area(&mesh), 10.0);
    // Every triangle has real area: re-inserting the corner must not
    // leave a sliver along the straight run.
    assert!(normals(&mesh).iter().all(|n| n.length() > 1e-9));
}

#[test]
fn a_run_of_corners_on_a_hole_edge_is_kept() {
    // A 4 x 4 square with a 2 x 2 hole whose bottom edge carries two extra
    // corners (8, 9). Every hole corner must be used, else a neighbour
    // sharing that edge (a reveal, a window lining) would crack.
    let mut positions = vec![
        p(0.0, 0.0, 0.0),
        p(4.0, 0.0, 0.0),
        p(4.0, 4.0, 0.0),
        p(0.0, 4.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(3.0, 1.0, 0.0),
        p(3.0, 3.0, 0.0),
        p(1.0, 3.0, 0.0),
    ];
    positions.push(p(1.5, 1.0, 0.0));
    positions.push(p(2.5, 1.0, 0.0));
    let mesh = compile(PolygonMesh {
        positions,
        faces: vec![PolygonFace {
            outer: vec![0, 1, 2, 3],
            holes: vec![vec![4, 7, 6, 5, 9, 8]],
        }],
    })
    .unwrap();
    for corner in 0..10 {
        assert!(mesh.indices.contains(&corner), "corner {corner} is used");
    }
    close(area(&mesh), 12.0);
    all_face(&mesh, Vec3::Z);
    assert!(normals(&mesh).iter().all(|n| n.length() > 1e-9));
}

#[test]
fn authored_triangles_keep_their_exact_corner_order() {
    // Mixed with a quad, a triangle still passes through untouched.
    let mesh = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(2.0, 0.0, 0.0),
        ],
        faces: vec![face(&[1, 4, 2]), face(&[0, 1, 2, 3])],
    })
    .unwrap();
    assert_eq!(&mesh.indices[..3], &[1, 4, 2]);
    assert_eq!(mesh.indices.len(), 9);
}

#[test]
fn a_repeated_closing_corner_is_not_a_triangle() {
    // Exporters often repeat the first corner at the end of a ring.
    let mesh = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(2.0, 0.0, 0.0),
            p(2.0, 2.0, 0.0),
            p(0.0, 2.0, 0.0),
        ],
        faces: vec![face(&[0, 1, 2, 3, 0])],
    })
    .unwrap();
    assert_eq!(mesh.indices.len(), 6);
    close(area(&mesh), 4.0);
    all_face(&mesh, Vec3::Z);
}

/// A wall lining from a real ArchiCAD export (`IfcPolygonalFaceSet`, 24
/// corners, 14 faces, authored `Closed`), with the positions exactly as
/// `ifc-geometry` lowers them (placement applied, metres). The two big
/// faces are an outer ring with a notch plus one window hole; the notch
/// head (corners 0, 7) and the window head (corners 9, 10) lie on the same
/// line z = 2.2.
///
/// earcut bridges the collinear heads with one edge running from the notch
/// corner 7 across to the window, over corners that the reveal quads still
/// split that line at. Without the edge split this mesh had four boundary
/// edges on the back face and was reported as a surface. The oracle is
/// closure from the corner indices alone, plus the volume: front face area
/// (outer rectangle minus notch minus window) times the 10 mm thickness.
#[test]
fn collinear_notch_and_window_heads_close_the_real_lining() {
    let positions: Vec<Point3> = [
        (10.10403333333333, 10.940216690615548, 2.2000000000000006),
        (10.10403333333333, 10.940216690615548, -0.17499999988571016),
        (6.2354499850988505, 10.940216690615548, -0.17499995568073604),
        (6.23544998509885, 10.940216690615548, 2.7500000320375286),
        (13.56276667411723, 10.940216690615548, 2.7500000320375286),
        (13.56276667411723, 10.940216690615548, -0.17499999926667445),
        (13.10403333333333, 10.940216690615548, -0.17499999998435126),
        (13.10403333333333, 10.940216690615548, 2.2000000000000006),
        (8.984750000000002, 10.940216690615548, 1.1999999999999995),
        (8.984750000000002, 10.940216690615548, 2.2000000000000006),
        (6.984750000000002, 10.940216690615548, 2.2000000000000006),
        (6.984750000000002, 10.940216690615548, 1.1999999999999995),
        (10.10403333333333, 10.930216681567819, 2.2000000000000006),
        (10.10403333333333, 10.93021668156782, -0.17499999988571016),
        (6.2354499850988505, 10.93021668156782, -0.17499995568073604),
        (6.23544998509885, 10.930216681567819, 2.7500000320375286),
        (13.56276667411723, 10.93021668156782, 2.7500000320375286),
        (13.56276667411723, 10.93021668156782, -0.17499999926667445),
        (13.10403333333333, 10.93021668156782, -0.17499999998435126),
        (13.10403333333333, 10.930216681567819, 2.2000000000000006),
        (8.984750000000002, 10.93021668156782, 1.1999999999999995),
        (8.984750000000002, 10.93021668156782, 2.2000000000000006),
        (6.984750000000002, 10.93021668156782, 2.2000000000000006),
        (6.984750000000002, 10.93021668156782, 1.1999999999999995),
    ]
    .iter()
    .map(|&(x, y, z)| p(x, y, z))
    .collect();
    let with_hole = |outer: &[u32], hole: &[u32]| PolygonFace {
        outer: outer.to_vec(),
        holes: vec![hole.to_vec()],
    };
    let mut faces = vec![with_hole(&[0, 1, 2, 3, 4, 5, 6, 7], &[8, 9, 10, 11])];
    for q in [
        [1, 0, 12, 13],
        [1, 13, 14, 2],
        [15, 3, 2, 14],
        [16, 4, 3, 15],
        [4, 16, 17, 5],
        [18, 6, 5, 17],
        [7, 6, 18, 19],
        [0, 7, 19, 12],
        [8, 20, 21, 9],
        [9, 21, 22, 10],
        [10, 22, 23, 11],
        [11, 23, 20, 8],
    ] {
        faces.push(face(&q));
    }
    faces.push(with_hole(
        &[19, 18, 17, 16, 15, 14, 13, 12],
        &[22, 21, 20, 23],
    ));
    let mesh = compile(PolygonMesh { positions, faces }).unwrap();
    let adjacency = axiolid_mesh::EdgeAdjacency::build(&mesh);
    assert!(
        adjacency.is_closed_two_manifold(),
        "boundary edges: {:?}",
        adjacency.boundary_edges().collect::<Vec<_>>()
    );
    // Front face area by hand: outer rectangle, minus the notch under the
    // head, minus the window. Thickness is the y offset.
    let (x0, x1, z0, z1) = (
        6.23544998509885,
        13.56276667411723,
        -0.17499999988571016,
        2.7500000320375286,
    );
    let head = 2.2000000000000006;
    let notch = (13.10403333333333 - 10.10403333333333) * (head - z0);
    let window = (8.984750000000002 - 6.984750000000002) * (head - 1.1999999999999995);
    let thickness = 10.940216690615548 - 10.93021668156782;
    let want = ((x1 - x0) * (z1 - z0) - notch - window) * thickness;
    let volume: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0;
    // The bottom corners differ in z by up to 4.4e-8 m, so the one-z0 hand
    // formula is exact only to that times the width and thickness.
    assert!(
        (volume.abs() - want).abs() < 1e-8,
        "volume {volume}, want {want}"
    );
}

#[test]
fn a_non_planar_quad_is_refused_by_face_index() {
    // Corner 2 lifted 0.1 off the plane of the other three: a non-planar
    // quad has two different triangulations, so neither may be picked.
    let error = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.1),
            p(0.0, 1.0, 0.0),
            p(5.0, 0.0, 0.0),
            p(6.0, 0.0, 0.0),
            p(6.0, 1.0, 0.0),
        ],
        faces: vec![face(&[4, 5, 6]), face(&[0, 1, 2, 3])],
    })
    .expect_err("non-planar");
    assert!(
        is_invalid_naming(&error, "face 1") && is_invalid_naming(&error, "not planar"),
        "{error:?}"
    );
}

#[test]
fn a_quad_off_its_plane_by_less_than_the_tolerance_is_triangulated() {
    // 1e-7 is below the 1e-6 metre tolerance: exporter rounding, not shape.
    let mesh = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 1e-7),
            p(0.0, 1.0, 0.0),
        ],
        faces: vec![face(&[0, 1, 2, 3])],
    })
    .unwrap();
    assert_eq!(mesh.indices.len(), 6);
}

#[test]
fn a_hole_outside_its_face_is_refused() {
    let error = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(5.0, 5.0, 0.0),
            p(6.0, 5.0, 0.0),
            p(6.0, 6.0, 0.0),
        ],
        faces: vec![PolygonFace {
            outer: vec![0, 1, 2, 3],
            holes: vec![vec![4, 6, 5]],
        }],
    })
    .expect_err("a hole must lie inside its face");
    assert!(is_invalid_naming(&error, "face 0"), "{error:?}");
}

#[test]
fn a_self_crossing_bowtie_is_refused() {
    // (0,0) (1,1) (1,0) (0,1): the two diagonals cross. Its shoelace area
    // is zero, so no triangulation can match it.
    let error = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
        ],
        faces: vec![face(&[0, 1, 2, 3])],
    })
    .expect_err("a bowtie has no triangulation");
    assert!(
        matches!(&error, GeomError::Degenerate(t) | GeomError::InvalidInput(t) if t.contains("face 0")),
        "{error:?}"
    );
}

#[test]
fn a_collinear_ring_is_degenerate() {
    let error = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(2.0, 0.0, 0.0),
            p(3.0, 0.0, 0.0),
        ],
        faces: vec![face(&[0, 1, 2, 3])],
    })
    .expect_err("no area");
    assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
}

#[test]
fn a_hole_index_out_of_range_is_invalid_input() {
    let error = compile(PolygonMesh {
        positions: vec![
            p(0.0, 0.0, 0.0),
            p(4.0, 0.0, 0.0),
            p(4.0, 4.0, 0.0),
            p(0.0, 4.0, 0.0),
        ],
        faces: vec![PolygonFace {
            outer: vec![0, 1, 2, 3],
            holes: vec![vec![1, 2, 9]],
        }],
    })
    .expect_err("index 9 does not exist");
    assert!(is_invalid_naming(&error, "index 9"), "{error:?}");
}

/// A wall layer from a real ArchiCAD export (`IfcPolygonalFaceSet`, 42
/// corners, 27 faces, two window holes), authored closed and manifold:
/// every authored edge is used by exactly two faces.
///
/// On the back face earcut runs a diagonal along the window head, collinear
/// with the head of the notch and with the hole's own head edge. Splitting
/// that invented edge must not duplicate an authored edge: the mesh stays
/// closed and two-manifold, with no edge used three or four times.
#[test]
fn a_diagonal_along_a_window_head_keeps_the_real_layer_manifold() {
    let positions = vec![
        p(1.179999997168804, 4.611966666666668, 0.07500000037252952),
        p(1.1799999971688038, 4.611966666666668, 2.2000000000000006),
        p(1.1799999971688042, 6.361966666666667, 2.2000000000000006),
        p(1.1799999971688042, 6.361966666666667, 0.07500000037252968),
        p(1.179999999999999, 10.326949991655312, 0.07500000037252952),
        p(1.1800000032119902, 10.326949991655312, 2.7500000320375286),
        p(1.1800000032119906, 1.0200000083446863, 2.7500000320375286),
        p(1.179999999999999, 1.0200000083446863, 0.07500000037252952),
        p(1.179999985947661, 8.983933333333333, 0.9499999999999995),
        p(1.179999985947661, 7.858933333333333, 0.9499999999999995),
        p(1.1799999658728868, 7.858933333333333, 2.2000000000000006),
        p(1.1799999658728868, 8.983933333333333, 2.2000000000000006),
        p(1.1800000120448917, 2.364066666666668, 2.2000000000000006),
        p(1.1800000120448917, 3.489066666666668, 2.2000000000000006),
        p(1.180000032119666, 3.489066666666668, 0.9499999999999995),
        p(1.180000032119666, 2.364066666666668, 0.9499999999999995),
        p(1.0200000037252945, 4.611966666666668, 0.07500000037252952),
        p(1.0199999695981772, 4.611966666666668, 2.2000000000000006),
        p(1.0199999695981772, 6.361966666666667, 2.2000000000000006),
        p(1.020000003725295, 6.361966666666667, 0.07500000037252952),
        p(1.0200000037252959, 10.326949991655312, 0.07500000037252952),
        p(1.0200000037253198, 10.326949991655312, 2.9500000227242937),
        p(1.179999999999999, 10.326949991655312, 2.9500000227242937),
        p(1.179999999999999, 1.0200000083446863, 2.9500000227242937),
        p(1.0200000037252936, 1.0200000083446863, 0.07500000037252952),
        p(1.0200000037253203, 1.0200000083446863, 2.9500000227242937),
        p(1.109999997168804, 7.858933333333333, 0.9499999999999995),
        p(1.1099999971688044, 8.983933333333333, 0.9499999999999995),
        p(1.109999997168804, 7.858933333333333, 2.2000000000000006),
        p(1.1099999971688042, 8.983933333333333, 2.2000000000000006),
        p(1.1099999971688037, 2.364066666666668, 2.2000000000000006),
        p(1.109999997168804, 3.489066666666668, 2.2000000000000006),
        p(1.1099999971688042, 3.489066666666668, 0.9499999999999995),
        p(1.1099999971688037, 2.364066666666668, 0.9499999999999995),
        p(1.0199999896729512, 7.858933333333333, 0.9499999999999995),
        p(1.0199999896729512, 8.983933333333333, 0.9499999999999995),
        p(1.019999969598177, 8.983933333333333, 2.2000000000000006),
        p(1.019999969598177, 7.858933333333333, 2.2000000000000006),
        p(1.020000015770185, 3.489066666666668, 2.2000000000000006),
        p(1.020000015770185, 2.364066666666668, 2.2000000000000006),
        p(1.0200000358449595, 2.364066666666668, 0.9499999999999995),
        p(1.0200000358449595, 3.489066666666668, 0.9499999999999995),
    ];
    let faces = vec![
        PolygonFace {
            outer: vec![0, 1, 2, 3, 4, 5, 6, 7],
            holes: vec![vec![8, 9, 10, 11], vec![12, 13, 14, 15]],
        },
        PolygonFace {
            outer: vec![1, 0, 16, 17],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![2, 1, 17, 18],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![3, 2, 18, 19],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![3, 19, 20, 4],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![4, 20, 21, 22, 5],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![6, 5, 22, 23],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![24, 7, 6, 23, 25],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![16, 0, 7, 24],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![26, 9, 8, 27],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![28, 10, 9, 26],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![29, 11, 10, 28],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![27, 8, 11, 29],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![12, 30, 31, 13],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![13, 31, 32, 14],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![14, 32, 33, 15],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![15, 33, 30, 12],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![17, 16, 24, 25, 21, 20, 19, 18],
            holes: vec![vec![34, 35, 36, 37], vec![38, 39, 40, 41]],
        },
        PolygonFace {
            outer: vec![23, 22, 21, 25],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![34, 26, 27, 35],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![37, 28, 26, 34],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![36, 29, 28, 37],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![35, 27, 29, 36],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![30, 39, 38, 31],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![31, 38, 41, 32],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![32, 41, 40, 33],
            holes: vec![],
        },
        PolygonFace {
            outer: vec![33, 40, 39, 30],
            holes: vec![],
        },
    ];
    let mesh = compile(PolygonMesh { positions, faces }).unwrap();
    let adjacency = axiolid_mesh::EdgeAdjacency::build(&mesh);
    assert!(
        adjacency.is_closed_two_manifold(),
        "boundary {:?}, non-manifold {:?}",
        adjacency.boundary_edges().collect::<Vec<_>>(),
        adjacency.non_manifold_edges().collect::<Vec<_>>()
    );
}

/// A stair stringer from a real ArchiCAD export (`arc_boe` #21009,
/// `IfcPolygonalFaceSet`, 38 corners, 21 faces), authored closed. Its two
/// sawtooth side faces have 19 corners each; the inner step corners lie on
/// one line up to 1e-8..5e-8 m of export noise, well under a millimetre.
///
/// earcut bridges those steps with long diagonals. A collinearity band
/// scaled to the face's size (2e-8 of ~20 m) put some step corners inside
/// it and others outside, so one side of a diagonal was split and the
/// other was not: three unpaired edges, reported as a surface. The band is
/// now tied to the linear tolerance, which the corners all fall inside.
#[test]
fn a_stair_stringer_with_noisy_collinear_steps_closes() {
    let positions = vec![
        p(19.055730261854272, -2.9397826126246684, -2.8),
        p(18.838372983807588, -2.526678826402099, -2.8),
        p(18.231554326511727, -1.3733745196303921, -1.906627452209237),
        p(18.231554326511727, -1.3733745196303921, -1.754117659085965),
        p(18.161708907836992, -1.2406280702631145, -1.754117659085965),
        p(18.161708907836992, -1.2406280702631145, -1.654117588496447),
        p(18.329337942590353, -1.5592196056363878, -1.654117588496447),
        p(18.329337942590353, -1.5592196056363878, -1.832352844473842),
        p(18.450403334222592, -1.789313449805038, -1.832352844473842),
        p(18.450403334222592, -1.789313449805038, -2.010588133233825),
        p(18.57146869532532, -2.019407235950066, -2.010588133233825),
        p(18.57146869532532, -2.019407235950066, -2.18882338921122),
        p(18.69253411748707, -2.2495011381423398, -2.18882338921122),
        p(18.69253411748707, -2.2495011381423398, -2.36705871075379),
        p(18.81359949385455, -2.4795949532991797, -2.36705871075379),
        p(18.81359949385455, -2.4795949532991797, -2.545294065078947),
        p(18.934664854957276, -2.7096887394442053, -2.545294065078947),
        p(18.934664854957276, -2.7096887394442053, -2.723529321056342),
        p(19.055730261854272, -2.9397826126246684, -2.723529321056342),
        p(17.76755152328832, -3.0900985908654626, -2.8),
        p(17.984908801335, -3.503202377088032, -2.8),
        p(17.160732865992454, -1.9367942840937555, -1.906627452209237),
        p(17.160732865992454, -1.9367942840937555, -1.754117659085965),
        p(17.09088744731772, -1.8040478347264775, -1.754117659085965),
        p(17.09088744731772, -1.8040478347264775, -1.654117588496447),
        p(17.25851648207108, -2.1226393700997517, -1.654117588496447),
        p(17.25851648207108, -2.1226393700997517, -1.832352844473842),
        p(17.37958187370332, -2.352733214268402, -1.832352844473842),
        p(17.37958187370332, -2.352733214268402, -2.010588133233825),
        p(17.500647234806046, -2.58282700041343, -2.010588133233825),
        p(17.500647234806046, -2.58282700041343, -2.18882338921122),
        p(17.621712656967798, -2.8129209026057023, -2.18882338921122),
        p(17.621712656967798, -2.8129209026057023, -2.36705871075379),
        p(17.74277803333528, -3.043014717762543, -2.36705871075379),
        p(17.74277803333528, -3.043014717762543, -2.545294065078947),
        p(17.863843394438007, -3.273108503907569, -2.545294065078947),
        p(17.863843394438007, -3.273108503907569, -2.723529321056342),
        p(17.984908801335, -3.503202377088032, -2.723529321056342),
    ];
    let mesh = compile_at(
        PolygonMesh {
            positions,
            faces: vec![
                face(&[
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                ]),
                face(&[19, 1, 0, 20]),
                face(&[1, 19, 21, 2]),
                face(&[2, 21, 22, 3]),
                face(&[3, 22, 23, 4]),
                face(&[4, 23, 24, 5]),
                face(&[5, 24, 25, 6]),
                face(&[6, 25, 26, 7]),
                face(&[7, 26, 27, 8]),
                face(&[8, 27, 28, 9]),
                face(&[9, 28, 29, 10]),
                face(&[10, 29, 30, 11]),
                face(&[11, 30, 31, 12]),
                face(&[12, 31, 32, 13]),
                face(&[13, 32, 33, 14]),
                face(&[14, 33, 34, 15]),
                face(&[15, 34, 35, 16]),
                face(&[16, 35, 36, 17]),
                face(&[17, 36, 37, 18]),
                face(&[20, 0, 18, 37]),
                face(&[
                    19, 20, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21,
                ]),
            ],
        },
        Tolerance::MILLIMETRE,
    )
    .unwrap();
    let adjacency = axiolid_mesh::EdgeAdjacency::build(&mesh);
    assert!(
        adjacency.is_closed_two_manifold(),
        "boundary {:?}, non-manifold {:?}",
        adjacency.boundary_edges().collect::<Vec<_>>(),
        adjacency.non_manifold_edges().collect::<Vec<_>>()
    );
}
