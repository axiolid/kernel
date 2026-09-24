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
    let mut b = GeometryGraphBuilder::new();
    let node = b.push(GeometryNode::PolygonMesh(mesh)).unwrap();
    let graph = b.finish(vec![node]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(
        &graph,
        node,
        &ExecutionOptions::new(Tolerance::METRE),
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
