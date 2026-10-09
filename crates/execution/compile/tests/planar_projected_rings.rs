//! Planar faces whose projected rings the certified clipper used to refuse.
//!
//! A face is projected onto in-plane axes before it is ear clipped. For
//! some planes those axes carry `-0.0` components, so a corner level with
//! the origin projects to `-0.0` next to the origin's `0.0`, and the
//! clipper read the ring's orientation at that twin (#269). Each shape
//! here is placed in every axis plane and every plane at 45 degrees
//! between two axes, facing either way, from every start corner, and must
//! compile, authored and as a faceted B-rep face, to its exact area (in a
//! diagonal plane, within rounding). So must a face whose one ring
//! joins its hole by a seam run both ways (a keyhole, #270), which the
//! clipper refused as overlapping itself.
//!
//! Coordinates are dyadic, so areas are exact sums.

use std::collections::HashMap;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Point2, Point3, Tolerance, Vec3};
use axiolid_mesh::{PolygonFace, PolygonMesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId};
use axiolid_topology::{BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Vertex};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compile_node(node: GeometryNode) -> Result<TriMesh, GeomError> {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(node).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, root, &options())
}

/// Twice the area the triangles cover, and whether each faces along
/// `normal`.
fn twice_area(mesh: &TriMesh, normal: Vec3) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            let doubled = (b - a).cross(c - a);
            assert!(doubled.dot(normal) > 0.0, "{t:?} faces away from {normal}");
            doubled.length()
        })
        .sum()
}

/// One face over the corners `0..n` (outer ring) and the given hole ranges.
fn authored(corners: &[Point3], outer: &[u32], holes: &[Vec<u32>]) -> PolygonMesh {
    PolygonMesh {
        positions: corners.to_vec(),
        faces: vec![PolygonFace {
            outer: outer.to_vec(),
            holes: holes.to_vec(),
        }],
    }
}

/// One face without a surface, as a surface model; a ring that runs along
/// an edge both ways uses one edge twice, once in each direction.
fn faceted(corners: &[Point3], outer: &[u32], holes: &[Vec<u32>]) -> BRep<NodeId> {
    let mut brep = BRep::default();
    // Corners at one position are one vertex, as an exporter welds them.
    let mut at: HashMap<(u64, u64, u64), usize> = HashMap::new();
    let mut verts = Vec::new();
    let mut vertex_of = Vec::new();
    for &position in corners {
        let key = (
            (position.x + 0.0).to_bits(),
            (position.y + 0.0).to_bits(),
            (position.z + 0.0).to_bits(),
        );
        let index = *at.entry(key).or_insert_with(|| {
            verts.push(brep.add_vertex(Vertex { position }));
            verts.len() - 1
        });
        vertex_of.push(index);
    }
    let mut edges = HashMap::new();
    let mut bounds = Vec::new();
    for (index, ring) in core::iter::once(outer)
        .chain(holes.iter().map(Vec::as_slice))
        .enumerate()
    {
        let mut uses = Vec::new();
        for i in 0..ring.len() {
            let (a, b) = (
                vertex_of[ring[i] as usize],
                vertex_of[ring[(i + 1) % ring.len()] as usize],
            );
            let key = (a.min(b), a.max(b));
            let edge = *edges.entry(key).or_insert_with(|| {
                brep.add_edge(Edge {
                    start: verts[key.0],
                    end: verts[key.1],
                    curve: None,
                })
            });
            uses.push(EdgeUse {
                edge,
                orientation: if a == key.0 {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
                pcurve: None,
            });
        }
        bounds.push(FaceBound {
            loop_id: brep.add_loop(Loop { edges: uses }),
            orientation: Orientation::Forward,
            outer: index == 0,
        });
    }
    let face = brep.add_face(Face {
        surface: None,
        bounds,
        orientation: Orientation::Forward,
    });
    brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep
}

/// In-plane frames `(e1, e2)`: the six axis planes, each either way round,
/// and the planes at 45 degrees between two axes, where `e1` is a diagonal
/// such as `(-1, -1, 0)`. The face's normal is `e1 x e2`. A faceted B-rep
/// face in such a diagonal plane, projected from its first corner, gets
/// `-0.0` twins on its first axis (found by searching every such frame for
/// the #269 shape); the axis planes keep the authored path's corners
/// covered, which are projected from the centroid.
fn frames() -> Vec<(Vec3, Vec3)> {
    let axes = [Vec3::X, Vec3::Y, Vec3::Z];
    let mut out = Vec::new();
    for i in 0..3 {
        let (a, b, c) = (axes[i], axes[(i + 1) % 3], axes[(i + 2) % 3]);
        out.push((b, c));
        out.push((c, b));
        out.push((-b, c));
        out.push((b, -c));
        for (sa, sb) in [(1.0, 1.0), (-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0)] {
            let diagonal = a * sa + b * sb;
            out.push((diagonal, c));
            out.push((diagonal, -c));
        }
    }
    out
}

/// Exactly `want` when it is a dyadic sum, else within rounding.
fn area(got: f64, want: f64, what: &str) {
    if want.fract() == 0.0 {
        assert_eq!(got, want, "{what}: area");
    } else {
        assert!(
            (got - want).abs() <= 1e-12 * want,
            "{what}: area {got}, want {want}"
        );
    }
}

/// The face over `outer` (counter-clockwise in its own plane) and `holes`
/// compiles, authored and as a faceted B-rep face, with `twice` times the
/// frame's cell area as twice its area, in every frame of [`frames`] and from every start corner.
fn check_face(outer: &[Point2], holes: &[Vec<Point2>], twice: f64, what: &str) {
    let n = outer.len();
    for (e1, e2) in frames() {
        let normal = e1.cross(e2);
        // Twice-areas scale by the frame's cell: 1 in an axis plane,
        // sqrt 2 in a diagonal one.
        let scale = normal.length();
        let place = |q: Point2| {
            let r = e1 * q.x + e2 * q.y;
            Point3::new(r.x, r.y, r.z)
        };
        for start in 0..n {
            let mut corners: Vec<Point3> = (0..n).map(|k| place(outer[(start + k) % n])).collect();
            let outer_ids: Vec<u32> = (0..n as u32).collect();
            let mut hole_ids = Vec::new();
            for hole in holes {
                let first = corners.len() as u32;
                corners.extend(hole.iter().map(|&q| place(q)));
                hole_ids.push((first..corners.len() as u32).collect::<Vec<u32>>());
            }
            let what = format!("{what} (frame {e1} {e2}, start {start})");
            let mesh = compile_node(GeometryNode::PolygonMesh(authored(
                &corners, &outer_ids, &hole_ids,
            )))
            .unwrap_or_else(|e| panic!("{what}: authored face refused: {e:?}"));
            area(
                twice_area(&mesh, normal),
                twice * scale,
                &format!("{what}: authored"),
            );
            let mesh = compile_node(GeometryNode::BRep(faceted(&corners, &outer_ids, &hole_ids)))
                .unwrap_or_else(|e| panic!("{what}: B-rep face refused: {e:?}"));
            area(
                twice_area(&mesh, normal),
                twice * scale,
                &format!("{what}: B-rep"),
            );
        }
    }
}

#[test]
fn a_rectangle_with_a_straight_corner_compiles_in_every_plane() {
    // #269: a 1 x 3 rectangle with a straight corner on its left side.
    // Facing down an axis, it projects with that side at -0.0 from some
    // start corners, which the clipper read as clockwise.
    let ring = [
        p(0.0, 0.0),
        p(1.0, 0.0),
        p(1.0, 3.0),
        p(0.0, 3.0),
        p(0.0, 2.0),
    ];
    check_face(&ring, &[], 6.0, "rectangle with a straight corner");
    // The corpus shape: straight corners on more than one side.
    let ring = [
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 2.0),
        p(2.0, 2.0),
        p(0.0, 2.0),
        p(0.0, 1.0),
    ];
    check_face(&ring, &[], 16.0, "rectangle with straight corners");
    // The same shape as a hole.
    let hole = vec![
        p(1.0, 1.0),
        p(1.0, 2.0),
        p(1.0, 4.0),
        p(2.0, 4.0),
        p(2.0, 1.0),
    ];
    let outer = [p(0.0, 0.0), p(3.0, 0.0), p(3.0, 5.0), p(0.0, 5.0)];
    check_face(&outer, &[hole], 24.0, "hole with a straight corner");
}

#[test]
fn a_keyhole_face_compiles_in_every_plane() {
    // #270: a 4 x 4 face less a 2 x 2 hole as one ring, the hole joined to
    // the outer boundary by a seam run both ways, as `IfcPolygonalFaceSet`
    // exporters write it.
    let ring = [
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 4.0),
        p(0.0, 4.0),
        p(0.0, 0.0),
        p(1.0, 1.0),
        p(1.0, 3.0),
        p(3.0, 3.0),
        p(3.0, 1.0),
        p(1.0, 1.0),
    ];
    check_face(&ring, &[], 24.0, "keyhole");
    // The seam leaving the outer ring inside its bottom edge for the
    // corner of a diamond.
    let ring = [
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(2.0, 1.0),
        p(1.0, 2.0),
        p(2.0, 3.0),
        p(3.0, 2.0),
        p(2.0, 1.0),
        p(2.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 4.0),
        p(0.0, 4.0),
    ];
    check_face(&ring, &[], 28.0, "seam from an edge");
    // The corpus shape: a rim and its inner outline 2^-13 apart.
    let e = 1.0 / 8192.0;
    let ring = [
        p(0.0, 0.0),
        p(0.5, 0.0),
        p(0.5, 0.375),
        p(0.0, 0.375),
        p(0.0, 0.0),
        p(e, e),
        p(e, 0.375 - e),
        p(0.5 - e, 0.375 - e),
        p(0.5 - e, e),
        p(e, e),
    ];
    let twice = 2.0 * (0.5 * 0.375 - (0.5 - 2.0 * e) * (0.375 - 2.0 * e));
    check_face(&ring, &[], twice, "thin rim");
}

#[test]
fn a_keyhole_face_naming_each_seam_corner_once_compiles() {
    // The exporter's other spelling: the seam's ends are one corner index
    // each, named twice in the ring.
    let positions = [
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (1.0, 1.0),
        (1.0, 3.0),
        (3.0, 3.0),
        (3.0, 1.0),
    ]
    .map(|(x, y)| Point3::new(x, y, 0.0))
    .to_vec();
    let outer = [0, 1, 2, 3, 0, 4, 5, 6, 7, 4];
    let mesh = compile_node(GeometryNode::PolygonMesh(authored(&positions, &outer, &[])))
        .expect("a keyhole is a surface patch");
    assert_eq!(twice_area(&mesh, Vec3::Z), 24.0);
    let mesh = compile_node(GeometryNode::BRep(faceted(&positions, &outer, &[])))
        .expect("a keyhole B-rep face");
    assert_eq!(twice_area(&mesh, Vec3::Z), 24.0);
}

#[test]
fn a_seam_crossing_the_hole_is_refused_by_name() {
    // The hole loop comes back across its own seam: not a keyhole.
    let positions = [
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 4.0),
        (0.0, 4.0),
        (2.0, 2.0),
        (2.0, 3.0),
        (3.0, 3.0),
        (3.0, 0.5),
        (0.5, 1.5),
    ]
    .map(|(x, y)| Point3::new(x, y, 0.0))
    .to_vec();
    let outer = [0, 1, 2, 3, 0, 4, 5, 6, 7, 8, 4];
    let error = compile_node(GeometryNode::PolygonMesh(authored(&positions, &outer, &[])))
        .expect_err("a seam crossing the ring");
    assert!(error.to_string().contains("has a seam from"), "{error:?}");
}
