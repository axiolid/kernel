//! A keyhole face of a closed shell never has its seam as a triangle edge
//! (#283).
//!
//! A keyhole face carries its hole in one ring: the outer boundary, a seam
//! to the hole, the hole, and the seam back (#270). The ear clipper may
//! make the seam a triangle edge, and a thin rim seamed at its corner
//! forces it. Inside one face that is a diagonal, but other faces can run
//! along the seam too (a rim's two end faces at its cut, here a fin of two
//! coincident triangles hanging from it), and the segment then had four
//! triangle uses: a declared B-rep solid read as `OpenSolid` and authored
//! polygons as a surface. A face of a shell now takes the seam out by a
//! flip or a Steiner point inside the face, so these boxes, whose top face
//! is a keyhole into a pyramid pocket, compile as closed solids with their
//! exact volume, as B-reps and as authored polygons. The volume oracle is
//! closed form: the box less the pocket's `area * depth / 3`.

use std::collections::HashMap;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_mesh::{EdgeAdjacency, PolygonFace, PolygonMesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::{CompileOutcome, MeshClosure};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation};
use axiolid_topology::{
    BRep, Edge, EdgeId, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

/// The detail under which a face reports Steiner points off an axis plane.
const STEINER: &str = "point added inside a keyhole face, off an axis plane";

/// The detail of a boolean measured against its operands flattened (#266).
const FLATTENED: &str =
    "measured against the exact boolean of operands flattened within their slab width, plus it";

/// A polyhedron as corners and faces, each face one ring of corner
/// indices, wound so the material is on the left seen from outside.
struct Polyhedron {
    corners: Vec<[f64; 3]>,
    faces: Vec<Vec<usize>>,
    /// The top face's seam, by corner index.
    seam: (usize, usize),
}

impl Polyhedron {
    fn polygons(&self) -> PolygonMesh {
        PolygonMesh {
            positions: self
                .corners
                .iter()
                .map(|c| Point3::new(c[0], c[1], c[2]))
                .collect(),
            faces: self
                .faces
                .iter()
                .map(|ring| PolygonFace {
                    outer: ring.iter().map(|&i| i as u32).collect(),
                    holes: Vec::new(),
                })
                .collect(),
        }
    }

    /// The same faces as a faceted B-rep solid, every edge one topological
    /// edge between its two corners, the seam's used both ways by the top
    /// face.
    fn brep(&self) -> BRep<NodeId> {
        let mut brep = BRep::default();
        let vertices: Vec<_> = self
            .corners
            .iter()
            .map(|c| {
                brep.add_vertex(Vertex {
                    position: Vec3::new(c[0], c[1], c[2]),
                })
            })
            .collect();
        let mut edges: HashMap<(usize, usize), EdgeId> = HashMap::new();
        let mut faces = Vec::new();
        for ring in &self.faces {
            let mut uses = Vec::new();
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                let key = (a.min(b), a.max(b));
                let edge = *edges.entry(key).or_insert_with(|| {
                    brep.add_edge(Edge {
                        start: vertices[key.0],
                        end: vertices[key.1],
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
            let loop_id = brep.add_loop(Loop { edges: uses });
            faces.push(brep.add_face(Face {
                surface: None,
                bounds: vec![FaceBound {
                    loop_id,
                    orientation: Orientation::Forward,
                    outer: true,
                }],
                orientation: Orientation::Forward,
            }));
        }
        let shell = brep.add_shell(Shell {
            faces: faces.iter().map(|&f| (f, Orientation::Forward)).collect(),
            closed: true,
        });
        brep.add_solid(Solid {
            outer: shell,
            voids: Vec::new(),
        });
        brep
    }

    /// The same body with a fin hanging from the seam: two coincident
    /// triangles down to a point inside the material, as a rim's two end
    /// faces at its cut meet along the seam. Every edge is still used
    /// twice, in opposite directions, once the top face keeps the seam out
    /// of its triangles; the volume is unchanged.
    fn with_fin(mut self) -> Self {
        let (a, b) = self.seam;
        let (pa, pb) = (self.corners[a], self.corners[b]);
        let tip = self.corners.len();
        self.corners
            .push([(pa[0] + pb[0]) * 0.5, (pa[1] + pb[1]) * 0.5, pa[2] * 0.5]);
        self.faces.push(vec![a, b, tip]);
        self.faces.push(vec![b, a, tip]);
        self
    }

    /// Every corner mapped by `map`.
    fn mapped(mut self, map: impl Fn([f64; 3]) -> [f64; 3]) -> Self {
        for corner in &mut self.corners {
            *corner = map(*corner);
        }
        self
    }
}

/// A `w x d x h` box whose top face is a keyhole: its outer boundary from
/// the corner (0, 0, h), the seam to the rim's first corner, the rim
/// (clockwise seen from above) and the seam back, the rim bounding a
/// pyramid pocket down to `apex`.
fn keyhole_box(w: f64, d: f64, h: f64, rim: &[[f64; 2]], apex: [f64; 3]) -> Polyhedron {
    let mut corners = vec![
        [0.0, 0.0, 0.0],
        [w, 0.0, 0.0],
        [w, d, 0.0],
        [0.0, d, 0.0],
        [0.0, 0.0, h],
        [w, 0.0, h],
        [w, d, h],
        [0.0, d, h],
    ];
    let mut faces = vec![
        vec![0, 3, 2, 1],
        vec![0, 1, 5, 4],
        vec![1, 2, 6, 5],
        vec![2, 3, 7, 6],
        vec![3, 0, 4, 7],
    ];
    let start = corners.len();
    corners.extend(rim.iter().map(|p| [p[0], p[1], h]));
    let tip = corners.len();
    corners.push(apex);
    let ring: Vec<usize> = (start..tip).collect();
    for i in 0..ring.len() {
        faces.push(vec![ring[(i + 1) % ring.len()], ring[i], tip]);
    }
    let mut top = vec![4, 5, 6, 7, 4];
    top.extend(&ring);
    top.push(start);
    faces.push(top);
    Polyhedron {
        corners,
        faces,
        seam: (4, start),
    }
}

/// The 4 x 4 x 4 box with a 2 x 2 rim pocketed 1 deep, and its volume.
fn wide() -> (Polyhedron, f64) {
    let rim = [[1.0, 1.0], [1.0, 3.0], [3.0, 3.0], [3.0, 1.0]];
    (
        keyhole_box(4.0, 4.0, 4.0, &rim, [2.0, 2.0, 3.0]),
        64.0 - 4.0 / 3.0,
    )
}

/// The corpus shape, in metres: a 0.5 x 0.375 rim 2^-13 wide (0.12 mm),
/// seamed at its corner,
/// over a pocket 0.5 deep, and its volume. Any fan at the outer corner
/// that skips the inner one covers the pocket, so the seam is forced on
/// the face's own vertices.
fn thin_rim() -> (Polyhedron, f64) {
    let e = 1.0 / 8192.0;
    let rim = [[e, e], [e, 0.375 - e], [0.5 - e, 0.375 - e], [0.5 - e, e]];
    let pocket = (0.5 - 2.0 * e) * (0.375 - 2.0 * e);
    (
        keyhole_box(0.5, 0.375, 1.0, &rim, [0.25, 0.1875, 0.5]),
        0.5 * 0.375 - pocket * 0.5 / 3.0,
    )
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn compile(node: GeometryNode) -> (CompileOutcome, DeviationReport) {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(node).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options())
        .expect("the keyhole box tessellates")
}

fn volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

/// How many triangles use the segment between two corners of `body`,
/// found by position.
fn uses_of(mesh: &TriMesh, body: &Polyhedron, (a, b): (usize, usize)) -> usize {
    let at = |i: usize| {
        let c = body.corners[i];
        Vec3::new(c[0], c[1], c[2])
    };
    let index = |p: Vec3| {
        mesh.positions
            .iter()
            .position(|&q| q == p)
            .expect("every corner is a mesh position") as u32
    };
    let (a, b) = (index(at(a)), index(at(b)));
    mesh.indices
        .chunks_exact(3)
        .filter(|t| {
            (0..3).any(|k| (t[k], t[(k + 1) % 3]) == (a, b) || (t[k], t[(k + 1) % 3]) == (b, a))
        })
        .count()
}

/// Closed, a solid, the volume `want` within rounding, and the seam used
/// only by the fin (`fin_uses` triangles, none without one). Returns how
/// many Steiner points the mesh added (positions beyond the body's
/// corners), the deviation report and the mesh.
fn assert_closed(
    body: &Polyhedron,
    node: GeometryNode,
    want: f64,
    fin_uses: usize,
    what: &str,
) -> (usize, DeviationReport, TriMesh) {
    let (outcome, report) = compile(node);
    let mesh = outcome.mesh.clone();
    let adjacency = EdgeAdjacency::build(&mesh);
    let open: Vec<_> = adjacency
        .edges()
        .filter(|(_, uses)| uses.len() != 2)
        .map(|(key, uses)| (key, uses.len()))
        .collect();
    assert!(
        adjacency.is_closed_two_manifold(),
        "{what}: edges not used twice: {open:?}"
    );
    assert_eq!(outcome.closure, MeshClosure::Solid, "{what}");
    let solid = outcome.solid_mesh().expect("a closed solid");
    let got = volume(solid);
    assert!(
        (got - want).abs() <= 1e-12 * want,
        "{what}: volume {got}, want {want}"
    );
    // The oracle consumers use agrees. Its audit drops a fin's coincident
    // pair of triangles as a degenerate sliver, so it measures only bodies
    // without one.
    if fin_uses == 0 {
        let measured = axiolid_measure::mesh::volume_properties(solid, Tolerance::METRE)
            .expect("a closed two-manifold has a volume")
            .signed_volume;
        assert!(
            (measured - want).abs() <= 1e-12 * want,
            "{what}: {measured}"
        );
    }
    assert_eq!(
        uses_of(&mesh, body, body.seam),
        fin_uses,
        "{what}: the seam is a triangle edge"
    );
    let added = mesh.positions.len() - body.corners.len();
    (added, report, mesh)
}

#[test]
fn a_box_with_a_keyhole_top_face_is_a_closed_solid() {
    for (name, (body, want)) in [("wide", wide()), ("thin rim", thin_rim())] {
        for fin in [false, true] {
            let body = if fin {
                body_clone(&body).with_fin()
            } else {
                body_clone(&body)
            };
            let uses = if fin { 2 } else { 0 };
            let what = format!("{name}, fin {fin}");
            let (brep_added, report, mesh) = assert_closed(
                &body,
                GeometryNode::BRep(body.brep()),
                want,
                uses,
                &format!("{what}: B-rep"),
            );
            assert!(
                report.contributions.iter().all(|c| c.detail != STEINER),
                "{what}: an axis-plane face reports no Steiner offset: {report:?}"
            );
            let (polygon_added, report, _) = assert_closed(
                &body,
                GeometryNode::PolygonMesh(body.polygons()),
                want,
                uses,
                &format!("{what}: polygons"),
            );
            assert!(
                report.contributions.iter().all(|c| c.detail != STEINER),
                "{what}: {report:?}"
            );
            assert_eq!(brep_added, polygon_added, "{what}: the same face");
            if name == "thin rim" {
                // Forced: a Steiner point, exactly in the top face's plane.
                assert!(brep_added > 0, "{what}: no Steiner point");
                for p in &mesh.positions {
                    let corner = body
                        .corners
                        .iter()
                        .any(|c| Vec3::new(c[0], c[1], c[2]) == *p);
                    assert!(
                        corner || p.z == 1.0,
                        "{what}: Steiner point {p} off the plane"
                    );
                }
            }
        }
    }
}

fn body_clone(body: &Polyhedron) -> Polyhedron {
    Polyhedron {
        corners: body.corners.clone(),
        faces: body.faces.clone(),
        seam: body.seam,
    }
}

#[test]
fn a_tilted_thin_rim_reports_its_steiner_point_within_rounding() {
    // (x, y, z) -> (x, y - z, y + z): the top face lies in the plane
    // z' - y' = 2, in no axis plane, and every volume doubles. The
    // coordinates stay dyadic, so the corners lie exactly on their planes,
    // but a face off the axis planes reports its Steiner point by name.
    let (body, want) = thin_rim();
    let body = body.with_fin().mapped(|[x, y, z]| [x, y - z, y + z]);
    for (node, path, what) in [
        (
            GeometryNode::BRep(body.brep()),
            DeviationPath::BRepFace,
            "B-rep",
        ),
        (
            GeometryNode::PolygonMesh(body.polygons()),
            DeviationPath::AuthoredMesh,
            "polygons",
        ),
    ] {
        let (added, report, _) = assert_closed(&body, node, 2.0 * want, 2, what);
        assert!(added > 0, "{what}: no Steiner point");
        let steiner = report
            .contributions
            .iter()
            .find(|c| c.path == path && c.detail == STEINER)
            .unwrap_or_else(|| panic!("{what}: the Steiner point is not reported: {report:?}"));
        match steiner.bound {
            DeviationBound::Certified(bound) => {
                assert!(bound > 0.0 && bound < 1e-12, "{what}: {bound}");
            }
            other => panic!("{what}: {other:?}"),
        }
    }
}

/// A cube from `lo` to `hi`, outward, as authored quads, every corner
/// mapped by `map`.
fn cube(lo: [f64; 3], hi: [f64; 3], map: impl Fn([f64; 3]) -> [f64; 3]) -> Polyhedron {
    let corners = (0..8)
        .map(|k| {
            map([
                if k & 1 == 0 { lo[0] } else { hi[0] },
                if k & 2 == 0 { lo[1] } else { hi[1] },
                if k & 4 == 0 { lo[2] } else { hi[2] },
            ])
        })
        .collect();
    Polyhedron {
        corners,
        faces: vec![
            vec![0, 2, 3, 1],
            vec![4, 5, 7, 6],
            vec![0, 1, 5, 4],
            vec![2, 6, 7, 3],
            vec![0, 4, 6, 2],
            vec![1, 3, 7, 5],
        ],
        seam: (0, 1),
    }
}

#[test]
fn a_boolean_of_a_tilted_keyhole_box_keeps_its_steiner_report() {
    // The operand's Steiner offset is a slab width like a warp's (#266):
    // it is carried through the boolean by name, and the operand is
    // measured flattened within it.
    let tilt = |[x, y, z]: [f64; 3]| [x, y - z, y + z];
    let (body, _) = thin_rim();
    let body = body.mapped(tilt);
    let tool = cube([0.125, 0.0625, -0.25], [0.25, 0.125, 0.25], tilt);
    for (left, right, path, what) in [
        (
            GeometryNode::PolygonMesh(body.polygons()),
            GeometryNode::PolygonMesh(tool.polygons()),
            DeviationPath::AuthoredMesh,
            "polygons",
        ),
        (
            GeometryNode::BRep(body.brep()),
            GeometryNode::BRep(tool.brep()),
            DeviationPath::BRepFace,
            "B-rep",
        ),
    ] {
        let mut g = GeometryGraphBuilder::new();
        let left = g.push(left).unwrap();
        let right = g.push(right).unwrap();
        let root = g
            .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left,
                right,
                operator: axiolid_core::BooleanOperator::Union,
            }))
            .unwrap();
        let graph = g.finish(vec![root]).unwrap();
        let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_with_deviation(&graph, root, &options())
            .expect("the union compiles");
        let bound = |path: DeviationPath, detail: &str| {
            report
                .contributions
                .iter()
                .find(|c| c.path == path && c.detail == detail)
                .map(|c| c.bound)
        };
        let Some(DeviationBound::Certified(steiner)) = bound(path, STEINER) else {
            panic!("{what}: the Steiner offset is dropped by the boolean: {report:?}");
        };
        // A B-rep with Steiner points is still a polyhedron up to a slab
        // width, so the boolean is measured flattened, not refused.
        let Some(DeviationBound::Certified(measured)) = bound(DeviationPath::Boolean, FLATTENED)
        else {
            panic!("{what}: measured against the flattened boolean: {report:?}");
        };
        assert!(measured >= steiner, "{what}: {measured} below {steiner}");
    }
}
