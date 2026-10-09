//! A ring vertex inside a shared edge welds through to the neighbour (#265).
//!
//! The certified clipper accepts a vertex lying inside another ring's edge
//! of the same face and inserts it there (#262, #260). In a closed shell
//! the neighbouring face shares that edge, so it has to be split at the
//! same vertex, or the mesh has a T-junction: two boundary edges on one
//! side, one on the other. These boxes have pockets whose rim touches an
//! edge of the box or of another pocket at an interior point of that edge.
//! Every oracle is closed form, from the corners: the box volume less each
//! pyramid pocket's `area * depth / 3`.

use std::collections::HashMap;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_mesh::{EdgeAdjacency, PolygonFace, PolygonMesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{CompileOutcome, MeshClosure, MeshCompiler};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId};
use axiolid_topology::{
    BRep, Edge, EdgeId, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

/// A polyhedron as corners and faces, each face an outer ring and holes of
/// corner indices, wound so the material is on the left seen from outside.
struct Polyhedron {
    corners: Vec<[f64; 3]>,
    faces: Vec<(Vec<usize>, Vec<Vec<usize>>)>,
}

impl Polyhedron {
    /// The same faces as authored polygons over shared positions.
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
                .map(|(outer, holes)| PolygonFace {
                    outer: outer.iter().map(|&i| i as u32).collect(),
                    holes: holes
                        .iter()
                        .map(|h| h.iter().map(|&i| i as u32).collect())
                        .collect(),
                })
                .collect(),
        }
    }

    /// The same faces as a faceted B-rep solid, every edge shared between
    /// the two faces that use it, the way lowering builds it.
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
        for (outer, holes) in &self.faces {
            let mut bounds = Vec::new();
            for (k, ring) in std::iter::once(outer).chain(holes).enumerate() {
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
                bounds.push(FaceBound {
                    loop_id: brep.add_loop(Loop { edges: uses }),
                    orientation: Orientation::Forward,
                    outer: k == 0,
                });
            }
            faces.push(brep.add_face(Face {
                surface: None,
                bounds,
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
}

/// A 4 x 4 x 4 box, corners 0..8, whose top face (z = 4) has the given
/// holes, each the rim of a pyramid pocket down to its apex. A hole is
/// given clockwise seen from above, as the top face's inner ring; each
/// pocket wall is a triangle over one rim edge, traversed the other way.
fn box_with_pockets(rims: &[(&[[f64; 2]], [f64; 3])]) -> Polyhedron {
    let mut corners = vec![
        [0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.0, 4.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 4.0],
        [4.0, 0.0, 4.0],
        [4.0, 4.0, 4.0],
        [0.0, 4.0, 4.0],
    ];
    let mut faces: Vec<(Vec<usize>, Vec<Vec<usize>>)> = vec![
        (vec![0, 3, 2, 1], vec![]),
        (vec![0, 1, 5, 4], vec![]),
        (vec![1, 2, 6, 5], vec![]),
        (vec![2, 3, 7, 6], vec![]),
        (vec![3, 0, 4, 7], vec![]),
    ];
    let mut holes = Vec::new();
    for &(rim, apex) in rims {
        let start = corners.len();
        corners.extend(rim.iter().map(|p| [p[0], p[1], 4.0]));
        let top = corners.len();
        corners.push(apex);
        let ring: Vec<usize> = (start..top).collect();
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            faces.push((vec![b, a, top], vec![]));
        }
        holes.push(ring);
    }
    faces.push((vec![4, 5, 6, 7], holes));
    Polyhedron { corners, faces }
}

/// A diamond pocket whose rim corner (2, 0) lies inside the top face's
/// front edge (0, 0)-(4, 0), which the front face shares. Volume
/// `64 - 2 * 1 / 3`: the diamond has area 2 and the apex is 1 deep.
fn pocket_touching_the_box_edge() -> (Polyhedron, f64) {
    let rim: &[[f64; 2]] = &[[2.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 1.0]];
    (
        box_with_pockets(&[(rim, [2.0, 1.0, 3.0])]),
        64.0 - 2.0 / 3.0,
    )
}

/// Two pockets: a diamond of area 2, and a unit square whose corner
/// (2, 2.5) lies inside the diamond's rim edge (2.5, 2)-(1.5, 3), which the
/// diamond's pocket wall shares. Both 1 deep: volume `64 - 2/3 - 1/3`.
fn pocket_touching_another_pocket() -> (Polyhedron, f64) {
    let diamond: &[[f64; 2]] = &[[1.5, 1.0], [0.5, 2.0], [1.5, 3.0], [2.5, 2.0]];
    let square: &[[f64; 2]] = &[[2.0, 2.5], [2.0, 3.5], [3.0, 3.5], [3.0, 2.5]];
    (
        box_with_pockets(&[(diamond, [1.5, 2.0, 3.0]), (square, [2.5, 3.0, 3.0])]),
        63.0,
    )
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compile(node: GeometryNode) -> CompileOutcome {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(node).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_reported(&graph, root, &options())
        .expect("the pocketed box tessellates")
}

/// Every edge of the index connectivity used exactly twice, in opposite
/// directions: no boundary edge, no T-junction.
fn closed_two_manifold(mesh: &TriMesh) -> bool {
    let adjacency = EdgeAdjacency::build(mesh);
    adjacency.edge_count() > 0 && adjacency.is_closed_two_manifold()
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

fn assert_welded(outcome: &CompileOutcome, want: f64) {
    let mesh = &outcome.mesh;
    let boundary = EdgeAdjacency::build(mesh).boundary_edges().count();
    assert!(
        closed_two_manifold(mesh),
        "{boundary} boundary edges: {:?}",
        mesh.indices
    );
    assert_eq!(outcome.closure, MeshClosure::Solid);
    let solid = outcome.solid_mesh().expect("a welded solid is a solid");
    let got = volume(solid);
    assert!(
        (got - want).abs() <= 1e-12 * want,
        "volume {got}, want {want}"
    );
    // The oracle that consumers use agrees.
    let measured = axiolid_measure::mesh::volume_properties(solid, Tolerance::MILLIMETRE)
        .expect("a closed two-manifold has a volume")
        .signed_volume;
    assert!((measured - want).abs() <= 1e-12 * want, "{measured}");
}

#[test]
fn a_pocket_touching_a_box_edge_welds_the_front_face_as_a_brep() {
    let (body, want) = pocket_touching_the_box_edge();
    assert_welded(&compile(GeometryNode::BRep(body.brep())), want);
}

#[test]
fn a_pocket_touching_a_box_edge_welds_the_front_face_as_polygons() {
    let (body, want) = pocket_touching_the_box_edge();
    assert_welded(&compile(GeometryNode::PolygonMesh(body.polygons())), want);
}

#[test]
fn a_pocket_touching_another_pockets_rim_welds_its_wall_as_a_brep() {
    let (body, want) = pocket_touching_another_pocket();
    assert_welded(&compile(GeometryNode::BRep(body.brep())), want);
}

#[test]
fn a_pocket_touching_another_pockets_rim_welds_its_wall_as_polygons() {
    // The wall is an authored triangle: it gains the touching corner as a
    // fourth, collinear one and is triangulated as a polygon.
    let (body, want) = pocket_touching_another_pocket();
    assert_welded(&compile(GeometryNode::PolygonMesh(body.polygons())), want);
}

/// Two pockets whose rim corners (1, 0) and (3, 0) both lie inside the
/// top face's front edge: the front face takes both, in order along its
/// edge from (4, 0) to (0, 0). Each diamond has area 1/2 and is 1 deep.
fn two_pockets_touching_one_edge() -> (Polyhedron, f64) {
    let left: &[[f64; 2]] = &[[1.0, 0.0], [0.5, 0.5], [1.0, 1.0], [1.5, 0.5]];
    let right: &[[f64; 2]] = &[[3.0, 0.0], [2.5, 0.5], [3.0, 1.0], [3.5, 0.5]];
    (
        box_with_pockets(&[(left, [1.0, 0.5, 3.0]), (right, [3.0, 0.5, 3.0])]),
        64.0 - 1.0 / 3.0,
    )
}

#[test]
fn two_corners_inside_one_shared_edge_weld_in_order() {
    let (body, want) = two_pockets_touching_one_edge();
    assert_welded(&compile(GeometryNode::BRep(body.brep())), want);
    assert_welded(&compile(GeometryNode::PolygonMesh(body.polygons())), want);
}

#[test]
fn a_touching_face_with_an_exporters_closing_corner_still_welds() {
    // The top face's outer ring repeats its first corner at the end, as
    // exporters write it; the clipper drops it, so the touching corner's
    // index in what the clipper sees is one less than in the face.
    let (mut body, want) = pocket_touching_the_box_edge();
    let top = body.faces.last_mut().expect("the top face");
    let first = top.0[0];
    top.0.push(first);
    assert_welded(&compile(GeometryNode::PolygonMesh(body.polygons())), want);
}

/// The pocketed box with its front face (the neighbour that would take the
/// corner) left out: no face can be welded to close the shell.
fn without_front_face(mut body: Polyhedron) -> Polyhedron {
    body.faces.remove(1);
    body
}

#[test]
fn a_declared_solid_whose_mesh_is_open_is_not_reported_solid() {
    let (body, _) = pocket_touching_the_box_edge();
    let outcome = compile(GeometryNode::BRep(without_front_face(body).brep()));
    assert!(!closed_two_manifold(&outcome.mesh));
    assert_eq!(outcome.closure, MeshClosure::OpenSolid);
    let error = outcome
        .solid_mesh()
        .expect_err("an open mesh has no volume");
    assert!(error.to_string().contains("not closed"), "{error}");
}

#[test]
fn an_authored_mesh_that_is_open_is_a_surface_not_a_solid() {
    // Polygons declare no solid: an open face set is a surface (#161).
    let (body, _) = pocket_touching_the_box_edge();
    let outcome = compile(GeometryNode::PolygonMesh(
        without_front_face(body).polygons(),
    ));
    assert_eq!(outcome.closure, MeshClosure::Surface);
}

#[test]
fn a_declared_solid_with_one_face_turned_round_is_not_reported_solid() {
    // Every edge is still used twice, but two of them in the same
    // direction: the shell is not consistently oriented, so it bounds no
    // volume.
    let (mut body, _) = pocket_touching_the_box_edge();
    body.faces[2].0.reverse();
    let outcome = compile(GeometryNode::BRep(body.brep()));
    let adjacency = EdgeAdjacency::build(&outcome.mesh);
    assert!(adjacency.edges().all(|(_, uses)| uses.len() == 2));
    assert_eq!(outcome.closure, MeshClosure::OpenSolid);
}

#[test]
fn an_open_solid_keeps_its_name_in_a_collection_and_is_refused_by_a_boolean() {
    let (body, _) = pocket_touching_the_box_edge();
    let (open, _) = pocket_touching_the_box_edge();
    let mut b = GeometryGraphBuilder::new();
    let solid = b.push(GeometryNode::BRep(body.brep())).unwrap();
    let open = b
        .push(GeometryNode::BRep(without_front_face(open).brep()))
        .unwrap();
    let both = b.push(GeometryNode::Collection(vec![solid, open])).unwrap();
    let cut = b
        .push(GeometryNode::SolidOperation(
            axiolid_model::SolidOperation::Boolean {
                left: solid,
                right: open,
                operator: axiolid_core::BooleanOperator::Difference,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![both, cut]).unwrap();
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let outcome = compiler
        .compile_mesh_reported(&graph, both, &options())
        .unwrap();
    assert_eq!(outcome.closure, MeshClosure::OpenSolid);
    let error = compiler
        .compile_mesh(&graph, cut, &options())
        .expect_err("an open solid has no volume to combine");
    assert!(
        error.to_string().contains("tool") && error.to_string().contains("not closed"),
        "{error}"
    );
}

/// A boolean provider that drops the last triangle of the subject: its
/// result is open, whatever it was given.
#[derive(Debug)]
struct Opener;

impl axiolid_contracts::Backend for Opener {
    fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
        axiolid_contracts::BackendDescriptor::new(
            axiolid_contracts::BackendId::new("opener"),
            axiolid_contracts::ExecutionTarget::PortableCpu,
        )
    }
}

impl axiolid_mesh_boolean_contract::MeshBoolean for Opener {
    fn scratch_requirement(&self) -> axiolid_contracts::ScratchRequirement {
        axiolid_contracts::ScratchRequirement::None
    }

    fn boolean(
        &self,
        subject: &TriMesh,
        _tool: &TriMesh,
        _operation: axiolid_core::BooleanOperator,
        _options: &ExecutionOptions,
    ) -> axiolid_contracts::GeomResult<axiolid_mesh_boolean_contract::BooleanOutcome> {
        let mut open = subject.clone();
        open.indices.truncate(open.indices.len() - 3);
        let evidence =
            axiolid_mesh_boolean_contract::BooleanEvidence::record(open.triangle_count(), 0, 0, 0);
        Ok(axiolid_mesh_boolean_contract::BooleanOutcome::new(
            open, evidence,
        ))
    }
}

#[test]
fn a_boolean_result_that_is_open_is_not_reported_solid() {
    // The provider is trusted for its triangles, not for their closure.
    let (body, _) = pocket_touching_the_box_edge();
    let mut b = GeometryGraphBuilder::new();
    let solid = b.push(GeometryNode::BRep(body.brep())).unwrap();
    let cut = b
        .push(GeometryNode::SolidOperation(
            axiolid_model::SolidOperation::Boolean {
                left: solid,
                right: solid,
                operator: axiolid_core::BooleanOperator::Union,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![cut]).unwrap();
    let outcome = ReferenceMeshCompiler::new(Opener)
        .compile_mesh_reported(&graph, cut, &options())
        .unwrap();
    assert_eq!(outcome.closure, MeshClosure::OpenSolid);
}
