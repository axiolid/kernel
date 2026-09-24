//! Surface models: B-reps with shells but no solid (#161).
//!
//! IFC `IfcShellBasedSurfaceModel` / `IfcFaceBasedSurfaceModel` lower to
//! shell-only B-reps: the file authors a surface, never a volume. They must
//! tessellate, and the result must say it is a surface so no consumer reads
//! a volume from it -- including when the shell happens to be closed, where
//! a divergence sum would return a plausible, meaningless number.
//!
//! Oracles are closed-form areas from the input coordinates.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{CompileOutcome, MeshClosure, MeshCompiler};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation};
use axiolid_topology::{
    BRep, Edge, EdgeId, EdgeUse, Face, FaceBound, FaceId, Loop, Orientation, Shell, Solid, Vertex,
};
use std::collections::HashMap;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

/// Faces from corner-index rings over shared corners, edges shared between
/// faces the way a lowered shell shares them.
struct Faces {
    brep: BRep<NodeId>,
    corners: Vec<axiolid_topology::VertexId>,
    edges: HashMap<(usize, usize), EdgeId>,
}

impl Faces {
    fn new(points: &[[f64; 3]]) -> Self {
        let mut brep = BRep::default();
        let corners = points
            .iter()
            .map(|p| {
                brep.add_vertex(Vertex {
                    position: Vec3::new(p[0], p[1], p[2]),
                })
            })
            .collect();
        Self {
            brep,
            corners,
            edges: HashMap::new(),
        }
    }

    fn face(&mut self, ring: &[usize]) -> FaceId {
        let mut uses = Vec::new();
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let key = (a.min(b), a.max(b));
            let (brep, corners) = (&mut self.brep, &self.corners);
            let edge = *self.edges.entry(key).or_insert_with(|| {
                brep.add_edge(Edge {
                    start: corners[key.0],
                    end: corners[key.1],
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
        let wire = self.brep.add_loop(Loop { edges: uses });
        self.brep.add_face(Face {
            surface: None,
            bounds: vec![FaceBound {
                loop_id: wire,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        })
    }

    fn shell(&mut self, faces: &[FaceId], closed: bool) -> axiolid_topology::ShellId {
        self.brep.add_shell(Shell {
            faces: faces.iter().map(|&f| (f, Orientation::Forward)).collect(),
            closed,
        })
    }
}

const CUBE: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [2.0, 2.0, 0.0],
    [0.0, 2.0, 0.0],
    [0.0, 0.0, 2.0],
    [2.0, 0.0, 2.0],
    [2.0, 2.0, 2.0],
    [0.0, 2.0, 2.0],
];
const CUBE_FACES: [[usize; 4]; 6] = [
    [0, 3, 2, 1],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 0, 4, 7],
];

/// A closed 2 x 2 x 2 cube shell; a solid only if `as_solid`.
fn cube(as_solid: bool) -> BRep<NodeId> {
    let mut f = Faces::new(&CUBE);
    let faces: Vec<FaceId> = CUBE_FACES.iter().map(|q| f.face(q)).collect();
    let shell = f.shell(&faces, true);
    if as_solid {
        f.brep.add_solid(Solid {
            outer: shell,
            voids: Vec::new(),
        });
    }
    f.brep
}

fn area(mesh: &TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            (b - a).cross(c - a).length() / 2.0
        })
        .sum()
}

fn compile(brep: BRep<NodeId>) -> CompileOutcome {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(GeometryNode::BRep(brep)).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    compiler()
        .compile_mesh_reported(&graph, root, &options())
        .expect("a surface model tessellates")
}

#[test]
fn an_open_single_face_shell_tessellates_as_a_surface() {
    // One 3 x 2 rectangle, no solid: area 6, no volume.
    let mut f = Faces::new(&[
        [0.0, 0.0, 1.0],
        [3.0, 0.0, 1.0],
        [3.0, 2.0, 1.0],
        [0.0, 2.0, 1.0],
    ]);
    let face = f.face(&[0, 1, 2, 3]);
    f.shell(&[face], false);
    let outcome = compile(f.brep);
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!((area(&outcome.mesh) - 6.0).abs() < 1e-12);
    assert_eq!(outcome.mesh.indices.len(), 6, "two triangles");
    assert!(matches!(
        outcome.solid_mesh(),
        Err(GeomError::InvalidInput(message)) if message.contains("surface model")
    ));
}

#[test]
fn a_closed_shell_without_a_solid_is_still_a_surface() {
    // The mesh is watertight and outward-wound, so a divergence sum would
    // give 8 -- exactly the number the flag exists to stop being read.
    let outcome = compile(cube(false));
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!((area(&outcome.mesh) - 24.0).abs() < 1e-12);
    assert!(
        outcome.solid_mesh().is_err(),
        "volume must be refused, not zero"
    );
}

#[test]
fn the_same_shell_declared_a_solid_is_a_solid() {
    let outcome = compile(cube(true));
    assert_eq!(outcome.closure, MeshClosure::Solid);
    let mesh = outcome.solid_mesh().expect("a declared solid");
    assert!((area(mesh) - 24.0).abs() < 1e-12);
}

#[test]
fn every_shell_of_a_surface_model_is_tessellated() {
    // Two disjoint unit squares in separate shells: area 2, not 1.
    let mut f = Faces::new(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [5.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
        [6.0, 1.0, 0.0],
        [5.0, 1.0, 0.0],
    ]);
    let a = f.face(&[0, 1, 2, 3]);
    let b = f.face(&[4, 5, 6, 7]);
    f.shell(&[a], false);
    f.shell(&[b], false);
    let outcome = compile(f.brep);
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!((area(&outcome.mesh) - 2.0).abs() < 1e-12);
}

#[test]
fn a_brep_with_neither_solid_nor_shell_is_refused() {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(GeometryNode::BRep(BRep::default())).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let error = compiler()
        .compile_mesh(&graph, root, &options())
        .expect_err("nothing to tessellate");
    assert!(
        matches!(&error, GeomError::InvalidInput(m) if m.contains("neither a solid nor a shell")),
        "{error:?}"
    );
}

#[test]
fn a_collection_with_one_surface_member_is_a_surface() {
    // A solid next to a surface model: the whole is not a solid, because
    // a volume summed over it would include the surface's divergence sum.
    let mut b = GeometryGraphBuilder::new();
    let solid = b.push(GeometryNode::BRep(cube(true))).unwrap();
    let surface = b.push(GeometryNode::BRep(cube(false))).unwrap();
    let both = b
        .push(GeometryNode::Collection(vec![solid, surface]))
        .unwrap();
    let only_solids = b.push(GeometryNode::Collection(vec![solid])).unwrap();
    let graph = b.finish(vec![both, only_solids]).unwrap();
    let outcome = |root| {
        compiler()
            .compile_mesh_reported(&graph, root, &options())
            .unwrap()
            .closure
    };
    assert_eq!(outcome(both), MeshClosure::Surface);
    assert_eq!(outcome(only_solids), MeshClosure::Solid);
}

#[test]
fn a_placed_surface_model_stays_a_surface() {
    // Placement moves triangles; it does not make a surface a solid.
    let mut b = GeometryGraphBuilder::new();
    let surface = b.push(GeometryNode::BRep(cube(false))).unwrap();
    let placed = b
        .push(GeometryNode::Instance(axiolid_model::Instance {
            source: surface,
            transform: axiolid_core::Transform3::from_translation(Vec3::new(10.0, 0.0, 0.0)),
        }))
        .unwrap();
    let graph = b.finish(vec![placed]).unwrap();
    let outcome = compiler()
        .compile_mesh_reported(&graph, placed, &options())
        .unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!((area(&outcome.mesh) - 24.0).abs() < 1e-12);
}

#[test]
fn a_boolean_with_a_surface_operand_is_refused() {
    let mut b = GeometryGraphBuilder::new();
    let solid = b.push(GeometryNode::BRep(cube(true))).unwrap();
    let surface = b.push(GeometryNode::BRep(cube(false))).unwrap();
    let as_tool = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: solid,
            right: surface,
            operator: BooleanOperator::Difference,
        }))
        .unwrap();
    let as_subject = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: surface,
            right: solid,
            operator: BooleanOperator::Union,
        }))
        .unwrap();
    let graph = b.finish(vec![as_tool, as_subject]).unwrap();
    for (root, role) in [(as_tool, "tool"), (as_subject, "subject")] {
        let error = compiler()
            .compile_mesh(&graph, root, &options())
            .expect_err("a surface has no volume to combine");
        assert!(
            matches!(&error, GeomError::InvalidInput(m) if m.contains(role) && m.contains("surface model")),
            "{error:?}"
        );
    }
}

#[test]
fn authored_meshes_report_their_closure_from_their_structure() {
    // A closed cube as a polygon mesh bounds a solid; the same with one
    // face missing is an open sheet.
    let positions: Vec<_> = CUBE
        .iter()
        .map(|p| axiolid_core::Point3::new(p[0], p[1], p[2]))
        .collect();
    let faces = |count: usize| {
        CUBE_FACES[..count]
            .iter()
            .map(|q| axiolid_mesh::PolygonFace {
                outer: q.iter().map(|&i| i as u32).collect(),
                holes: Vec::new(),
            })
            .collect()
    };
    let mut b = GeometryGraphBuilder::new();
    let closed = b
        .push(GeometryNode::PolygonMesh(axiolid_mesh::PolygonMesh {
            positions: positions.clone(),
            faces: faces(6),
        }))
        .unwrap();
    let open = b
        .push(GeometryNode::PolygonMesh(axiolid_mesh::PolygonMesh {
            positions,
            faces: faces(5),
        }))
        .unwrap();
    let graph = b.finish(vec![closed, open]).unwrap();
    let closure = |root| {
        compiler()
            .compile_mesh_reported(&graph, root, &options())
            .unwrap()
            .closure
    };
    assert_eq!(closure(closed), MeshClosure::Solid);
    assert_eq!(closure(open), MeshClosure::Surface);
}

#[test]
fn an_untracked_outcome_is_not_a_solid() {
    let outcome = CompileOutcome::untracked(TriMesh::default());
    assert_eq!(outcome.closure, MeshClosure::Unknown);
    assert!(outcome.solid_mesh().is_err());
}
