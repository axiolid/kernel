//! Faceted B-rep tessellation gates.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::Vec3;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode};
use axiolid_topology::audit_brep;
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};
use std::collections::HashMap;

/// Build a unit cube as topology, the way lowering produces it.
fn cube() -> BRep<axiolid_model::NodeId> {
    cube_with_sense(Orientation::Forward)
}

/// Same cube, but every face registered with the given shell sense.
fn cube_with_sense(sense: Orientation) -> BRep<axiolid_model::NodeId> {
    let (mut brep, shell) = cube_shell(sense);
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    brep
}

/// The cube's faces and shell, leaving solid registration to the caller.
fn cube_shell(sense: Orientation) -> (BRep<axiolid_model::NodeId>, axiolid_topology::ShellId) {
    cube_shell_with_duplicate_outer(sense, false)
}

fn cube_shell_with_duplicate_outer(
    sense: Orientation,
    duplicate_outer: bool,
) -> (BRep<axiolid_model::NodeId>, axiolid_topology::ShellId) {
    let mut brep = BRep::default();
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let verts: Vec<_> = corners
        .iter()
        .map(|c| {
            brep.add_vertex(Vertex {
                position: Vec3::new(c[0], c[1], c[2]),
            })
        })
        .collect();
    // Outward-facing quads (CCW seen from outside).
    let quads = [
        [0usize, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ];
    let mut edges: HashMap<(usize, usize), axiolid_topology::EdgeId> = HashMap::new();
    let mut faces = Vec::new();
    for (face_index, quad) in quads.into_iter().enumerate() {
        let mut uses = Vec::new();
        for i in 0..4 {
            let (a, b) = (quad[i], quad[(i + 1) % 4]);
            let key = if a < b { (a, b) } else { (b, a) };
            let id = *edges.entry(key).or_insert_with(|| {
                brep.add_edge(Edge {
                    start: verts[key.0],
                    end: verts[key.1],
                    curve: None,
                })
            });
            let orientation = if a == key.0 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            uses.push(EdgeUse {
                edge: id,
                orientation,
                pcurve: None,
            });
        }
        let wire = brep.add_loop(Loop { edges: uses });
        let outer = FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        };
        let mut bounds = vec![outer];
        if duplicate_outer && face_index == 0 {
            bounds.push(outer);
        }
        faces.push(brep.add_face(Face {
            surface: None,
            bounds,
            orientation: Orientation::Forward,
        }));
    }
    let shell = brep.add_shell(Shell {
        faces: faces.iter().map(|f| (*f, sense)).collect(),
        closed: true,
    });
    (brep, shell)
}

fn planar_face(points: &[Vec3]) -> BRep<axiolid_model::NodeId> {
    planar_face_with_holes(points, &[])
}

fn planar_face_with_holes(outer: &[Vec3], holes: &[&[Vec3]]) -> BRep<axiolid_model::NodeId> {
    planar_face_with_bound_roles(outer, holes, true)
}

fn planar_face_without_outer(points: &[Vec3]) -> BRep<axiolid_model::NodeId> {
    planar_face_with_bound_roles(points, &[], false)
}

fn planar_face_with_bound_roles(
    outer: &[Vec3],
    holes: &[&[Vec3]],
    mark_outer: bool,
) -> BRep<axiolid_model::NodeId> {
    let mut brep = BRep::default();
    let mut bounds = Vec::new();
    for (ring_index, points) in core::iter::once(outer)
        .chain(holes.iter().copied())
        .enumerate()
    {
        let vertices: Vec<_> = points
            .iter()
            .map(|&position| brep.add_vertex(Vertex { position }))
            .collect();
        let edges: Vec<_> = (0..vertices.len())
            .map(|index| {
                brep.add_edge(Edge {
                    start: vertices[index],
                    end: vertices[(index + 1) % vertices.len()],
                    curve: None,
                })
            })
            .collect();
        let wire = brep.add_loop(Loop {
            edges: edges
                .into_iter()
                .map(|edge| EdgeUse {
                    edge,
                    orientation: Orientation::Forward,
                    pcurve: None,
                })
                .collect(),
        });
        bounds.push(FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: mark_outer && ring_index == 0,
        });
    }
    let face = brep.add_face(Face {
        surface: None,
        bounds,
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    brep
}

fn compile_result(
    brep: BRep<axiolid_model::NodeId>,
) -> axiolid_contracts::GeomResult<axiolid_mesh::TriMesh> {
    let mut builder = GeometryGraphBuilder::new();
    let root = builder.push(GeometryNode::BRep(brep)).expect("push");
    let graph = builder.finish(vec![root]).expect("finish");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    compiler.compile_mesh(
        &graph,
        root,
        &ExecutionOptions::new(axiolid_core::Tolerance::METRE),
    )
}

fn compile(brep: BRep<axiolid_model::NodeId>) -> axiolid_mesh::TriMesh {
    compile_result(brep).expect("cube compiles")
}

/// A cube tessellates to a welded, closed, outward mesh.
///
/// Per-face vertex copies would still render correctly but leave every edge
/// unshared, so the weld is asserted through edge parity rather than by
/// counting positions alone.
#[test]
fn a_solid_with_an_empty_outer_shell_is_rejected() {
    let mut brep = BRep::default();
    let outer = brep.add_shell(Shell {
        faces: Vec::new(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer,
        voids: Vec::new(),
    });

    let error = compile_result(brep).expect_err("a solid cannot have an empty outer shell");
    assert!(
        matches!(error, axiolid_contracts::GeomError::InvalidInput(_)),
        "expected malformed outer shell, got {error:?}"
    );
}

#[test]
fn a_face_without_an_outer_bound_is_rejected_before_tessellation() {
    let points = [Vec3::ZERO, Vec3::X, Vec3::new(1.0, 1.0, 0.0), Vec3::Y];
    let error = compile_result(planar_face_without_outer(&points))
        .expect_err("a face must designate one outer bound");
    match error {
        axiolid_contracts::GeomError::InvalidInput(message) => assert!(
            message.contains("faces_without_outer_bound: 1"),
            "diagnostic must name the missing outer bound: {message}"
        ),
        other => panic!("expected invalid topology, got {other:?}"),
    }
}

#[test]
fn a_face_with_an_empty_outer_loop_is_rejected_before_tessellation() {
    let mut brep = BRep::default();
    let empty = brep.add_loop(Loop { edges: Vec::new() });
    let face = brep.add_face(Face {
        surface: None,
        bounds: vec![FaceBound {
            loop_id: empty,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let outer = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer,
        voids: Vec::new(),
    });

    let error = compile_result(brep).expect_err("an empty loop cannot bound a face");
    match error {
        axiolid_contracts::GeomError::InvalidInput(message) => assert!(
            message.contains("empty_loops: 1"),
            "diagnostic must name the empty loop: {message}"
        ),
        other => panic!("expected invalid topology, got {other:?}"),
    }
}

#[test]
fn non_finite_planar_bounds_are_rejected() {
    for (label, coordinate) in [("NaN", f64::NAN), ("infinity", f64::INFINITY)] {
        let points = [Vec3::ZERO, Vec3::X, Vec3::new(0.0, coordinate, 0.0)];
        let error = compile_result(planar_face(&points))
            .expect_err("non-finite planar input cannot be tessellated");
        match error {
            axiolid_contracts::GeomError::Degenerate(message) => assert!(
                message.contains("non-finite"),
                "{label} diagnostic must identify non-finite area: {message}"
            ),
            other => panic!("expected degenerate {label} bound, got {other:?}"),
        }
    }
}

#[test]
fn multiple_outer_bounds_are_rejected_before_tessellation() {
    let (mut brep, outer) = cube_shell_with_duplicate_outer(Orientation::Forward, true);
    brep.add_solid(Solid {
        outer,
        voids: Vec::new(),
    });

    let error = compile_result(brep).expect_err("multiple outer bounds must fail closed");
    match error {
        axiolid_contracts::GeomError::InvalidInput(message) => assert!(
            message.contains("faces_with_multiple_outer_bounds: 1"),
            "diagnostic must name the malformed face: {message}"
        ),
        other => panic!("expected invalid topology, got {other:?}"),
    }
}

#[test]
fn a_two_edge_planar_bound_is_rejected() {
    let brep = planar_face(&[Vec3::ZERO, Vec3::X]);
    let error = compile_result(brep).expect_err("two edges cannot bound a planar face");
    assert!(
        matches!(error, axiolid_contracts::GeomError::Degenerate(_)),
        "expected degenerate planar bound, got {error:?}"
    );
}

#[test]
fn a_zero_area_planar_bound_is_rejected() {
    let brep = planar_face(&[Vec3::ZERO, Vec3::X, Vec3::X * 2.0]);
    let error = compile_result(brep).expect_err("collinear vertices define no planar face");
    assert!(
        matches!(error, axiolid_contracts::GeomError::Degenerate(_)),
        "expected zero-area planar bound, got {error:?}"
    );
}

#[test]
fn a_two_edge_inner_bound_is_rejected() {
    let outer = [
        Vec3::ZERO,
        Vec3::X * 4.0,
        Vec3::new(4.0, 4.0, 0.0),
        Vec3::Y * 4.0,
    ];
    let hole = [Vec3::new(1.0, 1.0, 0.0), Vec3::new(2.0, 1.0, 0.0)];
    let error = compile_result(planar_face_with_holes(&outer, &[&hole]))
        .expect_err("two edges cannot bound a planar hole");
    assert!(
        matches!(error, axiolid_contracts::GeomError::Degenerate(_)),
        "expected undersized inner bound, got {error:?}"
    );
}

#[test]
fn a_zero_area_inner_bound_is_rejected() {
    let outer = [
        Vec3::ZERO,
        Vec3::X * 4.0,
        Vec3::new(4.0, 4.0, 0.0),
        Vec3::Y * 4.0,
    ];
    let hole = [
        Vec3::new(1.0, 1.0, 0.0),
        Vec3::new(2.0, 1.0, 0.0),
        Vec3::new(3.0, 1.0, 0.0),
    ];
    let error = compile_result(planar_face_with_holes(&outer, &[&hole]))
        .expect_err("a collinear inner bound cannot define a hole");
    assert!(
        matches!(error, axiolid_contracts::GeomError::Degenerate(_)),
        "expected zero-area inner bound, got {error:?}"
    );
}

#[test]
fn a_cube_tessellates_to_a_closed_welded_mesh() {
    let mesh = compile(cube());

    assert_eq!(
        mesh.positions.len(),
        8,
        "one position per topological vertex"
    );
    assert_eq!(
        mesh.indices.len(),
        12 * 3,
        "six quads become twelve triangles"
    );

    let mut directed: HashMap<(u32, u32), i32> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for &(a, b) in &[(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((a, b)).or_default() += 1;
        }
    }
    for (&(a, b), &count) in &directed {
        assert_eq!(count, 1, "directed edge {a}->{b} repeats");
        assert_eq!(
            directed.get(&(b, a)).copied().unwrap_or(0),
            1,
            "edge {a}-{b} unpaired"
        );
    }
}

/// Signed volume is positive: the winding survived tessellation.
///
/// A flipped face is invisible to a triangle count and to edge parity, but
/// inverts the divergence integral.
#[test]
fn tessellated_winding_stays_outward() {
    let mesh = compile(cube());
    let mut volume = 0.0;
    for t in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        volume += a.dot(b.cross(c)) / 6.0;
    }
    assert!(
        (volume - 1.0).abs() < 1e-9,
        "unit cube volume, got {volume}"
    );
}

/// Add an axis-aligned box's six quads as a shell of `brep`, wound outward
/// from the box when `outward`, into it otherwise.
fn add_box_shell(
    brep: &mut BRep<axiolid_model::NodeId>,
    lo: f64,
    hi: f64,
    outward: bool,
) -> axiolid_topology::ShellId {
    let verts: Vec<_> = (0..8)
        .map(|i| {
            let pick = |bit: usize| if i >> bit & 1 == 1 { hi } else { lo };
            brep.add_vertex(Vertex {
                position: Vec3::new(pick(0), pick(1), pick(2)),
            })
        })
        .collect();
    // Corner index bits are (x, y, z); these quads run CCW seen from outside.
    let quads = [
        [0usize, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [1, 3, 7, 5],
        [3, 2, 6, 7],
        [2, 0, 4, 6],
    ];
    let mut edges: HashMap<(usize, usize), axiolid_topology::EdgeId> = HashMap::new();
    let mut faces = Vec::new();
    for mut quad in quads {
        if !outward {
            quad.reverse();
        }
        let mut uses = Vec::new();
        for i in 0..4 {
            let (a, b) = (quad[i], quad[(i + 1) % 4]);
            let key = if a < b { (a, b) } else { (b, a) };
            let id = *edges.entry(key).or_insert_with(|| {
                brep.add_edge(Edge {
                    start: verts[key.0],
                    end: verts[key.1],
                    curve: None,
                })
            });
            uses.push(EdgeUse {
                edge: id,
                orientation: if a == key.0 {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
                pcurve: None,
            });
        }
        let wire = brep.add_loop(Loop { edges: uses });
        faces.push((
            brep.add_face(Face {
                surface: None,
                bounds: vec![FaceBound {
                    loop_id: wire,
                    orientation: Orientation::Forward,
                    outer: true,
                }],
                orientation: Orientation::Forward,
            }),
            Orientation::Forward,
        ));
    }
    brep.add_shell(Shell {
        faces,
        closed: true,
    })
}

fn signed_volume(mesh: &axiolid_mesh::TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// A void shell is a cavity: it is tessellated, facing into the cavity, so
/// the mesh encloses the outer volume less the void (#120). It used to be
/// dropped as "boolean intent", which silently filled every cavity.
#[test]
fn a_void_shell_is_tessellated_as_a_cavity() {
    let (mut brep, outer) = cube_shell(Orientation::Forward);
    // Wound into the cavity, away from the material: the B-rep sense.
    let void = add_box_shell(&mut brep, 0.25, 0.75, false);
    brep.add_solid(Solid {
        outer,
        voids: vec![void],
    });

    let mesh = compile(brep);
    assert_eq!(mesh.positions.len(), 16, "both shells' corners");
    assert_eq!(mesh.indices.len(), 72, "both shells' triangles");
    let volume = signed_volume(&mesh);
    assert!(
        (volume - (1.0 - 0.125)).abs() < 1e-12,
        "cube less its cavity, got {volume}"
    );
}

/// A void authored facing out of the cavity (the STEP convention, reversed
/// on use) still removes material: a cavity cannot add volume.
#[test]
fn a_void_authored_outward_still_subtracts() {
    let (mut brep, outer) = cube_shell(Orientation::Forward);
    let void = add_box_shell(&mut brep, 0.25, 0.75, true);
    brep.add_solid(Solid {
        outer,
        voids: vec![void],
    });
    let volume = signed_volume(&compile(brep));
    assert!(
        (volume - (1.0 - 0.125)).abs() < 1e-12,
        "cube less its cavity, got {volume}"
    );
}

/// A void shell that does not close would leave a hole in the mesh.
#[test]
fn an_open_void_shell_is_refused() {
    let (mut brep, outer) = cube_shell(Orientation::Forward);
    let inner: Vec<_> = [
        [0.25, 0.25, 0.25],
        [0.75, 0.25, 0.25],
        [0.75, 0.75, 0.25],
        [0.25, 0.75, 0.25],
    ]
    .iter()
    .map(|c| {
        brep.add_vertex(Vertex {
            position: Vec3::new(c[0], c[1], c[2]),
        })
    })
    .collect();
    let mut uses = Vec::new();
    for index in 0..4 {
        let edge = brep.add_edge(Edge {
            start: inner[index],
            end: inner[(index + 1) % 4],
            curve: None,
        });
        uses.push(EdgeUse {
            edge,
            orientation: Orientation::Forward,
            pcurve: None,
        });
    }
    let wire = brep.add_loop(Loop { edges: uses });
    let face = brep.add_face(Face {
        surface: None,
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let void = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: true,
    });
    brep.add_solid(Solid {
        outer,
        voids: vec![void],
    });

    let error = compile_result(brep).expect_err("one face bounds no cavity");
    assert!(
        format!("{error:?}").contains("void shell is not closed"),
        "got {error:?}"
    );
}

/// A reversed shell sense flips the emitted winding.
///
/// The cube is symmetric enough that a missed flip still yields a closed
/// mesh, so orientation is measured by signed volume: inverting every face
/// negates it.
#[test]
fn a_reversed_shell_inverts_the_signed_volume() {
    let mesh = compile(cube_with_sense(Orientation::Reversed));
    let mut volume = 0.0;
    for t in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        volume += a.dot(b.cross(c)) / 6.0;
    }
    assert!(
        (volume + 1.0).abs() < 1e-9,
        "reversed shell must invert volume, got {volume}"
    );
}
/// A concave face triangulates by area, not by a fan from vertex zero.
///
/// An L-shape has a reflex corner and its first three vertices are collinear-
/// adjacent, so both a naive fan and a first-two-edge normal fail here while
/// looking fine on a convex quad.
#[test]
fn a_concave_face_keeps_its_true_area() {
    let mut brep = BRep::default();
    // L-shape in the z=0 plane, CCW, area 3 of a 2x2 square.
    // A collinear run at the start (0,0)->(1,0)->(2,0) makes a first-two-edge
    // normal degenerate, which is exactly the case Newell must survive.
    let ring = [
        [0.0, 0.0],
        [1.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let verts: Vec<_> = ring
        .iter()
        .map(|c| {
            brep.add_vertex(Vertex {
                position: Vec3::new(c[0], c[1], 0.0),
            })
        })
        .collect();
    let mut uses = Vec::new();
    for index in 0..verts.len() {
        let edge = brep.add_edge(Edge {
            start: verts[index],
            end: verts[(index + 1) % verts.len()],
            curve: None,
        });
        uses.push(EdgeUse {
            edge,
            orientation: Orientation::Forward,
            pcurve: None,
        });
    }
    let wire = brep.add_loop(Loop { edges: uses });
    let face = brep.add_face(Face {
        surface: None,
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let mesh = compile(brep);
    let mut area = 0.0;
    for t in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        area += (b - a).cross(c - a).length() / 2.0;
    }
    assert!((area - 3.0).abs() < 1e-9, "L-shape area is 3, got {area}");
    assert_eq!(mesh.indices.len(), 5 * 3, "a 7-gon becomes 5 triangles");
}
/// A reversed bound flips that loop, which flips the face normal.
///
/// IfcFaceBound.Orientation of .F. means traverse the loop backwards. Ignoring
/// it silently mirrors the facet.
#[test]
fn a_reversed_bound_reverses_the_loop() {
    fn triangle(bound_sense: Orientation) -> axiolid_mesh::TriMesh {
        let mut brep = BRep::default();
        let pts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let verts: Vec<_> = pts
            .iter()
            .map(|c| {
                brep.add_vertex(Vertex {
                    position: Vec3::new(c[0], c[1], 0.0),
                })
            })
            .collect();
        let mut uses = Vec::new();
        for index in 0..3 {
            let edge = brep.add_edge(Edge {
                start: verts[index],
                end: verts[(index + 1) % 3],
                curve: None,
            });
            uses.push(EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: None,
            });
        }
        let wire = brep.add_loop(Loop { edges: uses });
        let face = brep.add_face(Face {
            surface: None,
            bounds: vec![FaceBound {
                loop_id: wire,
                orientation: bound_sense,
                outer: true,
            }],
            orientation: Orientation::Forward,
        });
        let shell = brep.add_shell(Shell {
            faces: vec![(face, Orientation::Forward)],
            closed: false,
        });
        brep.add_solid(Solid {
            outer: shell,
            voids: Vec::new(),
        });
        compile(brep)
    }

    fn normal_z(mesh: &axiolid_mesh::TriMesh) -> f64 {
        let t = &mesh.indices[0..3];
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        (b - a).cross(c - a).z
    }

    let forward = normal_z(&triangle(Orientation::Forward));
    let reversed = normal_z(&triangle(Orientation::Reversed));
    assert!(forward > 0.0, "forward bound faces +Z, got {forward}");
    assert!(
        reversed < 0.0,
        "reversed bound must flip the normal, got {reversed}"
    );
}

// --- curved support surfaces ------------------------------------------------

/// A face with a curved support must be refused, not silently faceted.
///
/// Projecting a cylindrical face onto the plane of its boundary yields a mesh
/// that looks valid and is wrong everywhere between the vertices.
#[test]
fn a_curved_face_is_refused_not_flattened() {
    let brep = cube();
    let mut builder = GeometryGraphBuilder::new();
    let cyl = builder
        .push(GeometryNode::Surface(axiolid_surface::Surface::Cylinder(
            axiolid_surface::Cylinder {
                frame: axiolid_core::Frame3 {
                    origin: axiolid_core::Point3::ZERO,
                    x: axiolid_core::Vec3::X,
                    y: axiolid_core::Vec3::Y,
                    z: axiolid_core::Vec3::Z,
                },
                radius: 1.0,
            },
        )))
        .expect("push surface");
    // Attach the curved support to the first face.
    let mut faces: Vec<_> = brep.faces().to_vec();
    faces[0].surface = Some(cyl);
    let mut rebuilt = BRep::default();
    for v in brep.vertices() {
        rebuilt.add_vertex(*v);
    }
    for e in brep.edges() {
        rebuilt.add_edge(e.clone());
    }
    for l in brep.loops() {
        rebuilt.add_loop(l.clone());
    }
    for f in faces {
        rebuilt.add_face(f);
    }
    for s in brep.shells() {
        rebuilt.add_shell(s.clone());
    }
    for s in brep.solids() {
        rebuilt.add_solid(s.clone());
    }
    let root = builder
        .push(GeometryNode::BRep(rebuilt))
        .expect("push brep");
    let graph = builder.finish(vec![root]).expect("finish");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let outcome = compiler.compile_mesh(
        &graph,
        root,
        &ExecutionOptions::new(axiolid_core::Tolerance::METRE),
    );
    assert!(
        outcome.is_err(),
        "a cylindrical face must be refused, not flattened: {outcome:?}"
    );
}

/// The cube the tessellator actually consumes must audit as a closed
/// manifold. If it does not, either the audit or the fixture is wrong, and
/// every downstream volume claim rests on the answer.
#[test]
fn the_tessellation_cube_is_a_closed_manifold() {
    let health = audit_brep(&cube());
    assert!(
        health.is_closed_manifold(),
        "the cube used for tessellation must be sound: {health:?}"
    );
}

/// A curved face WITH a pcurve is sampled on its surface, not flattened.
///
/// The trim states the face boundary in `(u, v)`, so the sampler knows which
/// part of the cylinder to emit. Radius is the check that it really sampled
/// the surface: every vertex must sit `r` from the axis, which a planar
/// projection through the boundary could not achieve.
#[test]
fn a_curved_face_with_a_pcurve_is_sampled_on_its_surface() {
    let radius = 2.0;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(axiolid_surface::Surface::Cylinder(
            axiolid_surface::Cylinder {
                frame: axiolid_core::Frame3 {
                    origin: axiolid_core::Point3::ZERO,
                    x: axiolid_core::Vec3::X,
                    y: axiolid_core::Vec3::Y,
                    z: axiolid_core::Vec3::Z,
                },
                radius,
            },
        )))
        .expect("surface");
    // A rectangle in parameter space: a quarter turn, two metres tall.
    let quarter = std::f64::consts::FRAC_PI_2;
    let trim = builder
        .push(GeometryNode::Curve2(axiolid_curve::Curve2::Polyline(
            axiolid_curve::Polyline2 {
                points: vec![
                    axiolid_core::Point2::new(0.0, 0.0),
                    axiolid_core::Point2::new(quarter, 0.0),
                    axiolid_core::Point2::new(quarter, 2.0),
                    axiolid_core::Point2::new(0.0, 2.0),
                ],
                closed: true,
            },
        )))
        .expect("trim");

    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let v: Vec<_> = (0..4)
        .map(|i| {
            let a = quarter * f64::from(i % 2);
            brep.add_vertex(Vertex {
                position: axiolid_core::Point3::new(radius * a.cos(), radius * a.sin(), 0.0),
            })
        })
        .collect();
    let e: Vec<_> = (0..4)
        .map(|i| {
            brep.add_edge(Edge {
                start: v[i],
                end: v[(i + 1) % 4],
                curve: None,
            })
        })
        .collect();
    let wire = brep.add_loop(Loop {
        edges: e
            .iter()
            .map(|&edge| EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(trim),
            })
            .collect(),
    });
    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).expect("brep");
    let graph = builder.finish(vec![root]).expect("finish");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let mesh = compiler
        .compile_mesh(
            &graph,
            root,
            &ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
        )
        .expect("a curved face with a trim must tessellate");

    assert!(
        mesh.positions.len() > 8,
        "a sampled cylinder needs more than the boundary vertices, got {}",
        mesh.positions.len()
    );
    for p in &mesh.positions {
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (r - radius).abs() < 1e-9,
            "vertex {p:?} is {r} from the axis, not {radius}: the face was \
             flattened rather than sampled"
        );
    }
}

/// Broken topology is refused before any triangle is emitted.
///
/// An open loop cannot bound a face, but the planar path would happily
/// project its points and return a plausible-looking mesh. The audit runs
/// first precisely so that never reaches geometry.
#[test]
fn unsound_topology_is_refused_before_tessellation() {
    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let v: Vec<_> = [
        axiolid_core::Point3::new(0.0, 0.0, 0.0),
        axiolid_core::Point3::new(1.0, 0.0, 0.0),
        axiolid_core::Point3::new(1.0, 1.0, 0.0),
        axiolid_core::Point3::new(0.0, 1.0, 0.0),
    ]
    .into_iter()
    .map(|position| brep.add_vertex(Vertex { position }))
    .collect();
    // Three disjoint edges: the loop cannot close.
    let e0 = brep.add_edge(Edge {
        start: v[0],
        end: v[1],
        curve: None,
    });
    let e1 = brep.add_edge(Edge {
        start: v[2],
        end: v[3],
        curve: None,
    });
    let wire = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: e0,
                orientation: Orientation::Forward,
                pcurve: None,
            },
            EdgeUse {
                edge: e1,
                orientation: Orientation::Forward,
                pcurve: None,
            },
        ],
    });
    let face = brep.add_face(Face {
        surface: None,
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let mut builder = GeometryGraphBuilder::new();
    let root = builder.push(GeometryNode::BRep(brep)).expect("push");
    let graph = builder.finish(vec![root]).expect("finish");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let result = compiler.compile_mesh(
        &graph,
        root,
        &ExecutionOptions::new(axiolid_core::Tolerance::METRE),
    );
    assert!(
        result.is_err(),
        "an open loop must be refused, not tessellated into a plausible mesh"
    );
}

/// Two curved faces meeting at a seam must share its vertices.
///
/// Each face owns its own pcurve. Evaluating both independently lands on the
/// same 3D curve at different parameters, so the seam looks closed and is
/// not. Sampling each edge once and reusing the indices is what closes it.
#[test]
fn two_curved_faces_share_their_seam_vertices() {
    use axiolid_core::{Point2, Point3, Scalar, Vec3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Cylinder, Surface};

    let mut builder = GeometryGraphBuilder::new();
    let radius = 2.0;
    let frame = axiolid_core::Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let surface = builder
        .push(GeometryNode::Surface(Surface::Cylinder(Cylinder {
            frame,
            radius,
        })))
        .expect("surface");

    // Half-cylinder split at u = PI. Both faces use the same seam EDGE,
    // each with its own pcurve node.
    let pi = core::f64::consts::PI;
    let line = |b: &mut GeometryGraphBuilder, a: Point2, c: Point2| {
        b.push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![a, c],
            closed: false,
        })))
        .expect("pcurve")
    };

    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let p = |u: Scalar, v: Scalar| Point3::new(radius * u.cos(), radius * u.sin(), v);
    let v00 = brep.add_vertex(Vertex {
        position: p(0.0, 0.0),
    });
    let v01 = brep.add_vertex(Vertex {
        position: p(0.0, 3.0),
    });
    let v10 = brep.add_vertex(Vertex {
        position: p(pi, 0.0),
    });
    let v11 = brep.add_vertex(Vertex {
        position: p(pi, 3.0),
    });

    // The shared seam edge, plus the three others bounding face A.
    let seam = brep.add_edge(Edge {
        start: v10,
        end: v11,
        curve: None,
    });
    let bottom = brep.add_edge(Edge {
        start: v00,
        end: v10,
        curve: None,
    });
    let top = brep.add_edge(Edge {
        start: v01,
        end: v11,
        curve: None,
    });
    let start = brep.add_edge(Edge {
        start: v00,
        end: v01,
        curve: None,
    });

    // Face A trims: u from 0 to PI.
    let a_bottom = line(&mut builder, Point2::new(0.0, 0.0), Point2::new(pi, 0.0));
    let a_seam = line(&mut builder, Point2::new(pi, 0.0), Point2::new(pi, 3.0));
    let a_top = line(&mut builder, Point2::new(pi, 3.0), Point2::new(0.0, 3.0));
    let a_start = line(&mut builder, Point2::new(0.0, 3.0), Point2::new(0.0, 0.0));

    let loop_a = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: bottom,
                orientation: Orientation::Forward,
                pcurve: Some(a_bottom),
            },
            EdgeUse {
                edge: seam,
                orientation: Orientation::Forward,
                pcurve: Some(a_seam),
            },
            EdgeUse {
                edge: top,
                orientation: Orientation::Reversed,
                pcurve: Some(a_top),
            },
            EdgeUse {
                edge: start,
                orientation: Orientation::Reversed,
                pcurve: Some(a_start),
            },
        ],
    });
    let face_a = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: loop_a,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });

    // Face B covers u from PI to TAU and walks the SAME seam edge backwards.
    let tau = core::f64::consts::TAU;
    let v20 = brep.add_vertex(Vertex {
        position: p(tau, 0.0),
    });
    let v21 = brep.add_vertex(Vertex {
        position: p(tau, 3.0),
    });
    let b_bottom_e = brep.add_edge(Edge {
        start: v10,
        end: v20,
        curve: None,
    });
    let b_top_e = brep.add_edge(Edge {
        start: v11,
        end: v21,
        curve: None,
    });
    let b_end_e = brep.add_edge(Edge {
        start: v20,
        end: v21,
        curve: None,
    });

    let b_seam = line(&mut builder, Point2::new(pi, 3.0), Point2::new(pi, 0.0));
    let b_bottom = line(&mut builder, Point2::new(pi, 0.0), Point2::new(tau, 0.0));
    let b_end = line(&mut builder, Point2::new(tau, 0.0), Point2::new(tau, 3.0));
    let b_top = line(&mut builder, Point2::new(tau, 3.0), Point2::new(pi, 3.0));

    let loop_b = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: seam,
                orientation: Orientation::Reversed,
                pcurve: Some(b_seam),
            },
            EdgeUse {
                edge: b_bottom_e,
                orientation: Orientation::Forward,
                pcurve: Some(b_bottom),
            },
            EdgeUse {
                edge: b_end_e,
                orientation: Orientation::Forward,
                pcurve: Some(b_end),
            },
            EdgeUse {
                edge: b_top_e,
                orientation: Orientation::Reversed,
                pcurve: Some(b_top),
            },
        ],
    });
    let face_b = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: loop_b,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });

    let shell = brep.add_shell(Shell {
        faces: vec![
            (face_a, Orientation::Forward),
            (face_b, Orientation::Forward),
        ],
        closed: false,
    });
    let _ = brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).expect("push");
    let graph = builder.finish(vec![root]).expect("finish");
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    let mesh = compiler
        .compile_mesh(
            &graph,
            root,
            &ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
        )
        .expect("two curved faces compile");

    // Every vertex is on the cylinder.
    for q in &mesh.positions {
        let r = (q.x * q.x + q.y * q.y).sqrt();
        assert!((r - radius).abs() < 1e-9, "vertex off the cylinder: {r}");
    }

    // The seam is shared: the interior edge along u = PI must be used by
    // exactly two triangles, one from each face. Duplicated seam vertices
    // would make it two separate boundary edges instead.
    let mut edges: std::collections::HashMap<(u32, u32), i32> = std::collections::HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *edges.entry(key).or_default() += if a < b { 1 } else { -1 };
        }
    }
    let interior = edges.values().filter(|v| **v == 0).count();
    assert!(
        interior > 0,
        "no interior edge is shared: the seam was duplicated"
    );
}

/// A full cylinder closes across its periodic seam.
///
/// The seam edge appears TWICE in one loop: once at u = 0 and once at
/// u = TAU. Both uses name the same EdgeId, so both must resolve to the
/// same 3D vertices or the tube is split down its length.
#[test]
fn a_periodic_seam_closes_the_tube() {
    use axiolid_core::{Point2, Vec3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Cylinder, Surface};
    use core::f64::consts::TAU;

    let radius = 2.0;
    let height = 3.0;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(Surface::Cylinder(Cylinder {
            frame: axiolid_core::Frame3 {
                origin: axiolid_core::Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius,
        })))
        .expect("surface");
    // Trims in (u, v). The loop walks: bottom rim u 0->TAU, seam up at
    // u = TAU, top rim back TAU->0, seam down at u = 0. The two seam uses
    // name ONE edge, so the tube closes only if both resolve alike.
    let bottom = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(0.0, 0.0), Point2::new(TAU, 0.0)],
            closed: false,
        })))
        .expect("bottom");
    let top = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(TAU, height), Point2::new(0.0, height)],
            closed: false,
        })))
        .expect("top");
    // The seam's two uses have DIFFERENT pcurves -- u = TAU going up and
    // u = 0 coming down -- but name the same topological edge.
    let seam_up = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(TAU, 0.0), Point2::new(TAU, height)],
            closed: false,
        })))
        .expect("seam up");
    let seam_down = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(0.0, height), Point2::new(0.0, 0.0)],
            closed: false,
        })))
        .expect("seam down");
    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let p = |u: f64, v: f64| axiolid_core::Point3::new(radius * u.cos(), radius * u.sin(), v);
    let v00 = brep.add_vertex(Vertex {
        position: p(0.0, 0.0),
    });
    let v01 = brep.add_vertex(Vertex {
        position: p(0.0, height),
    });
    // Rims start and end at the seam vertices: a periodic loop has no
    // other corners.
    let e_bottom = brep.add_edge(Edge {
        start: v00,
        end: v00,
        curve: None,
    });
    let e_seam = brep.add_edge(Edge {
        start: v00,
        end: v01,
        curve: None,
    });
    let e_top = brep.add_edge(Edge {
        start: v01,
        end: v01,
        curve: None,
    });
    let wire = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: e_bottom,
                orientation: Orientation::Forward,
                pcurve: Some(bottom),
            },
            EdgeUse {
                edge: e_seam,
                orientation: Orientation::Forward,
                pcurve: Some(seam_up),
            },
            EdgeUse {
                edge: e_top,
                orientation: Orientation::Forward,
                pcurve: Some(top),
            },
            EdgeUse {
                edge: e_seam,
                orientation: Orientation::Reversed,
                pcurve: Some(seam_down),
            },
        ],
    });
    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).expect("brep");
    let graph = builder.finish(vec![root]).expect("finish");
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(
            &graph,
            root,
            &ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
        )
        .expect("periodic cylinder tessellates");

    // Every vertex sits on the cylinder.
    for p in &mesh.positions {
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (r - radius).abs() < 1e-6,
            "vertex off the cylinder: r = {r}"
        );
    }
    // The seam must be interior: no vertical run of boundary edges at u = 0.
    let mut edges: std::collections::HashMap<(u32, u32), i32> = Default::default();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *edges.entry(key).or_default() += if a < b { 1 } else { -1 };
        }
    }
    let seam_open = edges
        .iter()
        .filter(|(_, &bal)| bal != 0)
        .filter(|((a, b), _)| {
            let (pa, pb) = (mesh.positions[*a as usize], mesh.positions[*b as usize]);
            (pa.z - pb.z).abs() > 1e-9 && pa.y.abs() < 1e-6 && pb.y.abs() < 1e-6 && pa.x > 0.0
        })
        .count();
    // The seam is now closed outright. Its two uses resolve to the same
    // vertices AND the reversed use starts from the correct endpoint, so
    // every edge along the join pairs.
    assert_eq!(
        seam_open, 0,
        "periodic seam must be closed, found {seam_open} open edges at u = 0"
    );
    assert!(
        mesh.positions.len() > 100,
        "the tube must actually be sampled, got {} vertices",
        mesh.positions.len()
    );

    // Closure must not be bought by deleting the seam. Classify EVERY open
    // edge: a capless tube is allowed exactly two rims of open boundary,
    // and nothing else. An implementation that dropped the seam triangles
    // would satisfy `seam_open == 0` above while failing here.
    let mut rim = 0usize;
    let mut other = 0usize;
    for ((a, b), _) in edges.iter().filter(|(_, &bal)| bal != 0) {
        let pa = mesh.positions[*a as usize];
        let pb = mesh.positions[*b as usize];
        if (pa.z - pb.z).abs() < 1e-9 {
            rim += 1;
        } else {
            other += 1;
        }
    }
    assert_eq!(
        other, 0,
        "every open edge must lie on a rim; {other} span the height"
    );
    assert_eq!(
        rim, 256,
        "both rims must remain fully sampled, got {rim} open rim edges"
    );

    // The lateral area of a closed cylinder is TAU * r * h. Chord sampling
    // inscribes the tube, so the mesh area is slightly under; assert it is
    // close, which fails outright if triangles went missing at the seam.
    let area: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            (b - a).cross(c - a).length() * 0.5
        })
        .sum();
    let exact = core::f64::consts::TAU * radius * height;
    assert!(
        (exact - area) / exact < 1e-3 && area <= exact,
        "lateral area {area} should inscribe {exact}"
    );
}

/// A half cylinder must NOT be welded shut.
///
/// The grid mesher wraps a patch whose first and last u columns are the
/// same mesh vertices. A sector's columns are DIFFERENT vertices that
/// merely sit some distance apart, so wrapping it would fabricate a
/// surface across the opening and silently turn an open sheet into a
/// closed tube. This is the failure mode that a coordinate-comparison
/// wrap test would eventually hit; the topological test must refuse it.
#[test]
fn a_half_cylinder_is_not_welded_shut() {
    use axiolid_core::{Point2, Vec3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Cylinder, Surface};
    use core::f64::consts::PI;

    let radius = 2.0;
    let height = 3.0;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(Surface::Cylinder(Cylinder {
            frame: axiolid_core::Frame3 {
                origin: axiolid_core::Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius,
        })))
        .expect("surface");
    let bottom = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(0.0, 0.0), Point2::new(PI, 0.0)],
            closed: false,
        })))
        .expect("bottom");
    let right = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(PI, 0.0), Point2::new(PI, height)],
            closed: false,
        })))
        .expect("right");
    let top = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(PI, height), Point2::new(0.0, height)],
            closed: false,
        })))
        .expect("top");
    let left = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(0.0, height), Point2::new(0.0, 0.0)],
            closed: false,
        })))
        .expect("left");

    let p = |u: f64, v: f64| axiolid_core::Point3::new(radius * u.cos(), radius * u.sin(), v);
    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let v00 = brep.add_vertex(Vertex {
        position: p(0.0, 0.0),
    });
    let v10 = brep.add_vertex(Vertex {
        position: p(PI, 0.0),
    });
    let v11 = brep.add_vertex(Vertex {
        position: p(PI, height),
    });
    let v01 = brep.add_vertex(Vertex {
        position: p(0.0, height),
    });
    let e_b = brep.add_edge(Edge {
        start: v00,
        end: v10,
        curve: None,
    });
    let e_r = brep.add_edge(Edge {
        start: v10,
        end: v11,
        curve: None,
    });
    let e_t = brep.add_edge(Edge {
        start: v11,
        end: v01,
        curve: None,
    });
    let e_l = brep.add_edge(Edge {
        start: v01,
        end: v00,
        curve: None,
    });
    let wire = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: e_b,
                orientation: Orientation::Forward,
                pcurve: Some(bottom),
            },
            EdgeUse {
                edge: e_r,
                orientation: Orientation::Forward,
                pcurve: Some(right),
            },
            EdgeUse {
                edge: e_t,
                orientation: Orientation::Forward,
                pcurve: Some(top),
            },
            EdgeUse {
                edge: e_l,
                orientation: Orientation::Forward,
                pcurve: Some(left),
            },
        ],
    });
    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    let root = builder.push(GeometryNode::BRep(brep)).expect("brep");
    let graph = builder.finish(vec![root]).expect("finish");
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(
            &graph,
            root,
            &ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
        )
        .expect("half cylinder tessellates");

    // Half the lateral area of a full cylinder, inscribed.
    let area: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            (b - a).cross(c - a).length() * 0.5
        })
        .sum();
    let exact = PI * radius * height;
    assert!(
        (exact - area) / exact < 1e-3 && area <= exact,
        "half-cylinder area {area} should inscribe {exact}"
    );

    // The opening must survive: a sheet has boundary on all four sides.
    let mut balance: std::collections::HashMap<(u32, u32), i32> = Default::default();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *balance.entry(key).or_default() += if a < b { 1 } else { -1 };
        }
    }
    let vertical: usize = balance
        .iter()
        .filter(|(_, &b)| b != 0)
        .filter(|((a, b), _)| {
            let (pa, pb) = (mesh.positions[*a as usize], mesh.positions[*b as usize]);
            (pa.z - pb.z).abs() > 1e-9
        })
        .count();
    assert!(
        vertical > 0,
        "the sector's straight sides must remain open, found none"
    );
}

/// A periodic patch WITH interior rows must close across its seam.
///
/// The tube fixture cannot prove the wrap: its seam needs no interior
/// samples, so every grid vertex is also a boundary vertex and the
/// boundary walk supplies the seam column whether the grid wraps or not.
/// A torus is curved in BOTH parameters, so the grid has interior rows
/// whose seam vertices exist only if the grid wraps. Without the wrap
/// those rows get two separate columns at u = 0 and u = TAU and the
/// surface cracks along its whole length.
#[test]
fn a_periodic_patch_with_interior_rows_closes_across_the_seam() {
    use axiolid_core::{Point2, Vec3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Surface, Torus};
    use core::f64::consts::TAU;

    let major = 5.0;
    let minor = 1.0;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(Surface::Torus(Torus {
            frame: axiolid_core::Frame3 {
                origin: axiolid_core::Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            major_radius: major,
            minor_radius: minor,
        })))
        .expect("surface");

    // A band around the tube: u full turn, v a quarter of the minor circle
    // so the patch needs interior rows in v.
    let v_hi = TAU / 4.0;
    let bottom = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(0.0, 0.0), Point2::new(TAU, 0.0)],
            closed: false,
        })))
        .expect("bottom");
    let top = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(TAU, v_hi), Point2::new(0.0, v_hi)],
            closed: false,
        })))
        .expect("top");
    let seam_up = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(TAU, 0.0), Point2::new(TAU, v_hi)],
            closed: false,
        })))
        .expect("seam up");
    let seam_down = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(0.0, v_hi), Point2::new(0.0, 0.0)],
            closed: false,
        })))
        .expect("seam down");

    let p = |u: f64, v: f64| {
        let r = major + minor * v.cos();
        axiolid_core::Point3::new(r * u.cos(), r * u.sin(), minor * v.sin())
    };
    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let v00 = brep.add_vertex(Vertex {
        position: p(0.0, 0.0),
    });
    let v01 = brep.add_vertex(Vertex {
        position: p(0.0, v_hi),
    });
    let e_bottom = brep.add_edge(Edge {
        start: v00,
        end: v00,
        curve: None,
    });
    let e_seam = brep.add_edge(Edge {
        start: v00,
        end: v01,
        curve: None,
    });
    let e_top = brep.add_edge(Edge {
        start: v01,
        end: v01,
        curve: None,
    });
    let wire = brep.add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: e_bottom,
                orientation: Orientation::Forward,
                pcurve: Some(bottom),
            },
            EdgeUse {
                edge: e_seam,
                orientation: Orientation::Forward,
                pcurve: Some(seam_up),
            },
            EdgeUse {
                edge: e_top,
                orientation: Orientation::Forward,
                pcurve: Some(top),
            },
            EdgeUse {
                edge: e_seam,
                orientation: Orientation::Reversed,
                pcurve: Some(seam_down),
            },
        ],
    });
    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).expect("brep");
    let graph = builder.finish(vec![root]).expect("finish");
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(
            &graph,
            root,
            &ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
        )
        .expect("torus band tessellates");

    // No vertex may be duplicated: a seam that failed to wrap emits a
    // second copy of every interior row's seam vertex.
    let mut seen: std::collections::HashMap<(i64, i64, i64), usize> = Default::default();
    for q in &mesh.positions {
        let key = (
            (q.x * 1e9).round() as i64,
            (q.y * 1e9).round() as i64,
            (q.z * 1e9).round() as i64,
        );
        *seen.entry(key).or_default() += 1;
    }
    let duplicated = seen.values().filter(|&&c| c > 1).count();
    assert_eq!(
        duplicated, 0,
        "the seam must be one column of vertices, found {duplicated} duplicated positions"
    );

    // And the only open edges may be the two rims, never a vertical crack.
    let mut balance: std::collections::HashMap<(u32, u32), i32> = Default::default();
    for t in mesh.indices.chunks_exact(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            *balance.entry(key).or_default() += if a < b { 1 } else { -1 };
        }
    }
    let cracks = balance
        .iter()
        .filter(|(_, &b)| b != 0)
        .filter(|((a, b), _)| {
            let (pa, pb) = (mesh.positions[*a as usize], mesh.positions[*b as usize]);
            // A rim edge runs around the major circle at constant minor
            // angle, so its endpoints share |z|. Anything else is a crack.
            (pa.z - pb.z).abs() > 1e-9
        })
        .count();
    assert_eq!(
        cracks, 0,
        "the seam must be closed, found {cracks} open edges spanning v"
    );
}

/// A curved face with a hole must fall back to the polygon triangulator.
///
/// The grid mesher only handles rectangles in parameter space. A face with
/// an inner loop is not one, and forcing a grid onto it would pave straight
/// over the hole. `recognise_grid` has to decline, which is a capability
/// boundary worth pinning: the guard is invisible until something relies
/// on it.
#[test]
fn a_curved_face_with_a_hole_is_not_gridded() {
    use axiolid_core::{Point2, Vec3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Cylinder, Surface};
    use core::f64::consts::PI;

    let radius = 2.0;
    let height = 4.0;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(Surface::Cylinder(Cylinder {
            frame: axiolid_core::Frame3 {
                origin: axiolid_core::Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius,
        })))
        .expect("surface");

    let p = |u: f64, v: f64| axiolid_core::Point3::new(radius * u.cos(), radius * u.sin(), v);
    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();

    // Outer: a rectangular patch u in [0, PI], v in [0, height].
    let outer_pts = [(0.0, 0.0), (PI, 0.0), (PI, height), (0.0, height)];
    let mut outer_uses = Vec::new();
    let mut outer_vertices = Vec::new();
    for (u, v) in outer_pts {
        outer_vertices.push(brep.add_vertex(Vertex { position: p(u, v) }));
    }
    for i in 0..4 {
        let (u0, v0) = outer_pts[i];
        let (u1, v1) = outer_pts[(i + 1) % 4];
        let pc = builder
            .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
                points: vec![Point2::new(u0, v0), Point2::new(u1, v1)],
                closed: false,
            })))
            .expect("outer pcurve");
        let edge = brep.add_edge(Edge {
            start: outer_vertices[i],
            end: outer_vertices[(i + 1) % 4],
            curve: None,
        });
        outer_uses.push(EdgeUse {
            edge,
            orientation: Orientation::Forward,
            pcurve: Some(pc),
        });
    }
    let outer_loop = brep.add_loop(Loop { edges: outer_uses });

    // Inner: a small window, wound the other way.
    let inner_pts = [(1.0, 1.0), (1.0, 3.0), (2.0, 3.0), (2.0, 1.0)];
    let mut inner_uses = Vec::new();
    let mut inner_vertices = Vec::new();
    for (u, v) in inner_pts {
        inner_vertices.push(brep.add_vertex(Vertex { position: p(u, v) }));
    }
    for i in 0..4 {
        let (u0, v0) = inner_pts[i];
        let (u1, v1) = inner_pts[(i + 1) % 4];
        let pc = builder
            .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
                points: vec![Point2::new(u0, v0), Point2::new(u1, v1)],
                closed: false,
            })))
            .expect("inner pcurve");
        let edge = brep.add_edge(Edge {
            start: inner_vertices[i],
            end: inner_vertices[(i + 1) % 4],
            curve: None,
        });
        inner_uses.push(EdgeUse {
            edge,
            orientation: Orientation::Forward,
            pcurve: Some(pc),
        });
    }
    let inner_loop = brep.add_loop(Loop { edges: inner_uses });

    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![
            FaceBound {
                loop_id: outer_loop,
                orientation: Orientation::Forward,
                outer: true,
            },
            FaceBound {
                loop_id: inner_loop,
                orientation: Orientation::Forward,
                outer: false,
            },
        ],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).expect("brep");
    let graph = builder.finish(vec![root]).expect("finish");
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(
            &graph,
            root,
            &ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
        )
        .expect("holed face tessellates");

    // A gridded result would tile the whole rectangle, so its area would
    // reach the full outer patch and its vertex count would be the grid's
    // (nu+1)*(nv+1). Both must be strictly below that: the face went to the
    // polygon triangulator, which is what the hole guard exists to force.
    //
    // The hole must actually be subtracted, not merely dented. Refinement
    // splits the long UV chords earcut produced, so the mapped area now
    // converges on the analytic patch minus the window. A boundary-only
    // triangulation scored 24.18 here -- it kept three quarters of a hole
    // it claimed to have cut.
    let area: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            (b - a).cross(c - a).length() * 0.5
        })
        .sum();
    let outer_area = PI * radius * height;
    // u spans 1 radian and v spans 2, so the window is r * 1 * 2 = 4.
    let hole_area = radius * 1.0 * 2.0;
    let expected = outer_area - hole_area;
    let removed = outer_area - area;
    assert!(
        (removed - hole_area).abs() < 1.0e-2,
        "the window must be subtracted: removed {removed} of {hole_area} \
         (area {area}, expected about {expected})"
    );
    assert!(
        area < expected + 1.0e-6,
        "a chordal mesh cannot exceed the analytic patch: {area} vs {expected}"
    );
    // The hole's own boundary must appear in the mesh: a grid would never
    // place vertices on the window's rim.
    let on_hole_rim = mesh
        .positions
        .iter()
        .filter(|p| {
            let u = p.y.atan2(p.x);
            (1.0..=2.0).contains(&u) && (p.z - 1.0).abs() < 1e-9
        })
        .count();
    assert!(
        on_hole_rim > 0,
        "the hole's rim vertices must survive into the mesh"
    );
}

/// A slanted trim must not be gridded.
///
/// `recognise_grid` only claims a boundary whose every point lies on the
/// rectangle's border at uniform spacing. A diagonal trim leaves points in
/// the rectangle's INTERIOR, and gridding it would replace the slanted
/// edge with the bounding box, inventing surface the face does not have.
#[test]
fn a_slanted_trim_is_not_gridded() {
    use axiolid_core::{Point2, Vec3};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Cylinder, Surface};

    let radius = 2.0;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(Surface::Cylinder(Cylinder {
            frame: axiolid_core::Frame3 {
                origin: axiolid_core::Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius,
        })))
        .expect("surface");

    // A triangle in parameter space: (0,0) -> (1,0) -> (1,2) -> back.
    // The closing edge is diagonal, so its samples sit strictly inside the
    // bounding rectangle.
    let corners = [(0.0_f64, 0.0_f64), (1.0, 0.0), (1.0, 2.0)];
    let p = |u: f64, v: f64| axiolid_core::Point3::new(radius * u.cos(), radius * u.sin(), v);
    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let vertices: Vec<_> = corners
        .iter()
        .map(|&(u, v)| brep.add_vertex(Vertex { position: p(u, v) }))
        .collect();
    let mut uses = Vec::new();
    for i in 0..3 {
        let (u0, v0) = corners[i];
        let (u1, v1) = corners[(i + 1) % 3];
        let pc = builder
            .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
                points: vec![Point2::new(u0, v0), Point2::new(u1, v1)],
                closed: false,
            })))
            .expect("pcurve");
        let edge = brep.add_edge(Edge {
            start: vertices[i],
            end: vertices[(i + 1) % 3],
            curve: None,
        });
        uses.push(EdgeUse {
            edge,
            orientation: Orientation::Forward,
            pcurve: Some(pc),
        });
    }
    let wire = brep.add_loop(Loop { edges: uses });
    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Reversed,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).expect("brep");
    let graph = builder.finish(vec![root]).expect("finish");
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(
            &graph,
            root,
            &ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
        )
        .expect("slanted face tessellates");

    let area: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            (b - a).cross(c - a).length() * 0.5
        })
        .sum();
    // The triangle is half its bounding rectangle. A gridded result would
    // pave the whole rectangle and roughly double the area.
    let rectangle = radius * 1.0 * 2.0;
    assert!(
        area < rectangle * 0.75,
        "a slanted trim must not be squared off: area {area} vs rectangle {rectangle}"
    );
    assert!(area > 0.0, "the face must still be meshed");
    for triangle in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[triangle[0] as usize];
        let b = mesh.positions[triangle[1] as usize];
        let c = mesh.positions[triangle[2] as usize];
        let center = (a + b + c) / 3.0;
        let radial = Vec3::new(center.x, center.y, 0.0);
        assert!(
            (b - a).cross(c - a).dot(radial) < 0.0,
            "a reversed curved bound must reverse Earcut winding"
        );
    }
}

#[test]
fn a_trimmed_rational_bspline_face_refines_its_support_surface() {
    use axiolid_core::{Point2, Point3, Tolerance};
    use axiolid_curve::{Curve2, KnotSpec, Polyline2};
    use axiolid_surface::{BSplineSurface, Surface};

    let radius = 2.0;
    let height = 1.5;
    let rational_weight = core::f64::consts::FRAC_1_SQRT_2;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(Surface::BSpline(BSplineSurface {
            u_degree: 2,
            v_degree: 1,
            control_points: vec![
                vec![
                    Point3::new(radius, 0.0, 0.0),
                    Point3::new(radius, 0.0, height),
                ],
                vec![
                    Point3::new(radius, radius, 0.0),
                    Point3::new(radius, radius, height),
                ],
                vec![
                    Point3::new(0.0, radius, 0.0),
                    Point3::new(0.0, radius, height),
                ],
            ],
            u_knots: vec![0.0, 1.0],
            u_multiplicities: vec![3, 3],
            v_knots: vec![0.0, height],
            v_multiplicities: vec![2, 2],
            weights: Some(vec![
                vec![1.0, 1.0],
                vec![rational_weight, rational_weight],
                vec![1.0, 1.0],
            ]),
            u_closed: false,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: None,
        })))
        .expect("NURBS surface");

    let mut line = |a: Point2, b: Point2| {
        builder
            .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
                points: vec![a, b],
                closed: false,
            })))
            .expect("pcurve")
    };
    let pcurves = [
        line(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)),
        line(Point2::new(1.0, 0.0), Point2::new(1.0, height)),
        line(Point2::new(1.0, height), Point2::new(0.0, height)),
        line(Point2::new(0.0, height), Point2::new(0.0, 0.0)),
    ];

    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let vertices = [
        brep.add_vertex(Vertex {
            position: Point3::new(radius, 0.0, 0.0),
        }),
        brep.add_vertex(Vertex {
            position: Point3::new(0.0, radius, 0.0),
        }),
        brep.add_vertex(Vertex {
            position: Point3::new(0.0, radius, height),
        }),
        brep.add_vertex(Vertex {
            position: Point3::new(radius, 0.0, height),
        }),
    ];
    let edges = [
        brep.add_edge(Edge {
            start: vertices[0],
            end: vertices[1],
            curve: None,
        }),
        brep.add_edge(Edge {
            start: vertices[1],
            end: vertices[2],
            curve: None,
        }),
        brep.add_edge(Edge {
            start: vertices[2],
            end: vertices[3],
            curve: None,
        }),
        brep.add_edge(Edge {
            start: vertices[3],
            end: vertices[0],
            curve: None,
        }),
    ];
    let wire = brep.add_loop(Loop {
        edges: edges
            .into_iter()
            .zip(pcurves)
            .map(|(edge, pcurve)| EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(pcurve),
            })
            .collect(),
    });
    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Reversed,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).expect("B-rep");
    let graph = builder.finish(vec![root]).expect("graph");
    let tolerance = Tolerance::MILLIMETRE;
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, root, &ExecutionOptions::new(tolerance))
        .expect("trimmed rational B-spline face tessellates");

    assert!(!mesh.indices.is_empty());
    for point in &mesh.positions {
        let radial_error = (point.x.hypot(point.y) - radius).abs();
        assert!(radial_error <= 1e-9, "off-support NURBS vertex: {point:?}");
    }
    for triangle in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[triangle[0] as usize];
        let b = mesh.positions[triangle[1] as usize];
        let c = mesh.positions[triangle[2] as usize];
        let center = (a + b + c) / 3.0;
        let radial = Point3::new(center.x, center.y, 0.0);
        assert!(
            (b - a).cross(c - a).dot(radial) < 0.0,
            "a reversed curved bound must reverse structured-grid winding"
        );
    }
    let max_sagitta = mesh
        .indices
        .chunks_exact(3)
        .flat_map(|triangle| {
            [
                (triangle[0], triangle[1]),
                (triangle[1], triangle[2]),
                (triangle[2], triangle[0]),
            ]
        })
        .map(|(a, b)| {
            let midpoint = (mesh.positions[a as usize] + mesh.positions[b as usize]) * 0.5;
            radius - midpoint.x.hypot(midpoint.y)
        })
        .fold(0.0_f64, f64::max);
    assert!(
        max_sagitta <= tolerance.linear(),
        "NURBS interior sagitta {max_sagitta} exceeded {}",
        tolerance.linear()
    );
}

#[test]
fn a_curved_face_preserves_a_pcurve_hole_during_refinement() {
    use axiolid_core::{Point2, Point3, Tolerance};
    use axiolid_curve::{Curve2, Polyline2};
    use axiolid_surface::{Cylinder, Surface};

    fn add_ring(
        builder: &mut GeometryGraphBuilder,
        brep: &mut BRep<axiolid_model::NodeId>,
        radius: f64,
        points: [Point2; 4],
    ) -> axiolid_topology::LoopId {
        let vertices: Vec<_> = points
            .iter()
            .map(|point| {
                brep.add_vertex(Vertex {
                    position: Point3::new(radius * point.x.cos(), radius * point.x.sin(), point.y),
                })
            })
            .collect();
        let mut uses = Vec::new();
        for index in 0..4 {
            let next = (index + 1) % 4;
            let pcurve = builder
                .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
                    points: vec![points[index], points[next]],
                    closed: false,
                })))
                .unwrap();
            let edge = brep.add_edge(Edge {
                start: vertices[index],
                end: vertices[next],
                curve: None,
            });
            uses.push(EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(pcurve),
            });
        }
        brep.add_loop(Loop { edges: uses })
    }

    let radius = 2.0;
    let mut builder = GeometryGraphBuilder::new();
    let surface = builder
        .push(GeometryNode::Surface(Surface::Cylinder(Cylinder {
            frame: axiolid_core::Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius,
        })))
        .unwrap();
    let mut brep: BRep<axiolid_model::NodeId> = BRep::default();
    let outer = add_ring(
        &mut builder,
        &mut brep,
        radius,
        [
            Point2::new(0.0, 0.0),
            Point2::new(core::f64::consts::FRAC_PI_2, 0.0),
            Point2::new(core::f64::consts::FRAC_PI_2, 3.0),
            Point2::new(0.0, 3.0),
        ],
    );
    // Clockwise parameter-space ring.
    let hole = add_ring(
        &mut builder,
        &mut brep,
        radius,
        [
            Point2::new(0.4, 1.0),
            Point2::new(0.4, 2.0),
            Point2::new(1.0, 2.0),
            Point2::new(1.0, 1.0),
        ],
    );
    let face = brep.add_face(Face {
        surface: Some(surface),
        bounds: vec![
            FaceBound {
                loop_id: hole,
                orientation: Orientation::Forward,
                outer: false,
            },
            FaceBound {
                loop_id: outer,
                orientation: Orientation::Forward,
                outer: true,
            },
        ],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });

    let root = builder.push(GeometryNode::BRep(brep)).unwrap();
    let graph = builder.finish(vec![root]).unwrap();
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, root, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .expect("curved face with a hole tessellates");

    let parameter = |index: u32| {
        let point = mesh.positions[index as usize];
        Point2::new(point.y.atan2(point.x), point.z)
    };
    let mut parameter_area = 0.0;
    for triangle in mesh.indices.chunks_exact(3) {
        let [a, b, c] = [
            parameter(triangle[0]),
            parameter(triangle[1]),
            parameter(triangle[2]),
        ];
        let centroid = (a + b + c) / 3.0;
        assert!(
            !(centroid.x > 0.4 && centroid.x < 1.0 && centroid.y > 1.0 && centroid.y < 2.0),
            "triangle centroid entered the pcurve hole: {centroid:?}"
        );
        parameter_area += ((b - a).x * (c - a).y - (b - a).y * (c - a).x).abs() * 0.5;
    }
    let expected = core::f64::consts::FRAC_PI_2 * 3.0 - 0.6;
    assert!(
        (parameter_area - expected).abs() <= 1.0e-10,
        "parameter area {parameter_area} did not preserve hole area {expected}"
    );
}
