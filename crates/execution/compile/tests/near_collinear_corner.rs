//! A cap corner a fraction of a micrometre off its neighbours' line keeps
//! its triangles (#278).
//!
//! The pentagon below is the cap of a real building part (a prism over a
//! 51-gon) shrunk with the failure kept: `v` lies about 0.4 um off the
//! line `a-b`, and `a` itself about 1 um off the line `p1-b`. Both are
//! within the noise band of the planar faces' sliver split, which before
//! #278 dropped the clipper's slivers `(a, v, b)` and `(b, p1, a)` from the
//! top cap and then split `p1-b` at `a` only: `v` was in no triangle of the
//! cap while both side quads kept it, and the closed prism meshed open.
//!
//! Every oracle is closed form: a closed mesh of a prism with planar caps
//! has the cap's shoelace area times the height as its volume, whichever
//! way the caps are triangulated, so the mesh volume is compared with it
//! to rounding, measured about the prism's own origin.

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

/// The issue's counter-clockwise pentagon `p0, p1, a, v, b`, in metres.
const PENTAGON: [[f64; 2]; 5] = [
    [0.877258, 2.548603],
    [-0.586191, 1.488751],
    [0.0, 0.0],
    [0.054955, -0.13957],
    [0.183174, -0.465204],
];

const HEIGHT: f64 = 0.03;

/// The pentagon with `v` rebuilt at `t = 0.3` along `a-b`, `offset` metres
/// off the line (positive to the left of `a -> b`, which makes `v` reflex).
fn pentagon_with(offset: f64) -> [[f64; 2]; 5] {
    let mut cap = PENTAGON;
    let (a, b) = (cap[2], cap[4]);
    let d = [b[0] - a[0], b[1] - a[1]];
    let length = d[0].hypot(d[1]);
    let left = [-d[1] / length, d[0] / length];
    cap[3] = [
        a[0] + 0.3 * d[0] + offset * left[0],
        a[1] + 0.3 * d[1] + offset * left[1],
    ];
    cap
}

/// The cap mirrored in `x`, its ring reversed so it stays counter-clockwise.
fn mirrored(cap: [[f64; 2]; 5]) -> [[f64; 2]; 5] {
    let mut out = cap.map(|p| [-p[0], p[1]]);
    out.reverse();
    out
}

/// A prism of `HEIGHT` over a counter-clockwise `cap`, every corner moved
/// by `shift`: corners `0..n` at the bottom, `n..2n` at the top, faces
/// wound outward (the bottom cap reversed, then the side quads).
fn prism(cap: &[[f64; 2]], shift: [f64; 3]) -> (Vec<[f64; 3]>, Vec<Vec<usize>>) {
    let n = cap.len();
    let mut corners: Vec<[f64; 3]> = Vec::with_capacity(2 * n);
    for z in [0.0, HEIGHT] {
        corners.extend(
            cap.iter()
                .map(|p| [p[0] + shift[0], p[1] + shift[1], z + shift[2]]),
        );
    }
    let mut faces = vec![(0..n).rev().collect::<Vec<_>>(), (n..2 * n).collect()];
    for i in 0..n {
        let j = (i + 1) % n;
        faces.push(vec![i, j, n + j, n + i]);
    }
    (corners, faces)
}

/// The same faces as authored polygons over shared positions.
fn polygons(corners: &[[f64; 3]], faces: &[Vec<usize>]) -> PolygonMesh {
    PolygonMesh {
        positions: corners
            .iter()
            .map(|c| Point3::new(c[0], c[1], c[2]))
            .collect(),
        faces: faces
            .iter()
            .map(|ring| PolygonFace {
                outer: ring.iter().map(|&i| i as u32).collect(),
                holes: Vec::new(),
            })
            .collect(),
    }
}

/// The same faces as a faceted B-rep solid, every edge shared between the
/// two faces that use it, the way lowering builds an `IfcFacetedBrep`.
fn brep(corners: &[[f64; 3]], faces: &[Vec<usize>]) -> BRep<NodeId> {
    let mut brep = BRep::default();
    let vertices: Vec<_> = corners
        .iter()
        .map(|c| {
            brep.add_vertex(Vertex {
                position: Vec3::new(c[0], c[1], c[2]),
            })
        })
        .collect();
    let mut edges: HashMap<(usize, usize), EdgeId> = HashMap::new();
    let mut ids = Vec::new();
    for ring in faces {
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
        let bound = FaceBound {
            loop_id: brep.add_loop(Loop { edges: uses }),
            orientation: Orientation::Forward,
            outer: true,
        };
        ids.push(brep.add_face(Face {
            surface: None,
            bounds: vec![bound],
            orientation: Orientation::Forward,
        }));
    }
    let shell = brep.add_shell(Shell {
        faces: ids.iter().map(|&f| (f, Orientation::Forward)).collect(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    brep
}

fn compile(node: GeometryNode) -> CompileOutcome {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(node).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_reported(&graph, root, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .expect("the prism tessellates")
}

/// The prism's corners as stored, less `shift`: exact for every shift used
/// here (Sterbenz), so the oracle sees the stored corners.
fn local(corners: &[[f64; 3]], shift: [f64; 3]) -> Vec<[f64; 3]> {
    corners
        .iter()
        .map(|c| [c[0] - shift[0], c[1] - shift[1], c[2] - shift[2]])
        .collect()
}

/// Shoelace area of the stored bottom cap times the stored height.
fn closed_form_volume(corners: &[[f64; 3]], n: usize) -> f64 {
    let doubled: f64 = (0..n)
        .map(|i| {
            let (p, q) = (corners[i], corners[(i + 1) % n]);
            p[0] * q[1] - q[0] * p[1]
        })
        .sum();
    0.5 * doubled * (corners[n][2] - corners[0][2])
}

/// The mesh's signed volume about `shift`.
fn volume(mesh: &TriMesh, shift: [f64; 3]) -> f64 {
    let at = |i: u32| {
        let p = mesh.positions[i as usize];
        Vec3::new(p.x - shift[0], p.y - shift[1], p.z - shift[2])
    };
    mesh.indices
        .chunks_exact(3)
        .map(|t| at(t[0]).dot(at(t[1]).cross(at(t[2]))))
        .sum::<f64>()
        / 6.0
}

/// A closed, consistently wound two-manifold reported `Solid`, every
/// authored corner a mesh corner, with the closed-form volume.
fn assert_closed_prism(outcome: &CompileOutcome, cap: &[[f64; 2]], shift: [f64; 3], what: &str) {
    let (corners, _) = prism(cap, shift);
    let mesh = &outcome.mesh;
    let adjacency = EdgeAdjacency::build(mesh);
    let boundary = adjacency.boundary_edges().count();
    assert!(
        adjacency.edge_count() > 0 && adjacency.is_closed_two_manifold(),
        "{what}: {boundary} boundary edges, {} triangles: {:?}",
        mesh.triangle_count(),
        mesh.indices
    );
    assert_eq!(outcome.closure, MeshClosure::Solid, "{what}");
    // Each cap of n corners is n - 2 triangles, each side quad two.
    let n = cap.len();
    assert_eq!(mesh.triangle_count(), 2 * (n - 2) + 2 * n, "{what}");
    let mut used = vec![false; mesh.positions.len()];
    for &i in &mesh.indices {
        used[i as usize] = true;
    }
    assert!(
        used.iter().all(|&u| u),
        "{what}: a corner is in no triangle"
    );
    let want = closed_form_volume(&local(&corners, shift), n);
    let got = volume(mesh, shift);
    assert!(
        (got - want).abs() <= 1e-12 * want,
        "{what}: volume {got}, want {want}"
    );
}

/// The prism compiled as a faceted B-rep and as authored polygons.
fn both_nodes(cap: &[[f64; 2]], shift: [f64; 3]) -> [(CompileOutcome, &'static str); 2] {
    let (corners, faces) = prism(cap, shift);
    [
        (
            compile(GeometryNode::BRep(brep(&corners, &faces))),
            "faceted B-rep",
        ),
        (
            compile(GeometryNode::PolygonMesh(polygons(&corners, &faces))),
            "authored polygons",
        ),
    ]
}

#[test]
fn the_issues_prism_meshes_closed_as_a_brep_and_as_polygons() {
    for (outcome, what) in both_nodes(&PENTAGON, [0.0; 3]) {
        assert_closed_prism(&outcome, &PENTAGON, [0.0; 3], what);
        let solid = outcome.solid_mesh().expect("a closed prism is a solid");
        // The oracle that consumers use agrees.
        let measured = axiolid_measure::mesh::volume_properties(solid, Tolerance::MILLIMETRE)
            .expect("a closed two-manifold has a volume")
            .signed_volume;
        let want = closed_form_volume(&prism(&PENTAGON, [0.0; 3]).0, PENTAGON.len());
        assert!(
            (measured - want).abs() <= 1e-12 * want,
            "{what}: {measured}"
        );
    }
}

#[test]
fn the_issues_prism_with_mirrored_caps_meshes_closed() {
    let cap = mirrored(PENTAGON);
    for (outcome, what) in both_nodes(&cap, [0.0; 3]) {
        assert_closed_prism(&outcome, &cap, [0.0; 3], what);
    }
}

/// `v` from a nanometre to ten micrometres off `a-b`, on either side,
/// mirrored or not, at the origin and far from it: every corner stays a
/// corner of both caps and the prism closes.
#[test]
fn every_off_line_distance_and_placement_meshes_closed() {
    let offsets = [1e-9, 1e-8, 1e-7, 2e-7, 4e-7, 7e-7, 1e-6, 2e-6, 5e-6, 1e-5];
    let shifts = [
        [0.0; 3],
        [1.0e5, -1.0e5, 1.0e5],
        [5.0e6, 5.0e6, -5.0e6],
        [-5.0e6, 2.5e6, 125.0],
    ];
    for offset in offsets.into_iter().flat_map(|d| [d, -d]) {
        for mirror in [false, true] {
            let cap = pentagon_with(offset);
            let cap = if mirror { mirrored(cap) } else { cap };
            for shift in shifts {
                for (outcome, what) in both_nodes(&cap, shift) {
                    let what =
                        format!("{what}, offset {offset:e}, mirrored {mirror}, at {shift:?}");
                    assert_closed_prism(&outcome, &cap, shift, &what);
                }
            }
        }
    }
}
