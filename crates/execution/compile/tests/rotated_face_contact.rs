//! A column and a corbel whose faces are in contact, under a general
//! rotation about Z (#356).
//!
//! The column is `0.4 x 0.5 x 6` (x, y, z; `z` from 0), centred on x and y.
//! A corbel `0.25` deep, `0.2` wide and `0.5` tall stands on each of its
//! four side faces at `z` in `[3.5, 4]`, flush with the face in contact.
//! Both are placed by the same rotation `T` and a translation, so in
//! exact arithmetic the corbel's contact face lies on the column's face
//! plane. Placed, each vertex rounds on its own, so the two faces are a
//! rounding residue apart and the mesh boolean, which has no tolerance,
//! returned two disjoint bodies, a pinch, or a body with cracks where the
//! faces cross.
//!
//! The policy (ADR 0096): the union is one body, or is refused by name;
//! never two bodies within the tolerance.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Scalar, Tolerance, Transform3, Vec3};
use axiolid_mesh::{component_count, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{MeshClosure, MeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Profile, RectangleProfile};

/// The column and one corbel.
const VOLUME: Scalar = 1.2 + 0.025;

/// The four faces: the corbel's centre in the column's frame (x, y) and its
/// profile size along x and y.
const FACES: [([Scalar; 2], [Scalar; 2]); 4] = [
    ([-0.325, 0.0], [0.25, 0.2]),
    ([0.325, 0.0], [0.25, 0.2]),
    ([0.0, -0.375], [0.2, 0.25]),
    ([0.0, 0.375], [0.2, 0.25]),
];

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn rectangle(x: Scalar, y: Scalar) -> GeometryNode {
    GeometryNode::Profile(Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
}

/// The placement: a turn of `theta` about Z, then `translation`.
fn placement(theta: Scalar, translation: Vec3) -> Transform3 {
    Transform3::from_translation(translation) * Transform3::from_rotation_z(theta)
}

/// A column and the corbel on `face`, united, both placed by `transform`.
/// With `composed` the corbel is one `Instance` of `transform` and its
/// offset, rather than the offset inside the placement.
fn column_and_corbel(
    transform: Transform3,
    face: usize,
    composed: bool,
) -> (GeometryGraph, NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let (centre, size) = FACES[face];
    let mut graph = |node| builder.push(node).expect("a valid node");
    let profile = graph(rectangle(0.4, 0.5));
    let column = graph(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 6.0,
    }));
    let column = graph(GeometryNode::Instance(Instance {
        source: column,
        transform,
    }));
    let profile = graph(rectangle(size[0], size[1]));
    let corbel = graph(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 0.5,
    }));
    let offset = Transform3::from_translation(Vec3::new(centre[0], centre[1], 3.5));
    let corbel = if composed {
        graph(GeometryNode::Instance(Instance {
            source: corbel,
            transform: transform * offset,
        }))
    } else {
        let lifted = graph(GeometryNode::Instance(Instance {
            source: corbel,
            transform: offset,
        }));
        graph(GeometryNode::Instance(Instance {
            source: lifted,
            transform,
        }))
    };
    let root = graph(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: column,
        right: corbel,
        operator: BooleanOperator::Union,
    }));
    (builder.finish(vec![root]).expect("a valid graph"), root)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> Result<(TriMesh, MeshClosure), GeomError> {
    let options = ExecutionOptions::new(Tolerance::METRE)
        .with_chord_error(1e-3)
        .expect("a positive chord error");
    let outcome = compiler().compile_mesh_reported(graph, root, &options)?;
    Ok((outcome.mesh, outcome.closure))
}

/// The divergence volume. The certified `enclosed_volume` refuses a body
/// with ulp-thin cracks as self-intersecting; this reads it through.
fn volume(mesh: &TriMesh) -> Scalar {
    let origin = mesh.positions[0];
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [t[0], t[1], t[2]].map(|i| mesh.positions[i as usize] - origin);
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

#[test]
fn the_issue_sweep_is_one_solid_or_refused_by_name() {
    let translation = Vec3::new(1.0, 2.0, 0.0);
    for composed in [false, true] {
        for face in 0..FACES.len() {
            for k in 0..200 {
                let theta = 0.001 + k as Scalar * 0.031_415_9;
                let (graph, root) =
                    column_and_corbel(placement(theta, translation), face, composed);
                let what = format!("theta {theta}, face {face}, composed {composed}");
                match compile(&graph, root) {
                    Ok((mesh, closure)) => {
                        assert_eq!(component_count(&mesh), 1, "{what}: two bodies in contact");
                        assert_eq!(closure, MeshClosure::Solid, "{what}");
                        let measured = volume(&mesh);
                        assert!(
                            (measured - VOLUME).abs() <= 1e-9 * VOLUME,
                            "{what}: volume {measured}"
                        );
                    }
                    Err(GeomError::Degenerate(_)) => {}
                    Err(other) => panic!("{what}: {other}"),
                }
            }
        }
    }
}
