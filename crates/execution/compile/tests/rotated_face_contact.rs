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
use axiolid_mesh_compile::deviation::{REFEREED_UNION, REFEREED_UNION_EXACT};
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
    column_and_corbels(transform, &[face], composed)
}

/// The column and a corbel on each of `faces`, united one after another:
/// `Union(Union(column, first), second)` and so on.
fn column_and_corbels(
    transform: Transform3,
    faces: &[usize],
    composed: bool,
) -> (GeometryGraph, NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let mut graph = |node| builder.push(node).expect("a valid node");
    let profile = graph(rectangle(0.4, 0.5));
    let column = graph(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 6.0,
    }));
    let mut root = graph(GeometryNode::Instance(Instance {
        source: column,
        transform,
    }));
    for &face in faces {
        let (centre, size) = FACES[face];
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
        root = graph(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: root,
            right: corbel,
            operator: BooleanOperator::Union,
        }));
    }
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
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [t[0], t[1], t[2]].map(|i| mesh.positions[i as usize] - origin);
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// Compiled and checked: one `Solid` body of `expected` volume.
fn assert_one_solid(
    graph: &GeometryGraph,
    root: NodeId,
    expected: Scalar,
    slack: Scalar,
    what: &str,
) {
    let (mesh, closure) = compile(graph, root).unwrap_or_else(|e| panic!("{what}: {e}"));
    assert_eq!(component_count(&mesh), 1, "{what}: two bodies in contact");
    assert_eq!(closure, MeshClosure::Solid, "{what}");
    let measured = volume(&mesh);
    assert!(
        (measured - expected).abs() <= slack * expected,
        "{what}: volume {measured}"
    );
}

#[test]
fn the_issue_sweep_is_one_solid_on_every_face() {
    let translation = Vec3::new(1.0, 2.0, 0.0);
    for composed in [false, true] {
        for face in 0..FACES.len() {
            for k in 0..200 {
                let theta = 0.001 + k as Scalar * 0.031_415_9;
                let (graph, root) =
                    column_and_corbel(placement(theta, translation), face, composed);
                let what = format!("theta {theta}, face {face}, composed {composed}");
                assert_one_solid(&graph, root, VOLUME, 1e-9, &what);
            }
        }
    }
}

/// xorshift64*: a fixed stream, so the cases are the same on every run.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn unit(&mut self) -> Scalar {
        (self.next() >> 11) as Scalar / (1u64 << 53) as Scalar
    }
}

#[test]
fn the_result_stands_where_its_operands_are() {
    // Cut in the shared frame, placed by the subject's transform: the
    // column's corners are among the result's vertices, wherever the
    // corbel's own offset puts the tool's frame.
    for composed in [false, true] {
        let t = placement(0.5, Vec3::new(6.1e5, 5.6e6, 204.3));
        let (graph, root) = column_and_corbel(t, 0, composed);
        let (mesh, _) = compile(&graph, root).expect("one solid");
        for corner in [
            (-0.2, -0.25, 0.0),
            (0.2, 0.25, 0.0),
            (-0.2, 0.25, 6.0),
            (0.2, -0.25, 6.0),
        ] {
            let placed = t.transform_point3(Vec3::new(corner.0, corner.1, corner.2));
            assert!(
                mesh.positions
                    .iter()
                    .any(|p| (*p - placed).abs().max_element() < 1e-6),
                "composed {composed}: no vertex at {placed}"
            );
        }
    }
}

#[test]
fn random_turns_and_sites_are_one_solid() {
    let mut rng = Rng(356);
    for case in 0..200 {
        let theta = rng.unit() * std::f64::consts::TAU;
        let face = (rng.next() % FACES.len() as u64) as usize;
        let translation = if rng.next().is_multiple_of(2) {
            Vec3::new(1.0, 2.0, 0.0)
        } else {
            Vec3::new(6.1e5, 5.6e6, 204.3)
        };
        let composed = rng.next().is_multiple_of(2);
        let (graph, root) = column_and_corbel(placement(theta, translation), face, composed);
        let what = format!("case {case}: theta {theta}, face {face}, at {translation}");
        assert_one_solid(&graph, root, VOLUME, 1e-8, &what);
    }
}

#[test]
fn a_column_with_corbels_on_four_faces_is_one_solid() {
    // A chain of unions: each result stays in the column's frame for the
    // next (`Built::local`).
    for k in 0..20 {
        let theta = 0.001 + k as Scalar * 0.31;
        for composed in [false, true] {
            let (graph, root) = column_and_corbels(
                placement(theta, Vec3::new(1.0, 2.0, 0.0)),
                &[0, 1, 2, 3],
                composed,
            );
            let what = format!("theta {theta}, composed {composed}");
            assert_one_solid(&graph, root, 1.2 + 4.0 * 0.025, 1e-9, &what);
            // The chain is cut in the shared frame, not decided by the exact
            // referee (#391): a union that skipped the frame would still
            // fuse there, and this is what shows it did not.
            let options = ExecutionOptions::new(Tolerance::METRE)
                .with_chord_error(1e-3)
                .expect("a positive chord error");
            let (_, report) = compiler()
                .compile_mesh_with_deviation(&graph, root, &options)
                .unwrap_or_else(|e| panic!("{what}: {e}"));
            assert!(
                report
                    .contributions
                    .iter()
                    .all(|c| c.detail != REFEREED_UNION && c.detail != REFEREED_UNION_EXACT),
                "{what}: refereed"
            );
        }
    }
}

#[test]
fn operands_offset_in_their_shared_frame_fuse() {
    // The corbel under a transform with the column's linear part and a
    // translation offset by the rotated lift, computed in `f64`: the offset
    // between them in the column's frame rounds, and the snap closes it.
    for theta in [0.001, 0.5, 2.0, 5.0] {
        let t = placement(theta, Vec3::new(6.1e5, 5.6e6, 204.3));
        let face = 1;
        let (centre, size) = FACES[face];
        let mut builder = GeometryGraphBuilder::new();
        let mut graph = |node| builder.push(node).expect("a valid node");
        let profile = graph(rectangle(0.4, 0.5));
        let column = graph(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 6.0,
        }));
        let column = graph(GeometryNode::Instance(Instance {
            source: column,
            transform: t,
        }));
        let profile = graph(rectangle(size[0], size[1]));
        let corbel = graph(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 0.5,
        }));
        let lift = t.matrix3 * Vec3::new(centre[0], centre[1], 3.5);
        let corbel = graph(GeometryNode::Instance(Instance {
            source: corbel,
            transform: Transform3 {
                matrix3: t.matrix3,
                translation: t.translation + lift,
            },
        }));
        let root = graph(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: column,
            right: corbel,
            operator: BooleanOperator::Union,
        }));
        let graph = builder.finish(vec![root]).expect("a valid graph");
        assert_one_solid(&graph, root, VOLUME, 1e-8, &format!("theta {theta}"));
    }
}

#[test]
fn a_flush_opening_under_a_shared_turned_placement_cuts_through() {
    // A wall `4 x 0.25 x 3` less an opening exactly as deep as the wall,
    // both under one turned placement far from the origin: the opening's
    // faces are the wall's planes, which the world frame leaves a rounding
    // apart (ADR 0080's #276 limit).
    let t = placement(0.4, Vec3::new(6.0e5, 5.6e6, 0.0));
    let mut builder = GeometryGraphBuilder::new();
    let mut graph = |node| builder.push(node).expect("a valid node");
    let profile = graph(rectangle(4.0, 0.25));
    let wall = graph(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 3.0,
    }));
    let wall = graph(GeometryNode::Instance(Instance {
        source: wall,
        transform: t,
    }));
    let profile = graph(rectangle(1.0, 0.25));
    let opening = graph(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 2.0,
    }));
    let opening = graph(GeometryNode::Instance(Instance {
        source: opening,
        transform: t * Transform3::from_translation(Vec3::new(0.0, 0.0, 0.5)),
    }));
    let root = graph(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: wall,
        right: opening,
        operator: BooleanOperator::Difference,
    }));
    let graph = builder.finish(vec![root]).expect("a valid graph");
    assert_one_solid(&graph, root, 3.0 - 0.5, 1e-9, "flush opening");
}

#[test]
fn measured_unions_in_a_shared_frame_are_certified() {
    let options = ExecutionOptions::new(Tolerance::METRE)
        .with_chord_error(1e-3)
        .expect("a positive chord error");
    for theta in [0.001, 0.5, 2.0] {
        let (graph, root) = column_and_corbel(placement(theta, Vec3::new(1.0, 2.0, 0.0)), 0, false);
        let (outcome, report) = compiler()
            .compile_mesh_with_deviation(&graph, root, &options)
            .unwrap_or_else(|e| panic!("theta {theta}: {e}"));
        assert_eq!(component_count(&outcome.mesh), 1, "theta {theta}");
        let bound = report
            .bound
            .unwrap_or_else(|| panic!("theta {theta}: {report:?}"));
        assert!(bound <= 1e-3, "theta {theta}: {report:?}");
    }
}
