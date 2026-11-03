//! A union of solids that meet within the linear tolerance (#357): two unit
//! boxes `g` apart under `Tolerance::METRE` (`1e-6`).
//!
//! Policy (ADR 0096): within the tolerance, the operands of a union are in
//! contact, as the exact compiler already reads them (ADR 0080). The mesh
//! compiler merges them into one body or refuses the union by name
//! ([`UNMERGED_UNION_CONTACT`]); it never returns two bodies within the
//! tolerance of each other. Beyond the tolerance a gap is authored: two
//! bodies, as before.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Mat3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_inspect::enclosed_volume;
use axiolid_mesh::{component_count, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::SNAPPED_OPERANDS;
use axiolid_mesh_compile::{DeviationBound, ReferenceMeshCompiler, UNMERGED_UNION_CONTACT};
use axiolid_mesh_compile_contract::{MeshClosure, MeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Profile, RectangleProfile};

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

/// A unit box `[-0.5, 0.5]^2 x [0, 1]` under `first`, and another under
/// `second` moved `offset` in its own frame, united.
fn boxes(first: Transform3, second: Transform3, offset: Vec3) -> (GeometryGraph, NodeId) {
    boxes_by(BooleanOperator::Union, first, second, offset)
}

fn boxes_by(
    operator: BooleanOperator,
    first: Transform3,
    second: Transform3,
    offset: Vec3,
) -> (GeometryGraph, NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let mut push = |node| builder.push(node).expect("a valid node");
    let profile = push(GeometryNode::Profile(Profile::Rectangle(
        RectangleProfile {
            x: 1.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        },
    )));
    let solid = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 1.0,
    }));
    let a = push(GeometryNode::Instance(Instance {
        source: solid,
        transform: first,
    }));
    let b = push(GeometryNode::Instance(Instance {
        source: solid,
        transform: second * Transform3::from_translation(offset),
    }));
    let root = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: a,
        right: b,
        operator,
    }));
    (builder.finish(vec![root]).expect("a valid graph"), root)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> Result<(TriMesh, MeshClosure), GeomError> {
    compile_at(graph, root, Tolerance::METRE)
}

fn compile_at(
    graph: &GeometryGraph,
    root: NodeId,
    tolerance: Tolerance,
) -> Result<(TriMesh, MeshClosure), GeomError> {
    let options = ExecutionOptions::new(tolerance)
        .with_chord_error(1e-3)
        .expect("a positive chord error");
    let outcome = compiler().compile_mesh_reported(graph, root, &options)?;
    Ok((outcome.mesh, outcome.closure))
}

fn volume(mesh: &TriMesh) -> Scalar {
    let interval = enclosed_volume(mesh).expect("a measurable volume");
    0.5 * (interval.lower + interval.upper)
}

fn is_unmerged(error: &GeomError) -> bool {
    matches!(error, GeomError::Degenerate(message) if message.starts_with(UNMERGED_UNION_CONTACT))
}

fn is_pinch(error: &GeomError) -> bool {
    matches!(error, GeomError::Degenerate(message) if message.contains("touches itself"))
}

#[test]
fn a_gap_beyond_the_tolerance_stays_two_bodies() {
    for gap in [2e-6, 1e-3] {
        let (graph, root) = boxes(
            Transform3::IDENTITY,
            Transform3::IDENTITY,
            Vec3::new(1.0 + gap, 0.0, 0.0),
        );
        let (mesh, closure) = compile(&graph, root).unwrap_or_else(|e| panic!("gap {gap:e}: {e}"));
        assert_eq!(component_count(&mesh), 2, "gap {gap:e}");
        assert_eq!(closure, MeshClosure::Solid, "gap {gap:e}");
        assert!((volume(&mesh) - 2.0).abs() < 1e-12, "gap {gap:e}");
    }
}

#[test]
fn a_gap_within_the_tolerance_fuses() {
    for gap in [0.0, 1e-15, 1e-9, 0.5e-6, 0.9e-6] {
        let (graph, root) = boxes(
            Transform3::IDENTITY,
            Transform3::IDENTITY,
            Vec3::new(1.0 + gap, 0.0, 0.0),
        );
        let (mesh, closure) = compile(&graph, root).unwrap_or_else(|e| panic!("gap {gap:e}: {e}"));
        assert_eq!(component_count(&mesh), 1, "gap {gap:e}");
        assert_eq!(closure, MeshClosure::Solid, "gap {gap:e}");
        // The tool's near face lands on the host's: the union gains the gap.
        assert!((volume(&mesh) - (2.0 + gap)).abs() < 1e-12, "gap {gap:e}");

        let options = ExecutionOptions::new(Tolerance::METRE)
            .with_chord_error(1e-3)
            .expect("a positive chord error");
        let (_, report) = compiler()
            .compile_mesh_with_deviation(&graph, root, &options)
            .unwrap_or_else(|e| panic!("gap {gap:e}: {e}"));
        let snapped: Vec<Scalar> = report
            .contributions
            .iter()
            .filter(|c| c.detail == SNAPPED_OPERANDS)
            .map(|c| match c.bound {
                DeviationBound::Certified(moved) => moved,
                other => panic!("a snap is certified: {other:?}"),
            })
            .collect();
        if gap > 1e-12 {
            // Only the tool's near face moves, by the gap itself.
            assert_eq!(snapped.len(), 1, "gap {gap:e}");
            assert!((snapped[0] - gap).abs() < 1e-15, "gap {gap:e}");
        } else {
            // Sub-rounding gaps vanish in the operands' coordinates.
            assert!(snapped.iter().all(|m| *m <= 1e-14), "gap {gap:e}");
        }
    }
}

#[test]
fn a_gap_within_a_millimetre_tolerance_fuses() {
    let (graph, root) = boxes(
        Transform3::IDENTITY,
        Transform3::IDENTITY,
        Vec3::new(1.0 + 0.5e-3, 0.0, 0.0),
    );
    let (mesh, closure) = compile_at(&graph, root, Tolerance::MILLIMETRE).expect("fused");
    assert_eq!(component_count(&mesh), 1);
    assert_eq!(closure, MeshClosure::Solid);
    assert!((volume(&mesh) - (2.0 + 0.5e-3)).abs() < 1e-12);
}

#[test]
fn a_difference_keeps_a_skin_within_the_tolerance() {
    // #291 stands for the other operators: a skin thinner than the
    // tolerance is authored, not closed. The tool stops 0.5e-6 short of the
    // host's far face and leaves exactly that skin.
    let skin = 0.5e-6;
    let (graph, root) = boxes_by(
        BooleanOperator::Difference,
        Transform3::IDENTITY,
        Transform3::IDENTITY,
        Vec3::new(skin, 0.0, 0.0),
    );
    let (mesh, closure) = compile(&graph, root).expect("the skin is kept");
    assert_eq!(component_count(&mesh), 1);
    assert_eq!(closure, MeshClosure::Solid);
    assert!((volume(&mesh) - skin).abs() < 1e-12, "{}", volume(&mesh));
}

#[test]
fn a_gap_within_the_tolerance_under_a_turn_fuses() {
    // The #357 boxes under one turned placement at georeferenced
    // coordinates: cut in their shared frame, where the faces are
    // axis-aligned, and placed.
    let turn = Transform3::from_translation(Vec3::new(6.1e5, 5.6e6, 204.3))
        * Transform3::from_rotation_z(0.7);
    let gap = 0.5e-6;
    let (graph, root) = boxes(turn, turn, Vec3::new(1.0 + gap, 0.0, 0.0));
    let (mesh, closure) = compile(&graph, root).expect("fused");
    assert_eq!(component_count(&mesh), 1);
    assert_eq!(closure, MeshClosure::Solid);
    assert!(
        (volume(&mesh) - (2.0 + gap)).abs() < 1e-8,
        "{}",
        volume(&mesh)
    );
}

#[test]
fn differently_turned_operands_a_gap_within_the_tolerance_apart_are_refused_by_name() {
    // Two placements of the same turn, composed another way: their linear
    // parts differ in the last bits, so the faces that meet are not exactly
    // coplanar once placed, and are not cut in a shared frame (the limit
    // ADR 0096 records; TODO(#391)).
    let site = Vec3::new(6.0e5, 5.6e6, 200.0);
    let a = Transform3::from_translation(site) * Transform3::from_rotation_z(0.3);
    let b = Transform3::from_translation(site)
        * Transform3::from_rotation_z(0.15)
        * Transform3::from_rotation_z(0.15);
    assert_ne!(a.matrix3, b.matrix3, "the placements must differ in bits");
    let (graph, root) = boxes(a, b, Vec3::new(1.0 + 0.5e-6, 0.0, 0.0));
    let error = compile(&graph, root).expect_err("the contact is not merged");
    assert!(is_unmerged(&error), "{error}");
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

    fn pick<T: Copy>(&mut self, items: &[T]) -> T {
        items[(self.next() % items.len() as u64) as usize]
    }
}

#[test]
fn a_union_never_returns_bodies_in_contact() {
    let mut rng = Rng(357);
    let faces = [
        Vec3::X,
        Vec3::NEG_X,
        Vec3::Y,
        Vec3::NEG_Y,
        Vec3::Z,
        Vec3::NEG_Z,
    ];
    let gaps = [0.0, 1e-15, 0.5e-6, 0.9e-6];
    let sites = [Vec3::ZERO, Vec3::new(6.0e5, 5.6e6, 200.0)];
    for case in 0..200 {
        let axis = Vec3::new(rng.unit() - 0.5, rng.unit() - 0.5, rng.unit() - 0.5).normalize();
        let angle = rng.unit() * std::f64::consts::TAU;
        let site = rng.pick(&sites);
        let a = Transform3::from_mat3_translation(Mat3::from_axis_angle(axis, angle), site);
        let half = Mat3::from_axis_angle(axis, 0.5 * angle);
        let b = if rng.next().is_multiple_of(2) {
            a
        } else {
            Transform3::from_mat3_translation(half * half, site)
        };
        let gap = rng.pick(&gaps);
        // The box is a unit in every direction, so a face is a unit away; z
        // stands on 0 and runs to 1.
        let face = rng.pick(&faces);
        let offset = face * (1.0 + gap);
        let (graph, root) = boxes(a, b, offset);
        match compile(&graph, root) {
            Ok((mesh, closure)) => {
                assert_eq!(
                    component_count(&mesh),
                    1,
                    "case {case} (gap {gap:e}, face {face:?}): bodies left in contact"
                );
                assert_eq!(closure, MeshClosure::Solid, "case {case}");
            }
            Err(error) => assert!(
                is_unmerged(&error) || is_pinch(&error),
                "case {case} (gap {gap:e}, face {face:?}): {error}"
            ),
        }
    }
}
