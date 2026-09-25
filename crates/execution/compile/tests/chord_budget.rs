//! The chord budget is its own knob, separate from the linear tolerance (#165).
//!
//! The reference compiler used to flatten every curve with a chord error
//! equal to `Tolerance::linear()`. That is a coincidence tolerance: at
//! `Tolerance::MILLIMETRE` it leaves a 5 mm arc a handful of chords, and the
//! meshed area of small profiles is off by percent. `ExecutionOptions::
//! with_chord_error` sets the budget; without it the old behaviour holds.
//!
//! Every oracle is closed form: a chord polygon inscribed in a circle of
//! radius `r` loses at most `2 * perimeter * sagitta`-order area, and a
//! sagitta bound `s` caps the relative area error of a disc at about `2s/r`.

use std::f64::consts::PI;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Tolerance, Transform3, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::Primitive;
use axiolid_profile::{CircleProfile, Profile};

/// A 5 mm radius disc extruded 1 m: the issue's "small arc" case.
const RADIUS: f64 = 0.005;
const DEPTH: f64 = 1.0;

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
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

/// A graph with one extruded circle, and its root.
fn rod() -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let circle = b
        .push(GeometryNode::Profile(Profile::Circle(CircleProfile {
            radius: RADIUS,
            thickness: None,
        })))
        .unwrap();
    let rod = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: circle,
            direction: Vec3::Z,
            depth: DEPTH,
        }))
        .unwrap();
    (b.finish(vec![rod]).unwrap(), rod)
}

fn compile(graph: &GeometryGraph, root: NodeId, options: &ExecutionOptions) -> TriMesh {
    compiler()
        .compile_mesh(graph, root, options)
        .expect("a rod compiles")
}

fn relative_error(got: f64, exact: f64) -> f64 {
    (got - exact).abs() / exact
}

#[test]
fn without_a_chord_budget_the_linear_tolerance_still_governs() {
    // The default is unchanged: a 5 mm disc at a 1 mm budget is coarse,
    // percent-level short. This pins that nothing moved for callers that
    // did not ask.
    let (graph, root) = rod();
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    assert_eq!(options.chord_error(), None);
    let exact = PI * RADIUS * RADIUS * DEPTH;
    let coarse = volume(&compile(&graph, root, &options));
    assert!(
        relative_error(coarse, exact) > 1e-2,
        "default budget should still be the 1 mm linear tolerance: \
         volume {coarse}, exact {exact}"
    );
}

#[test]
fn a_chord_budget_tightens_a_small_arc_to_its_stated_bound() {
    // The issue's done-when: r = 5 mm at Tolerance::MILLIMETRE, area within
    // a stated bound. A sagitta of s = 1 um on r = 5 mm bounds the inscribed
    // polygon's area deficit by about 2s/r = 4e-4 relative.
    let (graph, root) = rod();
    let chord = 1e-6;
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(chord)
        .expect("a positive finite budget");
    let exact = PI * RADIUS * RADIUS * DEPTH;
    let fine = volume(&compile(&graph, root, &options));
    let bound = 2.0 * chord / RADIUS;
    assert!(
        fine <= exact * (1.0 + 1e-12),
        "an inscribed polygon cannot exceed the disc: {fine} > {exact}"
    );
    assert!(
        relative_error(fine, exact) <= bound,
        "volume {fine}, exact {exact}, relative error {} above bound {bound}",
        relative_error(fine, exact)
    );
}

#[test]
fn a_finer_budget_never_meshes_worse() {
    let (graph, root) = rod();
    let exact = PI * RADIUS * RADIUS * DEPTH;
    let mut previous = f64::INFINITY;
    for chord in [1e-3, 1e-4, 1e-5, 1e-6] {
        let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
            .with_chord_error(chord)
            .unwrap();
        let error = relative_error(volume(&compile(&graph, root, &options)), exact);
        assert!(
            error <= previous,
            "chord {chord}: error {error} worse than the coarser budget's {previous}"
        );
        previous = error;
    }
}

#[test]
fn an_invalid_chord_budget_is_refused_at_construction() {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    for bad in [0.0, -1e-6, f64::NAN, f64::INFINITY] {
        assert!(
            options.clone().with_chord_error(bad).is_none(),
            "{bad} must not be accepted as a chord budget"
        );
    }
}

#[test]
fn a_scaled_instance_keeps_the_world_chord_budget() {
    // A 1000x instance of a 5 um rod is the same 5 mm rod in world space.
    // The budget is a world distance, so the local budget must shrink by
    // the scale: the instanced volume must meet the same bound as a rod
    // authored at full size.
    let mut b = GeometryGraphBuilder::new();
    let circle = b
        .push(GeometryNode::Profile(Profile::Circle(CircleProfile {
            radius: RADIUS / 1000.0,
            thickness: None,
        })))
        .unwrap();
    let small = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: circle,
            direction: Vec3::Z,
            depth: DEPTH / 1000.0,
        }))
        .unwrap();
    let scaled = b
        .push(GeometryNode::Instance(Instance {
            source: small,
            transform: Transform3::from_scale(Vec3::splat(1000.0)),
        }))
        .unwrap();
    let graph = b.finish(vec![scaled]).unwrap();
    let chord = 1e-6;
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(chord)
        .unwrap();
    let exact = PI * RADIUS * RADIUS * DEPTH;
    let got = volume(&compile(&graph, scaled, &options));
    let bound = 2.0 * chord / RADIUS;
    assert!(
        relative_error(got, exact) <= bound,
        "instanced volume {got}, exact {exact}, relative error {} above {bound}",
        relative_error(got, exact)
    );
}

#[test]
fn a_primitive_cylinder_honours_the_chord_budget() {
    // CSG primitives take one tolerance; the budget must reach them too.
    let mut b = GeometryGraphBuilder::new();
    let cylinder = b
        .push(GeometryNode::Primitive(Primitive::Cylinder {
            radius: RADIUS,
            height: DEPTH,
        }))
        .unwrap();
    let graph = b.finish(vec![cylinder]).unwrap();
    let exact = PI * RADIUS * RADIUS * DEPTH;
    let coarse = volume(&compile(
        &graph,
        cylinder,
        &ExecutionOptions::new(Tolerance::MILLIMETRE),
    ));
    let chord = 1e-6;
    let fine = volume(&compile(
        &graph,
        cylinder,
        &ExecutionOptions::new(Tolerance::MILLIMETRE)
            .with_chord_error(chord)
            .unwrap(),
    ));
    assert!(
        relative_error(fine, exact) < relative_error(coarse, exact),
        "budget ignored: coarse {coarse}, fine {fine}, exact {exact}"
    );
    assert!(
        relative_error(fine, exact) <= 2.0 * chord / RADIUS,
        "cylinder volume {fine}, exact {exact}"
    );
}
