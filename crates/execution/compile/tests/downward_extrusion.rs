//! Exact compilation of extrusions against the profile normal (#275).
//!
//! An opening cut down from a slab's top (`ExtrudedDirection (0, 0, -1)`)
//! is the same solid as one cut up from its bottom. The exact compiler used
//! to refuse it as a "non-forward planar extrusion", so a boolean with such
//! an operand had no exact result and its mesh deviation stayed unbounded.
//! These cases pin the fix at the compiler: the exact result exists, is a
//! sound outward solid of the closed-form volume, and the mesh is certified
//! against it the same way whichever way the opening was authored.

use std::f64::consts::PI;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Plane3, Point3, Tolerance, Transform3, Vec3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationReport, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn push(b: &mut GeometryGraphBuilder, node: GeometryNode) -> NodeId {
    b.push(node).expect("a valid node")
}

fn rectangle(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn extrusion(
    b: &mut GeometryGraphBuilder,
    profile: Profile,
    direction: Vec3,
    depth: f64,
    at: Vec3,
) -> NodeId {
    let profile = push(b, GeometryNode::Profile(profile));
    let solid = push(
        b,
        GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction,
            depth,
        }),
    );
    push(
        b,
        GeometryNode::Instance(Instance {
            source: solid,
            transform: Transform3::from_translation(at),
        }),
    )
}

fn boolean(
    b: &mut GeometryGraphBuilder,
    left: NodeId,
    right: NodeId,
    operator: BooleanOperator,
) -> NodeId {
    push(
        b,
        GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }),
    )
}

fn exact(graph: &GeometryGraph, root: NodeId) -> ExactBRep {
    let solid = ReferenceExactCompiler::new()
        .compile_exact(graph, root, &options())
        .expect("the exact compiler builds it");
    let health = geometric_audit(&solid, options().tolerance());
    assert!(health.is_consistent(), "{:?}", health.defects());
    assert!(axiolid_topology::audit_brep(solid.topology()).is_closed_manifold());
    solid
}

fn volume(solid: &ExactBRep) -> f64 {
    exact_properties(solid, options().tolerance())
        .expect("measurable")
        .signed_volume
}

fn deviation(graph: &GeometryGraph, root: NodeId) -> DeviationReport {
    let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean)
        .compile_mesh_with_deviation(graph, root, &options())
        .expect("the mesh compiles");
    report
}

fn unbounded(report: &DeviationReport) -> Vec<&'static str> {
    report
        .contributions
        .iter()
        .filter_map(|c| match c.bound {
            DeviationBound::Unbounded(reason) => Some(reason),
            _ => None,
        })
        .collect()
}

/// The issue's slab: 6 x 9 x 0.25 less a round hole of radius 0.5 and a
/// 2 x 1 rectangle, cut up from the bottom or down from the top.
fn slab(down: bool) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let slab = extrusion(
        &mut b,
        rectangle(6.0, 9.0),
        Vec3::Z,
        0.25,
        Vec3::new(3.0, 4.5, 0.0),
    );
    let hole = extrusion(
        &mut b,
        Profile::Circle(CircleProfile {
            radius: 0.5,
            thickness: None,
        }),
        Vec3::Z,
        0.25,
        Vec3::new(2.0, 6.0, 0.0),
    );
    let (direction, z) = if down {
        (-Vec3::Z, 0.25)
    } else {
        (Vec3::Z, 0.0)
    };
    let opening = extrusion(
        &mut b,
        rectangle(2.0, 1.0),
        direction,
        0.25,
        Vec3::new(4.0, 2.0, z),
    );
    let less_hole = boolean(&mut b, slab, hole, BooleanOperator::Difference);
    let root = boolean(&mut b, less_hole, opening, BooleanOperator::Difference);
    (b.finish(vec![root]).expect("graph"), root)
}

/// The repro: the same slab certifies the same bound whichever way its
/// rectangular opening was extruded.
#[test]
fn a_slab_opening_cut_down_certifies_like_one_cut_up() {
    let (up_graph, up_root) = slab(false);
    let (down_graph, down_root) = slab(true);
    let up = deviation(&up_graph, up_root);
    let down = deviation(&down_graph, down_root);
    assert!(unbounded(&up).is_empty(), "{:?}", unbounded(&up));
    assert!(unbounded(&down).is_empty(), "{:?}", unbounded(&down));
    let (up_bound, down_bound) = (
        up.bound.expect("certified up"),
        down.bound.expect("certified down"),
    );
    assert!(
        down_bound <= options().tolerance().linear(),
        "bound {down_bound} above the requested chord error"
    );
    assert!(
        (up_bound - down_bound).abs() <= 1e-9,
        "bounds differ: up {up_bound}, down {down_bound}"
    );
}

/// The exact result of the downward cut is the same solid, of closed-form
/// volume `6 * 9 * 0.25 - (pi 0.25 + 2) * 0.25`.
#[test]
fn a_slab_opening_cut_down_is_exact_with_the_closed_form_volume() {
    let expected = (54.0 - PI * 0.25 - 2.0) * 0.25;
    let (up_graph, up_root) = slab(false);
    let (down_graph, down_root) = slab(true);
    let up = exact(&up_graph, up_root);
    let down = exact(&down_graph, down_root);
    for (what, solid) in [("up", &up), ("down", &down)] {
        let measured = volume(solid);
        assert!(
            (measured - expected).abs() <= 1e-9 * expected,
            "{what}: volume {measured}, expected {expected}"
        );
    }
    assert_eq!(
        up.topology().faces().len(),
        down.topology().faces().len(),
        "the two cuts build the same faces"
    );
}

/// A slab extruded down from its top, less an opening cut down too: every
/// operand points against its profile normal.
#[test]
fn a_downward_slab_less_a_downward_opening_is_exact() {
    let mut b = GeometryGraphBuilder::new();
    let slab = extrusion(
        &mut b,
        rectangle(6.0, 9.0),
        -Vec3::Z,
        0.25,
        Vec3::new(3.0, 4.5, 0.25),
    );
    let opening = extrusion(
        &mut b,
        rectangle(2.0, 1.0),
        -Vec3::Z,
        0.25,
        Vec3::new(4.0, 2.0, 0.25),
    );
    let root = boolean(&mut b, slab, opening, BooleanOperator::Difference);
    let graph = b.finish(vec![root]).expect("graph");
    let solid = exact(&graph, root);
    let expected = (54.0 - 2.0) * 0.25;
    let measured = volume(&solid);
    assert!(
        (measured - expected).abs() <= 1e-9 * expected,
        "volume {measured}, expected {expected}"
    );
    for vertex in solid.topology().vertices() {
        let z = vertex.position.z;
        assert!(z == 0.0 || z == 0.25, "vertex at z = {z}");
    }
    let report = deviation(&graph, root);
    assert!(unbounded(&report).is_empty(), "{:?}", unbounded(&report));
    assert!(report.bound.is_some());
}

/// A grille-like round extrusion along `-z`, clipped by an inclined
/// half-space: exact, of closed-form volume, and its mesh certified.
#[test]
fn a_round_extrusion_down_clipped_by_a_half_space_is_certified() {
    let (r, depth) = (0.5, 1.0);
    let mut b = GeometryGraphBuilder::new();
    let column = extrusion(
        &mut b,
        Profile::Circle(CircleProfile {
            radius: r,
            thickness: None,
        }),
        -Vec3::Z,
        depth,
        Vec3::ZERO,
    );
    // Above the plane z = -0.5 - 0.2 x, which crosses the whole column.
    let above = push(
        &mut b,
        GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: Point3::new(0.0, 0.0, -0.5),
                normal: Vec3::new(0.2, 0.0, 1.0),
            },
            agreement: true,
        }),
    );
    let root = boolean(&mut b, column, above, BooleanOperator::Difference);
    let graph = b.finish(vec![root]).expect("graph");

    // The kept height is linear in x, so the volume is the disc's area
    // times the height over its centre.
    let expected = PI * r * r * 0.5;
    let measured = volume(&exact(&graph, root));
    assert!(
        (measured - expected).abs() <= 1e-9 * expected,
        "volume {measured}, expected {expected}"
    );
    let report = deviation(&graph, root);
    assert!(unbounded(&report).is_empty(), "{:?}", unbounded(&report));
    assert!(report.bound.is_some(), "certified bound");
}
